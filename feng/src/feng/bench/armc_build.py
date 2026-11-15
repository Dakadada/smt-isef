"""Arm C bind->act gate: parse Holgersen S1/S4, lift the direct-bound positive
loci hg19->GRCh38, draw a bounded GENCODE v50 negative pool, and assemble a
SITE-LEVEL binding-feature table with bind-no-act labels.

REDESIGN (liftover + site-level features, 2026-09)
--------------------------------------------------
The S4 exon-interval coordinates are hg19/GRCh37, not GRCh38. Verified: of the
204 unique direct-bound positive intervals (numeric edit_dist<=5), exact/overlap
match to GENCODE v50 (GRCh38) exons is only ~10-17%, but after hg19->GRCh38
liftover 100% (204/204) lift cleanly. The earlier "positive == annotated GENCODE
exon_key" join therefore produced a near-empty positive class purely from an
assembly frame mismatch. The corrected pipeline:

1. POSITIVE CLASS = direct-bound-and-acted: S4 rows whose ``edit_dist`` is a
   NUMERIC value <= ``max_edit`` (a close antisense sequence match). Rows the
   paper marks "background" (NaN edit_dist) are INDIRECT effects -- not bound --
   and are EXCLUDED (counted separately, never labeled). Each qualifying interval
   is LIFTED hg19->GRCh38 (:mod:`feng.rbp.liftover`, chain from
   ``FENG_ARMC_CHAIN``); unmapped / multi-mapped intervals are dropped and
   counted. The POSITIVE TARGET is the lifted interval's OWN GRCh38 region
   sequence (fetched from the reference), NOT an annotated GENCODE exon -- the
   spliced off-target site is frequently novel / boundary-drifted relative to the
   v50 annotation, so requiring an exact annotated exon caps recovery at ~51%
   and discards real positives. The region sequence itself is the target.
2. NEGATIVE CLASS = bound-but-not-acted GENCODE v50 exons, sampled to stay
   tractable. A single deterministic RANDOM BACKGROUND POOL of ``pool_size``
   annotated exons is drawn once (fixed seed, label-blind, shared across all
   SBOs); per SBO, any pool exon overlapping one of that SBO's positive loci is
   excluded. A pool exon the SBO binds (edit<=5) but does not act on is
   actionable=0. Positives and negatives thus pass the identical binding gate.
3. FEATURES are SITE / DUPLEX level, computed by an IDENTICAL code path for both
   labels, to avoid a positives-are-longer-regions confound. For each
   (SBO, target sequence) the scan finds the SBO's BEST binding window (lowest
   edit distance) of ``revcomp(SBO)`` in that target; the feature vector is
   ``edit_distance``, ``exact_match_length``, ``duplex_dg`` (IntaRNA on the SBO x
   that best window), and ``binding_window_gc`` (GC of the best window). The
   whole-region ``exon_length`` and whole-region ``gc_fraction`` are DELIBERATELY
   NOT features: region-positives and exon-negatives differ systematically in
   length/composition, which is a substrate confound, not binding. A target where
   ``revcomp(SBO)`` has no hit at edit<=max_edit is 'not bound' and is excluded
   from BOTH classes (the gate is conditioned on binding). ``baseline_psi`` /
   ``case_psi`` / ``control_psi`` are never features (action-adjacent).
4. ``max_edit`` stays 5 (positives extend to edit 4-5). The low-complexity SBO
   exclusions are retained.

RESIDUAL LIMITATION (flagged, acceptable)
-----------------------------------------
Negatives are sampled from ANNOTATED GENCODE exons, some of which are not
expressed in the assayed lines (HEK293T / HepG2): such an exon cannot be acted
on for a trivial expression reason rather than a binding one, so a bound-but-not-
expressed exon is labeled 0 for a reason the binding features cannot capture.
Because every feature is binding-only, the model cannot exploit expression to
"cheat"; the effect is label noise biased toward NO-GO (a bound-and-would-act
exon that is simply silent in these lines is mislabeled negative). That is
CONSERVATIVE (it can only depress, never inflate, the model's apparent lift over
the baselines) and is recorded here rather than corrected.

Leakage discipline (read this before touching the label join)
---------------------------------------------------------------
Every binding feature is computed by the SAME code path for label=1 and label=0
rows -- the scan (:func:`scan_sbo`), the duplex-energy call (:func:`duplex_energy`
on the best window), and the window-GC derivation never see the label and run
identically for every bound target. In particular, S4's own ``rnaplex_dg`` /
``rnacofold_dg`` / ``case_psi`` / ``control_psi`` columns are read only inside the
label join to decide WHICH loci are positive (via ``id -> dg_id`` from S1 and
S4's ``edit_dist`` / interval) -- those paper-reported values are never copied
into ``duplex_dg`` or any other feature column; ``duplex_dg`` and
``edit_distance`` are recomputed here identically for both classes.

Discovered simplification vs. the original brief: pdfplumber's ``extract_tables()``
already reassembles S4's two-line-wrapped exon-interval cells into one logical row
(every row returns 11 cells, ``len(rows) == 724``, 75 unique IDs), so the raw-text
regex reconstruction the brief described is unnecessary.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Iterable, Iterator, Mapping

import numpy as np
import pandas as pd

from feng.bench.scores import (
    _fetch_transcript_sequences,
    exon_window_specs,
    find_intarna,
    intarna_min_energy,
    IntaRNAError,
    reverse_complement,
)
from feng.bench.strata import positive_interval_to_exon_key
from feng.rbp.liftover import ChainMap, _combine_strand, load_chain_file
from feng.scan.seed_extend import HybridizationHit, seed_and_extend

# --- default paths (repo-root-relative unless a scratch path is given) ------

S1_PDF_DEFAULT = Path("feng/sj-pdf-1-nat-10.1089_nat.2020.0921-2.pdf")
S4_PDF_DEFAULT = Path("feng/sj-pdf-4-nat-10.1089_nat.2020.0921.pdf")
SBO_TSV_DEFAULT = Path("feng/artifacts/phase0_arm_c_probe/sbo_sequences.tsv")
GTF_DEFAULT = Path("/home/dashii/orcd/scratch/feng/gencode/v50/gencode.v50.primary_assembly.annotation.gtf")
GENOME_FASTA_DEFAULT = Path("/home/dashii/orcd/scratch/feng/gencode/v50/GRCh38.primary_assembly.genome.fa")
FEATURES_OUT_DEFAULT = Path("feng/artifacts/phase0_arm_c_probe/armc_features.parquet")

# UCSC hg19->GRCh38 chain that lifts the hg19/GRCh37 S4 positive loci into the
# GRCh38 (GENCODE v50) frame the negatives/reference live in. Overridable via
# FENG_ARMC_CHAIN so the sbatch and tests can point at their staged copy.
CHAIN_DEFAULT = Path("/home/dashii/orcd/scratch/feng/liftover/hg19ToHg38.over.chain.gz")

DEFAULT_MAX_EDIT = 5

_CHAIN_CACHE: dict[str, ChainMap] = {}


def resolve_chain_path(path: str | Path | None = None) -> Path:
    """Chain path from an explicit arg, else ``FENG_ARMC_CHAIN``, else the default."""
    if path is not None:
        return Path(path)
    return Path(os.environ.get("FENG_ARMC_CHAIN", str(CHAIN_DEFAULT)))


def load_armc_chain(path: str | Path | None = None) -> ChainMap:
    """Load (and memoize) the hg19->GRCh38 :class:`ChainMap` used for S4 positives."""
    resolved = str(resolve_chain_path(path))
    if resolved not in _CHAIN_CACHE:
        _CHAIN_CACHE[resolved] = load_chain_file(resolved)
    return _CHAIN_CACHE[resolved]

# Bounded-pool design (see module docstring). The pool is drawn once with this
# fixed seed and shared across every SBO; positives are added on top so both
# label classes pass the identical binding gate. ~20k keeps a per-SBO scan near
# ~1 h at max_edit=5 (the canary scaled to ~3.5 h at 81k), while leaving a large
# bind-no-act negative class per SBO.
DEFAULT_POOL_SIZE = 20000
DEFAULT_POOL_SEED = 20260903

# Symmetric bp pad added to EVERY target region (positives AND pool negatives)
# before fetching its sequence. The S4 positive interval is the annotated
# off-target exon, but the 18-20 nt antisense site can straddle the exon/intron
# boundary or drift a few bp on liftover, so an unpadded fetch misses (or
# degrades the edit distance of) the true site for edge-case positives. A small
# pad guarantees containment; applied IDENTICALLY to both classes it is fetch
# geometry, not a feature (all features are computed from the best binding WINDOW
# inside the fetched sequence, which is pad-invariant).
DEFAULT_TARGET_PAD = 30

# --- S1 parsing ---------------------------------------------------------

_S1_COLUMNS = [
    "id", "dg_id", "cell", "cases", "controls", "sequence", "category",
    "splicing_events", "expr_events",
]


def parse_s1(pdf_path: Path = S1_PDF_DEFAULT) -> pd.DataFrame:
    """Parse Holgersen Supp Table S1 (SBO design table).

    Returns one row per (dg_id, cell-type/run) experiment, columns
    ``id, dg_id, cell, sequence, category``. ``id`` is the paper's unique
    per-experiment numeric key.
    """
    import pdfplumber

    rows: list[list] = []
    with pdfplumber.open(pdf_path) as pdf:
        for page in pdf.pages:
            for table in page.extract_tables():
                rows.extend(table)
    if not rows:
        raise ValueError(f"no tables extracted from {pdf_path}")
    header = rows[0]
    if header[0] != "ID":
        raise ValueError(f"unexpected S1 header: {header}")
    data = [r for r in rows[1:] if r and r[0] not in (None, "ID")]
    if not data:
        raise ValueError(f"no data rows parsed from S1 {pdf_path}")
    df = pd.DataFrame(data, columns=_S1_COLUMNS)
    df["id"] = df["id"].astype(int)
    df["dg_id"] = df["dg_id"].astype(str).str.strip()
    df["cell"] = df["cell"].astype(str).str.strip()
    df["sequence"] = df["sequence"].astype(str).str.replace("\n", "").str.strip().str.upper()
    df["category"] = df["category"].astype(str).str.replace("\n", " ").str.strip()
    if df["id"].duplicated().any():
        dupes = sorted(df.loc[df["id"].duplicated(), "id"].tolist())
        raise ValueError(f"S1 'ID' column is not unique: duplicate ids {dupes}")
    return df[["id", "dg_id", "cell", "sequence", "category"]].reset_index(drop=True)


# --- S4 parsing ---------------------------------------------------------

_S4_ROW_WIDTH = 11
_INTERVAL_RE = re.compile(r"^(chr\w+):(\d+)-(\d+):([+-])$")


def _parse_interval(raw: str | None) -> tuple[str, int, int, str]:
    text = (raw or "").replace("\n", "")
    match = _INTERVAL_RE.match(text.strip())
    if not match:
        raise ValueError(f"unparseable S4 exon interval: {raw!r}")
    chrom, start, end, strand = match.groups()
    return chrom, int(start), int(end), strand


def _to_float(raw) -> float | None:
    if raw is None:
        return None
    text = str(raw).strip()
    if text == "" or text.lower() == "background":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def parse_s4(pdf_path: Path = S4_PDF_DEFAULT) -> pd.DataFrame:
    """Parse Holgersen Supp Table S4 (action-positive spliced off-target sites).

    Returns columns ``id, chrom, start, end, strand, case_psi, control_psi,
    rnaplex_dg, rnacofold_dg, edit_dist, exact_match``. ``start``/``end`` are
    the paper's own 1-based inclusive coordinates, unconverted. ``edit_dist`` is
    ``None`` for rows the paper marks "background" (an INDIRECT off-target, not a
    close sequence match).
    """
    import pdfplumber

    rows: list[list] = []
    with pdfplumber.open(pdf_path) as pdf:
        for page in pdf.pages:
            for table in page.extract_tables():
                rows.extend(table)
    if not rows:
        raise ValueError(f"no tables extracted from {pdf_path}")
    header = rows[0]
    if header[0] != "ID":
        raise ValueError(f"unexpected S4 header: {header}")
    data = [r for r in rows[1:] if r and r[0] not in (None, "ID")]
    if not data:
        raise ValueError(f"no data rows parsed from S4 {pdf_path}")
    records = []
    for r in data:
        if len(r) != _S4_ROW_WIDTH:
            raise ValueError(f"unexpected S4 row width {len(r)} (want {_S4_ROW_WIDTH}): {r}")
        (raw_id, raw_interval, raw_case, raw_ctrl, edit_body, _mism_body,
         _edit_pad, _edit_intron, exact_match, rnaplex, rnacofold) = r
        chrom, start, end, strand = _parse_interval(raw_interval)
        records.append({
            "id": int(raw_id),
            "chrom": chrom,
            "start": start,
            "end": end,
            "strand": strand,
            "case_psi": _to_float(raw_case),
            "control_psi": _to_float(raw_ctrl),
            "rnaplex_dg": _to_float(rnaplex),
            "rnacofold_dg": _to_float(rnacofold),
            "edit_dist": _to_float(edit_body),
            "exact_match": _to_float(exact_match),
        })
    return pd.DataFrame.from_records(records, columns=[
        "id", "chrom", "start", "end", "strand", "case_psi", "control_psi",
        "rnaplex_dg", "rnacofold_dg", "edit_dist", "exact_match",
    ])


# --- GENCODE candidate universe (genome-wide) -------------------------------


def _exon_key(chrom: str, start_0base: int, end: int, strand: str) -> str:
    return f"{chrom}:{start_0base}-{end}:{strand}"


def _parse_exon_key(key: str) -> tuple[str, int, int, str]:
    """Inverse of :func:`_exon_key`. chrom carries no ':'; span is 'start-end'."""
    chrom, span, strand = key.split(":")
    start_0base, end = span.split("-")
    return chrom, int(start_0base), int(end), strand


def iter_gencode_exon_keys(gtf_path: Path) -> Iterator[str]:
    """Stream ``chrom:start_0base-end:strand`` keys for every ``exon`` feature.

    Reads the ONE GTF file line by line (memory-bounded; never recurses the
    filesystem). GTF is 1-based inclusive, so ``start_0base = start - 1`` and
    ``end`` is unchanged -- matching :func:`feng.bench.strata.positive_interval_to_exon_key`
    and the detectable-universe ``exon_key`` frame exactly. Not deduplicated;
    callers dedup (see :func:`load_gencode_exon_keys`).
    """
    with Path(gtf_path).open("rt") as fh:
        for line in fh:
            # Fast reject: the feature column is the only place '\texon\t'
            # occurs (attributes are space-delimited, tab-free).
            if "\texon\t" not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or fields[2] != "exon":
                continue
            chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            yield _exon_key(chrom, start - 1, end, strand)


def load_gencode_exon_keys(gtf_path: Path = GTF_DEFAULT) -> np.ndarray:
    """Return the sorted array of unique GENCODE exon keys (deduped intervals)."""
    keys = set(iter_gencode_exon_keys(gtf_path))
    if not keys:
        raise ValueError(f"no exon features parsed from {gtf_path}")
    return np.array(sorted(keys), dtype=object)


def exons_df_from_keys(keys: Iterable[str]) -> pd.DataFrame:
    """Materialize the exon-identity DataFrame (for windowing/fetch) from keys."""
    records = []
    for key in keys:
        chrom, start_0base, end, strand = _parse_exon_key(key)
        records.append({
            "exon_key": key, "chrom": chrom,
            "exon_start_0base": start_0base, "exon_end": end, "strand": strand,
        })
    return pd.DataFrame.from_records(records, columns=[
        "exon_key", "chrom", "exon_start_0base", "exon_end", "strand",
    ])


def sample_background_pool(
    all_keys: np.ndarray, *, pool_size: int = DEFAULT_POOL_SIZE,
    seed: int = DEFAULT_POOL_SEED,
) -> list[str]:
    """Draw ONE deterministic, label-blind background pool of exon keys.

    Sorted-then-seeded so the draw is reproducible for a given
    ``(GENCODE universe, pool_size, seed)`` regardless of GTF iteration order,
    and shared identically across every SBO. Returns all keys (sorted) when the
    universe is smaller than ``pool_size``.
    """
    keys = np.array(sorted(all_keys), dtype=object)
    if pool_size >= len(keys):
        return list(keys)
    rng = np.random.RandomState(seed)
    idx = rng.choice(len(keys), size=pool_size, replace=False)
    return sorted(keys[idx].tolist())


# --- positives (direct-bound-and-acted) -------------------------------------


def _iter_direct_positive_rows(
    s1: pd.DataFrame, s4: pd.DataFrame, max_edit: int,
) -> Iterator[tuple[str, tuple[str, int, int, str]]]:
    """Yield ``(dg_id, (chrom, start, end, strand))`` for every direct-bound S4 row.

    Direct-bound = a NUMERIC ``edit_dist`` <= ``max_edit`` (a close antisense
    match). Rows the paper marks "background" (NaN edit_dist) are INDIRECT and
    skipped. ``start``/``end`` are S4's own 1-based inclusive coordinates.
    """
    id_to_dg = s1.set_index("id")["dg_id"].to_dict()
    for row in s4.itertuples(index=False):
        edit = row.edit_dist
        if edit is None or pd.isna(edit) or float(edit) > max_edit:
            continue
        dg = id_to_dg.get(row.id)
        if dg is None:
            continue
        yield dg, (str(row.chrom), int(row.start), int(row.end), str(row.strand))


def _lift_positive_interval(
    interval: tuple[str, int, int, str], chain: ChainMap | None,
) -> tuple[str, str | None]:
    """Lift one 1-based-inclusive S4 positive interval to a GRCh38 ``exon_key``.

    Returns ``(status, key)`` where ``status`` is ``"lifted"``, ``"unmapped"`` or
    ``"multimapped"``. When ``chain is None`` the interval is taken to be already
    on the target frame and only the 1-based-inclusive -> 0-based-half-open
    conversion is applied (used by the synthetic unit tests). With a chain, the
    interval is converted to 0-based half-open (``start-1`` .. ``end``) and lifted;
    a strand-flipping (``-``) chain flips the feature strand
    (:func:`feng.rbp.liftover._combine_strand`), matching the eCLIP liftover path.
    """
    chrom, start, end, strand = interval
    if chain is None:
        try:
            return "lifted", positive_interval_to_exon_key(chrom, start, end, strand)
        except ValueError:
            return "unmapped", None
    results = chain.lift_interval(chrom, start - 1, end)
    if not results:
        return "unmapped", None
    if len(results) > 1:
        return "multimapped", None
    new_chrom, new_start, new_end, chain_strand = results[0]
    return "lifted", _exon_key(new_chrom, new_start, new_end, _combine_strand(strand, chain_strand))


def positive_exon_keys_by_dg(
    s1: pd.DataFrame, s4: pd.DataFrame, *, max_edit: int = DEFAULT_MAX_EDIT,
    chain: ChainMap | None = None,
) -> dict[str, set[str]]:
    """Map each ``dg_id`` to the GRCh38 positive-locus ``exon_key``\\ s it directly
    bound AND acted on.

    Join path: S4.id -> S1.id (unique) -> S1.dg_id, then S4's interval is LIFTED
    hg19->GRCh38 (:func:`_lift_positive_interval`) and keyed on its OWN lifted
    coordinates -- NOT snapped to an annotated GENCODE exon. ONLY S4 rows with a
    NUMERIC ``edit_dist`` <= ``max_edit`` qualify. Unmapped / multi-mapped
    intervals are dropped. A dg_id tested in multiple cell types contributes the
    union across its S1 ids (binding features are sequence-only).
    """
    out: dict[str, set[str]] = {}
    for dg, interval in _iter_direct_positive_rows(s1, s4, max_edit):
        _status, key = _lift_positive_interval(interval, chain)
        if key is None:
            continue
        out.setdefault(dg, set()).add(key)
    return out


def positive_lift_coverage(
    s1: pd.DataFrame, s4: pd.DataFrame, *, max_edit: int = DEFAULT_MAX_EDIT,
    chain: ChainMap | None = None,
) -> dict:
    """Report the lift recovery over the UNIQUE direct-positive intervals.

    ``lift_fraction = n_lifted / n_direct_positive_intervals``; the build-exonome
    tripwire refuses to proceed below 90% (a wrong / absent chain leaves the
    positive class empty). Counts are over deduplicated intervals so a locus hit
    by several SBOs is counted once.
    """
    intervals = {interval for _dg, interval in _iter_direct_positive_rows(s1, s4, max_edit)}
    n_lifted = n_unmapped = n_multimapped = 0
    for interval in intervals:
        status, key = _lift_positive_interval(interval, chain)
        if status == "lifted" and key is not None:
            n_lifted += 1
        elif status == "multimapped":
            n_multimapped += 1
        else:
            n_unmapped += 1
    n = len(intervals)
    return {
        "n_direct_positive_intervals": n,
        "n_lifted": n_lifted,
        "n_unmapped": n_unmapped,
        "n_multimapped": n_multimapped,
        "lift_fraction": (n_lifted / n) if n else float("nan"),
    }


def s4_direct_vs_indirect(
    s4: pd.DataFrame, *, max_edit: int = DEFAULT_MAX_EDIT,
) -> tuple[int, int]:
    """Count (direct-bound rows: numeric edit_dist<=max_edit, indirect rows: NaN)."""
    edit = pd.to_numeric(s4["edit_dist"], errors="coerce")
    n_direct = int(((edit <= max_edit) & edit.notna()).sum())
    n_indirect = int(edit.isna().sum())
    return n_direct, n_indirect


# --- exonome ------------------------------------------------------------


def _read_fasta(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    with Path(path).open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out[name] = "".join(chunks)
                name = line[1:]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            out[name] = "".join(chunks)
    return out


def build_exonome(
    exons: pd.DataFrame,
    genome_fasta: Path | None = None,
    cache_path: Path | None = None,
    *,
    pad: int = 0,
    samtools: Path | None = None,
) -> dict[str, str]:
    """Fetch (or load cached) transcript-strand sequence for each exon in
    ``exons`` (columns exon_key/chrom/exon_start_0base/exon_end/strand),
    optionally padded by ``pad`` bp on each side (clamped at the contig start).

    Guard: returns the cached FASTA if ``cache_path`` exists and is nonempty --
    callers running many parallel tasks MUST serialize the first build
    externally (e.g. flock); this function only checks-then-writes.
    """
    if cache_path is None:
        raise ValueError("cache_path is required")
    cache_path = Path(cache_path)
    if cache_path.is_file() and cache_path.stat().st_size > 0:
        return _read_fasta(cache_path)
    if genome_fasta is None:
        raise ValueError("genome_fasta is required to build a fresh exonome cache")

    required = {"exon_key", "chrom", "exon_start_0base", "exon_end", "strand"}
    missing = required - set(exons.columns)
    if missing:
        raise ValueError(f"exon table missing columns: {sorted(missing)}")
    exons = exons.drop_duplicates(subset="exon_key")

    specs = exon_window_specs(exons, window=pad)
    sequences = _fetch_transcript_sequences(specs, Path(genome_fasta), samtools=samtools)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    tmp = cache_path.with_name(cache_path.name + f".tmp.{os.getpid()}")
    with tmp.open("w") as fh:
        for key, seq in sequences.items():
            fh.write(f">{key}\n{seq}\n")
    tmp.replace(cache_path)
    return sequences


def build_target_exonome(
    *,
    gtf_path: Path,
    genome_fasta: Path,
    cache_path: Path,
    s1: pd.DataFrame,
    s4: pd.DataFrame,
    chain: ChainMap | None,
    pool_size: int = DEFAULT_POOL_SIZE,
    pool_seed: int = DEFAULT_POOL_SEED,
    max_edit: int = DEFAULT_MAX_EDIT,
    pad: int = DEFAULT_TARGET_PAD,
    samtools: Path | None = None,
) -> tuple[dict[str, str], dict]:
    """Build (or load) the shared scan target = background pool UNION the LIFTED
    direct-positive loci, and return ``(exonome, coverage_stats)``.

    The positive targets are the lifted GRCh38 regions' OWN sequences (fetched by
    their lifted coordinates), not annotated GENCODE exons. ``coverage_stats``
    records the lift recovery (>=90% tripwire): how many unique direct-bound S4
    intervals lift cleanly hg19->GRCh38.
    """
    positives_by_dg = positive_exon_keys_by_dg(s1, s4, max_edit=max_edit, chain=chain)
    positive_keys: set[str] = set().union(*positives_by_dg.values()) if positives_by_dg else set()

    coverage = positive_lift_coverage(s1, s4, max_edit=max_edit, chain=chain)

    all_keys = load_gencode_exon_keys(gtf_path)
    coverage["n_gencode_exons"] = len(all_keys)
    coverage["n_positive_loci"] = len(positive_keys)
    coverage["pool_size"] = pool_size
    coverage["pool_seed"] = pool_seed

    pool = sample_background_pool(all_keys, pool_size=pool_size, seed=pool_seed)
    # Scan target: background pool + every lifted positive locus (keyed on its own
    # GRCh38 coordinates), so both label classes are exposed to the identical
    # binding gate. Positive-locus sequences are fetched from the reference just
    # like the pool exons.
    target_keys = sorted(set(pool) | positive_keys)
    exons = exons_df_from_keys(target_keys)
    exonome = build_exonome(exons, genome_fasta=genome_fasta,
                            cache_path=cache_path, pad=pad, samtools=samtools)
    coverage["n_target_exons"] = len(exonome)
    coverage["target_pad"] = pad
    return exonome, coverage


# --- scan -----------------------------------------------------------------


def _longest_exact_match(a: str, b: str) -> int:
    """Longest exact contiguous match (longest common substring), O(len(a)*len(b))
    DP -- both inputs here are <=~40 nt so this is negligible per call."""
    if not a or not b:
        return 0
    prev = [0] * (len(b) + 1)
    best = 0
    for i in range(1, len(a) + 1):
        curr = [0] * (len(b) + 1)
        ai = a[i - 1]
        for j in range(1, len(b) + 1):
            if ai == b[j - 1]:
                v = prev[j - 1] + 1
                curr[j] = v
                if v > best:
                    best = v
        prev = curr
    return best


# --- low-complexity / promiscuity guard -------------------------------------

# Distinct-3mer floor below which an 18-20 nt SBO is treated as degenerate.
# Calibrated against the real 81-SBO panel: this flags exactly the repeat-motif
# / homopolymer designs (ATCATC..., CTGCTG..., poly-A, (AC)n) whose seeds at
# max_edit=5 (3 nt each) occur pervasively and would "bind" essentially every
# repeat-bearing exon -- spurious near-universal hits, almost all label=0, that
# blow the per-task wall clock and swamp the pooled class balance. They are
# EXCLUDED from the bind->act table with a recorded reason rather than scanned.
_MIN_DISTINCT_3MERS = 7
_MAX_HOMOPOLYMER_RUN = 8
_MAX_DINUCLEOTIDE_FRACTION = 0.66


def sbo_low_complexity_reason(sequence: str) -> str | None:
    """Return a short reason string if ``sequence`` is too low-complexity to
    scan specifically at ``max_edit=5``, else ``None``.

    Purely a function of the SBO sequence (never the label), so it is applied
    identically regardless of any exon's actionable status.
    """
    seq = sequence.strip().upper()
    if len(seq) < 3:
        return "too short"
    kmers = {seq[i : i + 3] for i in range(len(seq) - 2)}
    if len(kmers) < _MIN_DISTINCT_3MERS:
        return f"distinct_3mers={len(kmers)}<{_MIN_DISTINCT_3MERS}"
    longest = best = 1
    for i in range(1, len(seq)):
        best = best + 1 if seq[i] == seq[i - 1] else 1
        longest = max(longest, best)
    if longest >= _MAX_HOMOPOLYMER_RUN:
        return f"homopolymer_run={longest}>={_MAX_HOMOPOLYMER_RUN}"
    if len(seq) > 2:
        dinuc = sum(1 for i in range(2, len(seq)) if seq[i] == seq[i - 2]) / (len(seq) - 2)
        if dinuc >= _MAX_DINUCLEOTIDE_FRACTION:
            return f"dinucleotide_fraction={dinuc:.2f}>={_MAX_DINUCLEOTIDE_FRACTION}"
    return None


def scan_sbo(
    sequence: str,
    exonome: Mapping[str, str],
    *,
    max_edit: int = DEFAULT_MAX_EDIT,
) -> pd.DataFrame:
    """Scan one SBO against every exon in ``exonome``; return the best hit per
    bound exon.

    ``sequence`` is the SBO in its own 5'->3' orientation. SBOs are antisense
    to their target, and ``exonome`` sequences are transcript(sense)-strand,
    so the scanner query is ``reverse_complement(sequence)`` -- NOT
    ``sequence`` itself. (A bug here -- scanning the SBO directly, or
    reverse-complementing the exonome instead -- would silently search for the
    SBO's own sequence rather than its target site; see
    ``test_scan_orientation_matters``.)
    """
    query = reverse_complement(sequence.strip().upper())
    hits: tuple[HybridizationHit, ...] = seed_and_extend(query, exonome, max_edit=max_edit)
    best: dict[str, HybridizationHit] = {}
    for hit in hits:
        current = best.get(hit.reference_name)
        if current is None or (hit.edit_distance, hit.mismatches) < (current.edit_distance, current.mismatches):
            best[hit.reference_name] = hit
    records = []
    for exon_key, hit in best.items():
        window = exonome[exon_key][hit.start:hit.end]
        records.append({
            "exon_key": exon_key,
            "edit_distance": hit.edit_distance,
            "mismatches": hit.mismatches,
            "hit_start": hit.start,
            "hit_end": hit.end,
            "hit_length": hit.end - hit.start,
            "exact_match_length": _longest_exact_match(query, window),
        })
    return pd.DataFrame.from_records(records, columns=[
        "exon_key", "edit_distance", "mismatches", "hit_start", "hit_end",
        "hit_length", "exact_match_length",
    ])


# --- duplex energy ----------------------------------------------------------

# Proxy scale: a crude, deterministic complementarity-based dG stand-in, used
# ONLY if IntaRNA cannot be located. Not calibrated against real thermodynamics
# -- flagged per-row via `duplex_dg_is_proxy` so the gate can report exposure.
_PROXY_DG_PER_MATCHED_BASE = -0.4
_PROXY_DG_PER_GC_BASE = -0.2


def _gc_fraction(sequence: str) -> float:
    seq = sequence.upper()
    if not seq:
        return 0.0
    return sum(1 for b in seq if b in "GC") / len(seq)


def duplex_energy(
    sbo_sequence: str,
    target_sequence: str,
    *,
    intarna: Path | None = None,
) -> tuple[float, bool]:
    """Return ``(duplex_dg, is_proxy)`` for one (SBO, exon-target) pair.

    Runs IntaRNA (SBO as query, exon sequence as target -- the same
    query/target convention :mod:`feng.bench.scores` uses for nusinersen) when
    the executable resolves; otherwise falls back to a deterministic
    edit/GC-based proxy. Called identically for actionable and non-actionable
    rows -- see the module docstring's leakage-discipline note.
    """
    exe = find_intarna(intarna)
    if exe is not None:
        try:
            energy = intarna_min_energy(target_sequence, intarna=exe, query=sbo_sequence)
            return energy, False
        except IntaRNAError:
            pass
    rc = reverse_complement(sbo_sequence.strip().upper())
    exact = _longest_exact_match(rc, target_sequence.upper())
    gc = _gc_fraction(sbo_sequence)
    proxy = (
        _PROXY_DG_PER_MATCHED_BASE * exact
        + _PROXY_DG_PER_GC_BASE * gc * len(sbo_sequence)
    )
    return proxy, True


# --- label join + feature table ---------------------------------------------


def _pool_keys_overlapping_positives(
    pool_keys: Iterable[str], positive_keys: Iterable[str],
    *, margin: int = DEFAULT_TARGET_PAD,
) -> set[str]:
    """Pool exon keys whose interval (inflated by ``margin`` bp) overlaps any of
    ``positive_keys``.

    Overlap is coordinate-only (same chrom, half-open interval intersection),
    strand-agnostic and conservative: a background exon covering -- or, once both
    sequences are pad-fetched, adjacent to -- an SBO's spliced locus must not be
    handed back as a bind-no-act negative, so it is dropped from the negative pool
    for that SBO. ``margin`` inflates the positive interval to account for the
    symmetric fetch pad on both sides.
    """
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for key in positive_keys:
        chrom, start_0base, end, _strand = _parse_exon_key(key)
        by_chrom.setdefault(chrom, []).append((start_0base - margin, end + margin))
    drop: set[str] = set()
    for key in pool_keys:
        chrom, start_0base, end, _strand = _parse_exon_key(key)
        for pos_start, pos_end in by_chrom.get(chrom, ()):
            if start_0base < pos_end and pos_start < end:
                drop.add(key)
                break
    return drop


def build_features(
    *,
    s1: pd.DataFrame,
    s4: pd.DataFrame,
    sbo_sequences: pd.DataFrame,
    exonome: Mapping[str, str],
    max_edit: int = DEFAULT_MAX_EDIT,
    intarna: Path | None = None,
    positives_by_dg: dict[str, set[str]] | None = None,
    chain: ChainMap | None = None,
) -> pd.DataFrame:
    """Scan every SBO in ``sbo_sequences`` against its target set and assemble the
    labeled SITE-LEVEL binding-feature table (one row per bound (SBO, target) pair).

    Per SBO the target set is that SBO's lifted positive loci UNION the shared
    background pool, minus any pool exon overlapping one of the SBO's positive
    loci. ``actionable = 1`` iff the bound target is one of that SBO's lifted
    positive loci; every other bound target is a bind-no-act negative (``0``).
    Features are computed from the SBO's BEST binding window in each target by the
    SAME code path for both labels: ``edit_distance`` / ``exact_match_length`` from
    the scan, ``duplex_dg`` from IntaRNA on (SBO x best window), ``binding_window_gc``
    from that window. Whole-region length / GC are NOT features (substrate
    confound). No S4-derived value ever enters a feature column. Per-SBO
    positive/negative counts are reported on stderr; an SBO with 0 positives among
    its bound targets still contributes its negatives (the LOSO gate handles the
    degenerate held-out fold).
    """
    if positives_by_dg is None:
        positives_by_dg = positive_exon_keys_by_dg(s1, s4, max_edit=max_edit, chain=chain)
    all_positive_keys: set[str] = set().union(*positives_by_dg.values()) if positives_by_dg else set()
    pool_keys = set(exonome) - all_positive_keys

    rows = []
    for row in sbo_sequences.itertuples(index=False):
        dg_id, sbo_sequence = row.dg_id, row.sequence
        reason = sbo_low_complexity_reason(sbo_sequence)
        if reason is not None:
            print(f"SKIP low-complexity SBO {dg_id} ({sbo_sequence}): {reason}", file=sys.stderr)
            continue
        positive_keys = positives_by_dg.get(dg_id, set())
        # This SBO's target set: its own lifted positive loci + the pool, minus any
        # pool exon overlapping one of those loci. Restricting the scan to this
        # subset keeps a locus acted on by a DIFFERENT SBO out of this SBO's rows.
        drop = _pool_keys_overlapping_positives(pool_keys, positive_keys)
        target_keys = (pool_keys - drop) | positive_keys
        sub_exonome = {k: exonome[k] for k in target_keys if k in exonome}
        bound = scan_sbo(sbo_sequence, sub_exonome, max_edit=max_edit)
        n_pos = n_neg = 0
        for b in bound.itertuples(index=False):
            seq = sub_exonome[b.exon_key]
            window = seq[b.hit_start:b.hit_end]
            actionable = int(b.exon_key in positive_keys)
            n_pos += actionable
            n_neg += 1 - actionable
            chrom, start_0base, end, strand = _parse_exon_key(b.exon_key)
            dg, is_proxy = duplex_energy(sbo_sequence, window, intarna=intarna)
            rows.append({
                "sbo_group": dg_id,
                "dg_id": dg_id,
                "exon_key": b.exon_key,
                "chrom": chrom,
                "strand": strand,
                "exon_start_0base": start_0base,
                "exon_end": end,
                "edit_distance": b.edit_distance,
                "mismatches": b.mismatches,
                "hit_length": b.hit_length,
                "exact_match_length": b.exact_match_length,
                "duplex_dg": dg,
                "duplex_dg_is_proxy": is_proxy,
                "binding_window_gc": _gc_fraction(window),
                "actionable": actionable,
            })
        print(f"SBO {dg_id}: bound={n_pos + n_neg} pos={n_pos} neg={n_neg}"
              + ("  [0 positives among bound -> reported, excluded from evaluable folds]" if n_pos == 0 else ""),
              file=sys.stderr)
    return pd.DataFrame.from_records(rows, columns=[
        "sbo_group", "dg_id", "exon_key", "chrom", "strand", "exon_start_0base",
        "exon_end", "edit_distance", "mismatches", "hit_length",
        "exact_match_length", "duplex_dg", "duplex_dg_is_proxy",
        "binding_window_gc", "actionable",
    ])


# --- CLI --------------------------------------------------------------------


def _cmd_build_exonome(args: argparse.Namespace) -> int:
    s1 = parse_s1(args.s1_pdf)
    s4 = parse_s4(args.s4_pdf)
    chain = load_armc_chain(args.chain)
    exonome, cov = build_target_exonome(
        gtf_path=args.gtf, genome_fasta=args.genome_fasta,
        cache_path=args.cache_path, s1=s1, s4=s4, chain=chain,
        pool_size=args.pool_size, pool_seed=args.pool_seed,
        max_edit=args.max_edit, pad=args.pad, samtools=args.samtools,
    )
    n_direct, n_indirect = s4_direct_vs_indirect(s4, max_edit=args.max_edit)
    print(f"chain: {chain.n_chains} chains from {resolve_chain_path(args.chain)}")
    print(f"exonome: {len(exonome)} target exons ({cov['n_positive_loci']} positive "
          f"loci + pool) at {args.cache_path}")
    print(f"S4 rows: direct-bound (numeric edit<= {args.max_edit})={n_direct}  "
          f"indirect (background/NaN, not bound)={n_indirect}")
    print(f"positive-lift: {cov['n_lifted']}/{cov['n_direct_positive_intervals']} "
          f"unique direct-positive intervals lifted hg19->GRCh38 "
          f"({cov['lift_fraction']:.3%}); unmapped={cov['n_unmapped']} "
          f"multimapped={cov['n_multimapped']}")
    # Hard tripwire: refuse to proceed (and let the sbatch `set -e` abort the
    # whole pipeline) when the positive class is effectively empty, rather than
    # spend the 81-SBO scan producing a NO-GO-by-construction gate. A wrong /
    # absent chain would leave nothing lifted.
    if cov["n_direct_positive_intervals"] and cov["lift_fraction"] < 0.90:
        print(
            f"FAIL: positive lift {cov['lift_fraction']:.3%} < 90% -- S4 positive "
            "intervals did not lift to GRCh38 (wrong/absent chain; see module "
            "docstring). Refusing to build a 0-positive gate.",
            file=sys.stderr,
        )
        return 3
    return 0


def _cmd_scan(args: argparse.Namespace) -> int:
    s1 = parse_s1(args.s1_pdf)
    s4 = parse_s4(args.s4_pdf)
    exonome = _read_fasta(args.exonome_cache)
    if not exonome:
        raise SystemExit(f"empty/absent exonome cache: {args.exonome_cache}")
    sbo_sequences = pd.read_csv(args.sbo_tsv, sep="\t")
    row = sbo_sequences.loc[sbo_sequences["dg_id"] == args.dg_id]
    if row.empty:
        raise SystemExit(f"dg_id {args.dg_id!r} not found in {args.sbo_tsv}")
    chain = load_armc_chain(args.chain)
    out = build_features(
        s1=s1, s4=s4, sbo_sequences=row, exonome=exonome,
        max_edit=args.max_edit, intarna=args.intarna, chain=chain,
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_name(args.out.name + f".tmp.{os.getpid()}")
    out.to_parquet(tmp, index=False)
    tmp.replace(args.out)
    print(f"{args.dg_id}: {len(out)} bound exons -> {args.out}")
    return 0


def _cmd_merge(args: argparse.Namespace) -> int:
    paths = sorted(Path(args.hits_dir).glob("*.parquet"))
    if not paths:
        raise SystemExit(f"no per-SBO hits parquet files found in {args.hits_dir}")
    frames = [pd.read_parquet(p) for p in paths]
    merged = pd.concat(frames, ignore_index=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_name(args.out.name + f".tmp.{os.getpid()}")
    merged.to_parquet(tmp, index=False)
    tmp.replace(args.out)
    print(f"merged {len(paths)} shards -> {len(merged)} rows -> {args.out}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Arm C bind->act feature construction.")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_exonome = sub.add_parser("build-exonome", help="build (or verify cached) pool+positives target exonome")
    p_exonome.add_argument("--gtf", type=Path, default=GTF_DEFAULT)
    p_exonome.add_argument("--genome-fasta", type=Path, default=GENOME_FASTA_DEFAULT)
    p_exonome.add_argument("--cache-path", type=Path, required=True)
    p_exonome.add_argument("--s1-pdf", type=Path, default=S1_PDF_DEFAULT)
    p_exonome.add_argument("--s4-pdf", type=Path, default=S4_PDF_DEFAULT)
    p_exonome.add_argument("--chain", type=Path, default=None,
                           help="hg19->GRCh38 UCSC chain (default: $FENG_ARMC_CHAIN or the staged path)")
    p_exonome.add_argument("--pool-size", type=int, default=DEFAULT_POOL_SIZE)
    p_exonome.add_argument("--pool-seed", type=int, default=DEFAULT_POOL_SEED)
    p_exonome.add_argument("--max-edit", type=int, default=DEFAULT_MAX_EDIT)
    p_exonome.add_argument("--pad", type=int, default=DEFAULT_TARGET_PAD,
                           help="symmetric bp pad added to every target region before fetch")
    p_exonome.add_argument("--samtools", type=Path, default=None)
    p_exonome.set_defaults(func=_cmd_build_exonome)

    p_scan = sub.add_parser("scan", help="scan one SBO against the cached exonome, write its labeled hits")
    p_scan.add_argument("--dg-id", required=True)
    p_scan.add_argument("--sbo-tsv", type=Path, default=SBO_TSV_DEFAULT)
    p_scan.add_argument("--s1-pdf", type=Path, default=S1_PDF_DEFAULT)
    p_scan.add_argument("--s4-pdf", type=Path, default=S4_PDF_DEFAULT)
    p_scan.add_argument("--exonome-cache", type=Path, required=True)
    p_scan.add_argument("--chain", type=Path, default=None,
                        help="hg19->GRCh38 UCSC chain (default: $FENG_ARMC_CHAIN or the staged path)")
    p_scan.add_argument("--intarna", type=Path, default=None)
    p_scan.add_argument("--max-edit", type=int, default=DEFAULT_MAX_EDIT)
    p_scan.add_argument("--out", type=Path, required=True)
    p_scan.set_defaults(func=_cmd_scan)

    p_merge = sub.add_parser("merge", help="merge per-SBO hits parquets into armc_features.parquet")
    p_merge.add_argument("--hits-dir", type=Path, required=True)
    p_merge.add_argument("--out", type=Path, default=FEATURES_OUT_DEFAULT)
    p_merge.set_defaults(func=_cmd_merge)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
