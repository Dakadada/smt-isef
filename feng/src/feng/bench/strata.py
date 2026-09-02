"""Control-only detectable-exon universe for the Arm D locked benchmark.

This module builds the background universe of internal cassette exons that are
adequately measured in the GSE307431 baseline libraries, using ONLY the
scrambled and untreated control samples. It never reads a treatment (F18OMe /
F18MOE / F20PMO) alignment and never consumes a treatment-vs-control effect
table. The heavy rMATS-turbo enumeration runs in the companion batch script;
this module parses that control-only output and derives per-exon metrics.

The rMATS run is invoked with ``--statoff`` and with both of its sample groups
drawn from the control pool. ``--statoff`` disables statistical testing, but
rMATS still emits a dummy-group ``IncLevelDifference`` for some rows. Contrast
fields are ignored and explicitly forbidden from the published parquet. The two
groups exist only to satisfy the rMATS interface. Pooling the scrambled and
untreated controls is a baseline pool, not a chemistry pool -- the three ASO
chemistries are never combined anywhere in this arm.

Allowed evidence states for anything reported here are ``candidate`` and
``paper_supported`` only.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import argparse
import hashlib
import json
import math
import subprocess
import tempfile

import pandas as pd

# --- control selection contract (freeze guard) ------------------------------

# The only sample classes that may enter this universe. Anything with a
# targeted control_type or an ASO-chemistry condition is a treatment sample and
# must never be selected here.
CONTROL_TYPES = frozenset({"scrambled", "untreated"})
CONTROL_CONDITION = "control"
TREATMENT_CONDITIONS = frozenset({"F18OMe", "F18MOE", "F20PMO"})
EXON_IDENTITY_COLUMNS = (
    "chrom",
    "exon_start_0base",
    "exon_end",
    "strand",
)
CONTRAST_FIELDS = frozenset({"PValue", "FDR", "IncLevelDifference"})


class FreezeGuardError(RuntimeError):
    """Raised when a non-control sample is about to enter the universe build."""


def select_control_runs(conditions: pd.DataFrame) -> pd.DataFrame:
    """Return the control rows and refuse any treatment leakage.

    A row qualifies only when its ``condition`` is exactly ``control`` and its
    ``control_type`` is one of the allowed baseline classes. If any qualifying
    row still carries a treatment marker, or if the selection is empty, the
    build is aborted rather than silently narrowed.
    """
    required = {"run_accession", "condition", "control_type", "chemistry"}
    missing = required.difference(conditions.columns)
    if missing:
        raise ValueError(f"conditions table lacks columns: {sorted(missing)}")

    # Independently derive the set of run accessions that carry any treatment
    # marker anywhere in the table (an ASO-chemistry condition or a targeted
    # control_type). This is a cross-check against the selection predicate, so a
    # run double-listed as both control and treatment is caught rather than
    # trusted.
    treatment_mask = conditions["condition"].isin(TREATMENT_CONDITIONS) | (
        conditions["control_type"] == "targeted"
    )
    treatment_runs = set(conditions.loc[treatment_mask, "run_accession"])

    is_control = (conditions["condition"] == CONTROL_CONDITION) & conditions[
        "control_type"
    ].isin(CONTROL_TYPES)
    selected = conditions.loc[is_control].copy()

    duplicates = sorted(
        selected.loc[selected["run_accession"].duplicated(keep=False), "run_accession"]
        .astype(str)
        .unique()
    )
    if duplicates:
        raise FreezeGuardError(
            "duplicate control run accessions: " + ", ".join(duplicates)
        )
    leaked = sorted(set(selected["run_accession"]).intersection(treatment_runs))
    if leaked:
        raise FreezeGuardError(
            "run accessions listed as both control and treatment: "
            + ", ".join(leaked)
        )
    if selected.empty:
        raise FreezeGuardError("no control samples selected; refusing to build")
    return selected.reset_index(drop=True)


# --- detectability thresholds (explicit, documented) ------------------------


@dataclass(frozen=True, slots=True)
class DetectabilityThresholds:
    """Thresholds defining an adequately measured background exon.

    ``min_mean_total_junction_reads`` -- an exon must carry, on average across
    the pooled control samples, at least this many uniquely counted junction
    reads (inclusion junctions + skipping junction). This is the depth needed
    for a stable inclusion estimate.

    ``min_samples_with_psi`` -- the inclusion level must be estimable (not NA)
    in at least this many control samples, so a single well-covered library
    cannot alone qualify an exon.
    """

    min_mean_total_junction_reads: float = 20.0
    min_samples_with_psi: int = 6


DEFAULT_THRESHOLDS = DetectabilityThresholds()


# --- rMATS SE parsing --------------------------------------------------------


def _parse_int_list(cell: str) -> list[int]:
    cell = (cell or "").strip()
    if not cell:
        return []
    return [int(round(float(part))) for part in cell.split(",") if part != ""]


def _parse_psi_list(cell: str) -> list[float]:
    values: list[float] = []
    for part in (cell or "").split(","):
        part = part.strip()
        if part == "" or part.upper() == "NA":
            continue
        values.append(float(part))
    return values


def load_se_events(se_table: Path, *, expected_sample_count: int = 11) -> pd.DataFrame:
    """Read a control-only ``SE.MATS.JC.txt`` table into a tidy frame.

    Both rMATS sample groups are controls, so their junction counts and
    inclusion levels are pooled into a single baseline set per exon.
    """
    frame = pd.read_csv(se_table, sep="\t", dtype=str)
    records = []
    for _, row in frame.iterrows():
        ijc_groups = (
            _parse_int_list(row["IJC_SAMPLE_1"]),
            _parse_int_list(row["IJC_SAMPLE_2"]),
        )
        sjc_groups = (
            _parse_int_list(row["SJC_SAMPLE_1"]),
            _parse_int_list(row["SJC_SAMPLE_2"]),
        )
        event_id = row["ID"].strip().strip('"')
        for group_index, (ijc_group, sjc_group) in enumerate(
            zip(ijc_groups, sjc_groups), start=1
        ):
            if len(ijc_group) != len(sjc_group):
                raise ValueError(
                    f"event {event_id}: IJC/SJC sample-count mismatch in group "
                    f"{group_index} ({len(ijc_group)} != {len(sjc_group)})"
                )
        ijc = ijc_groups[0] + ijc_groups[1]
        sjc = sjc_groups[0] + sjc_groups[1]
        if len(ijc) != expected_sample_count:
            raise ValueError(
                f"event {event_id}: expected {expected_sample_count} pooled control "
                f"samples, found {len(ijc)}"
            )
        psi = _parse_psi_list(row["IncLevel1"]) + _parse_psi_list(row["IncLevel2"])
        totals = [i + s for i, s in zip(ijc, sjc)]
        gene_id = row["GeneID"].strip().strip('"')
        records.append(
            {
                "event_id": event_id,
                "gene_id": gene_id,
                "gene_id_base": gene_id.split(".")[0],
                "gene_symbol": row["geneSymbol"].strip().strip('"'),
                "chrom": row["chr"].strip(),
                "strand": row["strand"].strip(),
                "exon_start_0base": int(row["exonStart_0base"]),
                "exon_end": int(row["exonEnd"]),
                "upstream_ee": int(row["upstreamEE"]),
                "downstream_es": int(row["downstreamES"]),
                "mean_inclusion_junction_reads": _safe_mean(ijc),
                "mean_skipping_junction_reads": _safe_mean(sjc),
                "mean_total_junction_reads": _safe_mean(totals),
                "n_samples_total": len(totals),
                "n_samples_with_psi": len(psi),
                "baseline_psi": _safe_mean(psi),
            }
        )
    events = pd.DataFrame.from_records(records)
    # An SE record is an internal cassette exon by construction: it always has a
    # defined upstream and downstream flanking exon, so the 5'-most and 3'-most
    # exons of a transcript are already excluded.
    events["exon_length"] = events["exon_end"] - events["exon_start_0base"]
    events["exon_key"] = (
        events["chrom"].astype(str)
        + ":"
        + events["exon_start_0base"].astype(str)
        + "-"
        + events["exon_end"].astype(str)
        + ":"
        + events["strand"].astype(str)
    )
    assert_no_contrast_fields(events)
    return events


def _safe_mean(values: Iterable[float]) -> float:
    values = list(values)
    if not values:
        return math.nan
    return float(sum(values) / len(values))


# --- gene expression (control mean CPM) -------------------------------------


def load_control_expression(reads_per_gene_paths: list[Path]) -> pd.DataFrame:
    """Mean control CPM per gene from STAR ``ReadsPerGene.out.tab`` files.

    The reverse-strand column (index 3) is used because the libraries are
    fr-firststrand. The four leading summary rows are excluded from the library
    size and from the per-gene table.
    """
    per_sample_cpm: list[pd.Series] = []
    for path in reads_per_gene_paths:
        table = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["gene_id", "unstranded", "forward", "reverse"],
            dtype={"gene_id": str},
        )
        table = table[~table["gene_id"].str.startswith("N_")].copy()
        library_size = table["reverse"].sum()
        if library_size <= 0:
            continue
        cpm = table.set_index("gene_id")["reverse"] / library_size * 1e6
        per_sample_cpm.append(cpm)
    if not per_sample_cpm:
        return pd.DataFrame(columns=["gene_id_base", "control_mean_cpm"])
    mean_cpm = pd.concat(per_sample_cpm, axis=1).mean(axis=1)
    out = mean_cpm.reset_index()
    out.columns = ["gene_id", "control_mean_cpm"]
    out["gene_id_base"] = out["gene_id"].str.split(".").str[0]
    return out[["gene_id_base", "control_mean_cpm"]]


# --- GC content --------------------------------------------------------------


class GenomeFastaError(RuntimeError):
    """Raised when required reference sequence cannot be read completely."""


def _gc_fraction(sequence: str) -> float:
    sequence = sequence.upper()
    if not sequence:
        raise GenomeFastaError("reference FASTA returned an empty exon sequence")
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def _resolve_fasta_contig(chrom: str, available: set[str]) -> str:
    """Resolve UCSC/non-UCSC aliases without weakening missing-contig checks."""
    candidates = [chrom]
    if chrom.startswith("chr"):
        candidates.append(chrom[3:])
    else:
        candidates.append(f"chr{chrom}")
    if chrom in {"chrM", "M", "MT"}:
        candidates.extend(["chrM", "M", "MT"])
    for candidate in candidates:
        if candidate in available:
            return candidate
    raise GenomeFastaError(f"FASTA lacks contig {chrom!r}")


def _attach_gc_with_pysam(events: pd.DataFrame, genome_fasta: Path) -> list[float]:
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - depends on build environment
        raise GenomeFastaError("pysam is unavailable") from exc

    gc_values: list[float] = []
    try:
        with pysam.FastaFile(str(genome_fasta)) as fasta:
            available = set(fasta.references)
            for _, row in events.iterrows():
                chrom = str(row["chrom"])
                fasta_chrom = _resolve_fasta_contig(chrom, available)
                sequence = fasta.fetch(
                    fasta_chrom,
                    int(row["exon_start_0base"]),
                    int(row["exon_end"]),
                )
                if len(sequence) != int(row["exon_length"]):
                    raise GenomeFastaError(
                        f"FASTA returned {len(sequence)} bases for {row['exon_key']}; "
                        f"expected {row['exon_length']}"
                    )
                gc_values.append(_gc_fraction(sequence))
    except GenomeFastaError:
        raise
    except Exception as exc:
        raise GenomeFastaError(f"failed to open/read FASTA {genome_fasta}") from exc
    return gc_values


def _attach_gc_with_samtools(
    events: pd.DataFrame, genome_fasta: Path, samtools: Path
) -> list[float]:
    """Read all unique exon intervals in one fail-closed ``samtools faidx`` call."""
    with tempfile.TemporaryDirectory(prefix="feng-fasta-") as tmp:
        tmpdir = Path(tmp)
        regions = tmpdir / "regions.txt"
        fai = tmpdir / "reference.fai"
        try:
            subprocess.run(
                [str(samtools), "faidx", "--fai-idx", str(fai), str(genome_fasta)],
                check=True,
                capture_output=True,
                text=True,
            )
        except (OSError, subprocess.CalledProcessError) as exc:
            stderr = getattr(exc, "stderr", "")
            raise GenomeFastaError(
                f"failed to index FASTA {genome_fasta} with {samtools}: "
                f"{stderr.strip()}"
            ) from exc
        available = {
            line.split("\t", 1)[0]
            for line in fai.read_text(encoding="utf-8").splitlines()
            if line
        }
        regions.write_text(
            "".join(
                f"{_resolve_fasta_contig(str(row.chrom), available)}:"
                f"{int(row.exon_start_0base) + 1}-{int(row.exon_end)}\n"
                for row in events.itertuples(index=False)
            ),
            encoding="utf-8",
        )
        command = [
            str(samtools),
            "faidx",
            "--fai-idx",
            str(fai),
            "-r",
            str(regions),
            str(genome_fasta),
        ]
        try:
            result = subprocess.run(
                command, check=True, capture_output=True, text=True
            )
        except (OSError, subprocess.CalledProcessError) as exc:
            stderr = getattr(exc, "stderr", "")
            raise GenomeFastaError(
                f"failed to open/read FASTA {genome_fasta} with {samtools}: "
                f"{stderr.strip()}"
            ) from exc

    sequences: list[str] = []
    current: list[str] = []
    for line in result.stdout.splitlines():
        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
        else:
            current.append(line.strip())
    if current:
        sequences.append("".join(current))
    if len(sequences) != len(events):
        raise GenomeFastaError(
            f"FASTA returned {len(sequences)} exon sequences; expected {len(events)}"
        )
    values = [_gc_fraction(sequence) for sequence in sequences]
    for length, key, sequence in zip(
        events["exon_length"], events["exon_key"], sequences
    ):
        if len(sequence) != int(length):
            raise GenomeFastaError(
                f"FASTA returned {len(sequence)} bases for {key}; expected {length}"
            )
    return values


def attach_gc_content(
    events: pd.DataFrame,
    genome_fasta: Path,
    *,
    samtools: Path | None = None,
) -> pd.DataFrame:
    """Attach required exon GC fractions, aborting on any missing sequence.

    ``pysam`` is preferred. An explicitly available ``samtools`` executable is
    a supported random-access fallback for environments where the Python wheel
    is absent. If neither backend can read every exon exactly, no output is
    returned and the build fails.
    """
    genome_fasta = Path(genome_fasta)
    if not genome_fasta.is_file() or genome_fasta.stat().st_size == 0:
        raise GenomeFastaError(f"reference FASTA is missing or empty: {genome_fasta}")
    if (events["exon_length"] <= 0).any():
        raise GenomeFastaError("cannot compute GC for a non-positive exon length")

    events = events.copy()
    try:
        values = _attach_gc_with_pysam(events, genome_fasta)
    except GenomeFastaError as pysam_error:
        executable = Path(samtools) if samtools is not None else None
        if executable is None or not executable.is_file():
            raise GenomeFastaError(
                "pysam could not read the required FASTA and no usable samtools "
                "fallback was supplied"
            ) from pysam_error
        values = _attach_gc_with_samtools(events, genome_fasta, executable)
    events["gc_fraction"] = values
    if events["gc_fraction"].isna().any():
        raise GenomeFastaError("GC calculation produced null values")
    return events


# --- assembly + classification ----------------------------------------------


def classify_detectable(
    events: pd.DataFrame, thresholds: DetectabilityThresholds = DEFAULT_THRESHOLDS
) -> pd.DataFrame:
    """Flag exons that clear the detectability thresholds."""
    events = events.copy()
    events["detectable"] = (
        (events["mean_total_junction_reads"] >= thresholds.min_mean_total_junction_reads)
        & (events["n_samples_with_psi"] >= thresholds.min_samples_with_psi)
        & events["baseline_psi"].notna()
    )
    return events


def assert_no_contrast_fields(events: pd.DataFrame) -> None:
    """Refuse publication of any dummy-group or treatment contrast field."""
    leaked = sorted(CONTRAST_FIELDS.intersection(events.columns))
    if leaked:
        raise FreezeGuardError(
            "contrast fields are forbidden from the detectable universe: "
            + ", ".join(leaked)
        )


def collapse_to_unique_exons(events: pd.DataFrame) -> pd.DataFrame:
    """Collapse alternate-flank rMATS records to one genomic exon identity.

    rMATS may enumerate the same cassette exon repeatedly with different
    upstream/downstream transcript flanks. The representative record is chosen
    deterministically, preferring a detectable event and then the event with the
    strongest sample coverage and junction depth. All source event IDs and the
    multiplicity remain attached for auditability.
    """
    missing = set(EXON_IDENTITY_COLUMNS).difference(events.columns)
    if missing:
        raise ValueError(f"events lack exon identity columns: {sorted(missing)}")
    events = events.copy()
    if "exon_key" not in events:
        events["exon_key"] = (
            events["chrom"].astype(str)
            + ":"
            + events["exon_start_0base"].astype(str)
            + "-"
            + events["exon_end"].astype(str)
            + ":"
            + events["strand"].astype(str)
        )
    grouped = events.groupby("exon_key", sort=False)
    source_ids = (
        grouped["event_id"]
        .agg(lambda values: ",".join(sorted({str(value) for value in values})))
        .rename("rmats_event_ids")
    )
    multiplicity = grouped.size().rename("n_rmats_events")

    # One genomic interval can be annotated to multiple genes. Keep the
    # depth-ranked representative columns for backward compatibility, but attach
    # the complete gene provenance and an explicit ambiguity flag so downstream
    # expression matching cannot mistake the representative for a unique gene.
    gene_aggregates: list[pd.Series] = []
    for source_column, output_column in (
        ("gene_id", "all_gene_ids"),
        ("gene_id_base", "all_gene_id_bases"),
        ("gene_symbol", "all_gene_symbols"),
    ):
        if source_column in events.columns:
            gene_aggregates.append(
                grouped[source_column]
                .agg(
                    lambda values: ",".join(
                        sorted({str(value) for value in values if pd.notna(value)})
                    )
                )
                .rename(output_column)
            )

    expression_provenance: pd.Series | None = None
    if {"gene_id_base", "control_mean_cpm"}.issubset(events.columns):
        expression_provenance = grouped.apply(
            lambda rows: json.dumps(
                {
                    str(gene): (
                        None
                        if values.dropna().empty
                        else float(values.dropna().iloc[0])
                    )
                    for gene, values in rows.groupby("gene_id_base", sort=True)[
                        "control_mean_cpm"
                    ]
                },
                sort_keys=True,
                separators=(",", ":"),
            ),
            include_groups=False,
        ).rename("control_mean_cpm_by_gene")
    ranked = events.sort_values(
        [
            "detectable",
            "n_samples_with_psi",
            "mean_total_junction_reads",
            "event_id",
        ],
        ascending=[False, False, False, True],
        kind="mergesort",
    )
    collapsed = ranked.drop_duplicates("exon_key", keep="first").copy()
    collapsed = collapsed.merge(source_ids, on="exon_key", how="left", validate="one_to_one")
    collapsed = collapsed.merge(multiplicity, on="exon_key", how="left", validate="one_to_one")
    for aggregate in gene_aggregates:
        collapsed = collapsed.merge(
            aggregate, on="exon_key", how="left", validate="one_to_one"
        )
    if expression_provenance is not None:
        collapsed = collapsed.merge(
            expression_provenance, on="exon_key", how="left", validate="one_to_one"
        )
    if "all_gene_id_bases" in collapsed.columns:
        collapsed["n_gene_ids"] = collapsed["all_gene_id_bases"].map(
            lambda value: 0 if not value else len(str(value).split(","))
        )
        collapsed["multi_gene"] = collapsed["n_gene_ids"] > 1
    collapsed = collapsed.sort_values(
        list(EXON_IDENTITY_COLUMNS), kind="mergesort"
    ).reset_index(drop=True)
    if collapsed["exon_key"].duplicated().any():  # pragma: no cover - invariant
        raise AssertionError("exon collapse left duplicate identities")
    return collapsed


def positive_interval_to_exon_key(
    chrom: str, start_1based: int, end_1based: int, strand: str
) -> str:
    """Map a 1-based inclusive positive interval to the universe ``exon_key``.

    Positives (``feng/artifacts/bench/positives_v50.tsv``) record 1-based
    inclusive GTF coordinates, while the detectable universe keys exons in the
    rMATS 0-based half-open frame (``exon_start_0base``). Any coordinate join
    between the two MUST pass through this converter: the start loses exactly one
    base and the end is unchanged, so a positive of length ``end - start + 1``
    maps to a half-open interval of the same length. Skipping this conversion
    silently shifts every positive by one base and fails the >=50-background
    gate with zero matches.
    """
    exon_start_0base = int(start_1based) - 1
    exon_end = int(end_1based)
    if exon_end <= exon_start_0base:
        raise ValueError(
            f"non-positive interval for {chrom}:{start_1based}-{end_1based}"
        )
    return f"{chrom}:{exon_start_0base}-{exon_end}:{strand}"


def build_universe(
    se_table: Path,
    reads_per_gene_paths: list[Path],
    genome_fasta: Path,
    thresholds: DetectabilityThresholds = DEFAULT_THRESHOLDS,
    *,
    samtools: Path | None = None,
) -> pd.DataFrame:
    """Build one fail-closed row per unique internal cassette exon."""
    events = load_se_events(se_table)
    expression = load_control_expression(reads_per_gene_paths)
    events = events.merge(expression, on="gene_id_base", how="left")
    events = classify_detectable(events, thresholds)
    events = collapse_to_unique_exons(events)
    events = attach_gc_content(events, genome_fasta, samtools=samtools)
    assert_no_contrast_fields(events)
    return events


# --- provenance report -------------------------------------------------------


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def render_readme(
    events: pd.DataFrame,
    control_runs: pd.DataFrame,
    se_table: Path,
    thresholds: DetectabilityThresholds,
    genome_fasta: Path,
    provenance_files: tuple[tuple[str, Path], ...] = (),
) -> str:
    detectable = events[events["detectable"]]
    per_chrom = (
        detectable.groupby("chrom").size().sort_values(ascending=False)
    )
    psi = detectable["baseline_psi"].dropna()
    # Count the union of every gene a detectable interval is annotated to, not
    # just the depth-ranked representative gene, so multi-gene intervals are not
    # undercounted.
    if "all_gene_id_bases" in detectable.columns:
        gene_universe: set[str] = set()
        for value in detectable["all_gene_id_bases"].dropna():
            gene_universe.update(part for part in str(value).split(",") if part)
        n_distinct_genes = len(gene_universe)
    else:
        n_distinct_genes = int(detectable["gene_id_base"].nunique())
    lines = [
        "# Arm D control-only detectable-exon universe (GENCODE v50)",
        "",
        "> Claim boundary: descriptive background construction only. Allowed",
        "> result states are `candidate` and `paper_supported`. No treatment",
        "> effect is read or reported here.",
        "",
        "## What this is",
        "",
        "The background set of internal cassette exons (rMATS SE events) that are",
        "adequately measured in the GSE307431 baseline libraries, built from the",
        "scrambled and untreated control alignments only. It supplies the matched",
        "background pool for the locked nusinersen off-target ranking; it does not",
        "contain or depend on any ASO-treated measurement.",
        "",
        "## Control samples used [paper_supported]",
        "",
        "| run_accession | condition | chemistry | control_type |",
        "| --- | --- | --- | --- |",
    ]
    for _, row in control_runs.iterrows():
        lines.append(
            f"| {row['run_accession']} | {row['condition']} | "
            f"{row['chemistry']} | {row['control_type']} |"
        )
    lines += [
        "",
        f"Total control libraries pooled: {len(control_runs)} "
        f"({int((control_runs['control_type'] == 'scrambled').sum())} scrambled + "
        f"{int((control_runs['control_type'] == 'untreated').sum())} untreated).",
        "",
        "## Method [candidate]",
        "",
        "- rMATS-turbo SE events enumerated from GENCODE v50 with `--task both`",
        "  and `--statoff`; both rMATS groups are controls. `--statoff` disables",
        "  statistical testing, but rMATS still emits a dummy-control",
        "  `IncLevelDifference` in some source rows. The builder ignores and",
        "  explicitly rejects all contrast fields before parquet publication.",
        "  Junction counts and inclusion levels from both",
        "  groups are pooled into one baseline set per exon. This baseline pool",
        "  never combines the three ASO chemistries.",
        "- Alternate-flank rMATS events are collapsed by the genomic exon key",
        "  `chrom:exonStart_0base-exonEnd:strand`; each parquet row is one unique",
        "  exon interval. Source event IDs and multiplicity remain in audit columns.",
        "  If one interval is annotated to multiple genes, all gene IDs/symbols and",
        "  per-gene control CPM values are retained in aggregate provenance columns;",
        "  `multi_gene` makes that ambiguity explicit.",
        "- Inclusion level (PSI) is the rMATS junction-count inclusion level; an",
        "  exon's junction depth is the mean of inclusion-junction plus",
        "  skipping-junction uniquely counted reads across control samples.",
        "- Internal exons only: an SE event always has a defined flanking exon on",
        "  each side, so the 5'-most and 3'-most exons of a transcript are",
        "  excluded by construction.",
        "- Gene expression is mean control CPM from STAR reverse-strand gene",
        "  counts (fr-firststrand). GC fraction and exon length come from the",
        "  reference exon coordinates.",
        "",
        "## Detectability thresholds [candidate]",
        "",
        f"- mean total junction reads >= {thresholds.min_mean_total_junction_reads:g}",
        f"- inclusion level estimable in >= {thresholds.min_samples_with_psi} "
        "control samples",
        "",
        "## Counts [candidate]",
        "",
        f"- Unique internal exon intervals after alternate-flank collapse: {len(events)}",
        f"- Detectable background exons: {int(events['detectable'].sum())}",
        f"- Distinct genes among detectable exons: {n_distinct_genes}",
    ]
    if not psi.empty:
        lines.append(
            "- Detectable baseline PSI quartiles (Q1/median/Q3): "
            f"{psi.quantile(0.25):.3f} / {psi.median():.3f} / {psi.quantile(0.75):.3f}"
        )
    lines += [
        "",
        "### Detectable exons by chromosome (density) [candidate]",
        "",
        "| chrom | detectable_exons |",
        "| --- | --- |",
    ]
    for chrom, count in per_chrom.items():
        lines.append(f"| {chrom} | {int(count)} |")
    lines += [
        "",
        "## Provenance [candidate]",
        "",
        f"- SE source table: `{se_table}`",
        f"- SE source SHA-256: `{_sha256(se_table)}`",
    ]
    lines.append(f"- Reference FASTA: `{genome_fasta}`")
    for label, path in provenance_files:
        lines.append(f"- {label}: `{path}`")
        lines.append(f"- {label} SHA-256: `{_sha256(path)}`")
    lines += [
        "",
        "## Downstream gate note",
        "",
        "The locked benchmark requires at least 50 eligible background exons per",
        "positive off-target exon. That gate is evaluated once the positives",
        "table lands; this file reports the detectable count and per-chromosome",
        "density so the gate can be checked against the matched strata.",
        "",
    ]
    return "\n".join(lines)


# --- CLI ---------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--conditions", required=True, type=Path)
    parser.add_argument("--se-table", required=True, type=Path)
    parser.add_argument("--star-dir", required=True, type=Path)
    parser.add_argument("--genome-fasta", required=True, type=Path)
    parser.add_argument("--samtools", type=Path, default=None)
    parser.add_argument("--b1", type=Path, default=None)
    parser.add_argument("--b2", type=Path, default=None)
    parser.add_argument("--out-parquet", required=True, type=Path)
    parser.add_argument("--out-readme", required=True, type=Path)
    parser.add_argument(
        "--min-mean-junction-reads",
        type=float,
        default=DEFAULT_THRESHOLDS.min_mean_total_junction_reads,
    )
    parser.add_argument(
        "--min-samples-with-psi",
        type=int,
        default=DEFAULT_THRESHOLDS.min_samples_with_psi,
    )
    args = parser.parse_args(argv)

    conditions = pd.read_csv(args.conditions, sep="\t", dtype=str)
    control_runs = select_control_runs(conditions)

    reads_per_gene = []
    for run in control_runs["run_accession"]:
        path = args.star_dir / run / f"{run}.ReadsPerGene.out.tab"
        if path.exists():
            reads_per_gene.append(path)

    thresholds = DetectabilityThresholds(
        min_mean_total_junction_reads=args.min_mean_junction_reads,
        min_samples_with_psi=args.min_samples_with_psi,
    )
    events = build_universe(
        args.se_table,
        reads_per_gene,
        args.genome_fasta,
        thresholds,
        samtools=args.samtools,
    )

    args.out_parquet.parent.mkdir(parents=True, exist_ok=True)
    assert_no_contrast_fields(events)
    if events["exon_key"].duplicated().any():
        raise AssertionError("refusing to publish duplicate exon identities")
    if events["gc_fraction"].isna().any():
        raise GenomeFastaError("refusing to publish null GC fractions")
    events.to_parquet(args.out_parquet, index=False)
    provenance_files: list[tuple[str, Path]] = [("Conditions table", args.conditions)]
    if args.b1 is not None:
        provenance_files.append(("rMATS b1 control list", args.b1))
    if args.b2 is not None:
        provenance_files.append(("rMATS b2 control list", args.b2))
    readme = render_readme(
        events,
        control_runs,
        args.se_table,
        thresholds,
        args.genome_fasta,
        tuple(provenance_files),
    )
    args.out_readme.write_text(readme, encoding="utf-8")

    print(
        f"detectable_universe: {int(events['detectable'].sum())} detectable "
        f"of {len(events)} unique exon intervals; parquet -> {args.out_parquet}"
    )
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
