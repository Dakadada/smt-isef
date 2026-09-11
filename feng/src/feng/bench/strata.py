"""Control-only detectable-exon universe for the Arm D locked benchmark.

This module builds the background universe of annotated transcript-internal
exons that are adequately measured in the GSE307431 baseline libraries, using
ONLY the three scrambled MOE (ScrMOE) control samples. It never reads a treatment (F18OMe /
F18MOE / F20PMO) alignment and never consumes a treatment-vs-control effect
table. GENCODE supplies candidate identities and transcript flanks; ScrMOE STAR
junction tables supply the same JC depth/PSI estimator used by rMATS. The
control-only rMATS output remains attached as provenance when an SE event exists.

The rMATS run is invoked with ``--statoff`` and with both of its sample groups
drawn from the control pool. ``--statoff`` disables statistical testing, but
rMATS still emits a dummy-group ``IncLevelDifference`` for some rows. Contrast
fields are ignored and explicitly forbidden from the published parquet. The two
groups exist only to satisfy the rMATS interface.

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
import os
import re
import subprocess
import tempfile
import unicodedata

import numpy as np
import pandas as pd

# --- control selection contract (freeze guard) ------------------------------

# The only sample classes that may enter this universe. Anything with a
# targeted control_type or an ASO-chemistry condition is a treatment sample and
# must never be selected here.
CONTROL_TYPES = frozenset({"scrambled"})
CONTROL_CONDITION = "control"
CONTROL_CHEMISTRY = "MOE"
EXPECTED_CONTROL_COUNT = 3
CONTROL_RUN_IDS = ("SRR35298704", "SRR35298699", "SRR35298715")
CONTROL_ALLOWLIST_SHA256 = "295328ba55a16d2d009c0ff69657a2d4ab8547cbe254917bcc50b4344b4a412e"
CONTROL_RUN_IDS_CANONICAL = ",".join(CONTROL_RUN_IDS)
CONTROL_UNIVERSE_POLICY_ID = "armD-scrmoe-only-v1"
DETECTABILITY_POLICY_ID = "armD-detectability-scrmoe-v1"
# Frozen depth threshold; candidate enumeration must not alter this criterion.
FROZEN_MIN_MEAN_TOTAL_JUNCTION_READS = 20.0
# The ScrMOE migration has exactly three controls, so PSI must be estimable in
# at least two of three replicates (the prespecified strict-majority rule).
SCRMOE_MIN_SAMPLES_WITH_PSI = 2
DETECTABILITY_POLICY = json.dumps(
    {
        "policy_id": DETECTABILITY_POLICY_ID,
        "min_mean_total_junction_reads": FROZEN_MIN_MEAN_TOTAL_JUNCTION_READS,
        "min_samples_with_psi": SCRMOE_MIN_SAMPLES_WITH_PSI,
        "migration_from": "pooled-11-control minimum 6 of 11",
        "rationale": (
            "With exactly three prespecified ScrMOE replicates, require PSI in at "
            "least two independent replicates (a strict majority); retain the "
            "prespecified mean-junction-depth threshold of 20."
        ),
    },
    sort_keys=True,
    separators=(",", ":"),
)
UNIVERSE_PROVENANCE = {
    "control_universe_policy_id": CONTROL_UNIVERSE_POLICY_ID,
    "control_allowlist_sha256": CONTROL_ALLOWLIST_SHA256,
    "control_run_ids": CONTROL_RUN_IDS_CANONICAL,
    "detectability_policy_id": DETECTABILITY_POLICY_ID,
    "detectability_policy": DETECTABILITY_POLICY,
}
UNIVERSE_PROVENANCE_COLUMNS = tuple(UNIVERSE_PROVENANCE)
INPUT_PROVENANCE_COLUMNS = (
    "gencode_gtf_sha256",
    "scrmoe_sj_sha256",
    "genome_fasta_sha256",
    "scrmoe_reads_per_gene_sha256",
    "junction_read_length",
    "junction_read_length_sha256",
)
TREATMENT_CONDITIONS = frozenset({"F18OMe", "F18MOE", "F20PMO"})
EXON_IDENTITY_COLUMNS = (
    "chrom",
    "exon_start_0base",
    "exon_end",
    "strand",
)
CONTRAST_FIELDS = frozenset({"PValue", "FDR", "IncLevelDifference"})
# GENCODE transcript classes whose exon boundaries represent comparable,
# splice-competent coding transcripts. Retained-intron, processed-transcript,
# lncRNA, TEC, and other transcript classes are intentionally outside Arm D's
# internal-exon background.
ALLOWED_GENCODE_TRANSCRIPT_TYPES = frozenset(
    {"protein_coding", "nonsense_mediated_decay", "non_stop_decay"}
)
ALLOWED_GENCODE_FEATURE_TYPES = frozenset({"exon"})
EXCLUDED_GENCODE_TRANSCRIPT_TAG_PREFIXES = ("retained_intron",)


class FreezeGuardError(RuntimeError):
    """Raised when a non-control sample is about to enter the universe build."""


def select_control_runs(conditions: pd.DataFrame) -> pd.DataFrame:
    """Return the control rows and refuse any treatment leakage.

    A row qualifies only when its ``condition`` is exactly ``control``, its
    ``control_type`` is ``scrambled``, and its chemistry is ``MOE``. If any qualifying
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
    ].isin(CONTROL_TYPES) & (conditions["chemistry"] == CONTROL_CHEMISTRY)
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
    if len(selected) != EXPECTED_CONTROL_COUNT:
        raise FreezeGuardError(
            f"expected exactly {EXPECTED_CONTROL_COUNT} ScrMOE controls; "
            f"found {len(selected)}"
        )
    observed = set(selected["run_accession"].astype(str))
    if observed != set(CONTROL_RUN_IDS):
        raise FreezeGuardError(
            "ScrMOE controls differ from the pinned run IDs: "
            f"observed={sorted(observed)} expected={sorted(CONTROL_RUN_IDS)}"
        )
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

    min_mean_total_junction_reads: float = FROZEN_MIN_MEAN_TOTAL_JUNCTION_READS
    min_samples_with_psi: int = SCRMOE_MIN_SAMPLES_WITH_PSI


DEFAULT_THRESHOLDS = DetectabilityThresholds()


def require_frozen_thresholds(thresholds: DetectabilityThresholds) -> None:
    """Refuse publication under any criterion other than the frozen policy."""
    if thresholds != DEFAULT_THRESHOLDS:
        raise FreezeGuardError(
            "publishable universe requires frozen detectability thresholds "
            f"{FROZEN_MIN_MEAN_TOTAL_JUNCTION_READS:g}/{SCRMOE_MIN_SAMPLES_WITH_PSI}"
        )


# --- rMATS SE parsing --------------------------------------------------------


def _cell_or_empty(cell) -> str:
    """Coerce an rMATS cell to a stripped string.

    ``pd.read_csv(dtype=str)`` still leaves genuinely-missing cells as a float
    ``NaN`` rather than an empty string, so a missing IncLevel/count cell reaches
    the parsers as ``float('nan')``. Treat that (and ``None``) as an empty cell,
    which the callers already handle as "no values" — identical to an empty or
    ``NA`` cell. Non-missing values are stringified verbatim.
    """
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return ""
    return str(cell).strip()


def _parse_int_list(cell) -> list[int]:
    cell = _cell_or_empty(cell)
    if not cell:
        return []
    return [int(round(float(part))) for part in cell.split(",") if part != ""]


def _parse_psi_list(cell) -> list[float]:
    values: list[float] = []
    for part in _cell_or_empty(cell).split(","):
        part = part.strip()
        if part == "" or part.upper() == "NA":
            continue
        values.append(float(part))
    return values


def load_se_events(
    se_table: Path, *, expected_sample_count: int = EXPECTED_CONTROL_COUNT
) -> pd.DataFrame:
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


def _gtf_attributes(cell: str) -> dict[str, str]:
    """Parse the quoted scalar attributes needed from one GENCODE GTF row."""
    return dict(re.findall(r'(\w+) "([^"]*)";', cell))


def _gtf_attribute_values(cell: str, name: str) -> tuple[str, ...]:
    """Return every value for a possibly repeated GTF attribute (for tags)."""
    return tuple(re.findall(rf'{re.escape(name)} "([^"]*)";', cell))


def load_annotated_internal_exons(gtf_path: Path) -> pd.DataFrame:
    """Load transcript-internal GENCODE exons and their annotated flanks."""
    records: list[dict] = []
    transcript_id: str | None = None
    transcript_exons: list[dict] = []
    completed_transcripts: set[str] = set()

    def flush() -> None:
        if not transcript_exons:
            return
        transcript_types = {row["transcript_type"] for row in transcript_exons}
        if len(transcript_types) != 1:
            raise ValueError(
                f"GENCODE transcript has inconsistent transcript_type: {transcript_id}"
            )
        tags = {
            tag for row in transcript_exons for tag in row["transcript_tags"]
        }
        if not transcript_types.issubset(ALLOWED_GENCODE_TRANSCRIPT_TYPES) or any(
            tag.startswith(EXCLUDED_GENCODE_TRANSCRIPT_TAG_PREFIXES) for tag in tags
        ):
            return
        if len(transcript_exons) < 3:
            return
        ordered = sorted(transcript_exons, key=lambda row: (row["start"], row["end"]))
        for index in range(1, len(ordered) - 1):
            lower, exon, upper = ordered[index - 1 : index + 2]
            records.append(
                {
                    "event_id": f"annotation:{exon['transcript_id']}:{exon['exon_id']}",
                    "gene_id": exon["gene_id"],
                    "gene_id_base": exon["gene_id"].split(".")[0],
                    "gene_symbol": exon["gene_symbol"],
                    "chrom": exon["chrom"],
                    "strand": exon["strand"],
                    "exon_start_0base": exon["start"] - 1,
                    "exon_end": exon["end"],
                    "upstream_ee": lower["end"],
                    "downstream_es": upper["start"] - 1,
                    "exon_length": exon["end"] - exon["start"] + 1,
                }
            )

    with Path(gtf_path).open("rt", encoding="utf-8") as handle:
        for line in handle:
            if "\texon\t" not in line:
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] not in ALLOWED_GENCODE_FEATURE_TYPES:
                continue
            attrs = _gtf_attributes(fields[8])
            current = attrs.get("transcript_id")
            if not current:
                raise ValueError("GENCODE exon lacks transcript_id")
            if transcript_id is not None and current != transcript_id:
                flush()
                completed_transcripts.add(transcript_id)
                transcript_exons = []
            if current in completed_transcripts:
                raise ValueError(f"GENCODE transcript is not contiguous in GTF: {current}")
            transcript_id = current
            transcript_exons.append(
                {
                    "chrom": fields[0],
                    "start": int(fields[3]),
                    "end": int(fields[4]),
                    "strand": fields[6],
                    "gene_id": attrs.get("gene_id", ""),
                    "gene_symbol": attrs.get("gene_name", ""),
                    "transcript_id": current,
                    "exon_id": attrs.get("exon_id", f"exon{len(transcript_exons) + 1}"),
                    "transcript_type": attrs.get("transcript_type", ""),
                    "transcript_tags": _gtf_attribute_values(fields[8], "tag"),
                }
            )
        flush()
    if not records:
        raise ValueError(f"no transcript-internal exons parsed from {gtf_path}")
    events = pd.DataFrame.from_records(records)
    events["exon_key"] = (
        events["chrom"].astype(str)
        + ":"
        + events["exon_start_0base"].astype(str)
        + "-"
        + events["exon_end"].astype(str)
        + ":"
        + events["strand"].astype(str)
    )
    return events


def load_control_junctions(sj_paths: list[Path]) -> list[dict[tuple, int]]:
    """Load unique-junction counts from exactly the three pinned ScrMOE runs."""
    paths = [Path(path) for path in sj_paths]
    observed_runs = [path.parent.name for path in paths]
    if len(paths) != EXPECTED_CONTROL_COUNT or set(observed_runs) != set(CONTROL_RUN_IDS):
        raise FreezeGuardError(
            "STAR junction inputs must be exactly the three ScrMOE runs: "
            + ", ".join(CONTROL_RUN_IDS)
        )
    junctions: list[dict[tuple, int]] = []
    for path, run in zip(paths, observed_runs):
        if (
            path.name != f"{run}.SJ.out.tab"
            or not path.is_file()
            or path.stat().st_size == 0
        ):
            raise FileNotFoundError(
                f"expected STAR junction input is missing or invalid: {path}"
            )
        raw = pd.read_csv(path, sep="\t", header=None)
        if raw.shape[1] != 9 or raw.isna().any().any():
            raise ValueError(f"STAR junction input must have 9 complete columns: {path}")
        numeric = raw[[1, 2, 3, 6]].to_numpy(dtype=float)
        if (
            not np.isfinite(numeric).all()
            or (numeric < 0).any()
            or not np.equal(numeric, np.floor(numeric)).all()
        ):
            raise ValueError(f"STAR junction coordinates/counts must be nonnegative integers: {path}")
        strand = raw[3].map({1: "+", 2: "-"})
        usable = raw.loc[strand.notna()].copy()
        usable[3] = strand[strand.notna()]
        keys = list(zip(usable[0], usable[1], usable[2], usable[3]))
        if len(keys) != len(set(keys)):
            raise ValueError(f"STAR junction input contains duplicate junctions: {path}")
        junctions.append({key: int(count) for key, count in zip(keys, usable[6])})
    return junctions


def attach_annotated_junction_metrics(
    exons: pd.DataFrame,
    junctions_by_sample: list[dict[tuple, int]],
    *,
    read_length: int = 100,
) -> pd.DataFrame:
    """Apply the existing rMATS JC depth/PSI estimator to annotated exons."""
    if read_length <= 1:
        raise ValueError("read_length must be greater than one")
    records: list[dict] = []
    skip_form_len = read_length - 1
    for row in exons.itertuples(index=False):
        strand = str(row.strand)
        lower = (str(row.chrom), int(row.upstream_ee) + 1, int(row.exon_start_0base), strand)
        upper = (str(row.chrom), int(row.exon_end) + 1, int(row.downstream_es), strand)
        skipping = (str(row.chrom), int(row.upstream_ee) + 1, int(row.downstream_es), strand)
        inclusion = [sample.get(lower, 0) + sample.get(upper, 0) for sample in junctions_by_sample]
        skipped = [sample.get(skipping, 0) for sample in junctions_by_sample]
        totals = [inc + skip for inc, skip in zip(inclusion, skipped)]
        inc_form_len = min(2 * skip_form_len, skip_form_len + int(row.exon_length))
        psi = [
            (inc / inc_form_len) / ((inc / inc_form_len) + (skip / skip_form_len))
            for inc, skip in zip(inclusion, skipped)
            if inc + skip > 0
        ]
        records.append(
            {
                "mean_inclusion_junction_reads": _safe_mean(inclusion),
                "mean_skipping_junction_reads": _safe_mean(skipped),
                "mean_total_junction_reads": _safe_mean(totals),
                "n_samples_total": len(totals),
                "n_samples_with_psi": len(psi),
                "baseline_psi": _safe_mean(psi),
            }
        )
    return pd.concat([exons.reset_index(drop=True), pd.DataFrame.from_records(records)], axis=1)


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
    paths = [Path(path) for path in reads_per_gene_paths]
    observed_runs = [path.parent.name for path in paths]
    if len(paths) != EXPECTED_CONTROL_COUNT or set(observed_runs) != set(CONTROL_RUN_IDS):
        raise FreezeGuardError(
            "STAR expression inputs must be exactly the three ScrMOE runs: "
            + ", ".join(CONTROL_RUN_IDS)
        )
    if len(observed_runs) != len(set(observed_runs)):
        raise FreezeGuardError("duplicate ScrMOE STAR expression input")

    per_sample_cpm: list[pd.Series] = []
    for path, run in zip(paths, observed_runs):
        if path.name != f"{run}.ReadsPerGene.out.tab":
            raise FreezeGuardError(f"unexpected STAR count filename for {run}: {path.name}")
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"expected STAR count input is missing or empty: {path}")
        raw = pd.read_csv(
            path,
            sep="\t",
            header=None,
            dtype=str,
            keep_default_na=False,
        )
        if raw.shape[1] != 4:
            raise ValueError(f"STAR count input must have exactly 4 columns: {path}")
        raw.columns = ["gene_id", "unstranded", "forward", "reverse"]
        # Normalize compatibility forms, then fail closed on invisible/control
        # characters and anything outside STAR/GENCODE's ASCII identifier form.
        gene_ids = raw["gene_id"].astype(str).map(
            lambda value: unicodedata.normalize("NFKC", value)
        )
        has_forbidden_unicode = gene_ids.map(
            lambda value: any(unicodedata.category(char)[0] in {"C", "Z"} for char in value)
        )
        valid_grammar = gene_ids.str.fullmatch(r"[A-Za-z0-9_][A-Za-z0-9_.:-]*")
        if (
            gene_ids.eq("").any()
            or has_forbidden_unicode.any()
            or not valid_grammar.all()
            or gene_ids.duplicated().any()
        ):
            raise ValueError(
                f"STAR count input has invalid or duplicate gene IDs after normalization: {path}"
            )
        numeric = raw[["unstranded", "forward", "reverse"]].apply(
            pd.to_numeric, errors="coerce"
        )
        values = numeric.to_numpy(dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f"STAR count input has non-numeric or non-finite counts: {path}")
        if (values < 0).any() or not np.equal(values, np.floor(values)).all():
            raise ValueError(f"STAR counts must be finite nonnegative integers: {path}")
        table = pd.concat([gene_ids.rename("gene_id"), numeric.astype("int64")], axis=1)
        table = table[~table["gene_id"].str.startswith("N_")].copy()
        library_size = table["reverse"].sum()
        if library_size <= 0:
            raise ValueError(f"STAR count input has non-positive reverse library size: {path}")
        cpm = table.set_index("gene_id")["reverse"] / library_size * 1e6
        per_sample_cpm.append(cpm)
    mean_cpm = pd.concat(per_sample_cpm, axis=1).mean(axis=1)
    out = mean_cpm.reset_index()
    out.columns = ["gene_id", "control_mean_cpm"]
    out["gene_id_base"] = out["gene_id"].str.split(".").str[0]
    return out[["gene_id_base", "control_mean_cpm"]]


def attach_universe_provenance(events: pd.DataFrame) -> pd.DataFrame:
    """Attach the immutable ScrMOE cohort and detectability-policy contract."""
    out = events.copy()
    for column, expected in UNIVERSE_PROVENANCE.items():
        out[column] = expected
    return out


def attach_input_provenance(
    events: pd.DataFrame,
    gtf_path: Path,
    sj_paths: list[Path],
    genome_fasta: Path,
    reads_per_gene_paths: list[Path],
    read_length: int,
) -> pd.DataFrame:
    """Bind every identity and matching-covariate input to its exact content."""
    expected = _input_provenance_values(
        gtf_path, sj_paths, genome_fasta, reads_per_gene_paths, read_length
    )
    out = events.copy()
    for column, value in expected.items():
        out[column] = value
    return out


def _input_provenance_values(
    gtf_path: Path,
    sj_paths: list[Path],
    genome_fasta: Path,
    reads_per_gene_paths: list[Path],
    read_length: int,
) -> dict[str, str | int]:
    """Return the authoritative hashes for every universe/match input."""
    sj_hashes = {
        path.name.removesuffix(".SJ.out.tab"): _sha256(path)
        for path in sj_paths
    }
    if set(sj_hashes) != set(CONTROL_RUN_IDS):
        raise FreezeGuardError(
            "STAR junction provenance differs from the pinned ScrMOE run IDs"
        )
    count_hashes = {
        path.name.removesuffix(".ReadsPerGene.out.tab"): _sha256(path)
        for path in reads_per_gene_paths
    }
    if set(count_hashes) != set(CONTROL_RUN_IDS):
        raise FreezeGuardError(
            "STAR expression provenance differs from the pinned ScrMOE run IDs"
        )
    return {
        "gencode_gtf_sha256": _sha256(gtf_path),
        "scrmoe_sj_sha256": json.dumps(
            sj_hashes, sort_keys=True, separators=(",", ":")
        ),
        "genome_fasta_sha256": _sha256(genome_fasta),
        "scrmoe_reads_per_gene_sha256": json.dumps(
            count_hashes, sort_keys=True, separators=(",", ":")
        ),
        "junction_read_length": int(read_length),
        "junction_read_length_sha256": hashlib.sha256(
            f"{int(read_length)}\n".encode("ascii")
        ).hexdigest(),
    }


def validate_input_provenance(
    events: pd.DataFrame,
    gtf_path: Path,
    sj_paths: list[Path],
    genome_fasta: Path,
    reads_per_gene_paths: list[Path],
    read_length: int,
) -> dict[str, str | int]:
    """Fail if an artifact is not bound to the current authoritative inputs."""
    expected = _input_provenance_values(
        gtf_path, sj_paths, genome_fasta, reads_per_gene_paths, read_length
    )
    missing = sorted(set(expected).difference(events.columns))
    if missing:
        raise FreezeGuardError(
            f"artifact lacks authoritative input provenance columns: {missing}"
        )
    for column, value in expected.items():
        if not events[column].eq(value).all():
            raise FreezeGuardError(
                f"authoritative input provenance mismatch in {column}"
            )
    return expected


def assert_scrmoe_universe(events: pd.DataFrame) -> None:
    """Require every row to carry the exact ScrMOE-only provenance contract."""
    missing = set(UNIVERSE_PROVENANCE_COLUMNS).difference(events.columns)
    if missing:
        raise FreezeGuardError(
            "universe lacks ScrMOE provenance columns: " + ", ".join(sorted(missing))
        )
    if events.empty:
        raise FreezeGuardError("ScrMOE universe is empty")
    for column, expected in UNIVERSE_PROVENANCE.items():
        if not events[column].eq(expected).all():
            raise FreezeGuardError(
                f"universe has invalid ScrMOE provenance in {column!r}"
            )


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
    """Collapse alternate transcript-flank records to one genomic exon identity.

    GENCODE may annotate the same exon repeatedly with different transcript
    flanks. The representative record is chosen
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
        .rename("annotation_event_ids")
    )
    multiplicity = grouped.size().rename("n_annotation_events")

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
    sj_paths: list[Path],
    gtf_path: Path,
    genome_fasta: Path,
    thresholds: DetectabilityThresholds = DEFAULT_THRESHOLDS,
    *,
    samtools: Path | None = None,
    read_length: int = 100,
) -> pd.DataFrame:
    """Build one fail-closed row per annotated transcript-internal exon."""
    require_frozen_thresholds(thresholds)
    # rMATS remains a provenance input, but does not seed candidate identity:
    # constitutive annotated exons need not appear in an SE event table.
    rmats = load_se_events(se_table)
    events = load_annotated_internal_exons(gtf_path)
    events = attach_annotated_junction_metrics(
        events, load_control_junctions(sj_paths), read_length=read_length
    )
    expression = load_control_expression(reads_per_gene_paths)
    events = events.merge(expression, on="gene_id_base", how="left")
    events = classify_detectable(events, thresholds)
    events = collapse_to_unique_exons(events)
    rmats_ids = rmats.groupby("exon_key")["event_id"].agg(
        lambda values: ",".join(sorted({str(value) for value in values}))
    )
    rmats_counts = rmats.groupby("exon_key").size()
    events["rmats_event_ids"] = events["exon_key"].map(rmats_ids).fillna("")
    events["n_rmats_events"] = events["exon_key"].map(rmats_counts).fillna(0).astype(int)
    events = attach_gc_content(events, genome_fasta, samtools=samtools)
    events = attach_universe_provenance(events)
    events = attach_input_provenance(
        events,
        gtf_path,
        sj_paths,
        genome_fasta,
        reads_per_gene_paths,
        read_length,
    )
    assert_scrmoe_universe(events)
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
        "The background set of annotated transcript-internal exons that are",
        "adequately measured in the GSE307431 baseline libraries, built from the",
        "three scrambled MOE (ScrMOE) control alignments only. It supplies the matched",
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
        "(all ScrMOE).",
        "",
        "## Method [candidate]",
        "",
        "- Candidate identities and transcript flanks are seeded from every",
        "  transcript-internal exon in the GENCODE v50 annotation. Detectability",
        "  is therefore not gated by whether rMATS enumerates a skipped-exon event.",
        "- rMATS-turbo SE events are retained as source-event provenance. They were",
        "  enumerated from GENCODE v50 with `--task both`",
        "  and `--statoff`; both rMATS groups are controls. `--statoff` disables",
        "  statistical testing, but rMATS still emits a dummy-control",
        "  `IncLevelDifference` in some source rows. The builder ignores and",
        "  explicitly rejects all contrast fields before parquet publication.",
        "  Neither rMATS group supplies the published exon metrics; those metrics",
        "  are recomputed from the three ScrMOE STAR junction tables and never",
        "  combine the three ASO chemistries.",
        "- Alternate annotated transcript flanks are collapsed by the genomic exon key",
        "  `chrom:exonStart_0base-exonEnd:strand`; each parquet row is one unique",
        "  exon interval. Source annotation IDs and multiplicity remain in audit columns.",
        "  If one interval is annotated to multiple genes, all gene IDs/symbols and",
        "  per-gene control CPM values are retained in aggregate provenance columns;",
        "  `multi_gene` makes that ambiguity explicit.",
        "- Inclusion level (PSI) uses the rMATS JC estimator applied directly to",
        "  annotated exon flanks and the three ScrMOE STAR junction tables; an",
        "  exon's junction depth is the mean of inclusion-junction plus",
        "  skipping-junction uniquely counted reads across control samples.",
        "- Internal exons only: candidates come from allowed GENCODE transcript",
        "  classes with at least three annotated exons. Sorting each transcript",
        "  by genomic position and retaining only positions between its first and",
        "  last exon excludes the 5'-most and 3'-most exons by construction.",
        "- Gene expression is mean control CPM from STAR reverse-strand gene",
        "  counts (fr-firststrand). GC fraction and exon length come from the",
        "  reference exon coordinates.",
        f"- Frozen detectability policy `{DETECTABILITY_POLICY_ID}`: {DETECTABILITY_POLICY}",
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
    lines.append(f"- Reference FASTA SHA-256: `{_sha256(genome_fasta)}`")
    lines.append(f"- Junction read length: `{int(events['junction_read_length'].iloc[0])}`")
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


BUNDLE_DATA_FILES = (
    "detectable_universe_v50.parquet",
    "detectable_universe_v50_README.md",
)
BUNDLE_MANIFEST = "manifest.json"


def write_bundle_manifest(bundle_dir: Path) -> Path:
    """Hash the complete staged data bundle and write its manifest atomically."""
    bundle_dir = Path(bundle_dir)
    files: dict[str, dict[str, str | int]] = {}
    for name in BUNDLE_DATA_FILES:
        path = bundle_dir / name
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"staged artifact is missing or empty: {path}")
        files[name] = {"sha256": _sha256(path), "size": path.stat().st_size}
    payload = json.dumps(
        {"schema_version": 1, "files": files}, indent=2, sort_keys=True
    ) + "\n"
    handle = tempfile.NamedTemporaryFile(
        prefix=f".{BUNDLE_MANIFEST}.", dir=bundle_dir, delete=False
    )
    temporary = Path(handle.name)
    try:
        handle.write(payload.encode("utf-8"))
        handle.close()
        os.replace(temporary, bundle_dir / BUNDLE_MANIFEST)
    finally:
        if not handle.closed:
            handle.close()
        temporary.unlink(missing_ok=True)
    return bundle_dir / BUNDLE_MANIFEST


def validate_bundle_manifest(bundle_dir: Path) -> None:
    """Fail unless the manifest exactly authenticates both bundle data files."""
    bundle_dir = Path(bundle_dir)
    manifest_path = bundle_dir / BUNDLE_MANIFEST
    if not manifest_path.is_file() or manifest_path.stat().st_size == 0:
        raise FileNotFoundError(f"bundle manifest is missing or empty: {manifest_path}")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"could not read bundle manifest: {manifest_path}") from exc
    if manifest.get("schema_version") != 1:
        raise ValueError("bundle manifest has unsupported schema_version")
    files = manifest.get("files")
    if not isinstance(files, dict) or set(files) != set(BUNDLE_DATA_FILES):
        raise ValueError("bundle manifest does not name the exact published data files")
    for name in BUNDLE_DATA_FILES:
        path = bundle_dir / name
        record = files[name]
        if not path.is_file() or not isinstance(record, dict):
            raise ValueError(f"bundle manifest entry is invalid: {name}")
        if record.get("size") != path.stat().st_size:
            raise ValueError(f"bundle manifest size mismatch: {name}")
        if record.get("sha256") != _sha256(path):
            raise ValueError(f"bundle manifest SHA-256 mismatch: {name}")


def _write_atomic(path: Path, content: bytes) -> None:
    handle = tempfile.NamedTemporaryFile(
        prefix=f".{path.name}.", dir=path.parent, delete=False
    )
    temporary = Path(handle.name)
    try:
        handle.write(content)
        handle.close()
        os.replace(temporary, path)
    finally:
        if not handle.closed:
            handle.close()
        temporary.unlink(missing_ok=True)


def _temporary_symlink(link: Path, target: Path) -> Path:
    handle = tempfile.NamedTemporaryFile(
        prefix=f".{link.name}.", dir=link.parent, delete=False
    )
    temporary = Path(handle.name)
    handle.close()
    temporary.unlink()
    temporary.symlink_to(target, target_is_directory=True)
    return temporary


def _snapshot_entry(path: Path) -> tuple[str, bytes | str | None]:
    """Capture enough state to atomically restore a file, symlink, or absence."""
    if path.is_symlink():
        return "symlink", os.readlink(path)
    if not path.exists():
        return "missing", None
    if path.is_file():
        return "file", path.read_bytes()
    raise OSError(f"publish entry is neither a file nor symlink: {path}")


def _restore_entry(path: Path, snapshot: tuple[str, bytes | str | None]) -> None:
    """Restore a captured entry from its observed current filesystem state."""
    kind, value = snapshot
    if kind == "missing":
        if path.is_symlink() or path.exists():
            path.unlink()
        return
    if kind == "file":
        if not isinstance(value, bytes):
            raise TypeError("file snapshot lacks bytes")
        if not path.is_symlink() and path.is_file() and path.read_bytes() == value:
            return
        _write_atomic(path, value)
        return
    if kind == "symlink":
        if not isinstance(value, str):
            raise TypeError("symlink snapshot lacks a target")
        if path.is_symlink() and os.readlink(path) == value:
            return
        temporary = _temporary_symlink(path, Path(value))
        try:
            os.replace(temporary, path)
        finally:
            temporary.unlink(missing_ok=True)
        return
    raise ValueError(f"unknown publish-entry snapshot kind: {kind}")


def _resolved_pointer_target(pointer: Path) -> Path | None:
    """Return an existing pointer target, distinguishing absence from dangling."""
    if pointer.is_symlink():
        try:
            return pointer.resolve(strict=True)
        except OSError as exc:
            raise OSError(f"current bundle pointer is dangling: {pointer}") from exc
    if pointer.exists():
        raise OSError(f"current bundle pointer is not a symlink: {pointer}")
    return None


def _bundle_done_matches_manifest(bundle: Path) -> bool:
    done = bundle / ".done"
    if not done.is_file():
        return False
    return done.read_text(encoding="ascii").strip() == _sha256(
        bundle / BUNDLE_MANIFEST
    )


def publish_staged_bundle(
    staged_bundle: Path,
    versioned_bundle: Path,
    current_pointer: Path,
    legacy_done_file: Path,
) -> None:
    """Publish a validated bundle with one atomic ``current`` pointer switch."""
    staged_bundle = Path(staged_bundle)
    versioned_bundle = Path(versioned_bundle)
    current_pointer = Path(current_pointer)
    legacy_done_file = Path(legacy_done_file)
    versioned_bundle.parent.mkdir(parents=True, exist_ok=True)
    current_pointer.parent.mkdir(parents=True, exist_ok=True)
    if staged_bundle.parent.stat().st_dev != versioned_bundle.parent.stat().st_dev:
        raise OSError("staged and versioned bundles must share a filesystem")

    old_bundle = _resolved_pointer_target(current_pointer)
    versioned_exists = versioned_bundle.exists()
    staged_exists = staged_bundle.exists()
    if versioned_exists and old_bundle == versioned_bundle.resolve(strict=True):
        validate_bundle_manifest(versioned_bundle)
        if not _bundle_done_matches_manifest(versioned_bundle):
            raise ValueError("published bundle completion sentinel is invalid")
        if staged_exists:
            validate_bundle_manifest(staged_bundle)
            if (staged_bundle / BUNDLE_MANIFEST).read_bytes() != (
                versioned_bundle / BUNDLE_MANIFEST
            ).read_bytes():
                raise FileExistsError(
                    f"versioned bundle already exists with different content: {versioned_bundle}"
                )
        legacy_done_file.unlink(missing_ok=True)
        return
    if staged_exists and versioned_exists:
        raise FileExistsError(f"versioned bundle already exists: {versioned_bundle}")
    candidate = staged_bundle if staged_exists else versioned_bundle
    if not candidate.exists():
        raise FileNotFoundError(f"staged bundle is missing: {staged_bundle}")
    validate_bundle_manifest(candidate)

    manifest_digest = _sha256(candidate / BUNDLE_MANIFEST)
    _write_atomic(candidate / ".done", f"{manifest_digest}\n".encode("ascii"))
    old_done = old_bundle / ".done" if old_bundle is not None else None
    old_done_snapshot = _snapshot_entry(old_done) if old_done is not None else None
    legacy_done_snapshot = _snapshot_entry(legacy_done_file)
    compatibility_snapshots = {
        current_pointer.parent / name: _snapshot_entry(current_pointer.parent / name)
        for name in BUNDLE_DATA_FILES
    }
    pointer_target = Path(os.path.relpath(versioned_bundle, current_pointer.parent))
    temporary_pointer = _temporary_symlink(current_pointer, pointer_target)
    try:
        if staged_bundle.exists():
            os.replace(staged_bundle, versioned_bundle)
        elif not versioned_bundle.exists():
            raise FileNotFoundError("publish candidate disappeared before rename")
        if old_done is not None:
            old_done.unlink(missing_ok=True)
        legacy_done_file.unlink(missing_ok=True)
        for name in BUNDLE_DATA_FILES:
            compatibility = current_pointer.parent / name
            target = Path(current_pointer.name) / name
            temporary = _temporary_symlink(compatibility, target)
            try:
                os.replace(temporary, compatibility)
            finally:
                temporary.unlink(missing_ok=True)
        os.replace(temporary_pointer, current_pointer)
        if old_done is not None and old_done_snapshot is not None:
            _restore_entry(old_done, old_done_snapshot)
    except BaseException:
        try:
            pointer_is_published = (
                current_pointer.is_symlink()
                and current_pointer.resolve(strict=True)
                == versioned_bundle.resolve(strict=True)
            )
        except OSError:
            pointer_is_published = False
        if pointer_is_published:
            if old_done is not None and old_done_snapshot is not None:
                _restore_entry(old_done, old_done_snapshot)
        else:
            for compatibility, snapshot in compatibility_snapshots.items():
                _restore_entry(compatibility, snapshot)
            if old_done is not None and old_done_snapshot is not None:
                _restore_entry(old_done, old_done_snapshot)
            _restore_entry(legacy_done_file, legacy_done_snapshot)
            if versioned_bundle.exists() and not staged_bundle.exists():
                os.replace(versioned_bundle, staged_bundle)
        raise
    finally:
        temporary_pointer.unlink(missing_ok=True)


def publish_staged_pair(
    staged_first: Path,
    staged_second: Path,
    out_first: Path,
    out_second: Path,
    done_file: Path,
    done_contents: str = "",
) -> None:
    """Publish a legacy two-file pair with the validity sentinel removed."""
    staged = (Path(staged_first), Path(staged_second))
    final = (Path(out_first), Path(out_second))
    done_file = Path(done_file)
    for source, destination in zip(staged, final):
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(f"staged artifact is missing or empty: {source}")
        if source.parent.stat().st_dev != destination.parent.stat().st_dev:
            raise OSError(
                f"staged artifact must share a filesystem with destination: {source}"
            )
    backups: list[Path | None] = []
    for destination in final:
        if destination.exists():
            handle = tempfile.NamedTemporaryFile(
                prefix=f".{destination.name}.backup-",
                dir=destination.parent,
                delete=False,
            )
            backup = Path(handle.name)
            handle.close()
            backup.unlink()
            os.link(destination, backup)
            backups.append(backup)
        else:
            backups.append(None)
    published: list[Path] = []
    done_file.unlink(missing_ok=True)
    try:
        for source, destination in zip(staged, final):
            os.replace(source, destination)
            published.append(destination)
        _write_atomic(done_file, done_contents.encode("utf-8"))
    except BaseException:
        done_file.unlink(missing_ok=True)
        for destination, backup in reversed(list(zip(final, backups))):
            if backup is not None and backup.exists():
                os.replace(backup, destination)
            elif destination in published:
                destination.unlink(missing_ok=True)
        raise
    finally:
        for backup in backups:
            if backup is not None:
                backup.unlink(missing_ok=True)


# --- CLI ---------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--conditions", required=True, type=Path)
    parser.add_argument("--se-table", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--star-dir", required=True, type=Path)
    parser.add_argument("--genome-fasta", required=True, type=Path)
    parser.add_argument("--samtools", type=Path, default=None)
    parser.add_argument("--read-length", type=int, default=100)
    parser.add_argument("--b1", type=Path, default=None)
    parser.add_argument("--b2", type=Path, default=None)
    parser.add_argument("--out-parquet", required=True, type=Path)
    parser.add_argument("--out-readme", required=True, type=Path)
    parser.add_argument(
        "--min-mean-junction-reads",
        type=float,
        default=DEFAULT_THRESHOLDS.min_mean_total_junction_reads,
        choices=[FROZEN_MIN_MEAN_TOTAL_JUNCTION_READS],
        help="frozen publishable threshold; other values are refused",
    )
    parser.add_argument(
        "--min-samples-with-psi",
        type=int,
        default=DEFAULT_THRESHOLDS.min_samples_with_psi,
        choices=[SCRMOE_MIN_SAMPLES_WITH_PSI],
        help="frozen publishable threshold; other values are refused",
    )
    args = parser.parse_args(argv)

    conditions = pd.read_csv(args.conditions, sep="\t", dtype=str)
    control_runs = select_control_runs(conditions)

    reads_per_gene = []
    junctions = []
    for run in control_runs["run_accession"]:
        path = args.star_dir / run / f"{run}.ReadsPerGene.out.tab"
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(
                f"expected STAR count input for ScrMOE run {run} is missing or empty: {path}"
            )
        reads_per_gene.append(path)
        junction_path = args.star_dir / run / f"{run}.SJ.out.tab"
        if not junction_path.is_file() or junction_path.stat().st_size == 0:
            raise FileNotFoundError(
                f"expected STAR junction input for ScrMOE run {run} is missing or empty: {junction_path}"
            )
        junctions.append(junction_path)

    thresholds = DetectabilityThresholds(
        min_mean_total_junction_reads=args.min_mean_junction_reads,
        min_samples_with_psi=args.min_samples_with_psi,
    )
    events = build_universe(
        args.se_table,
        reads_per_gene,
        junctions,
        args.gtf,
        args.genome_fasta,
        thresholds,
        samtools=args.samtools,
        read_length=args.read_length,
    )

    args.out_parquet.parent.mkdir(parents=True, exist_ok=True)
    assert_no_contrast_fields(events)
    if events["exon_key"].duplicated().any():
        raise AssertionError("refusing to publish duplicate exon identities")
    if events["gc_fraction"].isna().any():
        raise GenomeFastaError("refusing to publish null GC fractions")
    events.to_parquet(args.out_parquet, index=False)
    provenance_files: list[tuple[str, Path]] = [
        ("GENCODE GTF", args.gtf),
        *[(f"ScrMOE STAR junction {path.name.removesuffix('.SJ.out.tab')}", path)
          for path in junctions],
        *[(
            f"ScrMOE STAR gene count {path.name.removesuffix('.ReadsPerGene.out.tab')}",
            path,
        ) for path in reads_per_gene],
        ("Conditions table", args.conditions),
    ]
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
