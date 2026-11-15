"""Physical-specificity metrics for the conference ISS-N1 ASO candidates.

The module consumes exact genomic hits for already IntaRNA-scored target
variants.  It never substitutes zero for an unavailable, failed, or saturated
measurement.
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
import math
from pathlib import Path
from typing import Iterable, Mapping, TextIO

from feng.titration.evaluate import rnaplfold_accessibility

from . import CANDIDATE
from .item1 import _write_tsv
from .offtarget import GeneIndex, parse_bowtie_line


GAS_CONSTANT_KCAL_PER_MOL_K = 1.987204e-3
TEMPERATURE_K = 310.15
CONCENTRATIONS_MOLAR = (1e-9, 1e-8, 1e-7, 1e-6)
ABSOLUTE_DG_THRESHOLDS = (-8.0, -10.0, -12.0, -14.0, -16.0)
ACCESSIBILITY_CUTOFF = 0.5
SMN2_GENE_ID = "ENSG00000205571.17"
ASO_ORDER = ("F18", "F14", "F10", "3UP8")


def fractional_occupancy(dg_kcal_mol: float, concentration_molar: float) -> float:
    """Return A/(A+Kd), with Kd=exp(dG/RT) in molar units."""
    if concentration_molar <= 0:
        raise ValueError("concentration must be positive")
    kd_molar = math.exp(float(dg_kcal_mol) / (GAS_CONSTANT_KCAL_PER_MOL_K * TEMPERATURE_K))
    return concentration_molar / (concentration_molar + kd_molar)


def _genes(row: Mapping) -> set[str]:
    return {gene for gene in str(row.get("compatible_gene_ids", "")).split(";") if gene}


def _ratio(numerator: float | None, denominator: float | None) -> float | None:
    if numerator is None or denominator is None or denominator == 0:
        return None
    return numerator / denominator


def occupancy_rows(
    loci: Iterable[Mapping],
    on_target_dg: Mapping[str, float],
    *,
    concentrations_molar: Iterable[float] = CONCENTRATIONS_MOLAR,
    accessibility_cutoff: float | None = None,
) -> list[dict]:
    """Summarize locus-multiplicity-weighted occupancy, optionally accessibility gated."""
    rows_by_aso = {aso: [] for aso in on_target_dg}
    saturated = set()
    for row in loci:
        aso = str(row["aso_id"])
        if row.get("site_class") == "saturated":
            saturated.add(aso)
        elif row.get("site_class") == "off_target" and row.get("occupancy_eligible_100x", True) in {True, "True", "true", "1", 1}:
            rows_by_aso.setdefault(aso, []).append(row)
    aso_ids = [aso for aso in ASO_ORDER if aso in rows_by_aso or aso in saturated]
    result = []
    for concentration in concentrations_molar:
        block = []
        for aso in aso_ids:
            base = {
                "state": CANDIDATE,
                "availability_state": "saturated" if aso in saturated else "available",
                "aso_id": aso,
                "concentration_molar": concentration,
                "concentration_label": {1e-9: "1nM", 1e-8: "10nM", 1e-7: "100nM", 1e-6: "1uM"}.get(concentration, str(concentration)),
                "temperature_k": TEMPERATURE_K,
                "gas_constant_kcal_per_mol_k": GAS_CONSTANT_KCAL_PER_MOL_K,
            }
            if aso in saturated:
                base.update({
                    "raw_off_target_locus_count": None,
                    "theta_on": None,
                    "sum_theta_off": None,
                    "theta_on_over_sum_theta_off": None,
                    "theta_on_over_theta_total": None,
                })
                block.append(base)
                continue
            all_off = rows_by_aso[aso]
            ungated_sum = sum(fractional_occupancy(float(row["dg"]), concentration) for row in all_off)
            selected = all_off
            if accessibility_cutoff is not None:
                selected = [
                    row for row in all_off
                    if row.get("accessibility") not in {None, ""}
                    and float(row["accessibility"]) >= accessibility_cutoff
                ]
                applicable = [row for row in all_off if _genes(row)]
                evaluated = sum(row.get("accessibility") not in {None, ""} for row in applicable)
                base.update({
                    "expression_state": "unavailable",
                    "expression_source": None,
                    "accessibility_cutoff": accessibility_cutoff,
                    "accessibility_applicable_locus_count": len(applicable),
                    "accessibility_evaluated_locus_count": evaluated,
                    "accessibility_failed_locus_count": len(applicable) - evaluated,
                    "accessible_off_target_locus_count": len(selected),
                    "accessible_distinct_gene_count": len(set().union(*(_genes(row) for row in selected))) if selected else 0,
                    "ungated_sum_theta_off": ungated_sum,
                })
            theta_on = fractional_occupancy(on_target_dg[aso], concentration)
            sum_off = sum(fractional_occupancy(float(row["dg"]), concentration) for row in selected)
            base.update({
                "raw_off_target_locus_count": len(all_off),
                "theta_on": theta_on,
                "sum_theta_off": sum_off,
                "theta_on_over_sum_theta_off": _ratio(theta_on, sum_off),
                "theta_on_over_theta_total": theta_on / (theta_on + sum_off),
            })
            block.append(base)
        lookup = {row["aso_id"]: row for row in block}
        f18 = lookup.get("F18", {})
        for row in block:
            row.update({
                "f10_over_f18_burden_ratio": None,
                "f14_over_f18_burden_ratio": None,
                "f10_over_f18_specificity_ratio": None,
                "f14_over_f18_specificity_ratio": None,
                "f10_over_f18_specificity_fraction_ratio": None,
                "f14_over_f18_specificity_fraction_ratio": None,
            })
        for comparison in ("F10", "F14"):
            if comparison not in lookup or not f18:
                continue
            comp = lookup[comparison]
            prefix = comparison.lower()
            values = {
                f"{prefix}_over_f18_burden_ratio": _ratio(comp.get("sum_theta_off"), f18.get("sum_theta_off")),
                f"{prefix}_over_f18_specificity_ratio": _ratio(comp.get("theta_on_over_sum_theta_off"), f18.get("theta_on_over_sum_theta_off")),
                f"{prefix}_over_f18_specificity_fraction_ratio": _ratio(comp.get("theta_on_over_theta_total"), f18.get("theta_on_over_theta_total")),
            }
            for row in block:
                row.update(values)
        result.extend(block)
    return result


def absolute_threshold_rows(
    loci: Iterable[Mapping], *, thresholds: Iterable[float] = ABSOLUTE_DG_THRESHOLDS
) -> list[dict]:
    """Count raw loci and distinct compatible genes under common dG bars."""
    grouped: dict[str, list[Mapping]] = {}
    saturated = set()
    for row in loci:
        aso = str(row["aso_id"])
        if row.get("site_class") == "saturated":
            saturated.add(aso)
        elif row.get("site_class") == "off_target":
            grouped.setdefault(aso, []).append(row)
    result = []
    for aso in ASO_ORDER:
        if aso not in grouped and aso not in saturated:
            continue
        for threshold in thresholds:
            if aso in saturated:
                count = genes = None
                availability = "saturated"
            else:
                selected = [row for row in grouped[aso] if float(row["dg"]) <= threshold]
                count = len(selected)
                genes = len(set().union(*(_genes(row) for row in selected))) if selected else 0
                availability = "available"
            result.append({
                "state": CANDIDATE,
                "availability_state": availability,
                "aso_id": aso,
                "absolute_dg_cutoff_kcal_mol": threshold,
                "raw_off_target_locus_count": count,
                "distinct_transcript_compatible_gene_count": genes,
            })
    return result


def load_query_scores(path: Path) -> dict[str, dict]:
    result = {}
    with Path(path).open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if str(row["intarna_no_duplex"]).lower() == "true":
                continue
            query_id = row["bowtie_query_id"]
            result[query_id] = {
                "aso_id": row["aso_id"],
                "dg": float(row["intarna_dg_kcal_mol"]),
                "is_exact": int(row["mismatch_count"]) == 0,
                "length": len(row["bowtie_query_rna_5to3"]),
                "target": row["variant_target_rna_5to3"],
                "query": row["bowtie_query_rna_5to3"],
                "viable_100x": str(row["viable_100x"]).lower() == "true",
            }
    return result


def write_replay_fasta(score_table: Path, output: Path) -> int:
    """Write the union needed by occupancy and the weakest absolute bar."""
    scores = load_query_scores(score_table)
    selected = [
        (query_id, row) for query_id, row in scores.items()
        if row["viable_100x"] or row["dg"] <= max(ABSOLUTE_DG_THRESHOLDS)
    ]
    temporary = Path(output).with_suffix(Path(output).suffix + ".tmp")
    with temporary.open("w") as handle:
        for query_id, row in selected:
            handle.write(f">{query_id}\n{row['query'].replace('U', 'T')}\n")
    temporary.replace(output)
    return len(selected)


def map_exact_loci(
    lines: Iterable[str], *, output: TextIO, gene_index: GeneIndex, query_scores: Mapping[str, Mapping]
) -> None:
    """Stream one record per exact genomic locus, retaining duplicates."""
    fieldnames = [
        "state", "aso_id", "query_id", "dg", "chrom", "start_1based",
        "end_1based", "transcript_strand", "target_5to3", "site_class",
        "compatible_gene_ids", "occupancy_eligible_100x",
    ]
    writer = csv.DictWriter(output, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    lengths = {query_id: int(row["length"]) for query_id, row in query_scores.items()}
    for line in lines:
        if not line.strip():
            continue
        hit = parse_bowtie_line(line, lengths)
        score = query_scores[hit.query_id]
        _any, compatible = gene_index.overlaps(hit)
        on_target = bool(score["is_exact"] and SMN2_GENE_ID in compatible)
        writer.writerow({
            "state": CANDIDATE,
            "aso_id": score["aso_id"],
            "query_id": hit.query_id,
            "dg": score["dg"],
            "chrom": hit.chrom,
            "start_1based": hit.start_0base + 1,
            "end_1based": hit.end_0base,
            "transcript_strand": "+" if hit.query_strand == "-" else "-",
            "target_5to3": score["target"],
            "site_class": "on_target" if on_target else "off_target",
            "compatible_gene_ids": ";".join(sorted(compatible)),
            "occupancy_eligible_100x": score["viable_100x"],
        })


def _read_tsv(path: Path) -> list[dict]:
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def add_accessibility(
    loci_path: Path, output: Path, *, rnaplfold: Path, genome_fasta: Path,
    samtools: Path, workers: int,
) -> None:
    rows = _read_tsv(loci_path)

    def evaluate(row: dict) -> tuple[str, float | None]:
        if (
            row["site_class"] != "off_target"
            or row.get("occupancy_eligible_100x", "False").lower() != "true"
            or not row["compatible_gene_ids"]
        ):
            return "not_applicable", None
        # rnaplfold_accessibility fail-closed-compares the rung target against a
        # GENOMIC (DNA-alphabet) fetch. Our targets are stored RNA (U), so they
        # must be transcribed back to DNA letters or the frame check rejects
        # every locus ("padded window does not contain the rung target").
        rung = {
            "chrom": row["chrom"], "start": row["start_1based"],
            "end": row["end_1based"], "strand": row["transcript_strand"],
            "target_5to3": row["target_5to3"].upper().replace("U", "T"),
        }
        try:
            return "available", rnaplfold_accessibility(
                rung, rnaplfold=rnaplfold, genome_fasta=genome_fasta, samtools=samtools
            )
        except Exception as exc:  # preserve the reason; never fabricate a value
            return f"failed:{type(exc).__name__}", None

    with ThreadPoolExecutor(max_workers=workers) as pool:
        evaluated = list(pool.map(evaluate, rows))
    for row, (state, value) in zip(rows, evaluated):
        row["accessibility_state"] = state
        row["accessibility"] = value
    _write_tsv(output, rows)


def finalize(loci_path: Path, item1_path: Path, output_dir: Path) -> None:
    loci = _read_tsv(loci_path)
    loci.append({"aso_id": "3UP8", "site_class": "saturated"})
    item1 = _read_tsv(item1_path)
    on_target = {
        row["aso_id"]: float(row["intarna_default_dg_kcal_mol"])
        for row in item1 if row["aso_id"] != "3UP8"
    }
    _write_tsv(output_dir / "item3_spec_occupancy.tsv", occupancy_rows(loci, on_target))
    _write_tsv(
        output_dir / "item3_spec_gated.tsv",
        occupancy_rows(loci, on_target, accessibility_cutoff=ACCESSIBILITY_CUTOFF),
    )
    _write_tsv(output_dir / "item3_spec_absolute.tsv", absolute_threshold_rows(loci))


def initialize_pending(output_dir: Path) -> None:
    """Publish honest placeholders until the required Slurm replay finishes."""
    occupancy = []
    gated = []
    for concentration in CONCENTRATIONS_MOLAR:
        for aso in ASO_ORDER:
            state = "saturated" if aso == "3UP8" else "pending_slurm"
            common = {
                "state": CANDIDATE, "availability_state": state, "aso_id": aso,
                "concentration_molar": concentration,
                "concentration_label": {1e-9: "1nM", 1e-8: "10nM", 1e-7: "100nM", 1e-6: "1uM"}[concentration],
                "temperature_k": TEMPERATURE_K,
                "gas_constant_kcal_per_mol_k": GAS_CONSTANT_KCAL_PER_MOL_K,
                "raw_off_target_locus_count": None, "theta_on": None,
                "sum_theta_off": None, "theta_on_over_sum_theta_off": None,
                "theta_on_over_theta_total": None,
                "f10_over_f18_burden_ratio": None, "f14_over_f18_burden_ratio": None,
                "f10_over_f18_specificity_ratio": None, "f14_over_f18_specificity_ratio": None,
                "f10_over_f18_specificity_fraction_ratio": None,
                "f14_over_f18_specificity_fraction_ratio": None,
            }
            occupancy.append(common)
            gated.append({
                **common, "expression_state": "unavailable", "expression_source": None,
                "accessibility_cutoff": ACCESSIBILITY_CUTOFF,
                "accessibility_applicable_locus_count": None,
                "accessibility_evaluated_locus_count": None,
                "accessibility_failed_locus_count": None,
                "accessible_off_target_locus_count": None,
                "accessible_distinct_gene_count": None, "ungated_sum_theta_off": None,
            })
    absolute = [
        {
            "state": CANDIDATE,
            "availability_state": "saturated" if aso == "3UP8" else "pending_slurm",
            "aso_id": aso, "absolute_dg_cutoff_kcal_mol": threshold,
            "raw_off_target_locus_count": None,
            "distinct_transcript_compatible_gene_count": None,
        }
        for aso in ASO_ORDER for threshold in ABSOLUTE_DG_THRESHOLDS
    ]
    _write_tsv(output_dir / "item3_spec_occupancy.tsv", occupancy)
    _write_tsv(output_dir / "item3_spec_gated.tsv", gated)
    _write_tsv(output_dir / "item3_spec_absolute.tsv", absolute)


def main() -> None:
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    prepare = commands.add_parser("prepare")
    prepare.add_argument("--scores", type=Path, required=True)
    prepare.add_argument("--output", type=Path, required=True)
    mapper = commands.add_parser("map")
    mapper.add_argument("--scores", type=Path, required=True)
    mapper.add_argument("--gene-table", type=Path, required=True)
    access = commands.add_parser("accessibility")
    access.add_argument("--loci", type=Path, required=True)
    access.add_argument("--output", type=Path, required=True)
    access.add_argument("--rnaplfold", type=Path, required=True)
    access.add_argument("--genome-fasta", type=Path, required=True)
    access.add_argument("--samtools", type=Path, required=True)
    access.add_argument("--workers", type=int, default=1)
    final = commands.add_parser("finalize")
    final.add_argument("--loci", type=Path, required=True)
    final.add_argument("--item1", type=Path, required=True)
    final.add_argument("--output-dir", type=Path, required=True)
    initialize = commands.add_parser("initialize")
    initialize.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        write_replay_fasta(args.scores, args.output)
    elif args.command == "map":
        import sys
        map_exact_loci(
            sys.stdin, output=sys.stdout, gene_index=GeneIndex.from_table(args.gene_table),
            query_scores=load_query_scores(args.scores),
        )
    elif args.command == "accessibility":
        add_accessibility(
            args.loci, args.output, rnaplfold=args.rnaplfold,
            genome_fasta=args.genome_fasta, samtools=args.samtools, workers=args.workers,
        )
    elif args.command == "finalize":
        finalize(args.loci, args.item1, args.output_dir)
    else:
        initialize_pending(args.output_dir)


if __name__ == "__main__":
    main()
