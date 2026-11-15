"""Energy-filtered Item-3 off-target analysis for the locked ISS-N1 ASOs."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
from itertools import combinations, product
import math
from pathlib import Path
from typing import Callable, Iterable, Mapping, TextIO

from . import CANDIDATE
from .item1 import _write_tsv
from .offtarget import GeneIndex, parse_bowtie_line
from .sequences import (
    ASOS,
    ISS_N1_CONTEXT,
    derive_target,
    reverse_complement_rna,
    to_rna,
)
from .thermo import intarna_many_targets


TEMPERATURE_C = 37.0
TEMPERATURE_K = 310.15
GAS_CONSTANT_KCAL_PER_MOL_K = 0.00198720425864083
THRESHOLD_FOLDS = (1, 3, 10, 30, 100)
PRIMARY_THRESHOLD_FOLD = 10
ASO_ORDER = ("F18", "F14", "F10")


@dataclass(frozen=True, slots=True)
class TargetVariant:
    sequence_rna_5to3: str
    mismatch_count: int


@dataclass(frozen=True, slots=True)
class EnergyQuery:
    aso_id: str
    length_nt: int
    viable_folds: tuple[int, ...]


def enumerate_target_variants(
    target_rna_5to3: str, *, max_mismatches: int = 3
):
    """Yield the complete unique RNA Hamming ball in deterministic order."""
    target = to_rna(target_rna_5to3)
    if max_mismatches < 0 or max_mismatches > len(target):
        raise ValueError("max_mismatches must be between zero and target length")
    alphabet = "ACGU"
    for mismatch_count in range(max_mismatches + 1):
        for positions in combinations(range(len(target)), mismatch_count):
            alternatives = [
                tuple(base for base in alphabet if base != target[position])
                for position in positions
            ]
            for replacements in product(*alternatives):
                sequence = list(target)
                for position, base in zip(positions, replacements):
                    sequence[position] = base
                yield TargetVariant("".join(sequence), mismatch_count)


def fold_margin_kcal_mol(max_fold: float) -> float:
    """Return RT ln(fold) at the predeclared 37 C temperature."""
    if max_fold < 1:
        raise ValueError("max_fold must be at least one")
    return GAS_CONSTANT_KCAL_PER_MOL_K * TEMPERATURE_K * math.log(max_fold)


def viable_at_fold(
    dg_kcal_mol: float,
    no_duplex: bool,
    *,
    on_target_dg: float,
    max_fold: float,
) -> bool:
    """Whether a variant is no more than ``max_fold`` weaker than on-target."""
    if no_duplex:
        return False
    return dg_kcal_mol <= on_target_dg + fold_margin_kcal_mol(max_fold)


def build_score_rows(
    aso_id: str,
    aso_rna_5to3: str,
    on_target_rna_5to3: str,
    variants: Iterable[TargetVariant],
    scores: Mapping[str, tuple[float, bool]],
) -> list[dict]:
    """Attach one apples-to-apples IntaRNA reference to an ASO's variants."""
    aso = to_rna(aso_rna_5to3)
    on_target = to_rna(on_target_rna_5to3)
    normalized_variants = list(variants)
    exact = [
        variant
        for variant in normalized_variants
        if variant.mismatch_count == 0 and variant.sequence_rna_5to3 == on_target
    ]
    if len(exact) != 1:
        raise ValueError("variants must contain the exact on-target exactly once")
    try:
        on_target_dg, on_target_no_duplex = scores[on_target]
    except KeyError as exc:
        raise ValueError("scores omit the exact on-target") from exc
    if on_target_no_duplex:
        raise ValueError(f"{aso_id} on-target has no favourable duplex")

    rt = GAS_CONSTANT_KCAL_PER_MOL_K * TEMPERATURE_K
    rows = []
    for variant in normalized_variants:
        sequence = to_rna(variant.sequence_rna_5to3)
        try:
            dg, no_duplex = scores[sequence]
        except KeyError as exc:
            raise ValueError(f"scores omit variant {sequence}") from exc
        delta = "" if no_duplex else dg - on_target_dg
        weakening = "" if no_duplex else math.exp((dg - on_target_dg) / rt)
        variant_id = (
            "EV_"
            + aso_id
            + "_"
            + hashlib.sha256(sequence.encode()).hexdigest()[:16]
        )
        row = {
            "state": CANDIDATE,
            "aso_id": aso_id,
            "aso_rna_5to3": aso,
            "on_target_rna_5to3": on_target,
            "variant_id": variant_id,
            "bowtie_query_id": variant_id + " [candidate]",
            "variant_target_rna_5to3": sequence,
            "bowtie_query_rna_5to3": reverse_complement_rna(sequence),
            "mismatch_count": variant.mismatch_count,
            "intarna_dg_kcal_mol": dg,
            "intarna_no_duplex": no_duplex,
            "on_target_dg_kcal_mol": on_target_dg,
            "delta_dg_vs_on_target_kcal_mol": delta,
            "equilibrium_weakening_fold": weakening,
        }
        for fold in THRESHOLD_FOLDS:
            row[f"viable_{fold}x"] = viable_at_fold(
                dg, no_duplex, on_target_dg=on_target_dg, max_fold=fold
            )
        rows.append(row)
    return rows


def score_variants(
    output_dir: Path,
    *,
    intarna: Path,
    threads: int = 1,
    aso_ids: Iterable[str] = ASO_ORDER,
    score_many: Callable = intarna_many_targets,
) -> list[dict]:
    """Enumerate and score every <=3-mismatch target variant for each ASO."""
    all_rows: list[dict] = []
    for aso_id in aso_ids:
        if aso_id not in ASO_ORDER:
            raise ValueError(f"energy-filtered analysis excludes ASO {aso_id!r}")
        aso = ASOS[aso_id]
        target = derive_target(aso.rna_sequence, ISS_N1_CONTEXT).sequence
        variants = list(enumerate_target_variants(target, max_mismatches=3))
        targets_by_id = {
            "EV_"
            + aso_id
            + "_"
            + hashlib.sha256(variant.sequence_rna_5to3.encode()).hexdigest()[:16]:
            variant.sequence_rna_5to3
            for variant in variants
        }
        id_scores = score_many(
            aso.rna_sequence,
            targets_by_id,
            intarna=Path(intarna),
            threads=threads,
        )
        scores_by_sequence = {
            sequence: id_scores[target_id]
            for target_id, sequence in targets_by_id.items()
        }
        all_rows.extend(
            build_score_rows(
                aso_id, aso.rna_sequence, target, variants, scores_by_sequence
            )
        )

    output_dir = Path(output_dir)
    _write_tsv(output_dir / "item3_energy_variant_scores.tsv", all_rows)
    fasta = output_dir / "item3_energy_viable_100x.fa"
    temporary = fasta.with_suffix(fasta.suffix + ".tmp")
    with temporary.open("w") as handle:
        for row in all_rows:
            if row["viable_100x"]:
                handle.write(
                    f">{row['bowtie_query_id']}\n"
                    f"{row['bowtie_query_rna_5to3'].replace('U', 'T')}\n"
                )
    temporary.replace(fasta)
    return all_rows


def _write_energy_raw_counts(
    path: Path, raw_counts: Mapping[tuple[str, int], int]
) -> None:
    temporary = Path(path).with_suffix(Path(path).suffix + ".tmp")
    with temporary.open("w") as handle:
        handle.write("aso_id\tthreshold_fold\texact_alignment_count\n")
        for (aso_id, fold), count in sorted(raw_counts.items()):
            handle.write(f"{aso_id}\t{fold}\t{count}\n")
    temporary.replace(path)


def map_energy_alignment_gene_keys(
    lines: Iterable[str],
    *,
    output: TextIO,
    raw_counts_path: Path,
    gene_index: GeneIndex,
    query_metadata: Mapping[str, EnergyQuery],
) -> None:
    """Map exact variant hits to threshold-level external-sort gene keys."""
    query_lengths = {
        query_id: metadata.length_nt for query_id, metadata in query_metadata.items()
    }
    groups = {
        (metadata.aso_id, fold)
        for metadata in query_metadata.values()
        for fold in metadata.viable_folds
    }
    raw_counts = {group: 0 for group in groups}
    for line in lines:
        if not line.strip():
            continue
        hit = parse_bowtie_line(line, query_lengths)
        metadata = query_metadata[hit.query_id]
        any_overlap, compatible_overlap = gene_index.overlaps(hit)
        for fold in metadata.viable_folds:
            raw_counts[(metadata.aso_id, fold)] += 1
            for gene_id in any_overlap:
                output.write(f"{metadata.aso_id}\t{fold}\tA\t{gene_id}\n")
            for gene_id in compatible_overlap:
                output.write(f"{metadata.aso_id}\t{fold}\tC\t{gene_id}\n")
    _write_energy_raw_counts(raw_counts_path, raw_counts)


def reduce_energy_gene_keys(
    lines: Iterable[str], *, raw_counts_path: Path
) -> dict[tuple[str, int], dict[str, int]]:
    """Reduce sorted unique keys to unioned gene counts for each threshold."""
    any_counts: dict[tuple[str, int], int] = {}
    compatible_counts: dict[tuple[str, int], int] = {}
    previous = None
    for line in lines:
        key = line.rstrip("\n")
        if not key or key == previous:
            continue
        previous = key
        fields = key.split("\t")
        if len(fields) != 4:
            raise ValueError(f"invalid energy gene key: {line!r}")
        aso_id, fold_text, overlap_class, _gene_id = fields
        group = (aso_id, int(fold_text))
        if overlap_class == "A":
            any_counts[group] = any_counts.get(group, 0) + 1
        elif overlap_class == "C":
            compatible_counts[group] = compatible_counts.get(group, 0) + 1
        else:
            raise ValueError(f"invalid overlap class: {overlap_class}")

    # The mapper publishes this atomic sidecar only after its stdout closes.
    # Reading it after stream exhaustion avoids a reducer-startup race.
    raw_counts: dict[tuple[str, int], int] = {}
    with Path(raw_counts_path).open() as handle:
        header = next(handle).rstrip("\n").split("\t")
        if header != ["aso_id", "threshold_fold", "exact_alignment_count"]:
            raise ValueError(f"unexpected energy raw-count header: {header}")
        for line in handle:
            aso_id, fold, count = line.rstrip("\n").split("\t")
            raw_counts[(aso_id, int(fold))] = int(count)
    unknown_groups = (set(any_counts) | set(compatible_counts)).difference(raw_counts)
    if unknown_groups:
        raise ValueError(f"energy gene keys contain unknown groups: {unknown_groups}")
    return {
        group: {
            "exact_alignment_count": raw_counts[group],
            "distinct_genes_any_strand": any_counts.get(group, 0),
            "distinct_genes_transcript_compatible": compatible_counts.get(group, 0),
        }
        for group in raw_counts
    }


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _as_bool(value: object) -> bool:
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise ValueError(f"invalid boolean value: {value!r}")


def load_energy_queries(score_table: Path) -> dict[str, EnergyQuery]:
    """Load only the loosest-cutoff variants and retain stricter memberships."""
    metadata = {}
    for row in _read_tsv(score_table):
        if not _as_bool(row["viable_100x"]):
            continue
        query_id = row["bowtie_query_id"]
        if query_id in metadata:
            raise ValueError(f"duplicate energy variant identifier: {query_id}")
        viable_folds = tuple(
            fold for fold in THRESHOLD_FOLDS if _as_bool(row[f"viable_{fold}x"])
        )
        metadata[query_id] = EnergyQuery(
            row["aso_id"], len(row["bowtie_query_rna_5to3"]), viable_folds
        )
    if not metadata:
        raise ValueError("score table contains no 100x-viable variants")
    return metadata


def write_energy_counts(
    counts: Mapping[tuple[str, int], Mapping[str, int]], output: Path
) -> None:
    order = {aso_id: index for index, aso_id in enumerate(ASO_ORDER)}
    rows = [
        {
            "state": CANDIDATE,
            "aso_id": aso_id,
            "threshold_fold": fold,
            **counts[(aso_id, fold)],
        }
        for aso_id, fold in sorted(
            counts, key=lambda group: (order.get(group[0], len(order)), group[1])
        )
    ]
    _write_tsv(Path(output), rows)


def finalize_energy(
    *,
    score_table: Path,
    exact_counts: Path,
    raw_gene_counts: Path,
    output_dir: Path,
) -> list[dict]:
    """Publish the predeclared sweep without mixing raw and derived counts."""
    scores = _read_tsv(score_table)
    exact_lookup = {
        (row["aso_id"], int(row["threshold_fold"])): row
        for row in _read_tsv(exact_counts)
    }
    raw_lookup = {
        row["aso_id"]: int(row["bowtie1_v3_distinct_genes"])
        for row in _read_tsv(raw_gene_counts)
    }
    sweep = []
    for aso_id in ASO_ORDER:
        aso_scores = [row for row in scores if row["aso_id"] == aso_id]
        if not aso_scores:
            raise ValueError(f"score table omits {aso_id}")
        on_target_values = {
            float(row["on_target_dg_kcal_mol"]) for row in aso_scores
        }
        if len(on_target_values) != 1:
            raise ValueError(f"inconsistent on-target energy for {aso_id}")
        if aso_id not in raw_lookup:
            raise ValueError(f"raw gene-count table omits {aso_id}")
        for fold in THRESHOLD_FOLDS:
            key = (aso_id, fold)
            if key not in exact_lookup:
                raise ValueError(f"exact gene-count table omits {key}")
            exact = exact_lookup[key]
            viable_count = sum(
                _as_bool(row[f"viable_{fold}x"]) for row in aso_scores
            )
            sweep.append(
                {
                    "state": CANDIDATE,
                    "aso_id": aso_id,
                    "threshold_fold": fold,
                    "threshold_margin_kcal_mol": fold_margin_kcal_mol(fold),
                    "threshold_is_primary": fold == PRIMARY_THRESHOLD_FOLD,
                    "on_target_dg_kcal_mol": next(iter(on_target_values)),
                    "total_variants_hamming_le_3": len(aso_scores),
                    "viable_variant_count": viable_count,
                    "viable_variant_fraction": viable_count / len(aso_scores),
                    "exact_alignment_count": int(exact["exact_alignment_count"]),
                    "energy_filtered_distinct_genes_any_strand": int(
                        exact["distinct_genes_any_strand"]
                    ),
                    "energy_filtered_distinct_genes": int(
                        exact["distinct_genes_transcript_compatible"]
                    ),
                    "raw_bowtie1_v3_distinct_genes": raw_lookup[aso_id],
                }
            )

    output_dir = Path(output_dir)
    _write_tsv(output_dir / "item3_energy_filtered_sweep.tsv", sweep)
    primary = {
        row["aso_id"]: row
        for row in sweep
        if row["threshold_fold"] == PRIMARY_THRESHOLD_FOLD
    }
    f18, f10 = primary["F18"], primary["F10"]
    confirmed = (
        f10["viable_variant_count"] < f18["viable_variant_count"]
        and f10["energy_filtered_distinct_genes"]
        < f18["energy_filtered_distinct_genes"]
    )
    headline = [
        {
            "state": CANDIDATE,
            "threshold_rule": "variant_dG <= on_target_dG + RT*ln(max_fold)",
            "temperature_c": TEMPERATURE_C,
            "primary_threshold_fold": PRIMARY_THRESHOLD_FOLD,
            "primary_threshold_margin_kcal_mol": fold_margin_kcal_mol(
                PRIMARY_THRESHOLD_FOLD
            ),
            "f18_raw_bowtie1_v3_distinct_genes": f18[
                "raw_bowtie1_v3_distinct_genes"
            ],
            "f10_raw_bowtie1_v3_distinct_genes": f10[
                "raw_bowtie1_v3_distinct_genes"
            ],
            "f18_viable_variant_count": f18["viable_variant_count"],
            "f10_viable_variant_count": f10["viable_variant_count"],
            "f18_energy_filtered_distinct_genes": f18[
                "energy_filtered_distinct_genes"
            ],
            "f10_energy_filtered_distinct_genes": f10[
                "energy_filtered_distinct_genes"
            ],
            "prediction": (
                "F10 has fewer viable variants and fewer energy-filtered "
                "transcript-compatible genes than F18 at 10x"
            ),
            "prediction_result": "confirmed" if confirmed else "not_confirmed",
        }
    ]
    _write_tsv(
        output_dir / "item3_energy_filtered_headline_f18_vs_f10.tsv", headline
    )
    return sweep


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    score = subparsers.add_parser("score")
    score.add_argument("--output-dir", type=Path, required=True)
    score.add_argument("--intarna", type=Path, required=True)
    score.add_argument("--threads", type=int, default=1)
    map_parser = subparsers.add_parser("map")
    map_parser.add_argument("--gene-table", type=Path, required=True)
    map_parser.add_argument("--score-table", type=Path, required=True)
    map_parser.add_argument("--raw-counts", type=Path, required=True)
    reduce_parser = subparsers.add_parser("reduce")
    reduce_parser.add_argument("--raw-counts", type=Path, required=True)
    reduce_parser.add_argument("--output", type=Path, required=True)
    final = subparsers.add_parser("finalize")
    final.add_argument("--score-table", type=Path, required=True)
    final.add_argument("--exact-counts", type=Path, required=True)
    final.add_argument("--raw-gene-counts", type=Path, required=True)
    final.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.command == "score":
        score_variants(
            args.output_dir, intarna=args.intarna, threads=args.threads
        )
    elif args.command == "map":
        import sys

        map_energy_alignment_gene_keys(
            sys.stdin,
            output=sys.stdout,
            raw_counts_path=args.raw_counts,
            gene_index=GeneIndex.from_table(args.gene_table),
            query_metadata=load_energy_queries(args.score_table),
        )
    elif args.command == "reduce":
        import sys

        write_energy_counts(
            reduce_energy_gene_keys(sys.stdin, raw_counts_path=args.raw_counts),
            args.output,
        )
    else:
        finalize_energy(
            score_table=args.score_table,
            exact_counts=args.exact_counts,
            raw_gene_counts=args.raw_gene_counts,
            output_dir=args.output_dir,
        )


if __name__ == "__main__":
    main()
