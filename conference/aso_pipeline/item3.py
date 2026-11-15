"""Item 3 preflight, streamed counting CLI, and Item 6 alignment summaries."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
from pathlib import Path

from . import CANDIDATE
from .item1 import _write_tsv
from .item6 import summarize_null
from .offtarget import (
    GeneIndex,
    expected_random_hits,
    map_alignment_gene_keys,
    reduce_sorted_gene_keys,
)
from .sequences import to_rna


GRCH38_PRIMARY_ASSEMBLY_BASES = 3_088_269_832
MISMATCHES = 3
SATURATION_EXPECTED_HITS = 10_000_000.0
MANDATORY_FAMILIES = {"F18", "F10"}


@dataclass(frozen=True, slots=True)
class QueryDecision:
    aso_id: str
    length_nt: int
    expected_hits: float
    status: str
    reason: str
    run_alignment: bool
    gene_count: None = None


def classify_query(
    aso_id: str,
    sequence: str,
    *,
    mandatory_family: bool,
    saturation_threshold: float = SATURATION_EXPECTED_HITS,
) -> QueryDecision:
    normalized = to_rna(sequence)
    estimate = expected_random_hits(
        len(normalized), MISMATCHES, GRCH38_PRIMARY_ASSEMBLY_BASES
    )
    if estimate >= saturation_threshold and not mandatory_family:
        return QueryDecision(
            aso_id,
            len(normalized),
            estimate,
            "saturated",
            "saturated / non-specific by construction - count uninformative",
            False,
        )
    if estimate >= saturation_threshold:
        return QueryDecision(
            aso_id,
            len(normalized),
            estimate,
            "scheduled_mandatory",
            "mandatory F18/F10 family override for headline and null comparison",
            True,
        )
    return QueryDecision(
        aso_id,
        len(normalized),
        estimate,
        "scheduled",
        "below predeclared expected-hit saturation threshold",
        True,
    )


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


def _identity_rows(mappings: list[dict]) -> list[dict]:
    return [
        {
            "state": CANDIDATE,
            "unique_query_id": row["unique_query_id"],
            "query_id": row["query_id"],
            "aso_id": row["aso_id"],
            "is_real": str(_as_bool(row["is_real"])).lower(),
            "scramble_index": row["draw_id"],
            "sequence_rna_5to3": row["sequence_rna_5to3"],
            "preflight_status": row["preflight_status"],
        }
        for row in mappings
    ]


def prepare_manifest(query_table: Path, output_dir: Path) -> list[dict]:
    source = _read_tsv(query_table)
    rows = []
    unique: dict[str, dict] = {}
    for row in source:
        sequence = to_rna(row["sequence_rna_5to3"])
        unique_id = "Q_" + hashlib.sha256(sequence.encode()).hexdigest()[:16]
        decision = classify_query(
            row["aso_id"],
            sequence,
            mandatory_family=row["aso_id"] in MANDATORY_FAMILIES,
        )
        mapped = {
            **row,
            "unique_query_id": unique_id,
            "length_nt": len(sequence),
            "expected_random_both_strand_hits_v3": decision.expected_hits,
            "preflight_status": decision.status,
            "preflight_reason": decision.reason,
            "run_alignment": decision.run_alignment,
        }
        rows.append(mapped)
        prior = unique.get(sequence)
        if prior is None:
            unique[sequence] = mapped
        elif prior["run_alignment"] != mapped["run_alignment"]:
            raise ValueError("identical sequence received conflicting preflight decisions")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_tsv(output_dir / "item3_query_manifest.tsv", rows)
    _write_tsv(output_dir / "item3_query_identity_map.tsv", _identity_rows(rows))
    def write_fasta(path: Path, entries: list[tuple[str, dict]]) -> None:
        temporary = path.with_suffix(path.suffix + ".tmp")
        with temporary.open("w") as handle:
            for sequence, entry in sorted(
                entries, key=lambda pair: pair[1]["unique_query_id"]
            ):
                handle.write(
                    f">{entry['unique_query_id']}\n{sequence.replace('U', 'T')}\n"
                )
        temporary.replace(path)

    eligible = [pair for pair in unique.items() if pair[1]["run_alignment"]]
    write_fasta(output_dir / "item3_unique_queries.fa", eligible)
    for family in ("F18", "F14", "F10"):
        family_entries = [pair for pair in eligible if pair[1]["aso_id"] == family]
        write_fasta(output_dir / f"item3_unique_queries_{family}.fa", family_entries)

    f10_dir = output_dir / "item3_unique_queries_F10"
    f10_dir.mkdir(exist_ok=True)
    for stale in f10_dir.glob("*.fa"):
        stale.unlink()
    for sequence, entry in [pair for pair in eligible if pair[1]["aso_id"] == "F10"]:
        write_fasta(f10_dir / f"{entry['unique_query_id']}.fa", [(sequence, entry)])
    return rows


def _read_fasta_lengths(path: Path) -> dict[str, int]:
    lengths = {}
    current = None
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                lengths[current] = 0
            elif line:
                if current is None:
                    raise ValueError("FASTA sequence before header")
                lengths[current] += len(line)
    return lengths


def map_alignments(
    *,
    aligner: str,
    gene_table: Path,
    query_fasta: Path,
    raw_counts: Path,
) -> None:
    import sys

    map_alignment_gene_keys(
        sys.stdin,
        output=sys.stdout,
        raw_counts_path=raw_counts,
        aligner=aligner,
        gene_index=GeneIndex.from_table(gene_table),
        query_lengths=_read_fasta_lengths(query_fasta),
    )


def reduce_gene_keys(
    *,
    aligner: str,
    query_fasta: Path,
    raw_counts: Path,
    output: Path,
) -> None:
    import sys

    lengths = _read_fasta_lengths(query_fasta)
    counts = reduce_sorted_gene_keys(
        sys.stdin, raw_counts_path=raw_counts, query_lengths=lengths
    )
    rows = [
        {"state": CANDIDATE, "unique_query_id": query_id, "aligner": aligner, **value}
        for query_id, value in sorted(counts.items())
    ]
    _write_tsv(output, rows)


def merge_count_tables(inputs: list[Path], output: Path) -> None:
    rows = []
    seen = set()
    for path in inputs:
        for row in _read_tsv(path):
            key = (row["unique_query_id"], row["aligner"])
            if key in seen:
                raise ValueError(f"duplicate count row while merging: {key}")
            seen.add(key)
            rows.append(row)
    rows.sort(key=lambda row: (row["aligner"], row["unique_query_id"]))
    _write_tsv(output, rows)


def finalize(
    *, manifest: Path, bowtie_counts: Path, blast_counts: Path, output_dir: Path
) -> None:
    mappings = _read_tsv(manifest)
    output_dir = Path(output_dir)
    _write_tsv(output_dir / "item3_query_identity_map.tsv", _identity_rows(mappings))
    count_lookup = {}
    for path in (bowtie_counts, blast_counts):
        for row in _read_tsv(path):
            count_lookup[(row["unique_query_id"], row["aligner"])] = row

    expanded = []
    for mapping in mappings:
        for aligner in ("bowtie1-v3", "blastn-short"):
            key = (mapping["unique_query_id"], aligner)
            count = count_lookup.get(key)
            saturated = not _as_bool(mapping["run_alignment"])
            if count is None and not saturated:
                raise ValueError(f"missing nonsaturated count: {key}")
            expanded.append(
                {
                    "state": CANDIDATE,
                    "query_id": mapping["query_id"],
                    "aso_id": mapping["aso_id"],
                    "draw_id": mapping["draw_id"],
                    "is_real": str(_as_bool(mapping["is_real"])).lower(),
                    "sequence_rna_5to3": mapping["sequence_rna_5to3"],
                    "aligner": aligner,
                    "preflight_status": mapping["preflight_status"],
                    "preflight_reason": mapping["preflight_reason"],
                    "expected_random_both_strand_hits_v3": mapping[
                        "expected_random_both_strand_hits_v3"
                    ],
                    "raw_alignment_count": "" if count is None else count["raw_alignment_count"],
                    "distinct_genes_any_strand": (
                        "" if count is None else count["distinct_genes_any_strand"]
                    ),
                    "distinct_genes_transcript_compatible": (
                        ""
                        if count is None
                        else count["distinct_genes_transcript_compatible"]
                    ),
                }
            )

    real_rows = [row for row in expanded if row["is_real"] == "true"]
    null_rows = [row for row in expanded if row["is_real"] == "false"]
    _write_tsv(output_dir / "item3_gene_counts.tsv", real_rows)
    _write_tsv(output_dir / "item6_alignment_null_raw.tsv", null_rows)

    real_wide_rows = []
    for aso_id in ("F18", "F14", "F10"):
        by_aligner = {
            row["aligner"]: row
            for row in real_rows
            if row["aso_id"] == aso_id
        }
        bowtie = by_aligner["bowtie1-v3"]
        blast = by_aligner["blastn-short"]
        real_wide_rows.append(
            {
                "state": CANDIDATE,
                "aso_id": aso_id,
                "is_real": "true",
                "sequence_rna_5to3": bowtie["sequence_rna_5to3"],
                "bowtie1_v3_raw_alignments": bowtie["raw_alignment_count"],
                "bowtie1_v3_distinct_genes_any_strand": bowtie[
                    "distinct_genes_any_strand"
                ],
                "bowtie1_v3_distinct_genes": bowtie[
                    "distinct_genes_transcript_compatible"
                ],
                "blastn_short_raw_alignments": blast["raw_alignment_count"],
                "blastn_short_distinct_genes_any_strand": blast[
                    "distinct_genes_any_strand"
                ],
                "blastn_short_distinct_genes": blast[
                    "distinct_genes_transcript_compatible"
                ],
            }
        )
    _write_tsv(output_dir / "item3_real_aso_gene_counts.tsv", real_wide_rows)

    real_lookup = {row["aso_id"]: row for row in real_wide_rows}
    _write_tsv(
        output_dir / "item3_headline_f18_vs_f10.tsv",
        [
            {
                "state": CANDIDATE,
                "metric": "distinct_genes_transcript_compatible",
                "f18_bowtie1_v3": real_lookup["F18"]["bowtie1_v3_distinct_genes"],
                "f10_bowtie1_v3": real_lookup["F10"]["bowtie1_v3_distinct_genes"],
                "f18_blastn_short": real_lookup["F18"]["blastn_short_distinct_genes"],
                "f10_blastn_short": real_lookup["F10"]["blastn_short_distinct_genes"],
            }
        ],
    )

    summaries = []
    for aso_id in ("F18", "F14", "F10", "3UP8"):
        for aligner in ("bowtie1-v3", "blastn-short"):
            real = next(
                row for row in real_rows if row["aso_id"] == aso_id and row["aligner"] == aligner
            )
            null = [
                row
                for row in null_rows
                if row["aso_id"] == aso_id and row["aligner"] == aligner
            ]
            values = [
                int(row["distinct_genes_transcript_compatible"])
                for row in null
                if row["distinct_genes_transcript_compatible"] != ""
            ]
            base = {
                "state": CANDIDATE,
                "aso_id": aso_id,
                "aligner": aligner,
                "real_gene_count": real["distinct_genes_transcript_compatible"],
                "null_n_requested": len(null),
                "null_n_numeric": len(values),
                "status": real["preflight_status"],
                "reason": real["preflight_reason"],
            }
            if real["distinct_genes_transcript_compatible"] != "" and values:
                metrics = summarize_null(
                    int(real["distinct_genes_transcript_compatible"]),
                    values,
                    lower_is_stronger=False,
                )
                metrics.pop("null_n_numeric")
                base.update(metrics)
            summaries.append(base)
    _write_tsv(output_dir / "item6_alignment_null_summary.tsv", summaries)

    clean_null_summaries = []
    for aso_id in ("F18", "F14", "F10"):
        for aligner in ("bowtie1-v3", "blastn-short"):
            values = [
                int(row["distinct_genes_transcript_compatible"])
                for row in null_rows
                if row["aso_id"] == aso_id
                and row["aligner"] == aligner
                and row["distinct_genes_transcript_compatible"] != ""
            ]
            real_value = next(
                int(row["distinct_genes_transcript_compatible"])
                for row in real_rows
                if row["aso_id"] == aso_id and row["aligner"] == aligner
            )
            metrics = summarize_null(real_value, values, lower_is_stronger=False)
            clean_null_summaries.append(
                {
                    "state": CANDIDATE,
                    "aso_id": aso_id,
                    "aligner": aligner,
                    "scramble_n": metrics["null_n_numeric"],
                    "scramble_mean": metrics["null_mean"],
                    "scramble_sd_population": metrics["null_sd_population"],
                    "scramble_q05": metrics["null_q05"],
                    "scramble_median": metrics["null_median"],
                    "scramble_q95": metrics["null_q95"],
                }
            )
    _write_tsv(
        output_dir / "item6_alignment_scramble_summary.tsv", clean_null_summaries
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)
    prepare = subparsers.add_parser("prepare")
    prepare.add_argument("--queries", type=Path, required=True)
    prepare.add_argument("--output-dir", type=Path, required=True)
    genes = subparsers.add_parser("genes")
    genes.add_argument("--gtf", type=Path, required=True)
    genes.add_argument("--output", type=Path, required=True)
    map_parser = subparsers.add_parser("map")
    map_parser.add_argument(
        "--aligner", choices=("bowtie1-v3", "blastn-short"), required=True
    )
    map_parser.add_argument("--gene-table", type=Path, required=True)
    map_parser.add_argument("--query-fasta", type=Path, required=True)
    map_parser.add_argument("--raw-counts", type=Path, required=True)
    reduce_parser = subparsers.add_parser("reduce")
    reduce_parser.add_argument(
        "--aligner", choices=("bowtie1-v3", "blastn-short"), required=True
    )
    reduce_parser.add_argument("--query-fasta", type=Path, required=True)
    reduce_parser.add_argument("--raw-counts", type=Path, required=True)
    reduce_parser.add_argument("--output", type=Path, required=True)
    merge = subparsers.add_parser("merge")
    merge.add_argument("--inputs", type=Path, nargs="+", required=True)
    merge.add_argument("--output", type=Path, required=True)
    final = subparsers.add_parser("finalize")
    final.add_argument("--manifest", type=Path, required=True)
    final.add_argument("--bowtie-counts", type=Path, required=True)
    final.add_argument("--blast-counts", type=Path, required=True)
    final.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.command == "prepare":
        prepare_manifest(args.queries, args.output_dir)
    elif args.command == "genes":
        GeneIndex.from_gtf(args.gtf).write_table(args.output)
    elif args.command == "map":
        map_alignments(
            aligner=args.aligner,
            gene_table=args.gene_table,
            query_fasta=args.query_fasta,
            raw_counts=args.raw_counts,
        )
    elif args.command == "reduce":
        reduce_gene_keys(
            aligner=args.aligner,
            query_fasta=args.query_fasta,
            raw_counts=args.raw_counts,
            output=args.output,
        )
    elif args.command == "merge":
        merge_count_tables(args.inputs, args.output)
    else:
        finalize(
            manifest=args.manifest,
            bowtie_counts=args.bowtie_counts,
            blast_counts=args.blast_counts,
            output_dir=args.output_dir,
        )


if __name__ == "__main__":
    main()
