from math import isclose
import csv
import io
from pathlib import Path

import pytest

from conference.aso_pipeline.item3_energy import (
    TargetVariant,
    EnergyQuery,
    build_score_rows,
    enumerate_target_variants,
    finalize_energy,
    fold_margin_kcal_mol,
    load_energy_queries,
    map_energy_alignment_gene_keys,
    reduce_energy_gene_keys,
    score_variants,
    viable_at_fold,
    write_energy_counts,
)
from conference.aso_pipeline.offtarget import GeneIndex


def _write_tsv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_enumeration_is_the_complete_unique_hamming_ball_through_three_mismatches():
    variants = list(enumerate_target_variants("ACGUACGUAC", max_mismatches=3))

    assert len(variants) == 3676
    assert len({variant.sequence_rna_5to3 for variant in variants}) == 3676
    assert sum(variant.mismatch_count == 0 for variant in variants) == 1
    assert sum(variant.mismatch_count == 1 for variant in variants) == 30
    assert sum(variant.mismatch_count == 2 for variant in variants) == 405
    assert sum(variant.mismatch_count == 3 for variant in variants) == 3240


def test_fold_thresholds_use_rt_log_ratio_at_37_celsius():
    assert isclose(fold_margin_kcal_mol(1), 0.0, abs_tol=1e-12)
    assert isclose(fold_margin_kcal_mol(10), 1.419, abs_tol=0.001)
    assert isclose(fold_margin_kcal_mol(100), 2.838, abs_tol=0.001)


def test_viability_keeps_more_negative_energy_and_rejects_no_duplex():
    assert viable_at_fold(-8.581, False, on_target_dg=-10.0, max_fold=10)
    assert not viable_at_fold(-8.580, False, on_target_dg=-10.0, max_fold=10)
    assert not viable_at_fold(0.0, True, on_target_dg=-10.0, max_fold=100)


def test_score_rows_reference_the_exact_target_and_orient_bowtie_queries():
    variants = [
        TargetVariant("CCAGCAUUAU", 0),
        TargetVariant("CCAGCAUUAA", 1),
    ]
    scores = {
        "CCAGCAUUAU": (-10.0, False),
        "CCAGCAUUAA": (-8.6, False),
    }

    rows = build_score_rows(
        "F10", "AUAAUGCUGG", "CCAGCAUUAU", variants, scores
    )

    exact, mismatch = rows
    assert exact["state"] == "[candidate]"
    assert exact["bowtie_query_id"] == exact["variant_id"] + " [candidate]"
    assert exact["intarna_dg_kcal_mol"] == exact["on_target_dg_kcal_mol"] == -10.0
    assert exact["bowtie_query_rna_5to3"] == "AUAAUGCUGG"
    assert mismatch["bowtie_query_rna_5to3"] == "UUAAUGCUGG"
    assert mismatch["delta_dg_vs_on_target_kcal_mol"] == pytest.approx(1.4)
    assert mismatch["viable_3x"] is False
    assert mismatch["viable_10x"] is True


def test_exact_target_no_duplex_fails_closed():
    with pytest.raises(ValueError, match="on-target.*no favourable duplex"):
        build_score_rows(
            "F10",
            "AUAAUGCUGG",
            "CCAGCAUUAU",
            [TargetVariant("CCAGCAUUAU", 0)],
            {"CCAGCAUUAU": (0.0, True)},
        )


def test_score_artifacts_report_viable_counts_and_align_only_loosest_set(
    tmp_path: Path,
):
    exact_target = "CCAGCAUUAU"

    def controlled_scorer(aso, targets, *, intarna, threads):
        return {
            target_id: (-10.0 if sequence == exact_target else -8.0, False)
            for target_id, sequence in targets.items()
        }

    rows = score_variants(
        tmp_path,
        intarna=Path("unused-in-controlled-test"),
        threads=2,
        aso_ids=("F10",),
        score_many=controlled_scorer,
    )

    assert len(rows) == 3676
    assert sum(row["viable_10x"] for row in rows) == 1
    assert sum(row["viable_100x"] for row in rows) == 3676
    with (tmp_path / "item3_energy_variant_scores.tsv").open() as handle:
        written = list(csv.DictReader(handle, delimiter="\t"))
    assert len(written) == 3676
    assert {row["state"] for row in written} == {"[candidate]"}
    fasta_headers = [
        line
        for line in (tmp_path / "item3_energy_viable_100x.fa").read_text().splitlines()
        if line.startswith(">")
    ]
    assert len(fasta_headers) == 3676
    metadata = load_energy_queries(tmp_path / "item3_energy_variant_scores.tsv")
    assert len(metadata) == 3676
    exact_id = next(
        row["variant_id"]
        for row in rows
        if row["variant_target_rna_5to3"] == exact_target
    )
    assert metadata[exact_id + " [candidate]"].viable_folds == (1, 3, 10, 30, 100)
    mismatch_id = next(
        row["variant_id"]
        for row in rows
        if row["variant_target_rna_5to3"] != exact_target
    )
    assert metadata[mismatch_id + " [candidate]"].viable_folds == (30, 100)


def test_energy_gene_counter_unions_genes_across_variants_and_keeps_thresholds(
    tmp_path: Path,
):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id "PLUS";\n'
    )
    metadata = {
        "q_exact": EnergyQuery("F10", 10, (1, 3, 10, 30, 100)),
        "q_mismatch": EnergyQuery("F10", 10, (100,)),
    }
    raw_counts = tmp_path / "raw.tsv"
    keys = io.StringIO()

    map_energy_alignment_gene_keys(
        [
            "q_exact\t-\tchr1\t10\n",
            "q_exact\t-\tchr1\t10\n",
            "q_mismatch\t-\tchr1\t20\n",
        ],
        output=keys,
        raw_counts_path=raw_counts,
        gene_index=GeneIndex.from_gtf(gtf),
        query_metadata=metadata,
    )
    counts = reduce_energy_gene_keys(
        sorted(set(keys.getvalue().splitlines())), raw_counts_path=raw_counts
    )

    assert counts[("F10", 10)] == {
        "exact_alignment_count": 2,
        "distinct_genes_any_strand": 1,
        "distinct_genes_transcript_compatible": 1,
    }
    assert counts[("F10", 100)] == {
        "exact_alignment_count": 3,
        "distinct_genes_any_strand": 1,
        "distinct_genes_transcript_compatible": 1,
    }
    output = tmp_path / "counts.tsv"
    write_energy_counts(counts, output)
    with output.open() as handle:
        written = list(csv.DictReader(handle, delimiter="\t"))
    assert written[0]["state"] == "[candidate]"
    assert written[0]["aso_id"] == "F10"


def test_energy_reducer_waits_for_key_stream_before_mapper_sidecar(tmp_path: Path):
    raw_counts = tmp_path / "raw.tsv"

    def keys_then_sidecar():
        raw_counts.write_text(
            "aso_id\tthreshold_fold\texact_alignment_count\nF10\t10\t2\n"
        )
        yield "F10\t10\tC\tGENE1\n"

    counts = reduce_energy_gene_keys(
        keys_then_sidecar(), raw_counts_path=raw_counts
    )

    assert counts[("F10", 10)]["exact_alignment_count"] == 2
    assert counts[("F10", 10)]["distinct_genes_transcript_compatible"] == 1


def test_finalize_keeps_raw_counts_separate_and_emits_predeclared_headline(
    tmp_path: Path,
):
    score_rows = []
    for aso_id, on_target_dg in (("F18", -20.0), ("F14", -15.0), ("F10", -10.0)):
        for index in range(2):
            row = {
                "state": "[candidate]",
                "aso_id": aso_id,
                "on_target_dg_kcal_mol": on_target_dg,
            }
            for fold in (1, 3, 10, 30, 100):
                row[f"viable_{fold}x"] = index == 0 or aso_id != "F10"
            score_rows.append(row)
    scores = tmp_path / "scores.tsv"
    _write_tsv(scores, score_rows)

    exact_rows = []
    for aso_id, primary_genes in (("F18", 9), ("F14", 5), ("F10", 2)):
        for fold in (1, 3, 10, 30, 100):
            exact_rows.append(
                {
                    "state": "[candidate]",
                    "aso_id": aso_id,
                    "threshold_fold": fold,
                    "exact_alignment_count": primary_genes * 2,
                    "distinct_genes_any_strand": primary_genes + 1,
                    "distinct_genes_transcript_compatible": primary_genes,
                }
            )
    exact = tmp_path / "exact.tsv"
    _write_tsv(exact, exact_rows)
    raw = tmp_path / "raw.tsv"
    _write_tsv(
        raw,
        [
            {"aso_id": "F18", "bowtie1_v3_distinct_genes": 1518},
            {"aso_id": "F14", "bowtie1_v3_distinct_genes": 37713},
            {"aso_id": "F10", "bowtie1_v3_distinct_genes": 71148},
        ],
    )

    finalize_energy(
        score_table=scores,
        exact_counts=exact,
        raw_gene_counts=raw,
        output_dir=tmp_path,
    )

    with (tmp_path / "item3_energy_filtered_sweep.tsv").open() as handle:
        sweep = list(csv.DictReader(handle, delimiter="\t"))
    assert len(sweep) == 15
    assert {row["state"] for row in sweep} == {"[candidate]"}
    f10_primary = next(
        row
        for row in sweep
        if row["aso_id"] == "F10" and row["threshold_fold"] == "10"
    )
    assert f10_primary["viable_variant_count"] == "1"
    assert f10_primary["raw_bowtie1_v3_distinct_genes"] == "71148"
    assert f10_primary["energy_filtered_distinct_genes"] == "2"

    with (tmp_path / "item3_energy_filtered_headline_f18_vs_f10.tsv").open() as handle:
        headline = list(csv.DictReader(handle, delimiter="\t"))
    assert headline[0]["primary_threshold_fold"] == "10"
    assert headline[0]["f18_energy_filtered_distinct_genes"] == "9"
    assert headline[0]["f10_energy_filtered_distinct_genes"] == "2"
    assert headline[0]["prediction_result"] == "confirmed"
