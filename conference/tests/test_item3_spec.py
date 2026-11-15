import csv
import io
from math import isclose
from pathlib import Path

import pytest

from conference.aso_pipeline.item3_spec import (
    GAS_CONSTANT_KCAL_PER_MOL_K,
    TEMPERATURE_K,
    absolute_threshold_rows,
    fractional_occupancy,
    map_exact_loci,
    occupancy_rows,
)
from conference.aso_pipeline.offtarget import GeneIndex


def test_fractional_occupancy_uses_molar_concentration_and_kd_from_dg():
    # dG = RT ln(1e-9) makes Kd exactly 1 nM, hence theta=1/2 at 1 nM.
    dg = GAS_CONSTANT_KCAL_PER_MOL_K * TEMPERATURE_K * __import__("math").log(1e-9)
    assert isclose(fractional_occupancy(dg, 1e-9), 0.5, abs_tol=1e-12)
    with pytest.raises(ValueError, match="positive"):
        fractional_occupancy(-10.0, 0.0)


def test_mapper_preserves_locus_multiplicity_and_marks_only_smn2_exact_hit_on_target(
    tmp_path: Path,
):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        'chr5\tt\tgene\t1\t100\t.\t+\t.\tgene_id "ENSG00000205571.17";\n'
        'chr5\tt\tgene\t201\t300\t.\t+\t.\tgene_id "OTHER";\n'
    )
    scores = {
        "q0": {"aso_id": "F10", "dg": -12.0, "is_exact": True, "length": 10, "target": "AAAAAAAAAA", "viable_100x": True},
        "q1": {"aso_id": "F10", "dg": -10.0, "is_exact": False, "length": 10, "target": "CAAAAAAAAA", "viable_100x": True},
    }
    output = io.StringIO()

    map_exact_loci(
        [
            "q0\t-\tchr5\t10\n",
            "q0\t-\tchr5\t210\n",
            "q1\t-\tchr5\t220\n",
            "q1\t-\tchr5\t220\n",
        ],
        output=output,
        gene_index=GeneIndex.from_gtf(gtf),
        query_scores=scores,
    )

    rows = list(csv.DictReader(io.StringIO(output.getvalue()), delimiter="\t"))
    assert len(rows) == 4
    assert [row["site_class"] for row in rows] == ["on_target", "off_target", "off_target", "off_target"]
    assert rows[2]["compatible_gene_ids"] == "OTHER"


def test_occupancy_rows_keep_raw_and_derived_values_and_cross_aso_ratios():
    loci = [
        {"aso_id": "F18", "site_class": "off_target", "dg": -9.0},
        {"aso_id": "F18", "site_class": "off_target", "dg": -9.0},
        {"aso_id": "F14", "site_class": "off_target", "dg": -9.0},
        {"aso_id": "F10", "site_class": "off_target", "dg": -9.0},
        {"aso_id": "F10", "site_class": "off_target", "dg": -9.0},
        {"aso_id": "F10", "site_class": "off_target", "dg": -9.0},
        {"aso_id": "3UP8", "site_class": "saturated", "dg": None},
    ]
    on_target = {"F18": -12.0, "F14": -11.0, "F10": -10.0}

    rows = occupancy_rows(loci, on_target, concentrations_molar=(1e-9,))

    f18 = next(row for row in rows if row["aso_id"] == "F18")
    f14 = next(row for row in rows if row["aso_id"] == "F14")
    f10 = next(row for row in rows if row["aso_id"] == "F10")
    saturated = next(row for row in rows if row["aso_id"] == "3UP8")
    assert f18["state"] == "[candidate]"
    assert f18["raw_off_target_locus_count"] == 2
    assert f10["sum_theta_off"] == pytest.approx(3 * f18["sum_theta_off"] / 2)
    assert f10["f10_over_f18_burden_ratio"] == pytest.approx(1.5)
    assert f14["f14_over_f18_burden_ratio"] == pytest.approx(0.5)
    assert f18["theta_on_over_sum_theta_off"] > 0
    assert f18["theta_on_over_theta_total"] > 0
    assert saturated["availability_state"] == "saturated"
    assert saturated["sum_theta_off"] is None


def test_occupancy_ignores_absolute_only_replay_loci():
    loci = [
        {"aso_id": "F10", "site_class": "off_target", "dg": -10.0, "occupancy_eligible_100x": True},
        {"aso_id": "F10", "site_class": "off_target", "dg": -8.0, "occupancy_eligible_100x": False},
    ]

    row = occupancy_rows(loci, {"F10": -12.0}, concentrations_molar=(1e-9,))[0]

    assert row["raw_off_target_locus_count"] == 1
    assert row["sum_theta_off"] == pytest.approx(fractional_occupancy(-10.0, 1e-9))


def test_accessibility_gate_excludes_failed_and_below_cutoff_without_inventing_expression():
    loci = [
        {"aso_id": "F10", "site_class": "off_target", "dg": -9.0, "compatible_gene_ids": "G1", "accessibility": 0.7},
        {"aso_id": "F10", "site_class": "off_target", "dg": -9.0, "compatible_gene_ids": "G2", "accessibility": 0.4},
        {"aso_id": "F10", "site_class": "off_target", "dg": -9.0, "compatible_gene_ids": "G3", "accessibility": None},
    ]

    row = next(
        row for row in occupancy_rows(
            loci,
            {"F10": -10.0},
            concentrations_molar=(1e-9,),
            accessibility_cutoff=0.5,
        ) if row["aso_id"] == "F10"
    )

    assert row["expression_state"] == "unavailable"
    assert row["accessibility_evaluated_locus_count"] == 2
    assert row["accessibility_failed_locus_count"] == 1
    assert row["accessible_off_target_locus_count"] == 1
    assert row["accessible_distinct_gene_count"] == 1


def test_absolute_threshold_counts_distinct_compatible_genes_and_raw_loci():
    loci = [
        {"aso_id": "F10", "site_class": "off_target", "dg": -12.0, "compatible_gene_ids": "G1;G2"},
        {"aso_id": "F10", "site_class": "off_target", "dg": -11.0, "compatible_gene_ids": "G2"},
        {"aso_id": "F10", "site_class": "off_target", "dg": -9.0, "compatible_gene_ids": ""},
    ]

    rows = absolute_threshold_rows(loci, thresholds=(-10.0, -12.0))

    minus10 = next(row for row in rows if row["absolute_dg_cutoff_kcal_mol"] == -10.0)
    minus12 = next(row for row in rows if row["absolute_dg_cutoff_kcal_mol"] == -12.0)
    assert minus10["raw_off_target_locus_count"] == 2
    assert minus10["distinct_transcript_compatible_gene_count"] == 2
    assert minus12["raw_off_target_locus_count"] == 1
    assert minus12["distinct_transcript_compatible_gene_count"] == 2
