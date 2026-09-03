"""Tests for the Arm E (StericRuler) length-titration report builder."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from feng.titration import report as R


# --- synthetic fixtures ------------------------------------------------------

def _base_row(length: int, **overrides) -> dict:
    """A window-shaped row (ladder + evaluate + compute_window keys)."""
    row = {
        "length": length,
        "chrom": "chr5",
        "start": 70076584 + (18 - length) // 2,
        "end": 70076601 - (18 - length) // 2,
        "strand": "+",
        "target_5to3": "A" * length,
        "aso_5to3": "T" * length,
        "state": "candidate",
        "intarna_dg": -1.0 * length,
        "intarna_no_duplex": False,
        "rnaduplex_dg": -1.1 * length,
        "rnacofold_dg": -1.1 * length,
        "rnaup_dg": -0.9 * length,
        "rnaplfold_accessibility": 0.4,
        "rnafold_mfe": 0.0,
        "off_target_hits": None,
        "potency": float(length),
        "potency_incomplete": False,
        "off_target_burden": None,
        "therapeutic_window": None,
        "window_incomplete": True,
    }
    row.update(overrides)
    return row


def _complete_rows() -> list[dict]:
    """A small window curve WITH off-target burden populated."""
    rows = []
    for length in (14, 16, 18, 20):
        burden = {14: 8, 16: 4, 18: 2, 20: 1}[length]
        pot = float(length)
        rows.append(
            _base_row(
                length,
                off_target_hits=burden,
                off_target_burden=burden,
                therapeutic_window=pot / (burden + 1),
                window_incomplete=False,
            )
        )
    return rows


@pytest.fixture
def synthetic_rows() -> list[dict]:
    return [_base_row(L) for L in (2, 8, 18, 24, 30)]


# --- artifact production -----------------------------------------------------

def test_build_writes_all_package_files(tmp_path, synthetic_rows):
    out = tmp_path / "arm_e_titration"
    R.build_arm_e_report(synthetic_rows, out_dir=out)
    for name in R.OUTPUT_NAMES:
        path = out / name
        assert path.exists(), f"missing artifact {name}"
        assert path.stat().st_size > 0, f"empty artifact {name}"


def test_figure_files_render_nonempty(tmp_path, synthetic_rows):
    out = tmp_path / "arm_e_titration"
    R.build_arm_e_report(synthetic_rows, out_dir=out)
    assert (out / R.FIG_SVG).stat().st_size > 0
    assert (out / R.FIG_PDF).stat().st_size > 0


def test_table_roundtrip_schema_and_state(tmp_path, synthetic_rows):
    out = tmp_path / "arm_e_titration"
    R.build_arm_e_report(synthetic_rows, out_dir=out)
    df = pd.read_parquet(out / R.TABLE_PARQUET)
    for col in R.TABLE_COLUMNS:
        assert col in df.columns, f"parquet missing column {col}"
    assert (df["state"] == "candidate").all()
    assert list(df["length"]) == sorted(r["length"] for r in synthetic_rows)
    # tsv reads back too
    tsv = pd.read_csv(out / R.TABLE_TSV, sep="\t")
    assert len(tsv) == len(synthetic_rows)


def test_incomplete_rows_render_and_are_flagged(tmp_path):
    # potency None (no usable on-target axis) + burden None
    rows = [
        _base_row(2, intarna_dg=None, intarna_no_duplex=None, rnaup_dg=None,
                  potency=None, potency_incomplete=True),
        _base_row(18),
    ]
    out = tmp_path / "arm_e_titration"
    R.build_arm_e_report(rows, out_dir=out)  # must not raise
    df = pd.read_parquet(out / R.TABLE_PARQUET)
    incomplete = df.loc[df["length"] == 2].iloc[0]
    assert bool(incomplete["potency_incomplete"]) is True
    assert bool(incomplete["window_incomplete"]) is True


def test_complete_curve_with_burden_plots(tmp_path):
    out = tmp_path / "arm_e_titration"
    R.build_arm_e_report(_complete_rows(), out_dir=out)
    df = pd.read_parquet(out / R.TABLE_PARQUET)
    assert df["therapeutic_window"].notna().all()
    assert (out / R.FIG_SVG).stat().st_size > 0


# --- honesty lint (load-bearing) ---------------------------------------------

def test_readme_has_claim_boundary_and_does_not_claim(tmp_path, synthetic_rows):
    out = tmp_path / "arm_e_titration"
    R.build_arm_e_report(synthetic_rows, out_dir=out)
    md = (out / R.REPORT_NAME).read_text()
    assert R.CLAIM_BOUNDARY_HEADING in md
    assert R.CLAIM_BOUNDARY_TEXT in md


def test_lint_passes_on_generated_report(tmp_path, synthetic_rows):
    out = tmp_path / "arm_e_titration"
    R.build_arm_e_report(synthetic_rows, out_dir=out)
    assert R.lint_scientific_outputs(out) == []


def test_lint_trips_on_forbidden_string(tmp_path, synthetic_rows):
    out = tmp_path / "arm_e_titration"
    R.build_arm_e_report(synthetic_rows, out_dir=out)
    # inject a forbidden phrase from feng.honesty.linter.FORBIDDEN_STRINGS
    md_path = out / R.REPORT_NAME
    md_path.write_text(md_path.read_text() + "\n\nThis result is generalizable.\n")
    violations = R.lint_scientific_outputs(out)
    assert any("generalizable" in v for v in violations)


def test_build_raises_if_report_would_contain_forbidden_language(tmp_path, synthetic_rows):
    # A forbidden hypothesis string must fail the build's internal lint.
    with pytest.raises(ValueError):
        R.build_arm_e_report(
            synthetic_rows,
            out_dir=tmp_path / "arm_e_titration",
            hypothesis="this ladder beats nusinersen",
        )


# --- provenance --------------------------------------------------------------

def test_provenance_records_draft_unsealed(tmp_path, synthetic_rows):
    out = tmp_path / "arm_e_titration"
    R.build_arm_e_report(synthetic_rows, out_dir=out)
    prov = json.loads((out / R.PROVENANCE).read_text())
    assert prov["prereg_locked"] is False
    assert "formula_version" in prov
    assert prov["state"] == "candidate"


def test_stale_files_rejected(tmp_path, synthetic_rows):
    out = tmp_path / "arm_e_titration"
    out.mkdir(parents=True)
    (out / "rogue_file.txt").write_text("not part of the package")
    with pytest.raises(ValueError):
        R.build_arm_e_report(synthetic_rows, out_dir=out)
