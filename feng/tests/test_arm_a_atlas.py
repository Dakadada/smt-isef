"""Acceptance and invariant tests for the descriptive Arm A atlas."""

from pathlib import Path

import pandas as pd
import pytest

from feng.atlas.arm_a import (
    BLOCKS,
    CLAIM_BOUNDARY_TEXT,
    OUTPUT_NAMES,
    STUDIES,
    expected_s4_geometry,
    assign_block,
    lint_scientific_outputs,
    load_inputs,
    parse_segments,
    prepare_atlas,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
GUIDES = Path("/home/dashii/orcd/scratch/feng/phase0/splicerush_dcas13d_field.csv")
ASOS = REPO_ROOT / "feng" / "artifacts" / "phase0_arm_a_overlap" / "harmonized_aso_positions.csv"


def test_corrected_coordinate_examples_and_artificial_g_segments() -> None:
    assert expected_s4_geometry("I6", "-60..-46") == (260, 274, "260-274")
    assert expected_s4_geometry("I7", "+10..+29") == (384, 403, "384-403")
    assert expected_s4_geometry("E7", "+26..+40") == (345, 360, "345-355;357-360")
    assert parse_segments("345-355;357-360") == ((345, 355), (357, 360))
    assert assign_block(329, 340) == "E7-2"  # midpoint 334.5, right-hand base bin


@pytest.mark.skipif(not GUIDES.is_file(), reason="read-only source field is not staged")
def test_loader_rejects_pre_correction_frame(tmp_path: Path) -> None:
    frame = pd.read_csv(ASOS, dtype={"aso_id": str})
    frame["transformed_minigene_start"] -= 9
    frame["transformed_minigene_end"] -= 9
    shifted = tmp_path / "old_frame.csv"
    frame.to_csv(shifted, index=False)
    with pytest.raises(ValueError, match="nine nucleotides lower"):
        load_inputs(GUIDES, shifted)


@pytest.mark.skipif(not GUIDES.is_file(), reason="read-only source field is not staged")
def test_actual_inputs_satisfy_atlas_acceptance_invariants() -> None:
    data = prepare_atlas(GUIDES, ASOS)
    assert len(data.guides) == 395
    assert len(data.windows) == 125
    assert data.windows["distinct_window_key"].nunique() == 115
    assert (data.windows["n_guides_containing_hull"] > 0).all()
    assert len(data.blocks) == len(BLOCKS) * len(STUDIES) == 33
    assert set(data.blocks.groupby("block_id")["study"].nunique()) == {3}
    assert data.blocks["total_aso_records"].sum() == 125
    assert set(data.windows["result_state"]) == {"candidate"}
    assert set(data.blocks["result_state"]) == {"candidate"}


def test_package_linter_allows_only_fixed_claim_boundary(tmp_path: Path) -> None:
    for name in OUTPUT_NAMES:
        path = tmp_path / name
        if path.suffix == ".md":
            path.write_text("# Report\n\nBounded result.\n\n" + CLAIM_BOUNDARY_TEXT, encoding="utf-8")
        elif path.suffix in {".tsv", ".svg"}:
            path.write_text("bounded descriptive output", encoding="utf-8")
        else:
            path.write_bytes(b"vector output")
    assert lint_scientific_outputs(tmp_path) == []
    (tmp_path / "arm_a_block_summary.tsv").write_text("held-out result", encoding="utf-8")
    assert lint_scientific_outputs(tmp_path) == ["arm_a_block_summary.tsv: held-out"]
