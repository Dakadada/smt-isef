"""Tests for the Arm E therapeutic-window aggregator + preregistration lock."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from feng.bench.scores import OLIGO_LENGTH
from feng.honesty.states import ResultState
from feng.titration.window import (
    AXIS_DEFINITIONS,
    OFF_TARGET_PSEUDOCOUNT,
    compute_window,
    directional_verdict,
    lock_prereg,
    prereg_hash,
    prereg_payload,
    verify_prereg,
)

_CANDIDATE = ResultState.CANDIDATE.value


def _row(length, **scores):
    """A minimal evaluate-style row with defaults for every score column."""
    base = {
        "length": length,
        "chrom": "chr5",
        "start": 70076584,
        "end": 70076601,
        "strand": "+",
        "target_5to3": "N" * length,
        "aso_5to3": "N" * length,
        "state": _CANDIDATE,
        "intarna_dg": None,
        "intarna_no_duplex": None,
        "rnaduplex_dg": None,
        "rnacofold_dg": None,
        "rnaup_dg": None,
        "rnaplfold_accessibility": None,
        "rnafold_mfe": None,
        "off_target_hits": None,
    }
    base.update(scores)
    return base


# --- potency: the three IntaRNA states ---------------------------------------


def test_potency_uses_real_intarna_and_rnaup():
    row = _row(18, intarna_dg=-27.0, intarna_no_duplex=False, rnaup_dg=-25.0)
    out = compute_window([row])[0]
    # potency = -mean(-27, -25) = +26
    assert out["potency"] == pytest.approx(26.0)
    assert out["potency_incomplete"] is False


def test_potency_ignores_intarna_when_no_duplex():
    # no_duplex 0.0 must NOT be credited as real binding; only rnaup counts.
    row = _row(4, intarna_dg=0.0, intarna_no_duplex=True, rnaup_dg=-3.0)
    out = compute_window([row])[0]
    assert out["potency"] == pytest.approx(3.0)  # -mean(-3) only
    assert out["potency_incomplete"] is False


def test_potency_ignores_intarna_when_none():
    # sub-seed IntaRNA error -> None; must not enter the mean.
    row = _row(2, intarna_dg=None, intarna_no_duplex=None, rnaup_dg=-1.0)
    out = compute_window([row])[0]
    assert out["potency"] == pytest.approx(1.0)


def test_potency_incomplete_when_no_backbone_axis():
    # no usable IntaRNA and no rnaup -> potency None, flagged incomplete.
    row = _row(2, intarna_dg=0.0, intarna_no_duplex=True, rnaup_dg=None)
    out = compute_window([row])[0]
    assert out["potency"] is None
    assert out["potency_incomplete"] is True


def test_no_duplex_mer_is_not_credited_spurious_binding():
    real = compute_window([_row(18, intarna_dg=-27.0, intarna_no_duplex=False)])[0]
    noduplex = compute_window([_row(4, intarna_dg=0.0, intarna_no_duplex=True, rnaup_dg=-2.0)])[0]
    assert real["potency"] > noduplex["potency"]


# --- off-target / window -----------------------------------------------------


def test_missing_off_target_flags_window_incomplete_not_zero_burden():
    row = _row(18, intarna_dg=-27.0, intarna_no_duplex=False, off_target_hits=None)
    out = compute_window([row])[0]
    assert out["off_target_burden"] is None
    assert out["therapeutic_window"] is None
    assert out["window_incomplete"] is True


def test_zero_off_target_burden_does_not_divide_by_zero():
    row = _row(18, intarna_dg=-27.0, intarna_no_duplex=False, off_target_hits=0)
    out = compute_window([row])[0]
    # window = potency / (0 + pseudocount) = potency / 1
    assert out["therapeutic_window"] == pytest.approx(27.0 / (0 + OFF_TARGET_PSEUDOCOUNT))
    assert out["window_incomplete"] is False


def test_window_penalises_off_target_burden():
    clean = compute_window([_row(18, intarna_dg=-27.0, intarna_no_duplex=False, off_target_hits=0)])[0]
    dirty = compute_window([_row(18, intarna_dg=-27.0, intarna_no_duplex=False, off_target_hits=50)])[0]
    assert clean["therapeutic_window"] > dirty["therapeutic_window"]


def test_window_incomplete_when_potency_missing_even_with_off_target():
    row = _row(2, intarna_dg=0.0, intarna_no_duplex=True, rnaup_dg=None, off_target_hits=3)
    out = compute_window([row])[0]
    assert out["therapeutic_window"] is None
    assert out["window_incomplete"] is True


def test_state_preserved_candidate():
    out = compute_window([_row(18, intarna_dg=-27.0, intarna_no_duplex=False)])[0]
    assert out["state"] == _CANDIDATE


# --- directional verdict -----------------------------------------------------


def test_directional_verdict_supports_when_peak_sub_18():
    curve = [(12, 5.0), (16, 9.0), (18, 8.0), (24, 6.0), (30, 4.0)]
    assert directional_verdict(curve) == "supports_shorter_better"


def test_directional_verdict_refutes_when_monotonic_increasing():
    curve = [(12, 1.0), (16, 2.0), (18, 3.0), (24, 4.0), (30, 5.0)]
    assert directional_verdict(curve) == "refutes_shorter_better"


def test_directional_verdict_indeterminate_peak_at_or_above_18_nonmonotonic():
    curve = [(12, 1.0), (18, 5.0), (24, 3.0), (30, 4.0)]
    assert directional_verdict(curve) == "indeterminate"


def test_directional_verdict_skips_none_values():
    curve = [(2, None), (16, 9.0), (18, 8.0), (30, None)]
    assert directional_verdict(curve) == "supports_shorter_better"


def test_directional_verdict_indeterminate_too_few_points():
    assert directional_verdict([(18, 5.0)]) == "indeterminate"
    assert directional_verdict([(2, None), (18, None)]) == "indeterminate"


def test_peak_boundary_at_18_is_not_shorter_better():
    # peak exactly at 18 (not sub-18) and monotonic up-to-18 then flat/down.
    curve = [(12, 1.0), (16, 2.0), (18, 3.0), (24, 3.0), (30, 3.0)]
    # monotonic non-decreasing -> refutes (longer never worse)
    assert directional_verdict(curve) == "refutes_shorter_better"


# --- MED-1: a sub-18 plateau top must NOT be read as shorter-better -----------


def test_sub18_plateau_is_not_supports_shorter_better():
    # The 16-mer only TIES the 18/24/30-mers (length-insensitive plateau).
    # A first-occurrence argmax at L16 must NOT manufacture a "supports" verdict.
    curve = [(12, 1.0), (16, 3.0), (18, 3.0), (24, 3.0), (30, 3.0)]
    assert directional_verdict(curve) != "supports_shorter_better"
    # non-decreasing -> refutes.
    assert directional_verdict(curve) == "refutes_shorter_better"


def test_strict_sub18_peak_still_supports():
    # A sub-18 length strictly above every longer rung -> genuine support.
    curve = [(12, 4.0), (16, 9.0), (18, 8.0), (24, 6.0), (30, 4.0)]
    assert directional_verdict(curve) == "supports_shorter_better"


def test_tie_at_global_max_spanning_sub18_and_over18_is_not_supports():
    # sub-18 peak (L12) ties a >=18 rung (L18) at the max -> not strict -> not supports.
    curve = [(12, 9.0), (16, 5.0), (18, 9.0), (24, 4.0)]
    assert directional_verdict(curve) != "supports_shorter_better"


def test_sub18_only_plateau_top_with_later_drop_is_not_supports():
    # Two sub-18 lengths tie at the top then it drops: non-monotonic, but the
    # peak only ties L16, so no STRICT sub-18 maximum -> indeterminate, not supports.
    curve = [(12, 9.0), (16, 9.0), (24, 3.0)]
    assert directional_verdict(curve) != "supports_shorter_better"


# --- LOW-1 / LOW-2: hardening -------------------------------------------------


def test_pseudocount_below_one_is_rejected():
    row = _row(18, intarna_dg=-27.0, intarna_no_duplex=False, off_target_hits=0)
    with pytest.raises(ValueError):
        compute_window([row], off_target_pseudocount=0)


def test_unfavourable_positive_energy_is_flagged_not_silent():
    # A usable axis with a POSITIVE (unfavourable) dg -> negative potency, flagged.
    row = _row(6, intarna_dg=2.0, intarna_no_duplex=False, rnaup_dg=1.0)
    out = compute_window([row])[0]
    assert out["potency"] < 0
    assert out["potency_unfavourable"] is True


def test_favourable_potency_not_flagged_unfavourable():
    out = compute_window([_row(18, intarna_dg=-27.0, intarna_no_duplex=False)])[0]
    assert out["potency_unfavourable"] is False


# --- preregistration lock ----------------------------------------------------


def _window_rows():
    rows = [
        _row(16, intarna_dg=-22.0, intarna_no_duplex=False, rnaup_dg=-22.0, off_target_hits=2),
        _row(18, intarna_dg=-27.0, intarna_no_duplex=False, rnaup_dg=-26.0, off_target_hits=5),
        _row(24, intarna_dg=-38.0, intarna_no_duplex=False, rnaup_dg=-38.0, off_target_hits=40),
    ]
    return compute_window(rows)


def test_prereg_payload_has_sorted_curve_and_axes():
    payload = prereg_payload(_window_rows(), hypothesis="shorter_is_better")
    lengths = [c["length"] for c in payload["curve"]]
    assert lengths == sorted(lengths)
    assert payload["hypothesis"] == "shorter_is_better"
    assert payload["axis_definitions"] == AXIS_DEFINITIONS
    assert payload["state"] == _CANDIDATE
    assert "directional_prediction" in payload


def test_prereg_hash_deterministic():
    p1 = prereg_payload(_window_rows(), hypothesis="shorter_is_better")
    p2 = prereg_payload(_window_rows(), hypothesis="shorter_is_better")
    assert prereg_hash(p1) == prereg_hash(p2)


def test_prereg_hash_changes_with_curve():
    base = prereg_payload(_window_rows(), hypothesis="shorter_is_better")
    perturbed_rows = _window_rows()
    perturbed_rows[0]["potency"] = perturbed_rows[0]["potency"] + 1.0
    perturbed = prereg_payload(perturbed_rows, hypothesis="shorter_is_better")
    assert prereg_hash(base) != prereg_hash(perturbed)


def test_prereg_hash_stable_across_dict_ordering():
    payload = prereg_payload(_window_rows(), hypothesis="shorter_is_better")
    reordered = dict(reversed(list(payload.items())))
    assert prereg_hash(payload) == prereg_hash(reordered)


def test_lock_prereg_defaults_to_draft(tmp_path):
    payload = prereg_payload(_window_rows(), hypothesis="shorter_is_better")
    path = tmp_path / "prereg.json"
    record = lock_prereg(payload, path)
    assert record["locked"] is False  # DRAFT by default; real seal is gated
    on_disk = json.loads(Path(path).read_text())
    assert on_disk["locked"] is False
    assert on_disk["hash"] == prereg_hash(payload)


def test_verify_prereg_round_trips(tmp_path):
    payload = prereg_payload(_window_rows(), hypothesis="shorter_is_better")
    path = tmp_path / "prereg.json"
    lock_prereg(payload, path)
    record = verify_prereg(path)
    assert record["hash"] == prereg_hash(payload)


def test_verify_prereg_detects_tamper(tmp_path):
    payload = prereg_payload(_window_rows(), hypothesis="shorter_is_better")
    path = tmp_path / "prereg.json"
    lock_prereg(payload, path)
    record = json.loads(path.read_text())
    record["payload"]["curve"][0]["potency"] = 999.0  # tamper, hash unchanged
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError):
        verify_prereg(path)


def test_prereg_prediction_falls_back_to_potency_without_off_target():
    # no off-target -> window incomplete for all; prediction uses potency axis.
    rows = compute_window([
        _row(16, intarna_dg=-22.0, intarna_no_duplex=False, rnaup_dg=-22.0),
        _row(18, intarna_dg=-27.0, intarna_no_duplex=False, rnaup_dg=-26.0),
        _row(24, intarna_dg=-38.0, intarna_no_duplex=False, rnaup_dg=-38.0),
    ])
    payload = prereg_payload(rows, hypothesis="shorter_is_better")
    assert payload["directional_prediction_axis"] == "potency"


def test_oligo_length_constant_is_18():
    # guards the sub-18 boundary the verdict logic depends on.
    assert OLIGO_LENGTH == 18
