"""Tests for the Arm E therapeutic-window aggregator + preregistration lock."""

from __future__ import annotations

import hashlib
import json
import math
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

from feng.bench.scores import OLIGO_LENGTH
from feng.honesty.states import ResultState
from feng.titration import arm_e_hoff, arm_e_scan
from feng.titration.window import (
    AXIS_DEFINITIONS,
    EXPECTED_LADDER_LENGTHS,
    OFF_TARGET_PSEUDOCOUNT,
    compute_window,
    assemble_off_target_windows,
    directional_verdict,
    lock_prereg,
    prereg_hash,
    prereg_payload,
    verify_prereg,
)

_CANDIDATE = ResultState.CANDIDATE.value

_BOOLEAN_ENCODERS = (
    pytest.param(lambda value: bool(value), id="python-bool"),
    pytest.param(lambda value: np.bool_(value), id="numpy-bool-scalar"),
    pytest.param(lambda value: np.array(value, dtype=bool), id="numpy-bool-0d"),
    pytest.param(
        lambda value: np.array([value, not value], dtype=bool),
        id="numpy-bool-array",
    ),
)

_MALFORMED_NUMERICS = (
    pytest.param("2.0", id="string"),
    pytest.param(b"2.0", id="bytes"),
    pytest.param(Decimal("2.0"), id="decimal"),
    pytest.param(Fraction(2, 1), id="fraction"),
    pytest.param(None, id="none"),
    pytest.param(2 + 0j, id="complex"),
)


class _MyInt(int):
    pass


class _MyFloat(float):
    pass


class _MyNumpyFloat(np.float64):
    pass


_SUBCLASS_NUMERICS = (
    pytest.param(_MyInt(2), id="int-subclass"),
    pytest.param(_MyFloat(2.0), id="float-subclass"),
    pytest.param(_MyNumpyFloat(2.0), id="numpy-real-subclass"),
)


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


def _full_curve(overrides=None):
    values = {length: float(length) for length in EXPECTED_LADDER_LENGTHS}
    values.update(overrides or {})
    return sorted(values.items())


# --- potency: the three IntaRNA states ---------------------------------------


def test_potency_uses_real_intarna_and_rnaup():
    row = _row(18, intarna_dg=-27.0, intarna_no_duplex=False, rnaup_dg=-25.0)
    out = compute_window([row])[0]
    # potency = -mean(-27, -25) = +26
    assert out["potency"] == pytest.approx(26.0)
    assert out["potency_incomplete"] is False


def test_secondary_on_target_axes_never_change_potency_or_window():
    backbone = _row(18, intarna_dg=-27.0, intarna_no_duplex=False,
                    rnaup_dg=-25.0, off_target_hits=4)
    secondary = dict(backbone, rnaduplex_dg=-999.0, rnacofold_dg=-888.0,
                     rnaplfold_accessibility=0.001)
    a, b = compute_window([backbone, secondary])
    assert b["potency"] == a["potency"]
    assert b["therapeutic_window"] == a["therapeutic_window"]


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


@pytest.mark.parametrize("bad_energy", [float("nan"), float("inf"), float("-inf")])
def test_nonfinite_backbone_energy_forces_incomplete_rung(bad_energy):
    out = compute_window([
        _row(16, intarna_dg=bad_energy, intarna_no_duplex=False,
             rnaup_dg=-10.0, off_target_hits=0)
    ])[0]
    assert out["potency"] is None
    assert out["potency_incomplete"] is True
    assert out["therapeutic_window"] is None
    assert out["window_incomplete"] is True


@pytest.mark.parametrize("bad_burden", [float("nan"), float("inf"), -1, 1.5, True])
def test_invalid_burden_forces_incomplete_rung(bad_burden):
    out = compute_window([
        _row(18, intarna_dg=-20.0, intarna_no_duplex=False,
             rnaup_dg=-10.0, off_target_hits=bad_burden)
    ])[0]
    assert out["off_target_burden"] is None
    assert out["therapeutic_window"] is None
    assert out["window_incomplete"] is True


@pytest.mark.parametrize("encode_bool", _BOOLEAN_ENCODERS)
def test_compute_window_rejects_boolean_forms_for_energy_and_burden(encode_bool):
    out = compute_window([
        _row(
            18,
            intarna_dg=encode_bool(True),
            intarna_no_duplex=False,
            off_target_hits=encode_bool(False),
        )
    ])[0]
    assert out["potency"] is None
    assert out["off_target_burden"] is None
    assert out["therapeutic_window"] is None
    assert out["window_incomplete"] is True


@pytest.mark.parametrize("malformed", _MALFORMED_NUMERICS)
def test_compute_window_rejects_non_allowlisted_numeric_like_values(malformed):
    out = compute_window([
        _row(
            18,
            intarna_dg=malformed,
            intarna_no_duplex=False,
            off_target_hits=malformed,
        )
    ])[0]
    assert out["potency"] is None
    assert out["off_target_burden"] is None
    assert out["therapeutic_window"] is None
    assert out["window_incomplete"] is True


@pytest.mark.parametrize("subclass_value", _SUBCLASS_NUMERICS)
def test_compute_window_rejects_real_number_subclasses(subclass_value):
    out = compute_window([
        _row(
            18,
            intarna_dg=subclass_value,
            intarna_no_duplex=False,
            off_target_hits=subclass_value,
        )
    ])[0]
    assert out["potency"] is None
    assert out["off_target_burden"] is None
    assert out["therapeutic_window"] is None
    assert out["window_incomplete"] is True


def test_compute_window_accepts_numpy_real_scalars():
    out = compute_window([
        _row(
            np.int64(18),
            intarna_dg=np.float64(-20.0),
            intarna_no_duplex=False,
            off_target_hits=np.int64(1),
        )
    ], off_target_pseudocount=np.float64(1.0))[0]
    assert out["length"] == 18
    assert out["potency"] == pytest.approx(20.0)
    assert out["therapeutic_window"] == pytest.approx(10.0)
    assert out["window_incomplete"] is False


def test_state_preserved_candidate():
    out = compute_window([_row(18, intarna_dg=-27.0, intarna_no_duplex=False)])[0]
    assert out["state"] == _CANDIDATE


# --- directional verdict -----------------------------------------------------


def test_directional_verdict_supports_when_peak_sub_18():
    curve = _full_curve({length: 20.0 - abs(length - 16) for length in EXPECTED_LADDER_LENGTHS})
    assert directional_verdict(curve) == "supports_shorter_better"


def test_directional_verdict_refutes_when_monotonic_increasing():
    curve = _full_curve()
    assert directional_verdict(curve) == "refutes_shorter_better"


def test_directional_verdict_accepts_finite_float_scalar_lengths_and_values():
    curve = [(float(length), float(length)) for length in EXPECTED_LADDER_LENGTHS]
    assert directional_verdict(curve) == "refutes_shorter_better"


def test_directional_verdict_accepts_numpy_real_scalars():
    curve = [
        (np.int64(length), np.float64(length))
        for length in EXPECTED_LADDER_LENGTHS
    ]
    assert directional_verdict(curve) == "refutes_shorter_better"


@pytest.mark.parametrize("encode_bool", _BOOLEAN_ENCODERS)
def test_directional_verdict_rejects_boolean_forms_in_values_and_lengths(encode_bool):
    boolean_values = [
        (length, encode_bool(False)) for length in EXPECTED_LADDER_LENGTHS
    ]
    assert directional_verdict(boolean_values) == "indeterminate"
    boolean_length = list(_full_curve())
    boolean_length[0] = (encode_bool(True), boolean_length[0][1])
    assert directional_verdict(boolean_length) == "indeterminate"


@pytest.mark.parametrize("malformed", _MALFORMED_NUMERICS)
def test_directional_verdict_rejects_non_allowlisted_numeric_like_values(malformed):
    curve = [(length, malformed) for length in EXPECTED_LADDER_LENGTHS]
    assert directional_verdict(curve) == "indeterminate"


@pytest.mark.parametrize("subclass_value", _SUBCLASS_NUMERICS)
def test_directional_verdict_rejects_real_number_subclasses(subclass_value):
    curve = [
        (length, subclass_value) for length in EXPECTED_LADDER_LENGTHS
    ]
    assert directional_verdict(curve) == "indeterminate"


def test_directional_verdict_indeterminate_peak_at_or_above_18_nonmonotonic():
    curve = _full_curve({18: 40.0, 20: 3.0, 22: 4.0, 24: 5.0, 26: 6.0, 28: 7.0, 30: 8.0})
    assert directional_verdict(curve) == "indeterminate"


def test_directional_verdict_skips_none_values():
    curve = _full_curve({2: None})
    assert directional_verdict(curve) == "refutes_shorter_better"


def test_directional_verdict_indeterminate_too_few_points():
    assert directional_verdict([(18, 5.0)]) == "indeterminate"
    assert directional_verdict([(2, None), (18, None)]) == "indeterminate"


def test_peak_boundary_at_18_is_not_shorter_better():
    # peak exactly at 18 (not sub-18) and monotonic up-to-18 then flat/down.
    curve = _full_curve({length: min(float(length), 18.0) for length in EXPECTED_LADDER_LENGTHS})
    # monotonic non-decreasing -> refutes (longer never worse)
    assert directional_verdict(curve) == "refutes_shorter_better"


# --- MED-1: a sub-18 plateau top must NOT be read as shorter-better -----------


def test_sub18_plateau_is_not_supports_shorter_better():
    # The 16-mer only TIES the 18/24/30-mers (length-insensitive plateau).
    # A first-occurrence argmax at L16 must NOT manufacture a "supports" verdict.
    curve = _full_curve({length: min(float(length), 16.0) for length in EXPECTED_LADDER_LENGTHS})
    assert directional_verdict(curve) != "supports_shorter_better"
    # non-decreasing -> refutes.
    assert directional_verdict(curve) == "refutes_shorter_better"


def test_strict_sub18_peak_still_supports():
    # A sub-18 length strictly above every longer rung -> genuine support.
    curve = _full_curve({length: 20.0 - abs(length - 16) for length in EXPECTED_LADDER_LENGTHS})
    assert directional_verdict(curve) == "supports_shorter_better"


def test_tie_at_global_max_spanning_sub18_and_over18_is_not_supports():
    # sub-18 peak (L12) ties a >=18 rung (L18) at the max -> not strict -> not supports.
    curve = _full_curve({12: 99.0, 18: 99.0, 20: 4.0, 22: 3.0, 24: 2.0, 26: 1.0, 28: 0.0, 30: -1.0})
    assert directional_verdict(curve) != "supports_shorter_better"


def test_sub18_only_plateau_top_with_later_drop_is_not_supports():
    # Two sub-18 lengths tie at the top then it drops: non-monotonic, but the
    # peak only ties L16, so no STRICT sub-18 maximum -> indeterminate, not supports.
    curve = _full_curve({12: 99.0, 16: 99.0, 18: 8.0, 20: 7.0, 22: 6.0, 24: 5.0, 26: 4.0, 28: 3.0, 30: 2.0})
    assert directional_verdict(curve) != "supports_shorter_better"


# --- LOW-1 / LOW-2: hardening -------------------------------------------------


@pytest.mark.parametrize("noncanonical", [0, 2, True])
def test_noncanonical_pseudocount_is_rejected_deterministically(noncanonical):
    row = _row(18, intarna_dg=-27.0, intarna_no_duplex=False, off_target_hits=0)
    with pytest.raises(ValueError, match="canonical value"):
        compute_window([row], off_target_pseudocount=noncanonical)


@pytest.mark.parametrize("encode_bool", _BOOLEAN_ENCODERS)
def test_all_boolean_pseudocount_forms_are_rejected_before_equality(encode_bool):
    row = _row(18, intarna_dg=-27.0, intarna_no_duplex=False, off_target_hits=0)
    with pytest.raises(ValueError, match="canonical value"):
        compute_window([row], off_target_pseudocount=encode_bool(True))


@pytest.mark.parametrize("malformed", _MALFORMED_NUMERICS)
def test_non_allowlisted_pseudocounts_are_rejected_before_equality(malformed):
    row = _row(18, intarna_dg=-27.0, intarna_no_duplex=False, off_target_hits=0)
    with pytest.raises(ValueError, match="canonical value"):
        compute_window([row], off_target_pseudocount=malformed)


@pytest.mark.parametrize("subclass_value", _SUBCLASS_NUMERICS)
def test_real_number_subclass_pseudocounts_are_rejected(subclass_value):
    row = _row(18, intarna_dg=-27.0, intarna_no_duplex=False, off_target_hits=0)
    with pytest.raises(ValueError, match="canonical value"):
        compute_window([row], off_target_pseudocount=subclass_value)


def test_canonical_pseudocount_complete_ladder_is_accepted():
    rows = compute_window([
        _row(length, intarna_dg=-float(length), intarna_no_duplex=False,
             rnaup_dg=-float(length), off_target_hits=0)
        for length in EXPECTED_LADDER_LENGTHS
    ], off_target_pseudocount=OFF_TARGET_PSEUDOCOUNT)
    payload = prereg_payload(rows, hypothesis="shorter_is_better")
    assert payload["curve_complete"] is True
    assert payload["directional_prediction"] == "refutes_shorter_better"


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


def test_prereg_prediction_is_none_without_complete_window_curve():
    # Sparse/missing window values cannot emit an endpoint verdict.
    rows = compute_window([
        _row(16, intarna_dg=-22.0, intarna_no_duplex=False, rnaup_dg=-22.0),
        _row(18, intarna_dg=-27.0, intarna_no_duplex=False, rnaup_dg=-26.0),
        _row(24, intarna_dg=-38.0, intarna_no_duplex=False, rnaup_dg=-38.0),
    ])
    payload = prereg_payload(rows, hypothesis="shorter_is_better")
    assert payload["directional_prediction_axis"] == "therapeutic_window"
    assert payload["directional_prediction"] is None
    assert payload["curve_complete"] is False


def test_prereg_prediction_requires_every_ladder_rung():
    rows = compute_window([
        _row(length, intarna_dg=-float(length), intarna_no_duplex=False,
             rnaup_dg=-float(length), off_target_hits=1)
        for length in EXPECTED_LADDER_LENGTHS
    ])
    complete = prereg_payload(rows, hypothesis="shorter_is_better")
    assert complete["curve_complete"] is True
    assert complete["directional_prediction"] == "refutes_shorter_better"
    missing = prereg_payload(rows[:-1], hypothesis="shorter_is_better")
    assert missing["curve_complete"] is False
    assert missing["directional_prediction"] is None


def test_nonfinite_complete_curve_cannot_emit_directional_prediction():
    rows = compute_window([
        _row(length, intarna_dg=(-float("inf") if length == 16 else -float(length)),
             intarna_no_duplex=False, rnaup_dg=-float(length), off_target_hits=0)
        for length in EXPECTED_LADDER_LENGTHS
    ])
    payload = prereg_payload(rows, hypothesis="shorter_is_better")
    rung = next(item for item in payload["curve"] if item["length"] == 16)
    assert rung["potency"] is None and rung["therapeutic_window"] is None
    assert payload["curve_complete"] is False
    assert payload["directional_prediction"] is None
    assert not any(
        isinstance(value, float) and not math.isfinite(value)
        for item in payload["curve"]
        for value in item.values()
    )


@pytest.mark.parametrize("field", ["potency", "off_target_burden", "therapeutic_window"])
def test_prereg_public_rows_fail_closed_when_any_window_component_nonfinite(field):
    rows = [
        {
            "length": length,
            "potency": float(length),
            "potency_incomplete": False,
            "off_target_burden": 0,
            "therapeutic_window": float(length),
            "window_incomplete": False,
        }
        for length in EXPECTED_LADDER_LENGTHS
    ]
    rows[3][field] = float("inf")
    payload = prereg_payload(rows, hypothesis="shorter_is_better")
    bad = payload["curve"][3]
    assert bad["therapeutic_window"] is None
    assert bad["window_incomplete"] is True
    assert payload["curve_complete"] is False
    assert payload["directional_prediction"] is None


def test_prereg_rejects_formula_mismatch_across_rounding_boundary():
    potency = 50.00000049
    rows = [
        {
            "length": length,
            "potency": potency,
            "potency_incomplete": False,
            "off_target_burden": 0,
            "therapeutic_window": (
                50.00000051 if length == 16 else potency
            ),
            "window_incomplete": False,
        }
        for length in EXPECTED_LADDER_LENGTHS
    ]
    payload = prereg_payload(rows, hypothesis="shorter_is_better")
    rung = next(row for row in payload["curve"] if row["length"] == 16)
    assert rung["therapeutic_window"] is None
    assert rung["window_incomplete"] is True
    assert payload["curve_complete"] is False
    assert payload["directional_prediction"] is None


@pytest.mark.parametrize("encode_bool", _BOOLEAN_ENCODERS)
def test_prereg_boolean_quantity_tuple_cannot_emit_verdict(encode_bool):
    rows = [
        {
            "length": length,
            "potency": encode_bool(length == 16),
            "potency_incomplete": False,
            "off_target_burden": encode_bool(False),
            "therapeutic_window": encode_bool(length == 16),
            "window_incomplete": False,
        }
        for length in EXPECTED_LADDER_LENGTHS
    ]
    payload = prereg_payload(rows, hypothesis="shorter_is_better")
    assert payload["curve_complete"] is False
    assert payload["directional_prediction"] is None
    assert all(row["therapeutic_window"] is None for row in payload["curve"])


@pytest.mark.parametrize("malformed", _MALFORMED_NUMERICS)
def test_prereg_non_allowlisted_numeric_like_tuple_cannot_emit_verdict(malformed):
    rows = [
        {
            "length": length,
            "potency": malformed,
            "potency_incomplete": False,
            "off_target_burden": malformed,
            "therapeutic_window": malformed,
            "window_incomplete": False,
        }
        for length in EXPECTED_LADDER_LENGTHS
    ]
    payload = prereg_payload(rows, hypothesis="shorter_is_better")
    assert payload["curve_complete"] is False
    assert payload["directional_prediction"] is None
    assert all(row["therapeutic_window"] is None for row in payload["curve"])


@pytest.mark.parametrize("subclass_value", _SUBCLASS_NUMERICS)
def test_prereg_real_number_subclass_tuple_cannot_emit_verdict(subclass_value):
    rows = [
        {
            "length": length,
            "potency": subclass_value,
            "potency_incomplete": False,
            "off_target_burden": subclass_value,
            "therapeutic_window": subclass_value,
            "window_incomplete": False,
        }
        for length in EXPECTED_LADDER_LENGTHS
    ]
    payload = prereg_payload(rows, hypothesis="shorter_is_better")
    assert payload["curve_complete"] is False
    assert payload["directional_prediction"] is None
    assert all(row["therapeutic_window"] is None for row in payload["curve"])


@pytest.mark.parametrize("encode_bool", _BOOLEAN_ENCODERS)
def test_prereg_rejects_all_boolean_length_forms(encode_bool):
    row = {
        "length": encode_bool(True),
        "potency": 1.0,
        "potency_incomplete": False,
        "off_target_burden": 0,
        "therapeutic_window": 1.0,
        "window_incomplete": False,
    }
    with pytest.raises(ValueError, match="length"):
        prereg_payload([row], hypothesis="shorter_is_better")


@pytest.mark.parametrize("bad_length", [True, "2.0", 2.9, float("inf"), 32])
def test_prereg_rejects_noncanonical_ladder_lengths_before_conversion(bad_length):
    rows = [{
        "length": bad_length,
        "potency": 1.0,
        "potency_incomplete": False,
        "off_target_burden": 0,
        "therapeutic_window": 1.0,
        "window_incomplete": False,
    }]
    with pytest.raises(ValueError, match="length"):
        prereg_payload(rows, hypothesis="shorter_is_better")


def test_locked_prereg_cannot_be_relocked_or_overwritten(tmp_path):
    payload = prereg_payload(_window_rows(), hypothesis="shorter_is_better")
    path = tmp_path / "prereg.json"
    lock_prereg(payload, path, locked=True)
    with pytest.raises(FileExistsError, match="durably sealed"):
        lock_prereg(payload, path, locked=True)
    changed = dict(payload, hypothesis="changed")
    with pytest.raises(FileExistsError, match="durably sealed"):
        lock_prereg(changed, path)


def test_verify_prereg_detects_locked_flag_mutation(tmp_path):
    payload = prereg_payload(_window_rows(), hypothesis="shorter_is_better")
    path = tmp_path / "prereg.json"
    lock_prereg(payload, path, locked=True)
    record = json.loads(path.read_text())
    record["locked"] = False
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="seal metadata"):
        verify_prereg(path)


def test_lock_prereg_migrates_only_hash_valid_legacy_unlocked_draft(tmp_path):
    old_payload = prereg_payload(_window_rows(), hypothesis="old")
    path = tmp_path / "prereg.json"
    legacy = {"payload": old_payload, "hash": prereg_hash(old_payload), "locked": False}
    path.write_text(json.dumps(legacy))

    new_payload = prereg_payload(_window_rows(), hypothesis="new")
    migrated = lock_prereg(new_payload, path)
    assert migrated["locked"] is False
    assert "seal_hash" in migrated
    assert verify_prereg(path)["payload"]["hypothesis"] == "new"


def test_lock_prereg_rejects_legacy_locked_record_without_seal(tmp_path):
    payload = prereg_payload(_window_rows(), hypothesis="old")
    path = tmp_path / "prereg.json"
    path.write_text(json.dumps({
        "payload": payload, "hash": prereg_hash(payload), "locked": True
    }))
    with pytest.raises(ValueError, match="legacy locked"):
        lock_prereg(payload, path)


def test_lock_prereg_uses_interprocess_lock_and_unique_tempfiles(tmp_path, monkeypatch):
    payload = prereg_payload(_window_rows(), hypothesis="shorter_is_better")
    path = tmp_path / "prereg.json"
    names = []
    real_mkstemp = __import__("tempfile").mkstemp

    def recording_mkstemp(*args, **kwargs):
        fd, name = real_mkstemp(*args, **kwargs)
        names.append(name)
        return fd, name

    monkeypatch.setattr("feng.titration.window.tempfile.mkstemp", recording_mkstemp)
    lock_prereg(payload, path)
    lock_prereg(payload, path)
    assert len(names) == 2 and len(set(names)) == 2
    assert path.with_name("prereg.json.lock").exists()


def test_oligo_length_constant_is_18():
    # guards the sub-18 boundary the verdict logic depends on.
    assert OLIGO_LENGTH == 18


def _write_terminal_off_target_artifacts(scan_dir: Path, hoff_dir: Path) -> list[dict]:
    rows = []
    scan_dir.mkdir()
    hoff_dir.mkdir()
    for length in EXPECTED_LADDER_LENGTHS:
        row = _row(
            length,
            intarna_dg=-float(length),
            intarna_no_duplex=False,
            rnaup_dg=-float(length),
        )
        rows.append(row)
        manifest = scan_dir / f"sites_L{length:02d}.tsv"
        informative = length >= arm_e_scan.MINIMUM_INFORMATIVE_LENGTH
        manifest.write_text(
            "\t".join(arm_e_scan.SITE_COLUMNS) + "\n" +
            ("".join(
                f"{length}\t{'N' * length}\tchr1\t{i + 1}\t{i + length}\t+\t"
                f"{'N' * length}\tr{i}\t0\t0\n"
                for i in range(length)
            ) if informative else "")
        )
        manifest_sha = hashlib.sha256(manifest.read_bytes()).hexdigest()
        (scan_dir / f"burden_L{length:02d}.json").write_text(json.dumps({
            "schema_version": arm_e_scan.SCAN_ARTIFACT_SCHEMA_VERSION,
            "algorithm_version": arm_e_scan.SCAN_ALGORITHM_VERSION,
            "length": length,
            "aso_5to3": "N" * length,
            "off_target_hits": length if informative else None,
            "off_target_burden": length if informative else None,
            "complete": informative,
            "terminal": True,
            "reason": None if informative else "length_below_minimum_informative_length",
            "manifest_sha256": manifest_sha,
            "input_digest": "0" * 64,
            "state": _CANDIDATE,
        }))
        (hoff_dir / f"pangolin_hoff_L{length:02d}.json").write_text(json.dumps({
            "schema_version": arm_e_hoff.HOFF_ARTIFACT_SCHEMA_VERSION,
            "algorithm_version": arm_e_hoff.HOFF_ALGORITHM_VERSION,
            "length": length,
            "aso_5to3": "N" * length,
            "pangolin_off_target_burden": float(length) if informative else None,
            "pangolin_off_target_incomplete": not informative,
            "reason": None if informative else "length_below_minimum_informative_length",
            "n_sites": length if informative else 0,
            "manifest_sha256": manifest_sha,
            "input_digest": "1" * 64,
            "terminal": True,
            "state": _CANDIDATE,
        }))
    return rows


def test_artifact_assembler_requires_all_15_and_feeds_scan_burdens(tmp_path):
    scan_dir, hoff_dir = tmp_path / "scan", tmp_path / "hoff"
    rows = _write_terminal_off_target_artifacts(scan_dir, hoff_dir)
    assembled = assemble_off_target_windows(rows, scan_dir=scan_dir, hoff_dir=hoff_dir)
    assert [row["off_target_hits"] for row in assembled[:3]] == [None, None, None]
    assert [row["off_target_hits"] for row in assembled[3:]] == list(EXPECTED_LADDER_LENGTHS[3:])
    assert all(row["window_incomplete"] for row in assembled[:3])
    assert all(not row["window_incomplete"] for row in assembled[3:])
    assert all(row["off_target_artifacts_terminal"] for row in assembled)
    assert prereg_payload(assembled, hypothesis="shorter_is_better")["curve_complete"] is True

    (hoff_dir / "pangolin_hoff_L30.json").unlink()
    incomplete = assemble_off_target_windows(rows, scan_dir=scan_dir, hoff_dir=hoff_dir)
    assert next(row for row in incomplete if row["length"] == 30)["off_target_hits"] is None
    payload = prereg_payload(incomplete, hypothesis="shorter_is_better")
    assert payload["curve_complete"] is False
    assert payload["directional_prediction"] is None

    rows = _write_terminal_off_target_artifacts(tmp_path / "scan3", tmp_path / "hoff3")
    (tmp_path / "hoff3" / "pangolin_hoff_L04.json").unlink()
    missing_short = assemble_off_target_windows(
        rows, scan_dir=tmp_path / "scan3", hoff_dir=tmp_path / "hoff3"
    )
    payload = prereg_payload(missing_short, hypothesis="shorter_is_better")
    assert payload["curve_complete"] is False
    assert payload["directional_prediction"] is None


def test_artifact_assembler_rejects_manifest_hash_or_schema_tamper(tmp_path):
    scan_dir, hoff_dir = tmp_path / "scan", tmp_path / "hoff"
    rows = _write_terminal_off_target_artifacts(scan_dir, hoff_dir)
    (scan_dir / "sites_L18.tsv").write_text("tampered\n")
    with pytest.raises(ValueError, match="manifest SHA-256"):
        assemble_off_target_windows(rows, scan_dir=scan_dir, hoff_dir=hoff_dir)

    rows = _write_terminal_off_target_artifacts(tmp_path / "scan2", tmp_path / "hoff2")
    bad = tmp_path / "scan2" / "burden_L18.json"
    payload = json.loads(bad.read_text())
    payload["schema_version"] = -1
    bad.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="schema"):
        assemble_off_target_windows(rows, scan_dir=tmp_path / "scan2", hoff_dir=tmp_path / "hoff2")
