"""Therapeutic-window aggregation + preregistration lock for Arm E (StericRuler).

This layer collapses the per-mer CPU score panel from
:func:`feng.titration.evaluate.evaluate_ladder` into the arm's pre-declared
endpoint -- a per-length **therapeutic window** = on-target potency divided by
off-target burden -- and provides a deterministic preregistration lock so the
predicted window-vs-length curve (and the directional "shorter is better"
prediction) can be hashed and sealed *before* the wet-lab length series is read.

Every formula here is a **candidate**, pre-registration-eligible design choice;
they are kept transparent and untuned on purpose. The definitions:

* **Potency** (``potency``) -- ``-mean`` of the usable on-target binding
  energies. The backbone is IntaRNA on-target ΔG and RNAup
  (accessibility-corrected) ΔG, both negative-is-favourable, so negating the
  mean turns a deeper duplex into a larger positive potency. Two hard rules from
  the adversarial review are honoured:

  - IntaRNA's ``0.0`` "no favourable duplex" reading (``intarna_no_duplex``) and
    its sub-seed ``None`` are NEVER averaged as real binding.
  - ``rnaduplex_dg`` / ``rnacofold_dg`` (a perfect reverse-complement duplex by
    construction) and ``rnaplfold_accessibility`` (already folded into RNAup) are
    carried through as secondary/diagnostic columns but do NOT enter the potency
    backbone -- no double-counting of accessibility, no composition-driven axis
    dominating. A mer with no usable backbone axis has ``potency = None`` and
    ``potency_incomplete = True`` rather than a spurious zero.

* **Off-target burden** (``off_target_burden``) -- the near-complement hit count
  from the off-target scan. If the scan was not run (``off_target_hits is None``)
  the burden is ``None`` and the window is flagged incomplete; a missing scan is
  never silently treated as zero burden (that would fabricate perfect
  specificity).

* **Therapeutic window** (``therapeutic_window``) --
  ``potency / (burden + OFF_TARGET_PSEUDOCOUNT)``. The Laplace pseudocount both
  guards divide-by-zero and keeps a zero-off-target mer finite (window = potency,
  i.e. maximal specificity). ``None`` + ``window_incomplete = True`` whenever
  potency or burden is missing.

Freeze-safe: consumes only the candidate sequence-derived score panel; reads no
treatment alignment, contrast, rank, or sealed ΔPSI. The preregistration writer
defaults to an unlocked DRAFT -- cutting the real ``locked = true`` seal is a
separate, gated step pending user sign-off and the GPU potency axes (MD /
deep-splice), which this CPU-only curve does not yet include.
"""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Mapping, Sequence

from feng.bench.scores import OLIGO_LENGTH
from feng.honesty.states import ResultState

_CANDIDATE = ResultState.CANDIDATE.value

# Candidate, pre-registration-eligible constants.
OFF_TARGET_PSEUDOCOUNT = 1
FORMULA_VERSION = "armE-window-v1"
_HASH_PRECISION = 6  # decimals floats are rounded to before canonical hashing

# Human-readable record of every scoring choice, embedded verbatim in the
# preregistration payload so the sealed hash covers the formulae, not just the
# numbers.
AXIS_DEFINITIONS = {
    "potency": (
        "negative mean of the usable on-target binding energies "
        "(IntaRNA on-target dg and RNAup accessibility-corrected dg); "
        "IntaRNA no_duplex/None excluded; RNAduplex/RNAcofold/RNAplfold "
        "accessibility excluded from the backbone (secondary/diagnostic)"
    ),
    "off_target_burden": (
        "count of near-complement off-target hits; None (not zero) when the "
        "off-target scan was not run"
    ),
    "therapeutic_window": (
        "potency / (off_target_burden + OFF_TARGET_PSEUDOCOUNT); "
        f"OFF_TARGET_PSEUDOCOUNT={OFF_TARGET_PSEUDOCOUNT}"
    ),
    "formula_version": FORMULA_VERSION,
}

_WINDOW_KEYS = (
    "potency",
    "potency_incomplete",
    "potency_unfavourable",
    "off_target_burden",
    "therapeutic_window",
    "window_incomplete",
)


def _potency(row: Mapping) -> tuple[float | None, bool]:
    """``(potency, incomplete)`` from a row's usable on-target energies.

    Returns ``(None, True)`` when no backbone axis is usable.
    """
    energies: list[float] = []

    dg = row.get("intarna_dg")
    no_duplex = row.get("intarna_no_duplex")
    if dg is not None and not no_duplex:
        energies.append(float(dg))

    rnaup = row.get("rnaup_dg")
    if rnaup is not None:
        energies.append(float(rnaup))

    if not energies:
        return None, True
    # Favourable binding is a negative energy; negate so a deeper duplex is a
    # larger positive potency.
    return -sum(energies) / len(energies), False


def compute_window(
    evaluated_rows: Sequence[Mapping],
    *,
    off_target_pseudocount: int = OFF_TARGET_PSEUDOCOUNT,
) -> list[dict]:
    """Attach ``potency`` / ``off_target_burden`` / ``therapeutic_window`` columns.

    Each input row (from :func:`evaluate_ladder`) is copied and extended with the
    keys in ``_WINDOW_KEYS``. The honesty ``state`` is preserved as ``candidate``.
    """
    # The divide-by-zero guard relies on a >=1 pseudocount; enforce it rather than
    # let a caller override to 0 and crash on a zero-off-target mer.
    if off_target_pseudocount < 1:
        raise ValueError(
            f"off_target_pseudocount must be >= 1 (got {off_target_pseudocount})"
        )
    out: list[dict] = []
    for row in evaluated_rows:
        r = dict(row)
        r["state"] = r.get("state", _CANDIDATE)

        potency, incomplete = _potency(row)
        r["potency"] = potency
        r["potency_incomplete"] = incomplete
        # A negative potency means a usable axis reported an UNFAVOURABLE (positive)
        # binding energy; flag it explicitly so a sign-flipped window is never read
        # as a real result.
        r["potency_unfavourable"] = potency is not None and potency < 0

        burden = row.get("off_target_hits")
        r["off_target_burden"] = burden

        if potency is None or burden is None:
            r["therapeutic_window"] = None
            r["window_incomplete"] = True
        else:
            r["therapeutic_window"] = potency / (burden + off_target_pseudocount)
            r["window_incomplete"] = False

        out.append(r)
    return out


def directional_verdict(curve: Sequence[tuple[int, float | None]]) -> str:
    """Does ``curve`` support, refute, or leave open "shorter is better"?

    ``curve`` is a sequence of ``(length, value)`` pairs (``None`` values are
    dropped). Returns:

    * ``"refutes_shorter_better"`` -- values are monotonically NON-DECREASING with
      length (a longer mer is never worse -> the drug 18-mer or longer wins). This
      is checked FIRST, so a length-insensitive flat/plateau curve refutes rather
      than manufacturing a false confirmation.
    * ``"supports_shorter_better"`` -- a sub-18 length is the peak and STRICTLY
      exceeds every longer-length value (a mere tie with a longer rung is NOT
      support: a plateau top does not show shorter is better).
    * ``"indeterminate"`` -- fewer than two usable points, or a non-monotonic curve
      whose peak is at/above 18 or whose sub-18 peak only ties a longer rung.
    """
    pts = [(int(length), float(value)) for length, value in curve if value is not None]
    if len(pts) < 2:
        return "indeterminate"
    pts.sort()
    lengths = [ln for ln, _ in pts]
    values = [v for _, v in pts]

    # Refutes first: a non-decreasing curve (including any flat plateau) means a
    # longer mer is never worse -- shorter-is-better cannot be read off it.
    if all(values[i] <= values[i + 1] for i in range(len(values) - 1)):
        return "refutes_shorter_better"

    # Support requires a sub-18 peak that STRICTLY beats every longer-length rung.
    peak_idx = max(range(len(pts)), key=lambda i: values[i])  # first occurrence of max
    peak_length, peak_value = lengths[peak_idx], values[peak_idx]
    if peak_length < OLIGO_LENGTH:
        longer = [values[i] for i in range(len(pts)) if lengths[i] > peak_length]
        if longer and all(peak_value > v for v in longer):
            return "supports_shorter_better"
    return "indeterminate"


def prereg_payload(
    window_rows: Sequence[Mapping],
    *,
    hypothesis: str,
    formula_version: str = FORMULA_VERSION,
) -> dict:
    """Canonical, deterministic preregistration payload for the window curve.

    Captures the length-sorted window-vs-length curve, the axis definitions, the
    formula version, and the current directional prediction. The prediction is
    read off the therapeutic-window curve when off-target burden is available for
    at least one mer, else off the potency curve (CPU-only pass) -- the axis used
    is recorded in ``directional_prediction_axis``.
    """
    curve = []
    for row in sorted(window_rows, key=lambda r: int(r["length"])):
        curve.append(
            {
                "length": int(row["length"]),
                "potency": _round(row.get("potency")),
                "therapeutic_window": _round(row.get("therapeutic_window")),
                "window_incomplete": bool(row.get("window_incomplete", True)),
            }
        )

    window_curve = [(c["length"], c["therapeutic_window"]) for c in curve]
    potency_curve = [(c["length"], c["potency"]) for c in curve]
    if any(value is not None for _, value in window_curve):
        prediction_axis = "therapeutic_window"
        prediction = directional_verdict(window_curve)
    else:
        prediction_axis = "potency"
        prediction = directional_verdict(potency_curve)

    return {
        "formula_version": formula_version,
        "hypothesis": hypothesis,
        "axis_definitions": AXIS_DEFINITIONS,
        "curve": curve,
        "directional_prediction": prediction,
        "directional_prediction_axis": prediction_axis,
        "state": _CANDIDATE,
    }


def _round(value):
    """Round a float to the hashing precision; pass ``None`` through."""
    if value is None:
        return None
    return round(float(value), _HASH_PRECISION)


def _canonicalize(obj):
    """Recursively round floats so the canonical JSON hash is reproducible."""
    if isinstance(obj, float):
        return round(obj, _HASH_PRECISION)
    if isinstance(obj, dict):
        return {key: _canonicalize(val) for key, val in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_canonicalize(item) for item in obj]
    return obj


def prereg_hash(payload: Mapping) -> str:
    """Stable SHA-256 over the canonical JSON of ``payload`` (sorted keys)."""
    canonical = json.dumps(
        _canonicalize(payload),
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    )
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def lock_prereg(payload: Mapping, path: str | Path, *, locked: bool = False) -> dict:
    """Write a preregistration lock file for ``payload`` and return the record.

    Defaults to an unlocked **DRAFT** (``locked = False``). Cutting the real
    ``locked = True`` seal is a deliberate, gated step -- it must happen only
    after user sign-off and once the GPU potency axes (MD / deep-splice) are part
    of the curve -- so this function does not fire it by default; a caller passes
    ``locked=True`` only as that explicit gated action.
    """
    record = {
        "payload": dict(payload),
        "hash": prereg_hash(payload),
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "locked": bool(locked),
        "state": _CANDIDATE,
    }
    Path(path).write_text(json.dumps(record, indent=2, sort_keys=True))
    return record


def verify_prereg(path: str | Path) -> dict:
    """Read a lock file and re-verify its hash against the stored payload.

    Raises ``ValueError`` if the recomputed hash does not match -- i.e. the
    payload was edited after the lock, or the canonical format drifted.
    """
    record = json.loads(Path(path).read_text())
    recomputed = prereg_hash(record["payload"])
    if recomputed != record["hash"]:
        raise ValueError(
            "preregistration hash mismatch: payload tampered or format drift "
            f"(stored {record['hash']}, recomputed {recomputed})"
        )
    return record
