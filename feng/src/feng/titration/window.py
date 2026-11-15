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
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Mapping, Sequence

import numpy as np

from feng.bench.scores import OLIGO_LENGTH
from feng.honesty.states import ResultState
from feng.titration import arm_e_hoff, arm_e_scan

_CANDIDATE = ResultState.CANDIDATE.value

# Candidate, pre-registration-eligible constants.
OFF_TARGET_PSEUDOCOUNT = 1
FORMULA_VERSION = "armE-window-v2"
EXPECTED_LADDER_LENGTHS = tuple(range(2, 31, 2))
MINIMUM_INFORMATIVE_LENGTH = arm_e_scan.MINIMUM_INFORMATIVE_LENGTH
INFORMATIVE_LADDER_LENGTHS = tuple(
    length for length in EXPECTED_LADDER_LENGTHS
    if length >= MINIMUM_INFORMATIVE_LENGTH
)
_HASH_PRECISION = 6  # decimals floats are rounded to before canonical hashing
POTENCY_BACKBONE = ("intarna_dg", "rnaup_dg")
_ALLOWED_REAL_SCALAR_TYPES = frozenset({
    int,
    float,
    np.float64,
    np.float32,
    np.float16,
    np.int64,
    np.int32,
    np.int16,
    np.int8,
    np.uint64,
    np.uint32,
    np.uint16,
    np.uint8,
})

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
        "count of canonical near-complement genomic sites (one best alignment "
        "per stranded contiguous exonic-union locus); None (not zero) when the "
        "scan/H-off artifact set is missing or incomplete. Universe excludes "
        "exon-exon junction sequences and mature-transcript junction off-targets"
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

    for axis in POTENCY_BACKBONE:
        value = row.get(axis)
        if value is None:
            continue
        if axis == "intarna_dg" and row.get("intarna_no_duplex"):
            continue
        numeric = _strict_finite_scalar(value)
        if numeric is None:
            return None, True
        energies.append(numeric)

    if not energies:
        return None, True
    # Favourable binding is a negative energy; negate so a deeper duplex is a
    # larger positive potency.
    potency = -sum(energies) / len(energies)
    return (potency, False) if math.isfinite(potency) else (None, True)


def _strict_finite_scalar(value) -> float | None:
    """Return a finite scalar from the explicit real-number type allowlist.

    Only the exact Python ``int``/``float`` and enumerated NumPy scalar types are
    accepted. Everything else is rejected before ``float`` coercion, including
    subclasses, booleans, strings, bytes, Decimal/Fraction, complex values, and
    every ndarray.
    """
    if type(value) not in _ALLOWED_REAL_SCALAR_TYPES:
        return None
    numeric = float(value)
    return numeric if math.isfinite(numeric) else None


def _burden(value) -> int | None:
    """Normalize an off-target count without inventing a zero."""
    numeric = _strict_finite_scalar(value)
    if numeric is None or numeric < 0 or not numeric.is_integer():
        return None
    return int(numeric)


def _canonical_pseudocount(value) -> int:
    """Validate the single formula-bound pseudocount without coercive equality."""
    numeric = _strict_finite_scalar(value)
    if numeric is None or numeric != OFF_TARGET_PSEUDOCOUNT:
        raise ValueError(
            "off_target_pseudocount is fixed at the canonical value "
            f"{OFF_TARGET_PSEUDOCOUNT} (got {value!r})"
        )
    return OFF_TARGET_PSEUDOCOUNT


def normalize_window_row(
    row: Mapping, *, off_target_pseudocount: int = OFF_TARGET_PSEUDOCOUNT
) -> dict:
    """Fail closed unless potency, burden, and window are jointly valid.

    This is the single public-input validator used by computation,
    preregistration, and report generation. A caller cannot make a complete
    curve by pairing a finite supplied window with a non-finite potency or
    burden. The supplied window must also agree with the declared formula.
    """
    pseudocount = _canonical_pseudocount(off_target_pseudocount)
    out = dict(row)
    out["length"] = _canonical_ladder_length(row.get("length"))
    potency = _strict_finite_scalar(row.get("potency"))
    burden = _burden(row.get("off_target_burden"))
    window = _strict_finite_scalar(row.get("therapeutic_window"))
    expected = (
        potency / (burden + pseudocount)
        if potency is not None and burden is not None
        else None
    )
    valid = (
        expected is not None
        and math.isfinite(expected)
        and window is not None
        # Consistency is deliberately exact and precedes any canonical-output
        # rounding. Near-equal caller-supplied endpoints are not formula values.
        and window == expected
        and not bool(row.get("potency_incomplete", False))
        and not bool(row.get("window_incomplete", True))
    )
    out["potency"] = potency
    out["potency_incomplete"] = potency is None or bool(
        row.get("potency_incomplete", False)
    )
    out["potency_unfavourable"] = potency is not None and potency < 0
    out["off_target_burden"] = burden
    if "off_target_hits" in out:
        out["off_target_hits"] = burden
    # Publish the recomputed endpoint, never the caller's copy.
    out["therapeutic_window"] = expected if valid else None
    out["window_incomplete"] = not valid
    return out


def normalize_window_rows(window_rows: Sequence[Mapping]) -> list[dict]:
    """Normalize public/report window rows through the canonical validator."""
    return [normalize_window_row(row) for row in window_rows]


def _canonical_ladder_length(value) -> int:
    """Return an exact canonical ladder rung; reject coercive aliases."""
    numeric = _strict_finite_scalar(value)
    if numeric is None or not numeric.is_integer():
        raise ValueError(f"ladder length must be a finite non-boolean integer: {value!r}")
    length = int(numeric)
    if length not in EXPECTED_LADDER_LENGTHS:
        raise ValueError(f"length is not a canonical ladder rung: {value!r}")
    return length


def compute_window(
    evaluated_rows: Sequence[Mapping],
    *,
    off_target_pseudocount: int = OFF_TARGET_PSEUDOCOUNT,
) -> list[dict]:
    """Attach ``potency`` / ``off_target_burden`` / ``therapeutic_window`` columns.

    Each input row (from :func:`evaluate_ladder`) is copied and extended with the
    keys in ``_WINDOW_KEYS``. The honesty ``state`` is preserved as ``candidate``.
    """
    # The formula metadata binds one canonical pseudocount; accepting overrides
    # here would make later prereg/report normalization interpret the rows differently.
    pseudocount = _canonical_pseudocount(off_target_pseudocount)
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

        burden = _burden(row.get("off_target_hits"))
        r["off_target_burden"] = burden

        if potency is None or burden is None:
            r["therapeutic_window"] = None
            r["window_incomplete"] = True
        else:
            window = potency / (burden + pseudocount)
            if math.isfinite(window):
                r["therapeutic_window"] = window
                r["window_incomplete"] = False
            else:
                r["therapeutic_window"] = None
                r["window_incomplete"] = True

        out.append(normalize_window_row(r, off_target_pseudocount=pseudocount))
    return out


def _artifact_json(path: Path, *, kind: str, schema_version: int,
                   algorithm_version: str) -> dict | None:
    if not path.is_file() or path.stat().st_size == 0:
        return None
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"{kind} artifact is not valid JSON: {path}") from exc
    if payload.get("schema_version") != schema_version:
        raise ValueError(f"{kind} artifact schema mismatch: {path}")
    if payload.get("algorithm_version") != algorithm_version:
        raise ValueError(f"{kind} artifact algorithm version mismatch: {path}")
    digest = payload.get("input_digest")
    if not isinstance(digest, str) or len(digest) != 64 or any(
        char not in "0123456789abcdef" for char in digest
    ):
        raise ValueError(f"{kind} artifact input/config digest missing or malformed: {path}")
    return payload


def assemble_off_target_windows(
    evaluated_rows: Sequence[Mapping], *, scan_dir: str | Path,
    hoff_dir: str | Path,
) -> list[dict]:
    """Authenticate all 15 scan/H-off shards and compute the window curve.

    Contract: input rows contain each canonical even rung exactly once. For each
    rung, ``sites_Lxx.tsv``, ``burden_Lxx.json``, and
    ``pangolin_hoff_Lxx.json`` are required. JSON schema/algorithm/config-digest,
    rung/ASO identity, and both stored manifest hashes are checked against the
    actual canonical manifest. Corruption raises. A missing or honestly terminal-
    incomplete shard sets that rung's ``off_target_hits=None``; therefore the
    downstream curve remains incomplete and preregistration emits prediction
    ``None``. Only a complete matched pair contributes the scan's canonical-site
    count. H-off is a required completeness gate, not a second hit count.
    """
    scan_dir, hoff_dir = Path(scan_dir), Path(hoff_dir)
    rows_by_length: dict[int, Mapping] = {}
    for row in evaluated_rows:
        length = _canonical_ladder_length(row.get("length"))
        if length in rows_by_length:
            raise ValueError(f"duplicate evaluated ladder rung: {length}")
        rows_by_length[length] = row
    if sorted(rows_by_length) != list(EXPECTED_LADDER_LENGTHS):
        raise ValueError("evaluated rows must contain exactly all 15 canonical ladder rungs")

    assembled: list[dict] = []
    for length in EXPECTED_LADDER_LENGTHS:
        source = rows_by_length[length]
        aso = source.get("aso_5to3")
        manifest = scan_dir / f"sites_L{length:02d}.tsv"
        scan_path = scan_dir / f"burden_L{length:02d}.json"
        hoff_path = hoff_dir / f"pangolin_hoff_L{length:02d}.json"
        scan_payload = _artifact_json(
            scan_path, kind="scan", schema_version=arm_e_scan.SCAN_ARTIFACT_SCHEMA_VERSION,
            algorithm_version=arm_e_scan.SCAN_ALGORITHM_VERSION,
        )
        hoff_payload = _artifact_json(
            hoff_path, kind="H-off", schema_version=arm_e_hoff.HOFF_ARTIFACT_SCHEMA_VERSION,
            algorithm_version=arm_e_hoff.HOFF_ALGORITHM_VERSION,
        )
        row = dict(source)
        row["off_target_artifacts_terminal"] = False
        row["off_target_artifacts_complete"] = False
        row["off_target_artifact_reason"] = None
        if scan_payload is None or hoff_payload is None or not manifest.is_file():
            row["off_target_hits"] = None
            row["off_target_artifact_reason"] = "missing terminal scan/H-off artifact"
            assembled.append(row)
            continue
        row["off_target_artifacts_terminal"] = (
            scan_payload.get("terminal") is True
            and hoff_payload.get("terminal") is True
        )
        manifest_sha = hashlib.sha256(manifest.read_bytes()).hexdigest()
        if (scan_payload.get("manifest_sha256") != manifest_sha
                or hoff_payload.get("manifest_sha256") != manifest_sha):
            raise ValueError(f"manifest SHA-256 mismatch for L={length}")
        if scan_payload.get("length") != length or hoff_payload.get("length") != length:
            raise ValueError(f"artifact rung mismatch for L={length}")
        if scan_payload.get("aso_5to3") != aso or (
            hoff_payload.get("aso_5to3") not in {aso, None}
        ):
            raise ValueError(f"artifact ASO mismatch for L={length}")
        try:
            manifest_rows = list(arm_e_scan.read_site_manifest(manifest))
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"canonical manifest schema invalid for L={length}") from exc
        if any(
            item["length"] != length or item["aso_5to3"] != aso
            for item in manifest_rows
        ):
            raise ValueError(f"canonical manifest rung/ASO mismatch for L={length}")
        complete = (
            scan_payload.get("terminal") is True
            and scan_payload.get("complete") is True
            and hoff_payload.get("terminal") is True
            and hoff_payload.get("pangolin_off_target_incomplete") is False
        )
        burden = _burden(scan_payload.get("off_target_hits")) if complete else None
        if complete and burden is None:
            raise ValueError(f"complete scan has invalid burden for L={length}")
        hoff_burden = _strict_finite_scalar(hoff_payload.get("pangolin_off_target_burden"))
        if complete and (burden != len(manifest_rows)
                         or hoff_payload.get("n_sites") != len(manifest_rows)
                         or hoff_burden is None or hoff_burden < 0):
            raise ValueError(f"terminal artifact/site-count mismatch for L={length}")
        row["off_target_hits"] = burden
        row["off_target_artifacts_complete"] = complete
        if not complete:
            row["off_target_artifact_reason"] = (
                scan_payload.get("reason") or hoff_payload.get("reason")
                or "terminal artifact marked incomplete"
            )
        assembled.append(row)
    return compute_window(assembled)


def directional_verdict(curve: Sequence[tuple[int, float | None]]) -> str:
    """Does ``curve`` support, refute, or leave open "shorter is better"?

    ``curve`` must contain every preregistered even-length rung from 2 through
    30 exactly once. Predetermined-uninformative L2/L4/L6 may be null; every
    informative rung L>=8 must be finite. Missing, duplicate, extra, or null
    informative rungs return ``"indeterminate"`` before classification. Returns:

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
    try:
        pts = [
            (_canonical_ladder_length(length), _strict_finite_scalar(value))
            for length, value in curve
        ]
    except ValueError:
        return "indeterminate"
    raw_lengths = [length for length, _ in pts]
    if (
        len(raw_lengths) != len(set(raw_lengths))
        or sorted(raw_lengths) != list(EXPECTED_LADDER_LENGTHS)
        or any(
            value is None
            for length, value in pts
            if length >= MINIMUM_INFORMATIVE_LENGTH
        )
    ):
        return "indeterminate"
    pts = sorted(
        (length, value) for length, value in pts
        if length >= MINIMUM_INFORMATIVE_LENGTH
    )
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
    read off the therapeutic-window curve only when every preregistered ladder
    rung is present exactly once and has a usable window value. An incomplete
    curve has an explicit ``None`` prediction; potency is never substituted for
    the therapeutic-window endpoint.
    """
    curve = []
    normalized_rows = sorted(
        normalize_window_rows(window_rows), key=lambda r: int(r["length"])
    )
    for row in normalized_rows:
        curve.append(
            {
                "length": int(row["length"]),
                "potency": _round(row.get("potency")),
                "off_target_burden": row.get("off_target_burden"),
                "therapeutic_window": _round(row.get("therapeutic_window")),
                "window_incomplete": bool(row.get("window_incomplete", True)),
            }
        )

    lengths = [int(row["length"]) for row in normalized_rows]
    # Classify on exact normalized floats. Rounding is only a serialization
    # concern and must never create or erase a directional inequality.
    window_curve = [
        (int(row["length"]), row.get("therapeutic_window"))
        for row in normalized_rows
    ]
    terminal_artifacts_present = all(
        row.get("off_target_artifacts_terminal") is True
        for row in normalized_rows
    ) if any("off_target_artifacts_terminal" in row for row in normalized_rows) else True
    informative_curve = [
        (length, value) for length, value in window_curve
        if length >= MINIMUM_INFORMATIVE_LENGTH
    ]
    curve_complete = (
        lengths == list(EXPECTED_LADDER_LENGTHS)
        and terminal_artifacts_present
        and all(value is not None for _, value in informative_curve)
        and not any(
            bool(row.get("window_incomplete", True))
            for row in normalized_rows
            if int(row["length"]) >= MINIMUM_INFORMATIVE_LENGTH
        )
    )
    prediction_axis = "therapeutic_window"
    prediction = directional_verdict(window_curve) if curve_complete else None

    return {
        "formula_version": formula_version,
        "hypothesis": hypothesis,
        "axis_definitions": AXIS_DEFINITIONS,
        "curve": curve,
        "directional_prediction": prediction,
        "directional_prediction_axis": prediction_axis,
        "curve_complete": curve_complete,
        "state": _CANDIDATE,
    }


def _round(value):
    """Round a float to the hashing precision; pass ``None`` through."""
    if value is None:
        return None
    numeric = _strict_finite_scalar(value)
    return round(numeric, _HASH_PRECISION) if numeric is not None else None


def _canonicalize(obj):
    """Recursively round floats so the canonical JSON hash is reproducible."""
    obj_type = type(obj)
    if obj_type is float:
        return round(obj, _HASH_PRECISION)
    if obj_type is dict:
        return {key: _canonicalize(val) for key, val in obj.items()}
    if obj_type in {list, tuple}:
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
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lock_path = path.with_name(f"{path.name}.lock")
    with lock_path.open("a+", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle.fileno(), fcntl.LOCK_EX)
        try:
            # Recheck only after acquiring the interprocess lock: another writer
            # may have sealed the file while this caller was waiting.
            if path.exists():
                # A hash-valid historical unlocked draft may be migrated. A
                # historical locked record without the v2 seal fails closed.
                existing = verify_prereg(path, allow_legacy_draft=True)
                if existing.get("locked") is True:
                    raise FileExistsError(
                        f"preregistration is durably sealed and cannot be overwritten: {path}"
                    )
            payload_digest = prereg_hash(payload)
            record = {
                "payload": dict(payload),
                "hash": payload_digest,
                "timestamp": datetime.now(timezone.utc).isoformat(),
                "locked": bool(locked),
                "state": _CANDIDATE,
            }
            record["seal_hash"] = prereg_hash(
                {"payload_hash": payload_digest, "locked": record["locked"]}
            )
            fd, tmp_name = tempfile.mkstemp(
                prefix=f"{path.name}.tmp.", dir=path.parent
            )
            try:
                with os.fdopen(fd, "w", encoding="utf-8") as tmp_handle:
                    json.dump(record, tmp_handle, indent=2, sort_keys=True)
                    tmp_handle.flush()
                    os.fsync(tmp_handle.fileno())
                os.replace(tmp_name, path)
                dir_fd = os.open(path.parent, os.O_RDONLY)
                try:
                    os.fsync(dir_fd)
                finally:
                    os.close(dir_fd)
            finally:
                if os.path.exists(tmp_name):
                    os.unlink(tmp_name)
            return record
        finally:
            fcntl.flock(lock_handle.fileno(), fcntl.LOCK_UN)


def verify_prereg(path: str | Path, *, allow_legacy_draft: bool = False) -> dict:
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
    if "seal_hash" not in record:
        if allow_legacy_draft and record.get("locked") is False:
            return record
        if record.get("locked") is True:
            raise ValueError("legacy locked preregistration lacks tamper-evident seal")
        raise ValueError("preregistration seal metadata tampered or missing")
    expected_seal_hash = prereg_hash(
        {"payload_hash": record["hash"], "locked": record.get("locked")}
    )
    if record.get("seal_hash") != expected_seal_hash:
        raise ValueError("preregistration seal metadata tampered or missing")
    return record
