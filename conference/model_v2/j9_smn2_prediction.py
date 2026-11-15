"""J9 - first reproducible SMN2 ASO feature/ordering dataset (DESCRIPTIVE, not supervised).

What this is: per-ASO feature values for the F18 / F14 / F10 ISS-N1 ASOs, every
value carrying an explicit evidence class, plus within-feature ordinal ranks and
rank stability across *prespecified* alternative settings.

What this is NOT: a trained, fitted, selected or validated predictor. N = 3
ranked ASOs, no verified negatives, no outcome labels are read or joined. See
``STATISTICAL_STATUS`` and ``request_metric`` - classification metrics are
refused by construction, not by convention.

Baselines (nested, each adds one channel):
  M0  seeded permutation / random-order reference ranks (null ordering)
  M1  + complementarity descriptors (J3 complementarity channel, derived)
  M2  + duplex free energy (J3 duplex channel; unmodified-RNA model_assumption)
  M3  + accessibility (J3 RNAplfold whole-footprint channel; proxy)
  M4  + splice geometry relative to SMN2 exon 7 donor/acceptor (F1, derived)
  M5  + N-ablation Pangolin delta read from the J4 atlas (sequence_ablation_proxy)

Sources are positively allow-listed (``check_source_allowed``). No lab outcome
workbook, gel table or PSI product is in the allowlist, and nothing in this
module opens a path that is not allow-listed.

CLI: ``python -m conference.model_v2.j9_smn2_prediction {build,status,plan}``
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from conference.aso_pipeline import sequences as aso_sequences  # noqa: E402
from conference.model_v2 import j3_scoring  # noqa: E402
from conference.model_v2.trans_intervention import (  # noqa: E402
    ANNOTATION,
    ASSEMBLY,
    EVIDENCE_CLASSES,
    ASOS as TI_ASOS,
    DEFAULT_FASTA,
    FaiFastaReader,
    GenomeAccessError,
    LOCI,
    place,
    reverse_complement,
)

DATASET_ID = "j9_smn2_v1"
DATASET_VERSION = "1.0.0"

REFERENCE_LOCUS_ID = "SMN2_exon7"
ALTERNATE_LOCUS_IDS = ("SMN1_exon7",)

#: ASOs carried into the ranked set. 3UP8 is excluded on purpose (see EXCLUDED_ASOS).
RANKED_ASO_IDS = ("F18", "F14", "F10")
#: id -> reason the ASO is kept out of every rank table.
EXCLUDED_ASOS = {
    "3UP8": "placement_unconfirmed_needs_human_confirmation",
}

BASELINES = ("M0", "M1", "M2", "M3", "M4", "M5")

#: Prespecified fold alternatives, taken verbatim from J2 (jobs/j2_fold.py
#: ``windows=((80, 40), (120, 60), (200, 100))``). J9 choice: the J3 context
#: ``flank`` is set equal to the fold window W so that BOTH the duplex context
#: and the RNAplfold window move together across the alternative.
VIENNA_ALTERNATIVES: Tuple[Tuple[str, int, int, int], ...] = (
    ("w80s40", 80, 80, 40),
    ("w120s60", 120, 120, 60),
    ("w200s100", 200, 200, 100),
)
DEFAULT_VARIANT_ID = "w120s60"

#: Prespecified mismatch-tolerance alternatives for the site-universe sensitivity
#: (J6 ``max_mismatch``; 3 is the J6 default).
MAX_MISMATCH_ALTERNATIVES: Tuple[int, ...] = (0, 1, 2, 3)

M0_PERMUTATIONS = 200
M0_SEED_SALT = "j9-m0-null"

ATLAS_SENTINEL_DIR = "~/orcd/scratch/feng/model_v2/j4_sentinel_v2"
ATLAS_FULL_DIR = "~/orcd/scratch/feng/model_v2/j4_atlas_v1"
DEFAULT_OUT_DIR = "~/orcd/scratch/feng/model_v2/" + DATASET_ID
ATLAS_EVIDENCE_LABEL = "sequence_ablation_proxy"

# --------------------------------------------------------------------------- statistical status
VERIFIED_NEGATIVES_AVAILABLE = False

STATISTICAL_STATUS: Dict[str, Any] = {
    "design": "descriptive ordering only",
    "n_asos_ranked": len(RANKED_ASO_IDS),
    "n_asos_excluded": len(EXCLUDED_ASOS),
    "supervision": "none; no outcome, PSI, gel, responder or dose-response label is read, joined or fitted",
    "model_selection": "none; every alternative setting is prespecified in this module, nothing is chosen by fit quality",
    "verified_negatives": VERIFIED_NEGATIVES_AVAILABLE,
    "verified_negatives_note": (
        "no verified non-responder ASO exists for this locus in this project; "
        "fabricating negatives is prohibited and positive-unlabeled assumptions are not made here"
    ),
    "permitted_outputs": [
        "per-feature values, each with an evidence class",
        "within-feature ordinal ranks over the 3 ranked ASOs",
        "rank stability (min/max/range) across prespecified alternatives",
    ],
    "prohibited_outputs": [
        "p-values or any null-hypothesis test",
        "AUC / AUROC",
        "AUPRC / average precision",
        "FDR / q-values",
        "accuracy, precision, recall, F1, MCC",
        "calibration or Brier scores",
        "confidence intervals on an effect size",
    ],
    "why_prohibited": (
        "N=3 ranked ASOs with no verified negatives and no held-out outcome: every "
        "listed quantity would be either undefined or a restatement of the ordering"
    ),
    "july_gels": (
        "the July exon-7 gels are PREVIOUSLY SEEN data, not prospective validation; "
        "they are excluded from this dataset entirely and no gel or outcome product is read"
    ),
    "chemistry": (
        "ASO chemistry is not recorded in conference/aso_pipeline/sequences.py -> unknown; "
        "duplex energies assume unmodified RNA and are neither potency nor occupancy"
    ),
    "m5_note": (
        "M5 values are RAW Pangolin model deltas under N-masking (sequence_ablation_proxy), "
        "not delta-PSI and not occupancy"
    ),
    "paralog": "SMN1 is a distinct locus; placements are never joined across SMN1/SMN2",
}

FORBIDDEN_METRICS = frozenset({
    "accuracy", "auc", "auroc", "auprc", "average_precision", "balanced_accuracy",
    "brier", "calibration", "f1", "fdr", "fpr", "mcc", "npv", "p_value", "ppv",
    "precision", "q_value", "recall", "roc", "roc_auc", "sensitivity", "specificity",
    "tpr",
})


class ForbiddenMetricError(RuntimeError):
    """A supervised / inferential metric was requested; J9 emits ordering only."""


class ForbiddenSourceError(RuntimeError):
    """A path outside the positive source allowlist was requested."""


class AtlasCoverageError(RuntimeError):
    """The J4 atlas does not cover the requested footprint."""


def _metric_key(name: str) -> str:
    key = str(name).strip().lower().replace("-", "_").replace(" ", "_")
    for prefix in ("sklearn.metrics.", "metrics."):
        if key.startswith(prefix):
            key = key[len(prefix):]
    return key.replace("_score", "").replace("__", "_").strip("_")


def request_metric(name: str, *, negatives: Optional[Mapping[str, Any]] = None) -> Any:
    """Refuse every classification / inferential metric.

    Raises ``ForbiddenMetricError`` for AUPRC, FDR, accuracy and friends while no
    verified negatives exist, and ``NotImplementedError`` if a caller claims
    verified negatives (J9 still does not compute metrics at N=3).
    """
    key = _metric_key(name)
    if key not in FORBIDDEN_METRICS:
        raise ValueError(
            f"unknown metric {name!r}; J9 emits per-feature values and ordinal ranks only"
        )
    verified = bool(negatives) and bool(dict(negatives).get("verified"))
    if not verified:
        raise ForbiddenMetricError(
            f"{name}: refused - no verified negatives (verified_negatives="
            f"{VERIFIED_NEGATIVES_AVAILABLE}) and N={len(RANKED_ASO_IDS)} ranked ASOs; "
            "see STATISTICAL_STATUS['why_prohibited']"
        )
    raise NotImplementedError(
        f"{name}: J9 computes no classification metric even with negatives supplied; "
        f"N={len(RANKED_ASO_IDS)} supports ordering only"
    )


# --------------------------------------------------------------------------- source allowlist
def _allowed_roots() -> Tuple[Path, ...]:
    home = Path(os.path.expanduser("~")).resolve()
    scratch = home / "orcd" / "scratch" / "feng"
    roots = [
        REPO_ROOT / "conference" / "model_v2",
        REPO_ROOT / "conference" / "aso_pipeline" / "sequences.py",
        scratch / "gencode" / "v50",
        scratch / "model_v2" / "j4_sentinel_v2",
        scratch / "model_v2" / "j4_atlas_v1",
        scratch / "model_v2" / DATASET_ID,
    ]
    out = []
    for root in roots:
        try:
            out.append(root.resolve())
        except OSError:
            out.append(root)
    return tuple(out)


def check_source_allowed(path: Any) -> Path:
    """Positive allowlist on every path J9 may touch.

    Anything that is not the model_v2 tree, the locked ASO sequence module, the
    pinned GENCODE reference, a J4 atlas directory, or this dataset's own output
    directory is refused. Lab outcome products are not in the allowlist, so no
    code path in this module can reach one.
    """
    candidate = Path(os.path.expanduser(str(path)))
    try:
        resolved = candidate.resolve()
    except OSError:
        resolved = candidate.absolute()
    for root in _allowed_roots():
        if resolved == root or root in resolved.parents:
            return resolved
    raise ForbiddenSourceError(f"source not in J9 allowlist: {resolved}")


# --------------------------------------------------------------------------- small helpers
def _json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _digest(value: Any) -> str:
    return hashlib.sha256(_json(value).encode("utf-8")).hexdigest()


def _sha256(path: Any) -> str:
    check_source_allowed(path)
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def git_commit(repo_root: Any = REPO_ROOT) -> Dict[str, str]:
    try:
        head = subprocess.run(
            ["git", "-C", str(repo_root), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=20, check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", str(repo_root), "status", "--porcelain"],
            capture_output=True, text=True, timeout=30, check=True,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return {"commit": "unavailable", "worktree": "unknown"}
    return {"commit": head, "worktree": "dirty" if dirty.strip() else "clean"}


# --------------------------------------------------------------------------- ASO registry
def aso_registry() -> Dict[str, Dict[str, Any]]:
    """ASO ids/sequences from the locked ``conference/aso_pipeline/sequences.py``.

    Cross-checked against the trans_intervention registry; any divergence is a
    hard error rather than a silent preference.
    """
    check_source_allowed(Path(aso_sequences.__file__))
    out: Dict[str, Dict[str, Any]] = {}
    for aso_id, aso in sorted(aso_sequences.ASOS.items()):
        dna = aso.dna_sequence.upper()
        reference = TI_ASOS.get(aso_id)
        if reference is None or reference.sequence_5to3 != dna:
            raise ValueError(
                f"{aso_id}: sequences.py / trans_intervention registry disagree"
            )
        out[aso_id] = {
            "aso_id": aso_id,
            "sequence_5to3": dna,
            "length": len(dna),
            "target_sequence": reverse_complement(dna),
            "target_label": aso.target_label,
            "chemistry_family": reference.chemistry.family,
            "chemistry_evidence": reference.chemistry.evidence,
            "sequence_evidence": reference.sequence_evidence,
            "requires_human_confirmation": bool(reference.requires_human_confirmation),
        }
    return out


_SCRAMBLE_TOKENS = ("scr", "scramble", "scrambled", "mismatch_control", "control")


def scramble_controls(registry: Optional[Mapping[str, Any]] = None) -> Dict[str, Any]:
    """Matched scrambles are included ONLY if their sequences are in sequences.py."""
    registry = aso_registry() if registry is None else registry
    found = sorted(
        aso_id for aso_id in registry
        if any(token in aso_id.lower() for token in _SCRAMBLE_TOKENS)
    )
    return {
        "available": bool(found),
        "aso_ids": tuple(found),
        "reason": (
            "matched scramble/control sequences present in conference/aso_pipeline/sequences.py"
            if found else
            "no scramble or matched-control sequence is recorded in "
            "conference/aso_pipeline/sequences.py; none is invented here"
        ),
    }


# --------------------------------------------------------------------------- placement (F1)
def place_asos(genome: Any, registry: Optional[Mapping[str, Any]] = None) -> Dict[str, Dict[str, Any]]:
    """Exact reverse-complement placement on SMN2 exon 7 (+ SMN1 as a separate locus).

    ``genome is None`` degrades every placement to ``not_available`` instead of
    guessing coordinates. 3UP8 can never come back ``placed``: the F1 registry
    flags it ``requires_human_confirmation``.
    """
    registry = aso_registry() if registry is None else registry
    out: Dict[str, Dict[str, Any]] = {}
    for aso_id in sorted(registry):
        record: Dict[str, Any] = {
            "aso_id": aso_id,
            "locus_id": REFERENCE_LOCUS_ID,
            "status": "not_available",
            "reason": "indexed GRCh38 reference unavailable",
            "placement": None,
            "candidates": (),
            "alternates": {},
        }
        if genome is not None:
            aso = TI_ASOS[aso_id]
            result = place(aso, LOCI[REFERENCE_LOCUS_ID], genome)
            record.update(
                status=result.status,
                reason=result.reason,
                placement=result.placement,
                candidates=tuple(result.candidates),
            )
            record["alternates"] = {
                name: place(aso, LOCI[name], genome) for name in ALTERNATE_LOCUS_IDS
            }
        if aso_id in EXCLUDED_ASOS and record["status"] == "placed":
            raise AssertionError(f"{aso_id}: excluded ASO must never be uniquely placed")
        out[aso_id] = record
    return out


# --------------------------------------------------------------------------- backends
class StubBackend:
    """Deterministic ViennaRNA stand-in for unit tests - NEVER a scientific result.

    Rows produced through this backend are emitted with ``status='stub_not_a_result'``
    and ``evidence_class='not_available'``, and are excluded from every rank table.
    """

    kind = "stub"

    def __init__(self) -> None:
        self.identity = {
            "tool": "stub",
            "python": "stub",
            "note": "deterministic stand-in; not ViennaRNA; values are not results",
        }
        self.error = None
        self.calls: List[Dict[str, Any]] = []

    def call(self, method: str, **kwargs: Any) -> Dict[str, Any]:
        self.calls.append(dict(method=method, **kwargs))
        seed = int(_digest([method, sorted(kwargs.items(), key=lambda kv: kv[0])])[:12], 16)
        rng = random.Random(seed)
        if method == "duplex":
            return {"energy_kcal_mol": round(-30.0 + rng.random() * 20.0, 4),
                    "structure": "(((&)))"}
        if method == "accessibility":
            return {"p_unpaired": round(rng.random(), 6)}
        raise RuntimeError(f"stub backend: unknown method {method!r}")

    def close(self) -> None:
        return None

    def __enter__(self) -> "StubBackend":
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()


class NullBackend:
    """No RNA folding at all: M2/M3 come back ``not_available`` with a pointer.

    This is the login-node mode. It lets M0/M1/M4/M5 be built honestly without
    folding anything outside a batch job.
    """

    kind = "unavailable"
    reason = ("ViennaRNA is not run on the login node; submit "
              "conference/model_v2/jobs/j9_smn2.sbatch to fill M2/M3")

    def __init__(self) -> None:
        self.identity = {"tool": "none", "unavailable": self.reason}
        self.error = self.reason

    def call(self, method: str, **kwargs: Any) -> Dict[str, Any]:
        raise FileNotFoundError(f"{method}: {self.reason}")

    def close(self) -> None:
        return None

    def __enter__(self) -> "NullBackend":
        return self

    def __exit__(self, *args: Any) -> None:
        self.close()


def vienna_backend(python: str = j3_scoring.VIENNA_PYTHON, timeout: int = 60) -> Any:
    """Real J3 ViennaRNA worker. Only ever started from jobs/j9_smn2.sbatch."""
    backend = j3_scoring.ViennaBackend(python=python, timeout=timeout)
    setattr(backend, "kind", "viennarna")
    return backend


def _backend_kind(backend: Any) -> str:
    return str(getattr(backend, "kind", "unknown"))


# --------------------------------------------------------------------------- feature rows
_FEATURE_FIELDS = (
    "dataset_version", "aso_id", "aso_sequence_5to3", "aso_length", "locus_id",
    "site_id", "placement_status", "baseline", "feature_name", "variant_id",
    "variant_json", "value_double", "value_int", "value_string", "evidence_class",
    "evidence_label", "status", "reason", "is_hypothesis", "hypothesis_label",
    "backend_kind", "source",
)


def _row(
    *, aso: Mapping[str, Any], baseline: str, feature_name: str, evidence_class: str,
    status: str, reason: str, source: str, variant_id: str = "default",
    variant: Optional[Mapping[str, Any]] = None, value_double: Optional[float] = None,
    value_int: Optional[int] = None, value_string: Optional[str] = None,
    evidence_label: str = "", site_id: Optional[str] = None,
    placement_status: str = "not_available", is_hypothesis: bool = False,
    hypothesis_label: str = "", backend_kind: str = "none",
    locus_id: str = REFERENCE_LOCUS_ID,
) -> Dict[str, Any]:
    if evidence_class not in EVIDENCE_CLASSES:
        raise ValueError(f"bad evidence class {evidence_class!r}")
    if baseline not in BASELINES:
        raise ValueError(f"bad baseline {baseline!r}")
    has_value = value_double is not None or value_int is not None or value_string is not None
    if has_value and evidence_class == "not_available" and status == "ok":
        raise ValueError(f"{feature_name}: a present value cannot be evidence not_available and ok")
    if not has_value and status == "ok":
        raise ValueError(f"{feature_name}: status ok requires a value")
    if value_double is not None and not math.isfinite(float(value_double)):
        raise ValueError(f"{feature_name}: nonfinite value")
    return {
        "dataset_version": DATASET_VERSION,
        "aso_id": str(aso["aso_id"]),
        "aso_sequence_5to3": str(aso["sequence_5to3"]),
        "aso_length": int(aso["length"]),
        "locus_id": locus_id,
        "site_id": site_id,
        "placement_status": placement_status,
        "baseline": baseline,
        "feature_name": feature_name,
        "variant_id": variant_id,
        "variant_json": _json(dict(variant or {})),
        "value_double": None if value_double is None else float(value_double),
        "value_int": None if value_int is None else int(value_int),
        "value_string": None if value_string is None else str(value_string),
        "evidence_class": evidence_class,
        "evidence_label": evidence_label,
        "status": status,
        "reason": reason,
        "is_hypothesis": bool(is_hypothesis),
        "hypothesis_label": hypothesis_label,
        "backend_kind": backend_kind,
        "source": source,
    }


# ------------------------------------------------------------------ M0
def m0_rows(aso: Mapping[str, Any], *, order_rank: int, n_ordered: int) -> List[Dict[str, Any]]:
    """Seeded permutation / random-order reference. This is the null ordering."""
    target = aso["target_sequence"]
    bases = list(aso["sequence_5to3"])
    rng = random.Random(int(_digest([M0_SEED_SALT, DATASET_VERSION, aso["aso_id"]])[:16], 16))
    mismatches: List[int] = []
    for _ in range(M0_PERMUTATIONS):
        rng.shuffle(bases)
        permuted = "".join(bases)
        mismatches.append(sum(1 for a, b in zip(permuted, reverse_complement(target)) if a != b))
    mean = sum(mismatches) / len(mismatches)
    variance = sum((m - mean) ** 2 for m in mismatches) / (len(mismatches) - 1)
    observed = 0  # the ASO is the exact reverse complement of its own target
    at_or_below = sum(1 for m in mismatches if m <= observed)
    common = dict(
        aso=aso, baseline="M0", source="j9:seeded_permutation_null",
        evidence_label="seeded_permutation_null", status="ok", reason="computed",
        evidence_class="derived",
    )
    return [
        _row(feature_name="m0_permuted_mismatch_mean", value_double=mean, **common),
        _row(feature_name="m0_permuted_mismatch_sd", value_double=math.sqrt(variance), **common),
        _row(feature_name="m0_observed_mismatch_count", value_int=observed, **common),
        _row(feature_name="m0_observed_percentile_in_null",
             value_double=at_or_below / len(mismatches), **common),
        _row(feature_name="m0_n_permutations", value_int=M0_PERMUTATIONS, **common),
        _row(feature_name="m0_random_order_rank", value_int=order_rank,
             **dict(common, evidence_label="seeded_random_order",
                    source="j9:seeded_random_order")),
        _row(feature_name="m0_random_order_n", value_int=n_ordered, **common),
    ]


def m0_random_order(aso_ids: Sequence[str]) -> Dict[str, int]:
    ids = sorted(aso_ids)
    rng = random.Random(int(_digest([M0_SEED_SALT, DATASET_VERSION, "order", ids])[:16], 16))
    shuffled = list(ids)
    rng.shuffle(shuffled)
    return {aso_id: index + 1 for index, aso_id in enumerate(shuffled)}


# ------------------------------------------------------------------ M1
def m1_rows(aso: Mapping[str, Any], *, site_id: Optional[str], placement_status: str) -> List[Dict[str, Any]]:
    common = dict(
        aso=aso, baseline="M1", site_id=site_id, placement_status=placement_status,
        source="j9:j3_scoring.complementarity", evidence_label="j3_complementarity_channel",
    )
    try:
        values = j3_scoring.complementarity(aso["sequence_5to3"], aso["target_sequence"])
    except Exception as exc:  # pragma: no cover - equal lengths by construction
        return [_row(feature_name="m1_complementarity", evidence_class="not_available",
                     status="failed", reason=f"{type(exc).__name__}: {exc}", **common)]
    ok = dict(common, evidence_class="derived", status="ok", reason="computed")
    return [
        _row(feature_name="m1_wc_pairs", value_int=values["wc_pairs"], **ok),
        _row(feature_name="m1_gu_pairs", value_int=values["gu_pairs"], **ok),
        _row(feature_name="m1_mismatch_count", value_int=values["mismatch_count"], **ok),
        _row(feature_name="m1_mismatch_positions_1b",
             value_string=_json(values["mismatch_positions_1b"]), **ok),
        _row(feature_name="m1_aso_length", value_int=int(aso["length"]), **ok),
        _row(feature_name="m1_wc_fraction",
             value_double=values["wc_pairs"] / float(aso["length"]), **ok),
    ]


# ------------------------------------------------------------------ M2 / M3
_J3_EVIDENCE = {"duplex": "model_assumption", "accessibility": "proxy"}


def m23_rows(
    aso: Mapping[str, Any], placement: Any, *, scorer: Any, backend: Any,
    variant_id: str, variant: Mapping[str, Any],
) -> List[Dict[str, Any]]:
    """J3 duplex (M2) and accessibility (M3) for one ASO under one fold alternative."""
    kind = _backend_kind(backend)
    common = dict(
        aso=aso, variant_id=variant_id, variant=variant, backend_kind=kind,
        site_id=getattr(placement, "site_id", None),
        placement_status="placed" if placement is not None else "not_available",
    )
    specs = (
        ("M2", "m2_duplex_energy_kcal_mol", "duplex", "duplex_energy_kcal_mol", "double"),
        ("M2", "m2_duplex_structure", "duplex", "duplex_structure", "string"),
        ("M3", "m3_accessibility_p_unpaired", "accessibility", "accessibility_p_unpaired", "double"),
    )
    if placement is None:
        return [
            _row(baseline=baseline, feature_name=name, evidence_class="not_available",
                 status="not_available", reason="no confirmed placement; J3 channels not computed",
                 evidence_label=f"j3_{channel}_channel", source="j9:j3_scoring.Scorer", **common)
            for baseline, name, channel, _field, _t in specs
        ]
    row_in = {
        "site_id": placement.site_id, "aso_id": aso["aso_id"], "chrom": placement.chrom,
        "start1": int(placement.start), "end1": int(placement.end), "strand": placement.strand,
        "target_sequence": placement.target_sequence, "assembly": ASSEMBLY,
        "annotation": ANNOTATION,
    }
    try:
        scored = scorer.score(row_in)
    except Exception as exc:
        return [
            _row(baseline=baseline, feature_name=name, evidence_class="not_available",
                 status="failed", reason=f"{type(exc).__name__}: {exc}",
                 evidence_label=f"j3_{channel}_channel", source="j9:j3_scoring.Scorer", **common)
            for baseline, name, channel, _field, _t in specs
        ]
    rows: List[Dict[str, Any]] = []
    for baseline, name, channel, field, kindof in specs:
        status = scored[channel + "_status"]
        reason = scored[channel + "_reason"]
        value = scored.get(field)
        if kind == "stub":
            # A stub number is not a ViennaRNA result and must never be ranked.
            rows.append(_row(
                baseline=baseline, feature_name=name, evidence_class="not_available",
                status="stub_not_a_result",
                reason=f"stub backend ({reason}); real values require jobs/j9_smn2.sbatch",
                evidence_label=f"j3_{channel}_channel", source="j9:j3_scoring.Scorer:stub",
                **common))
            continue
        if status != "ok" or value is None:
            rows.append(_row(
                baseline=baseline, feature_name=name, evidence_class="not_available",
                status=status, reason=reason, evidence_label=f"j3_{channel}_channel",
                source="j9:j3_scoring.Scorer", **common))
            continue
        payload = {"value_string": str(value)} if kindof == "string" else {"value_double": float(value)}
        rows.append(_row(
            baseline=baseline, feature_name=name, evidence_class=_J3_EVIDENCE[channel],
            status="ok", reason=reason, evidence_label=f"j3_{channel}_channel",
            source="j9:j3_scoring.Scorer", **payload, **common))
    rows.append(_row(
        baseline="M2", feature_name="m2_model_assumption",
        value_string=j3_scoring.ASSUMPTION, evidence_class="model_assumption", status="ok",
        reason="J3 duplex assumption carried verbatim", evidence_label="j3_duplex_channel",
        source="j9:j3_scoring.ASSUMPTION", **common))
    return rows


# ------------------------------------------------------------------ M4
_M4_GEOMETRY_FIELDS = (
    "donor_dist_5p_from_last_exonic", "donor_dist_3p_from_last_exonic",
    "donor_dist_5p_from_first_intronic", "donor_dist_3p_from_first_intronic",
    "acceptor_dist_5p_from_first_exonic", "acceptor_dist_3p_from_first_exonic",
    "acceptor_dist_5p_from_last_intronic", "acceptor_dist_3p_from_last_intronic",
    "target_5p_pos", "target_3p_pos", "donor_last_exonic_base",
    "acceptor_first_exonic_base",
)


def m4_rows(
    aso: Mapping[str, Any], placement: Any, *, is_hypothesis: bool = False,
    hypothesis_label: str = "", placement_status: str = "placed",
) -> List[Dict[str, Any]]:
    common = dict(
        aso=aso, baseline="M4", source="j9:trans_intervention.splice_geometry",
        evidence_label="f1_splice_geometry", is_hypothesis=is_hypothesis,
        hypothesis_label=hypothesis_label,
    )
    if placement is None:
        return [_row(feature_name="m4_splice_geometry", evidence_class="not_available",
                     status="not_available",
                     reason="no confirmed placement; splice geometry undefined",
                     placement_status=placement_status, **common)]
    geometry = placement.geometry
    ok = dict(common, site_id=placement.site_id, placement_status=placement_status,
              evidence_class=geometry.evidence, status="ok", reason="derived from F1 placement")
    rows = [_row(feature_name="m4_" + field, value_int=int(getattr(geometry, field)), **ok)
            for field in _M4_GEOMETRY_FIELDS]
    rows.append(_row(feature_name="m4_abs_donor_dist_5p_from_last_exonic",
                     value_int=abs(int(geometry.donor_dist_5p_from_last_exonic)), **ok))
    rows.append(_row(feature_name="m4_abs_acceptor_dist_3p_from_first_exonic",
                     value_int=abs(int(geometry.acceptor_dist_3p_from_first_exonic)), **ok))
    rows.append(_row(feature_name="m4_region", value_string=geometry.region, **ok))
    rows.append(_row(feature_name="m4_region_5p", value_string=geometry.region_5p, **ok))
    rows.append(_row(feature_name="m4_region_3p", value_string=geometry.region_3p, **ok))
    rows.append(_row(feature_name="m4_target_span_nt",
                     value_int=int(placement.end - placement.start + 1), **ok))
    return rows


# ------------------------------------------------------------------ M5 (J4 atlas)
def atlas_inputs(atlas_dir: Any) -> List[Dict[str, Any]]:
    """Committed J4 ablation shards only: ``.done`` marker sha256 must match."""
    directory = check_source_allowed(atlas_dir)
    if not directory.is_dir():
        raise AtlasCoverageError(f"J4 atlas directory absent: {directory}")
    import pyarrow.parquet as pq

    identities: List[Dict[str, Any]] = []
    for done in sorted(directory.glob("*.done")):
        try:
            marker = json.loads(done.read_text())
            entry = marker["ablation"]
            path = check_source_allowed(entry["path"])
        except (OSError, ValueError, KeyError, TypeError) as exc:
            raise AtlasCoverageError(f"unusable J4 completion marker {done}: {exc}") from exc
        if not path.is_file():
            path = directory / "ablation_maps" / Path(str(entry["path"])).name
            check_source_allowed(path)
        sha = _sha256(path)
        if sha != entry.get("sha256"):
            raise AtlasCoverageError(f"J4 shard checksum mismatch: {path}")
        rows = pq.read_metadata(path).num_rows
        if rows != entry.get("rows"):
            raise AtlasCoverageError(f"J4 shard row-count mismatch: {path}")
        identities.append({"marker": str(done), "path": str(path), "sha256": sha, "rows": rows})
    if not identities:
        raise AtlasCoverageError(f"no committed J4 shard in {directory}")
    return identities


def read_atlas(atlas_dir: Any) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    identities = atlas_inputs(atlas_dir)
    import pyarrow.parquet as pq

    rows: List[Dict[str, Any]] = []
    for identity in identities:
        rows.extend(pq.read_table(identity["path"]).to_pylist())
    return rows, identities


def m5_rows(
    aso: Mapping[str, Any], placement: Any, atlas: Sequence[Mapping[str, Any]],
    *, atlas_source: str,
) -> List[Dict[str, Any]]:
    """Raw N-ablation model deltas, read ONLY where the atlas covers the footprint."""
    common = dict(aso=aso, baseline="M5", source=atlas_source,
                  evidence_label=ATLAS_EVIDENCE_LABEL)
    if placement is None:
        return [_row(feature_name="m5_sequence_ablation_proxy", evidence_class="not_available",
                     status="not_available",
                     reason="no confirmed placement; atlas footprint undefined", **common)]
    hits = [
        row for row in atlas
        if row.get("locus") == REFERENCE_LOCUS_ID
        and row.get("status") == "ok"
        and int(row.get("footprint_start", -1)) == int(placement.start)
        and int(row.get("footprint_end", -1)) == int(placement.end)
        and int(row.get("length", -1)) == int(aso["length"])
    ]
    base = dict(common, site_id=placement.site_id, placement_status="placed")
    if not hits:
        return [_row(
            feature_name="m5_sequence_ablation_proxy", evidence_class="not_available",
            status="not_available",
            reason=(f"J4 atlas does not cover footprint {placement.chrom}:"
                    f"{placement.start}-{placement.end} length={aso['length']}; "
                    "rerun with the full j4_atlas_v1 coverage"),
            **base)]
    ok = dict(base, evidence_class="proxy", status="ok",
              reason="raw Pangolin N-mask model delta; not delta-PSI")
    rows: List[Dict[str, Any]] = []
    deltas: List[float] = []
    for hit in sorted(hits, key=lambda r: int(r["tissue"])):
        tissue = int(hit["tissue"])
        delta = float(hit["raw_model_delta"])
        deltas.append(delta)
        rows.append(_row(feature_name=f"m5_raw_model_delta_tissue{tissue}", value_double=delta, **ok))
        rows.append(_row(feature_name=f"m5_reference_score_tissue{tissue}",
                         value_double=float(hit["reference_score"]), **ok))
        rows.append(_row(feature_name=f"m5_masked_score_tissue{tissue}",
                         value_double=float(hit["masked_score"]), **ok))
    rows.append(_row(feature_name="m5_raw_model_delta_mean",
                     value_double=sum(deltas) / len(deltas), **ok))
    rows.append(_row(feature_name="m5_raw_model_delta_max", value_double=max(deltas), **ok))
    rows.append(_row(feature_name="m5_n_tissue_channels", value_int=len(deltas), **ok))
    rows.append(_row(feature_name="m5_site_kind",
                     value_string=str(hits[0].get("site_kind", "")), **ok))
    rows.append(_row(feature_name="m5_interpretation",
                     value_string="N-ablation cis-feature sensitivity; not occupancy or delta PSI",
                     **dict(ok, evidence_class="model_assumption")))
    return rows


# ------------------------------------------------------------------ mismatch sensitivity
def mismatch_sensitivity_rows(
    aso: Mapping[str, Any], genome: Any, *, locus_id: str = REFERENCE_LOCUS_ID,
) -> List[Dict[str, Any]]:
    """How many candidate sites the locus window yields per prespecified max_mismatch."""
    locus = LOCI[locus_id]
    rows: List[Dict[str, Any]] = []
    common = dict(aso=aso, baseline="M1", evidence_label="j6_max_mismatch_alternative",
                  source="j9:feng.scan.seed_extend")
    if genome is None:
        return [_row(feature_name="m1_n_sites_in_locus_window", variant_id=f"mm{mm}",
                     variant={"max_mismatch": mm}, evidence_class="not_available",
                     status="not_available", reason="reference unavailable", **common)
                for mm in MAX_MISMATCH_ALTERNATIVES]
    try:
        from feng.scan.seed_extend import seed_and_extend
        reference = genome.fetch(locus.chrom, locus.search_start, locus.search_end).upper()
    except Exception as exc:
        return [_row(feature_name="m1_n_sites_in_locus_window", variant_id=f"mm{mm}",
                     variant={"max_mismatch": mm}, evidence_class="not_available",
                     status="failed", reason=f"{type(exc).__name__}: {exc}", **common)
                for mm in MAX_MISMATCH_ALTERNATIVES]
    for mm in MAX_MISMATCH_ALTERNATIVES:
        variant = {"max_mismatch": mm, "window": [locus.search_start, locus.search_end]}
        try:
            hits = seed_and_extend(aso["target_sequence"], reference, max_mismatch=mm)
            exact = sum(1 for hit in hits if hit.mismatches == 0)
            rows.append(_row(feature_name="m1_n_sites_in_locus_window", variant_id=f"mm{mm}",
                             variant=variant, value_int=len(hits), evidence_class="derived",
                             status="ok", reason="seed-and-extend over the F1 locus window",
                             **common))
            rows.append(_row(feature_name="m1_n_exact_sites_in_locus_window",
                             variant_id=f"mm{mm}", variant=variant, value_int=exact,
                             evidence_class="derived", status="ok",
                             reason="mismatch-free subset", **common))
        except Exception as exc:
            rows.append(_row(feature_name="m1_n_sites_in_locus_window", variant_id=f"mm{mm}",
                             variant=variant, evidence_class="not_available", status="failed",
                             reason=f"{type(exc).__name__}: {exc}", **common))
    return rows


# --------------------------------------------------------------------------- assembly
def _excluded_rows(aso: Mapping[str, Any], record: Mapping[str, Any]) -> List[Dict[str, Any]]:
    """One explicit not_available row per baseline, plus optional hypothesis rows."""
    aso_id = str(aso["aso_id"])
    reason = EXCLUDED_ASOS[aso_id]
    rows = [
        _row(aso=aso, baseline=baseline, feature_name=f"{baseline.lower()}_excluded",
             evidence_class="not_available", status="needs_human_confirmation",
             reason=(f"{aso_id} excluded from every baseline and rank table: {reason}; "
                     f"placement status={record['status']} ({record['reason']}), "
                     f"{len(record['candidates'])} candidate placement(s) enumerated, "
                     "none selected"),
             source="j9:excluded_aso", evidence_label="needs_human_confirmation",
             placement_status=record["status"])
        for baseline in BASELINES
    ]
    for index, candidate in enumerate(record["candidates"]):
        label = f"candidate_placement_{index}_needs_human_confirmation"
        rows.append(_row(
            aso=aso, baseline="M4", feature_name="m4_candidate_interval",
            value_string=f"{candidate.chrom}:{candidate.start}-{candidate.end}{candidate.strand}",
            evidence_class="unlabeled", status="needs_human_confirmation",
            reason="hypothesis only; placement NOT confirmed and never ranked",
            source="j9:excluded_aso_hypothesis", evidence_label="hypothesis",
            site_id=candidate.site_id, placement_status=record["status"],
            is_hypothesis=True, hypothesis_label=label))
        rows.extend(m4_rows(aso, candidate, is_hypothesis=True, hypothesis_label=label,
                            placement_status=record["status"]))
    return rows


def build_features(
    *, genome: Any, backend: Any, atlas_dir: Any = ATLAS_SENTINEL_DIR,
    cache_path: Any = None, registry: Optional[Mapping[str, Any]] = None,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Assemble the long-format feature table plus its provenance block."""
    registry = aso_registry() if registry is None else registry
    placements = place_asos(genome, registry)
    order = m0_random_order(RANKED_ASO_IDS)
    try:
        atlas, atlas_identities = read_atlas(atlas_dir)
        atlas_note = "committed"
    except (AtlasCoverageError, ImportError, OSError, ForbiddenSourceError) as exc:
        atlas, atlas_identities = [], []
        atlas_note = f"{type(exc).__name__}: {exc}"
    atlas_source = f"j4_atlas:{Path(os.path.expanduser(str(atlas_dir))).name}"

    rows: List[Dict[str, Any]] = []
    scorers: Dict[str, Any] = {}
    cache_root = Path(os.path.expanduser(str(cache_path))) if cache_path else None
    for variant_id, flank, window, span in VIENNA_ALTERNATIVES:
        scorers[variant_id] = j3_scoring.Scorer(
            backend=backend, genome=genome,
            cache_path=str((cache_root or Path(os.path.expanduser(DEFAULT_OUT_DIR)))
                           / "cache" / f"j9_cache_{variant_id}.sqlite"),
            asos={k: v["sequence_5to3"] for k, v in registry.items()},
            flank=flank, window=window, span=span)
    try:
        for aso_id in sorted(registry):
            aso = registry[aso_id]
            record = placements[aso_id]
            if aso_id in EXCLUDED_ASOS:
                rows.extend(_excluded_rows(aso, record))
                continue
            if aso_id not in RANKED_ASO_IDS:
                rows.append(_row(
                    aso=aso, baseline="M0", feature_name="m0_out_of_scope",
                    evidence_class="not_available", status="not_available",
                    reason=f"{aso_id} is not in the J9 ranked set {RANKED_ASO_IDS}",
                    source="j9:scope", placement_status=record["status"]))
                continue
            placement = record["placement"]
            site_id = getattr(placement, "site_id", None)
            rows.extend(m0_rows(aso, order_rank=order[aso_id], n_ordered=len(order)))
            rows.extend(m1_rows(aso, site_id=site_id, placement_status=record["status"]))
            rows.extend(mismatch_sensitivity_rows(aso, genome))
            for variant_id, flank, window, span in VIENNA_ALTERNATIVES:
                rows.extend(m23_rows(
                    aso, placement, scorer=scorers[variant_id], backend=backend,
                    variant_id=variant_id,
                    variant={"flank": flank, "plfold_window": window, "plfold_span": span}))
            rows.extend(m4_rows(aso, placement, placement_status=record["status"]))
            rows.extend(m5_rows(aso, placement, atlas, atlas_source=atlas_source))
            for name in sorted(record["alternates"]):
                alternate = record["alternates"][name]
                rows.append(_row(
                    aso=aso, baseline="M4", feature_name="m4_paralog_placement_status",
                    value_string=f"{alternate.status}:{alternate.reason}",
                    evidence_class="derived", status="ok",
                    reason="paralog locus reported separately; never joined to SMN2",
                    source="j9:trans_intervention.place", evidence_label="paralog_separate_locus",
                    locus_id=name, placement_status=record["status"]))
    finally:
        for scorer in scorers.values():
            scorer.__exit__(None, None, None)

    scrambles = scramble_controls(registry)
    if not scrambles["available"]:
        rows.append(_row(
            aso={"aso_id": "matched_scramble_controls", "sequence_5to3": "NA", "length": 0},
            baseline="M0", feature_name="matched_scramble_controls",
            evidence_class="not_available", status="not_available",
            reason=scrambles["reason"], source="j9:sequences.py",
            evidence_label="not_available"))

    provenance = {
        "atlas_dir": str(Path(os.path.expanduser(str(atlas_dir)))),
        "atlas_status": atlas_note,
        "atlas_inputs": atlas_identities,
        "atlas_rows_read": len(atlas),
        "backend_kind": _backend_kind(backend),
        "backend_identity": dict(getattr(backend, "identity", {})),
        "genome_available": genome is not None,
        "scorer_params": {
            variant_id: {k: v for k, v in scorers[variant_id].params.items() if k != "backend"}
            for variant_id, *_ in VIENNA_ALTERNATIVES
        },
        "scramble_controls": {k: list(v) if isinstance(v, tuple) else v
                              for k, v in scrambles.items()},
        "placements": {
            aso_id: {
                "status": record["status"], "reason": record["reason"],
                "site_id": getattr(record["placement"], "site_id", None),
                "interval": (None if record["placement"] is None else
                             f"{record['placement'].chrom}:{record['placement'].start}-"
                             f"{record['placement'].end}{record['placement'].strand}"),
                "n_candidates": len(record["candidates"]),
            }
            for aso_id, record in sorted(placements.items())
        },
    }
    rows.sort(key=lambda r: (r["baseline"], r["feature_name"], r["variant_id"],
                             r["aso_id"], r["hypothesis_label"], r["locus_id"]))
    for row in rows:
        if not row.get("evidence_class"):
            raise AssertionError(f"unlabeled value: {row}")
    return rows, provenance


# --------------------------------------------------------------------------- ranks
#: feature_name -> ordering direction. "asc" = smaller value ranks 1.
RANK_DIRECTIONS: Dict[str, str] = {
    "m0_random_order_rank": "asc",
    "m0_observed_percentile_in_null": "asc",
    "m1_wc_pairs": "desc",
    "m1_wc_fraction": "desc",
    "m1_mismatch_count": "asc",
    "m1_aso_length": "desc",
    "m1_n_sites_in_locus_window": "asc",
    "m1_n_exact_sites_in_locus_window": "asc",
    "m2_duplex_energy_kcal_mol": "asc",
    "m3_accessibility_p_unpaired": "desc",
    "m4_abs_donor_dist_5p_from_last_exonic": "asc",
    "m4_abs_acceptor_dist_3p_from_first_exonic": "asc",
    "m4_target_span_nt": "desc",
    "m5_raw_model_delta_mean": "desc",
    "m5_raw_model_delta_max": "desc",
}

RANK_METHOD = "competition_rank_ordering_only_no_significance"

_RANK_FIELDS = (
    "dataset_version", "section", "baseline", "feature_name", "variant_id", "aso_id",
    "direction", "value_double", "rank", "n_ranked", "tied_with", "min_rank", "max_rank",
    "rank_range", "n_variants", "stable", "rank_method", "evidence_class", "status",
    "reason",
)


def _rank_row(**kwargs: Any) -> Dict[str, Any]:
    row = {name: None for name in _RANK_FIELDS}
    row.update(dataset_version=DATASET_VERSION, rank_method=RANK_METHOD)
    row.update(kwargs)
    return row


def build_ranks(features: Sequence[Mapping[str, Any]]) -> List[Dict[str, Any]]:
    """Ordinal ranks per (baseline, feature, variant) plus cross-variant stability.

    Only ``status == 'ok'`` numeric values for the ranked ASOs enter; hypothesis
    rows, excluded ASOs and stub-backend rows are never ranked.
    """
    groups: Dict[Tuple[str, str, str], List[Tuple[str, float, str]]] = {}
    for row in features:
        if row["is_hypothesis"] or row["aso_id"] not in RANKED_ASO_IDS:
            continue
        if row["status"] != "ok" or row["feature_name"] not in RANK_DIRECTIONS:
            continue
        value = row["value_double"]
        if value is None:
            value = None if row["value_int"] is None else float(row["value_int"])
        if value is None:
            continue
        key = (row["baseline"], row["feature_name"], row["variant_id"])
        groups.setdefault(key, []).append((row["aso_id"], float(value), row["evidence_class"]))

    ranks: List[Dict[str, Any]] = []
    for (baseline, feature_name, variant_id), items in sorted(groups.items()):
        direction = RANK_DIRECTIONS[feature_name]
        sign = 1.0 if direction == "asc" else -1.0
        ordered = sorted(items, key=lambda item: (sign * item[1], item[0]))
        values = [sign * item[1] for item in ordered]
        for index, (aso_id, value, evidence) in enumerate(ordered):
            rank = 1 + sum(1 for other in values if other < sign * value)
            ties = sum(1 for other in values if other == sign * value) - 1
            ranks.append(_rank_row(
                section="rank", baseline=baseline, feature_name=feature_name,
                variant_id=variant_id, aso_id=aso_id, direction=direction,
                value_double=float(value), rank=rank, n_ranked=len(ordered),
                tied_with=ties, evidence_class=evidence, status="ok",
                reason="ordinal position only; no significance claim"))

    by_aso: Dict[Tuple[str, str, str], List[int]] = {}
    for row in ranks:
        by_aso.setdefault((row["baseline"], row["feature_name"], row["aso_id"]), []).append(row["rank"])
    for (baseline, feature_name, aso_id), observed in sorted(by_aso.items()):
        ranks.append(_rank_row(
            section="stability", baseline=baseline, feature_name=feature_name,
            variant_id="all", aso_id=aso_id, direction=RANK_DIRECTIONS[feature_name],
            min_rank=min(observed), max_rank=max(observed),
            rank_range=max(observed) - min(observed), n_variants=len(observed),
            stable=bool(max(observed) == min(observed)), evidence_class="derived",
            status="ok", reason="rank spread across prespecified alternatives"))
    ranks.sort(key=lambda r: (r["section"], r["baseline"], r["feature_name"],
                              r["variant_id"], r["aso_id"]))
    return ranks


# --------------------------------------------------------------------------- output
def _feature_schema():
    import pyarrow as pa
    types = {
        "aso_length": pa.int32(), "value_double": pa.float64(), "value_int": pa.int64(),
        "is_hypothesis": pa.bool_(),
    }
    return pa.schema([(name, types.get(name, pa.string())) for name in _FEATURE_FIELDS])


def _rank_schema():
    import pyarrow as pa
    ints = ("rank", "n_ranked", "tied_with", "min_rank", "max_rank", "rank_range", "n_variants")
    types = {name: pa.int32() for name in ints}
    types.update(value_double=pa.float64(), stable=pa.bool_())
    return pa.schema([(name, types.get(name, pa.string())) for name in _RANK_FIELDS])


def _write_parquet(path: Path, rows: Sequence[Mapping[str, Any]], schema: Any,
                   metadata: Mapping[bytes, bytes]) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    schema = schema.with_metadata(dict(metadata))
    table = pa.table({field.name: [row.get(field.name) for row in rows]
                      for field in schema}, schema=schema)
    pq.write_table(table, str(path), compression="none")


def _source_identities() -> List[Dict[str, Any]]:
    paths = [
        Path(__file__).resolve(),
        REPO_ROOT / "conference" / "model_v2" / "trans_intervention.py",
        REPO_ROOT / "conference" / "model_v2" / "j3_scoring.py",
        REPO_ROOT / "conference" / "model_v2" / "jobs" / "j3_vienna.py",
        REPO_ROOT / "conference" / "aso_pipeline" / "sequences.py",
        REPO_ROOT / "conference" / "model_v2" / "SCIENTIFIC_LOCK.md",
    ]
    out = []
    for path in paths:
        try:
            out.append({"path": str(path.relative_to(REPO_ROOT)), "sha256": _sha256(path)})
        except (OSError, ForbiddenSourceError, ValueError) as exc:
            out.append({"path": str(path), "sha256": f"unavailable: {type(exc).__name__}"})
    return sorted(out, key=lambda item: item["path"])


README_TEMPLATE = """# {dataset_id} - SMN2 ASO feature + ordering dataset

Descriptive feature tables for the ISS-N1 ASOs {ranked}. **Not a trained,
fitted or validated predictor.** Nothing here was selected against an outcome.

## Files
- `features.parquet` - long format, one row per (ASO, baseline, feature, alternative).
  Every row carries `evidence_class` ({classes}) plus an `evidence_label`,
  `status` and `reason`. Values are never present without a label.
- `ranks.parquet` - `section='rank'`: ordinal position per (baseline, feature,
  alternative). `section='stability'`: min/max/range of that rank across the
  prespecified alternatives.
- `manifest.json` - sha256 of every output, source file hashes, git commit,
  the statistical-status block and the prespecified alternative list.

## Baselines
M0 seeded permutation / random-order null; M1 complementarity (derived);
M2 duplex free energy (unmodified-RNA `model_assumption`, chemistry unknown);
M3 RNAplfold whole-footprint accessibility (`proxy`); M4 splice geometry to the
SMN2 exon-7 donor/acceptor (derived, F1); M5 N-ablation Pangolin delta read from
the J4 atlas and labelled `{atlas_label}` - a RAW model delta, not delta-PSI.

## Prespecified alternatives (rank stability)
Fold window/span {alternatives} (verbatim from J2 `jobs/j2_fold.py`; J9 sets the
J3 context flank equal to the window W so M2 and M3 both move with it) and
`max_mismatch` {mismatches} (J6 site-universe tolerance).

## Statistical status - read this before quoting any ordering
{status}

## Scope and exclusions
- 3UP8 is excluded: its placement is unconfirmed. It appears only as explicit
  `not_available` / `needs_human_confirmation` rows, plus candidate placements
  carried as `is_hypothesis=True` rows that enter no rank table.
- Matched scrambles: {scrambles}
- SMN1 is reported as a separate locus and is never joined to SMN2.
- No lab outcome workbook, gel table or PSI product is read; sources are
  positively allow-listed in `j9_smn2_prediction.check_source_allowed`.
- Requesting AUPRC / FDR / accuracy raises `ForbiddenMetricError`.

## Build state of this directory
- RNA-folding backend: `{backend_kind}`. `viennarna` = M2/M3 are real J3 channel
  values; `unavailable` = M2/M3 are explicit `not_available` rows awaiting
  `jobs/j9_smn2.sbatch`; `stub` = test stand-in, never a result.
- J4 atlas: `{atlas_dir}` ({atlas_status}), {atlas_rows} row(s) read.
- GRCh38 reference available: {genome_available}.

## Reproduce
```
{command}
```
"""


def write_dataset(
    out_dir: Any, features: Sequence[Mapping[str, Any]], ranks: Sequence[Mapping[str, Any]],
    provenance: Mapping[str, Any], *, command: str = "",
) -> Dict[str, Any]:
    directory = check_source_allowed(out_dir)
    directory.mkdir(parents=True, exist_ok=True)
    status_blob = _json(STATISTICAL_STATUS).encode("utf-8")
    metadata = {
        b"dataset_id": DATASET_ID.encode(),
        b"dataset_version": DATASET_VERSION.encode(),
        b"statistical_status": status_blob,
        b"coordinates": b"GRCh38; GENCODE v50; 1-based inclusive",
        b"supervision": b"none; descriptive ordering only; no verified negatives",
    }
    _write_parquet(directory / "features.parquet", features, _feature_schema(), metadata)
    _write_parquet(directory / "ranks.parquet", ranks, _rank_schema(), metadata)

    outputs = []
    for name, rows in (("features.parquet", features), ("ranks.parquet", ranks)):
        path = directory / name
        outputs.append({"name": name, "rows": len(rows), "bytes": path.stat().st_size,
                        "sha256": _sha256(path)})
    manifest = {
        "dataset_id": DATASET_ID,
        "dataset_version": DATASET_VERSION,
        "git": git_commit(),
        "assembly": ASSEMBLY,
        "annotation": ANNOTATION,
        "reference_locus": REFERENCE_LOCUS_ID,
        "alternate_loci": list(ALTERNATE_LOCUS_IDS),
        "ranked_aso_ids": list(RANKED_ASO_IDS),
        "excluded_asos": dict(EXCLUDED_ASOS),
        "baselines": list(BASELINES),
        "prespecified_alternatives": {
            "fold_window_span": [
                {"variant_id": v, "flank": f, "plfold_window": w, "plfold_span": s}
                for v, f, w, s in VIENNA_ALTERNATIVES
            ],
            "max_mismatch": list(MAX_MISMATCH_ALTERNATIVES),
        },
        "rank_directions": dict(sorted(RANK_DIRECTIONS.items())),
        "rank_method": RANK_METHOD,
        "statistical_status": STATISTICAL_STATUS,
        "evidence_classes": list(EVIDENCE_CLASSES),
        "outputs": outputs,
        "sources": _source_identities(),
        "provenance": json.loads(_json(provenance)),
        "reproduce_command": command,
        "note": "deterministic by construction: no timestamp is recorded in this manifest",
    }
    (directory / "manifest.json").write_text(_json(manifest) + "\n", encoding="utf-8")

    scrambles = provenance.get("scramble_controls", {})
    (directory / "README.md").write_text(README_TEMPLATE.format(
        dataset_id=DATASET_ID,
        ranked=", ".join(RANKED_ASO_IDS),
        classes=" | ".join(EVIDENCE_CLASSES),
        atlas_label=ATLAS_EVIDENCE_LABEL,
        alternatives=", ".join(f"({w}/{s})" for _v, _f, w, s in VIENNA_ALTERNATIVES),
        mismatches=", ".join(str(m) for m in MAX_MISMATCH_ALTERNATIVES),
        status="\n".join(f"- **{k}**: {v}" for k, v in sorted(STATISTICAL_STATUS.items())),
        scrambles=scrambles.get("reason", "unknown"),
        backend_kind=provenance.get("backend_kind", "unknown"),
        atlas_dir=provenance.get("atlas_dir", "unknown"),
        atlas_status=provenance.get("atlas_status", "unknown"),
        atlas_rows=provenance.get("atlas_rows_read", 0),
        genome_available=provenance.get("genome_available", False),
        command=command or "see conference/model_v2/jobs/j9_smn2.sbatch",
    ), encoding="utf-8")
    return manifest


def build(
    *, out_dir: Any = DEFAULT_OUT_DIR, atlas_dir: Any = ATLAS_SENTINEL_DIR,
    genome: Any = None, backend: Any = None, cache_path: Any = None,
    command: str = "",
) -> Dict[str, Any]:
    owns_backend = backend is None
    backend = StubBackend() if owns_backend else backend
    try:
        features, provenance = build_features(
            genome=genome, backend=backend, atlas_dir=atlas_dir,
            cache_path=cache_path if cache_path is not None else out_dir)
        ranks = build_ranks(features)
        return write_dataset(out_dir, features, ranks, provenance, command=command)
    finally:
        if owns_backend:
            backend.close()


# --------------------------------------------------------------------------- CLI
def _open_genome(path: Optional[str], enabled: bool) -> Any:
    if not enabled:
        return None
    try:
        return FaiFastaReader(check_source_allowed(path or DEFAULT_FASTA))
    except (GenomeAccessError, ForbiddenSourceError) as exc:
        print(f"reference unavailable: {exc}", file=sys.stderr)
        return None


def _cmd_build(args: argparse.Namespace) -> int:
    genome = _open_genome(args.fasta, not args.no_genome)
    backend: Any
    if args.stub_backend:
        backend = StubBackend()
    elif args.no_vienna:
        backend = NullBackend()
    else:
        backend = vienna_backend(args.vienna_python)
        if getattr(backend, "error", None):
            print(f"ViennaRNA worker unavailable: {backend.error}", file=sys.stderr)
    command = " ".join([Path(sys.executable).name, "-m",
                        "conference.model_v2.j9_smn2_prediction"] + list(sys.argv[1:]))
    try:
        manifest = build(out_dir=args.out_dir, atlas_dir=args.atlas_dir, genome=genome,
                         backend=backend, cache_path=args.cache or args.out_dir,
                         command=command)
    finally:
        backend.close()
    print(_json({"dataset_id": manifest["dataset_id"], "outputs": manifest["outputs"],
                 "git": manifest["git"],
                 "atlas_status": manifest["provenance"]["atlas_status"],
                 "backend_kind": manifest["provenance"]["backend_kind"]}))
    return 0


def _cmd_status(_args: argparse.Namespace) -> int:
    print(json.dumps(STATISTICAL_STATUS, indent=2, sort_keys=True))
    return 0


def _cmd_plan(_args: argparse.Namespace) -> int:
    print(json.dumps({
        "baselines": list(BASELINES),
        "ranked_aso_ids": list(RANKED_ASO_IDS),
        "excluded_asos": EXCLUDED_ASOS,
        "fold_alternatives": [{"variant_id": v, "flank": f, "window": w, "span": s}
                              for v, f, w, s in VIENNA_ALTERNATIVES],
        "max_mismatch_alternatives": list(MAX_MISMATCH_ALTERNATIVES),
        "rank_directions": RANK_DIRECTIONS,
        "allowed_source_roots": [str(p) for p in _allowed_roots()],
    }, indent=2, sort_keys=True))
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m conference.model_v2.j9_smn2_prediction")
    sub = parser.add_subparsers(dest="cmd", required=True)
    build_parser = sub.add_parser("build", help="write features/ranks/manifest/README")
    build_parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    build_parser.add_argument("--atlas-dir", default=ATLAS_SENTINEL_DIR)
    build_parser.add_argument("--fasta", default=DEFAULT_FASTA)
    build_parser.add_argument("--vienna-python", default=j3_scoring.VIENNA_PYTHON)
    build_parser.add_argument("--cache")
    build_parser.add_argument("--no-genome", action="store_true",
                              help="skip placement instead of guessing coordinates")
    build_parser.add_argument("--stub-backend", action="store_true",
                              help="no ViennaRNA; M2/M3 emitted as stub_not_a_result")
    build_parser.add_argument("--no-vienna", action="store_true",
                              help="login-node mode: M2/M3 emitted as not_available")
    build_parser.set_defaults(fn=_cmd_build)
    status_parser = sub.add_parser("status", help="print the statistical-status block")
    status_parser.set_defaults(fn=_cmd_status)
    plan_parser = sub.add_parser("plan", help="print the prespecified design")
    plan_parser.set_defaults(fn=_cmd_plan)
    args = parser.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
