"""Immutable, honesty-gated output for engine-internal score characterization."""
from __future__ import annotations

from dataclasses import dataclass, fields, is_dataclass
from datetime import datetime, timezone
from enum import Enum
import hashlib
import inspect
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping

from design_engine import PROSPECTIVE_BANNER, PROSPECTIVE_STATUS

from .approved import derive_characterization as _APPROVED_DERIVATION
from .policy import CLAIM_BOUNDARY, P1_STATUS


APPROVED_DERIVATION_SOURCE_SHA256 = (
    "8893b290964e76cb999297c996d874bd380da19b11bbae158cd0c2410c18bb23"
)
APPROVED_DERIVATION_MODULE_SHA256 = (
    "d3e58d4ec68f76cbbdfa38769c53fe79376f8a6273eccde7ebb393b189a3f3ff"
)


def _freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _freeze(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, set):
        return frozenset(_freeze(item) for item in value)
    return value


def _plain(value: Any) -> Any:
    if is_dataclass(value):
        return {field.name: _plain(getattr(value, field.name)) for field in fields(value)}
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list, set, frozenset)):
        return [_plain(item) for item in value]
    if isinstance(value, Enum):
        return value.value
    return value


@dataclass(frozen=True)
class CandidateScore:
    candidate_id: str
    effector_class: str
    chemistry: str
    start: int
    end: int
    length: int
    class_base: float
    chemistry_delta: float
    length_penalty: float
    engine_score: float
    grna_shadow_start: int | None
    grna_shadow_end: int | None

    def __post_init__(self) -> None:
        if self.length != self.end - self.start or not (0 <= self.start < self.end):
            raise ValueError("candidate length/coordinates are inconsistent")
        if (self.grna_shadow_start is None) != (self.grna_shadow_end is None):
            raise ValueError("gRNA shadow bounds must be both present or both absent")
        if self.effector_class == "dCas13_gRNA":
            if self.grna_shadow_start is None or not (
                    0 <= self.grna_shadow_start <= self.start < self.end
                    <= self.grna_shadow_end):
                raise ValueError("dCas13 candidates require a valid engine shadow interval")
        elif self.grna_shadow_start is not None:
            raise ValueError("only dCas13 candidates may carry a gRNA shadow interval")
        recomputed = round(max(0.0, min(1.0,
            self.class_base + self.chemistry_delta - self.length_penalty)), 6)
        if self.engine_score != recomputed:
            raise ValueError("engine score does not match disclosed additive components")


@dataclass(frozen=True)
class ScoreDistribution:
    effector_class: str
    chemistry: str
    length: int
    candidate_count: int
    observed_scores: tuple[float, ...]
    minimum: float | None
    maximum: float | None
    mean: float | None

    def __post_init__(self) -> None:
        if self.candidate_count != len(self.observed_scores):
            raise ValueError("distribution count does not match observed scores")
        expected = (None, None, None) if not self.observed_scores else (
            min(self.observed_scores), max(self.observed_scores),
            round(sum(self.observed_scores) / len(self.observed_scores), 6))
        if (self.minimum, self.maximum, self.mean) != expected:
            raise ValueError("distribution summary does not match observed scores")


@dataclass(frozen=True)
class EffectorClassSummary:
    effector_class: str
    configured_lengths: tuple[int, ...]
    candidate_count: int
    observed_scores: tuple[float, ...]
    minimum: float | None
    maximum: float | None
    mean: float | None
    availability: str
    limitation: str | None

    def __post_init__(self) -> None:
        if self.candidate_count != len(self.observed_scores):
            raise ValueError("class count does not match observed scores")
        expected = (None, None, None) if not self.observed_scores else (
            min(self.observed_scores), max(self.observed_scores),
            round(sum(self.observed_scores) / len(self.observed_scores), 6))
        if (self.minimum, self.maximum, self.mean) != expected:
            raise ValueError("class summary does not match observed scores")
        if self.availability not in {"observed-in-engine-output", "not-generated"}:
            raise ValueError("invalid class availability")
        if (self.candidate_count == 0) != (self.availability == "not-generated"):
            raise ValueError("class availability conflicts with candidate count")
        if self.candidate_count == 0 and not self.limitation:
            raise ValueError("an absent configured class requires a limitation")


@dataclass(frozen=True)
class EffectorScoreCharacterization:
    schema_version: str
    banner: str
    prediction_status: str
    claim_boundary: str
    p1_status: str
    generated_at: str
    target_id: str
    locus: str
    context_length: int
    score_name: str
    score_formula: str
    class_bases: Mapping[str, float]
    aso_chemistry_deltas: Mapping[str, float]
    length_penalty_threshold: int
    length_penalty_per_nt: float
    grna_shadow_nt: int
    candidates: tuple[CandidateScore, ...]
    distributions: tuple[ScoreDistribution, ...]
    class_summaries: tuple[EffectorClassSummary, ...]
    provenance: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "class_bases", _freeze(self.class_bases))
        object.__setattr__(self, "aso_chemistry_deltas", _freeze(self.aso_chemistry_deltas))
        object.__setattr__(self, "provenance", _freeze(self.provenance))
        self.validate_for_rendering()

    def validate_for_rendering(self, _derive=_APPROVED_DERIVATION) -> None:
        if self.schema_version != "engine-effector-score-characterization/1":
            raise ValueError("unsupported effector-score characterization schema")
        if self.banner != PROSPECTIVE_BANNER or self.prediction_status != PROSPECTIVE_STATUS:
            raise ValueError("characterization must remain exactly prospective")
        if self.claim_boundary != CLAIM_BOUNDARY or self.p1_status != P1_STATUS:
            raise ValueError("engine-internal claim boundary or P1 gate changed")
        if self.target_id != "SMN2_exon7_inclusion" or self.locus != "SMN2":
            raise ValueError("characterization is fail-closed outside committed SMN2")
        try:
            timestamp = datetime.fromisoformat(self.generated_at)
        except (TypeError, ValueError) as exc:
            raise ValueError("generated_at must be canonical UTC ISO-8601") from exc
        if (timestamp.tzinfo is None
                or timestamp.utcoffset() != timezone.utc.utcoffset(timestamp)
                or timestamp.isoformat() != self.generated_at):
            raise ValueError("generated_at must be canonical UTC ISO-8601")
        identifiers = tuple(row.candidate_id for row in self.candidates)
        if len(identifiers) != len(set(identifiers)):
            raise ValueError("characterized candidate IDs must be unique")
        derivation_source_sha256 = hashlib.sha256(
            inspect.getsource(_derive).encode("utf-8")).hexdigest()
        derivation_path = Path(inspect.getsourcefile(_derive) or "")
        derivation_module_sha256 = hashlib.sha256(derivation_path.read_bytes()).hexdigest()
        if (derivation_source_sha256 != APPROVED_DERIVATION_SOURCE_SHA256
                or derivation_module_sha256 != APPROVED_DERIVATION_MODULE_SHA256):
            raise ValueError("approved characterization derivation identity changed")
        if (self.provenance.get("approved_derivation_source_sha256")
                != APPROVED_DERIVATION_SOURCE_SHA256
                or self.provenance.get("approved_derivation_module_sha256")
                != APPROVED_DERIVATION_MODULE_SHA256):
            raise ValueError("characterization lacks the exact approved derivation identity")
        # Re-run from frozen files and compare every rendered field. A caller cannot
        # make an internally consistent but false aggregate pass by editing summaries.
        expected = _derive(timestamp)
        if _plain(self) != _plain(expected):
            raise ValueError(
                "characterization is not the deterministic authenticated SMN2 engine output")

    def to_dict(self) -> dict[str, Any]:
        self.validate_for_rendering()
        return _plain(self)
