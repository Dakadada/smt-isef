"""Public, serialisable types for prospective splice-effector design."""
from __future__ import annotations

from dataclasses import dataclass, fields, is_dataclass
from enum import Enum
import math
import hashlib
import json
from types import MappingProxyType
from typing import Any, Mapping

from physics_adapter import EventDefinition, TranscriptContext
from offtarget import Evidence, Transcript

PROSPECTIVE_BANNER = (
    "PROSPECTIVE PREDICTION — NOT EXPERIMENTALLY VALIDATED. "
    "Candidates require independent experimental and safety validation."
)
PROSPECTIVE_STATUS = "prospective-not-experimentally-validated"


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


def _manifest_hash(value: Any) -> str:
    payload = json.dumps(_plain(value), sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode()).hexdigest()


class ObjectiveDirection(str, Enum):
    INCREASE_INCLUSION = "inclusion"
    DECREASE_INCLUSION = "skipping"

    @property
    def sign(self) -> int:
        return 1 if self is ObjectiveDirection.INCREASE_INCLUSION else -1


@dataclass(frozen=True)
class TherapeuticObjective:
    """A signed objective bound to one named event and locus."""
    event_id: str
    locus: str
    direction: ObjectiveDirection
    sign_convention: str = "positive=increased-exon-inclusion"
    override_provenance: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.direction, ObjectiveDirection):
            raise TypeError("direction must be an ObjectiveDirection")
        if self.sign_convention != "positive=increased-exon-inclusion":
            raise ValueError("unsupported therapeutic sign convention")
        if not self.event_id.strip() or not self.locus.strip():
            raise ValueError("objective event and locus are required")
        builtins = {
            "SMN2_exon7_inclusion": ("SMN2", ObjectiveDirection.INCREASE_INCLUSION),
            "DMD_exon53_inclusion": ("DMD", ObjectiveDirection.DECREASE_INCLUSION),
        }
        expected = builtins.get(self.event_id)
        if expected is not None and (self.locus, self.direction) != expected:
            raise ValueError(f"built-in objective {self.event_id} has a fixed locus and direction")
        if expected is None and not (self.override_provenance and self.override_provenance.strip()):
            raise ValueError("custom objectives require explicit override provenance")


@dataclass(frozen=True)
class TargetDefinition:
    target_id: str
    locus: str
    objective: TherapeuticObjective
    context: TranscriptContext
    event: EventDefinition
    design_start: int
    design_end: int
    offtarget_transcripts: tuple[Transcript, ...] = ()
    offtarget_evidence: tuple[Evidence, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.objective, TherapeuticObjective):
            raise TypeError("target objective must be a TherapeuticObjective")
        if self.target_id != self.event.event_id or self.objective.event_id != self.event.event_id:
            raise ValueError("target, objective, and event identities must match")
        if self.locus != self.objective.locus:
            raise ValueError("target and objective loci must match")
        if not (0 <= self.design_start < self.design_end <= len(self.context.sequence)):
            raise ValueError("design interval lies outside transcript context")
        intervals = self.event.accessible_elements + self.event.silencers
        if any(iv.start < 0 or iv.end > len(self.context.sequence) for iv in intervals):
            raise ValueError("event interval lies outside transcript context")

    @property
    def therapeutic_polarity(self) -> str:
        return self.objective.direction.value


@dataclass(frozen=True)
class GenerationConfig:
    aso_lengths: tuple[int, ...] = (18, 20)
    grna_lengths: tuple[int, ...] = (22, 24)
    u7_lengths: tuple[int, ...] = (18, 20)
    aso_chemistries: tuple[str, ...] = ("2MOE", "OMe", "PMO", "LNA")
    stride: int = 4
    max_candidates: int = 160
    grna_shadow: int = 4

    def __post_init__(self) -> None:
        lengths = self.aso_lengths + self.grna_lengths + self.u7_lengths
        if not lengths or any(n < 8 for n in lengths):
            raise ValueError("candidate lengths must be at least 8 nt")
        if self.stride < 1 or self.max_candidates < 1 or self.grna_shadow < 0:
            raise ValueError("generation bounds must be positive")
        if not self.aso_chemistries:
            raise ValueError("at least one ASO chemistry is required")


@dataclass(frozen=True)
class GeneratedCandidate:
    candidate_id: str
    effector_class: str
    chemistry: str
    start: int
    end: int
    target_sequence: str
    effector_sequence: str


@dataclass(frozen=True)
class MagnitudeClaim:
    """A numeric claim carrying the exact evidence binding which permits it."""
    value: float
    unit: str
    evidence_provenance: str
    input_hash: str
    protocol_hash: str

    def __post_init__(self) -> None:
        if not math.isfinite(self.value):
            raise ValueError("magnitude value must be finite")
        if not self.unit.strip() or not self.evidence_provenance.strip():
            raise ValueError("magnitude unit and evidence provenance are required")
        for value in (self.input_hash, self.protocol_hash):
            if len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
                raise ValueError("magnitude evidence hashes must be SHA-256 digests")


@dataclass(frozen=True)
class ScoredCandidate:
    candidate: GeneratedCandidate
    predicted_inclusion_direction: str
    therapeutic_goal: ObjectiveDirection
    goal_alignment: str
    directional_confidence: float
    magnitude_claim: MagnitudeClaim | None
    prediction_input_hash: str
    prediction_protocol_hash: str
    abstained: bool
    uncertainty: tuple[str, ...]
    off_target_burden: float
    off_target_hit_count: int
    off_target_sign: str
    manufacturability: float
    span_robustness: float
    pareto: bool = False
    pareto_tier: int | None = None
    rank: int | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.therapeutic_goal, ObjectiveDirection):
            raise TypeError("therapeutic goal must be typed")
        if self.abstained and self.rank is not None:
            raise ValueError("abstained candidates cannot be ranked")
        if self.goal_alignment != "supports_goal" and self.rank is not None:
            raise ValueError("only goal-supporting candidates can be ranked")
        if self.magnitude_claim is not None and (
            self.magnitude_claim.input_hash != self.prediction_input_hash
            or self.magnitude_claim.protocol_hash != self.prediction_protocol_hash
        ):
            raise ValueError("magnitude claim evidence does not match prediction inputs/protocol")


@dataclass(frozen=True)
class DesignReport:
    schema_version: str
    banner: str
    prediction_status: str
    target_id: str
    locus: str
    objective: TherapeuticObjective
    candidates: tuple[ScoredCandidate, ...]
    pareto_candidate_ids: tuple[str, ...]
    provenance: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "provenance", _freeze(self.provenance))
        if self.banner != PROSPECTIVE_BANNER:
            raise ValueError("design report must carry the exact prospective honesty banner")
        if self.prediction_status != PROSPECTIVE_STATUS:
            raise ValueError("design report status must remain prospective and unvalidated")
        if self.provenance.get("prediction_status") != PROSPECTIVE_STATUS:
            raise ValueError("design report provenance must carry prospective status")
        required_hashes = ("run_manifest_sha256", "physics_adapter_identity_sha256",
                           "folding_backend_identity_sha256",
                           "offtarget_enumerator_identity_sha256")
        for key in required_hashes:
            value = self.provenance.get(key, "")
            if len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
                raise ValueError(f"design report provenance requires valid {key}")
        ids = [row.candidate.candidate_id for row in self.candidates]
        if len(ids) != len(set(ids)):
            raise ValueError("candidate IDs must be unique")
        if len(self.pareto_candidate_ids) != len(set(self.pareto_candidate_ids)):
            raise ValueError("Pareto identifiers must be unique")
        if set(self.pareto_candidate_ids) != {
            row.candidate.candidate_id for row in self.candidates if row.pareto
        }:
            raise ValueError("Pareto identifiers do not match candidate annotations")

    @property
    def therapeutic_polarity(self) -> str:
        return self.objective.direction.value

    def validate_for_rendering(self) -> None:
        self.__post_init__()
        manifest = self.provenance.get("run_manifest")
        if not isinstance(manifest, Mapping):
            raise ValueError("design report provenance requires run_manifest content")
        if _manifest_hash(manifest) != self.provenance.get("run_manifest_sha256"):
            raise ValueError("design report run manifest content does not match its digest")
        if (_plain(manifest.get("scored_candidates")) != _plain(self.candidates)
                or tuple(manifest.get("pareto_candidate_ids", ())) != self.pareto_candidate_ids):
            raise ValueError("design report candidates do not match run manifest")

    def to_dict(self) -> dict[str, Any]:
        self.validate_for_rendering()
        return _plain(self)
