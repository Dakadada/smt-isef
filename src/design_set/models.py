"""Immutable types for the SMN2 prospective predicted-superiority set."""
from __future__ import annotations

from dataclasses import dataclass, fields, is_dataclass
from datetime import datetime, timezone
from enum import Enum
from types import MappingProxyType
from typing import Any, Mapping

from design_engine import PROSPECTIVE_BANNER, PROSPECTIVE_STATUS


CLAIM_BOUNDARY = (
    "Claim boundary: predicted directional effect and reduced off-target overlap are "
    "prospective model outputs; predicted directional/specificity ranking only; "
    "NOT experimentally validated; NOT efficacy/potency."
)

EMPTY_OFF_TARGET_DISCLOSURE = (
    "Off-target burden was NOT evaluated because the committed off-target transcript "
    "corpus is empty (n=0). The off_target_burden axis is therefore non-discriminating, "
    "and the empty predicted-superior set must NOT be read as 'no better molecule exists.'"
)


def _freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _freeze(item) for key, item in value.items()})
    if isinstance(value, (tuple, list)):
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
class ReferencePlacement:
    label: str
    effector_class: str
    chemistry: str
    length: int
    coordinate_status: str
    start: int | None
    end: int | None
    candidate_id: str | None
    limitation: str | None

    def __post_init__(self) -> None:
        statuses = {"resolved-full-18mer", "reference-position-within-committed-window"}
        if self.coordinate_status not in statuses:
            raise ValueError("unsupported reference coordinate status")
        resolved = self.coordinate_status == "resolved-full-18mer"
        if resolved != all(value is not None for value in (self.start, self.end, self.candidate_id)):
            raise ValueError("reference coordinates and candidate identity must resolve together")
        if resolved and (self.end - self.start != self.length or self.limitation is not None):
            raise ValueError("resolved reference placement is internally inconsistent")
        if not resolved and (self.start is not None or self.end is not None
                             or self.candidate_id is not None or not self.limitation):
            raise ValueError("limited reference placement must not fabricate coordinates")


@dataclass(frozen=True)
class DesignSetCandidate:
    banner: str
    prediction_status: str
    candidate_id: str
    effector_class: str
    chemistry: str
    start: int
    end: int
    target_sequence: str
    effector_sequence: str
    directional_confidence: float
    off_target_burden: float
    span_robustness: float
    rank: int | None
    pareto: bool
    predicted_superior_to_reference: bool
    is_reference: bool

    def __post_init__(self) -> None:
        if self.banner != PROSPECTIVE_BANNER or self.prediction_status != PROSPECTIVE_STATUS:
            raise ValueError("every design-set candidate must remain exactly prospective")
        if not (0 <= self.start < self.end):
            raise ValueError("candidate coordinates must be zero-based half-open")
        if self.rank is None and self.pareto:
            raise ValueError("an unranked candidate cannot be on the Pareto front")
        if self.pareto != (self.rank == 1):
            raise ValueError("Pareto annotation must match the three-axis rank")
        if self.is_reference and self.predicted_superior_to_reference:
            raise ValueError("the reference cannot strictly dominate itself")


def _prediction_axis_dominates(a: DesignSetCandidate, b: DesignSetCandidate) -> bool:
    left = (a.directional_confidence, -a.off_target_burden, a.span_robustness)
    right = (b.directional_confidence, -b.off_target_burden, b.span_robustness)
    return all(x >= y for x, y in zip(left, right)) and any(
        x > y for x, y in zip(left, right))


@dataclass(frozen=True)
class ProspectiveDesignSet:
    schema_version: str
    banner: str
    prediction_status: str
    claim_boundary: str
    generated_at: str
    target_id: str
    locus: str
    coordinate_convention: str
    reference: ReferencePlacement
    candidates: tuple[DesignSetCandidate, ...]
    superior_candidate_ids: tuple[str, ...]
    provenance: Mapping[str, Any]

    def __post_init__(self) -> None:
        object.__setattr__(self, "provenance", _freeze(self.provenance))
        self.validate_for_rendering()

    def validate_for_rendering(self) -> None:
        if self.schema_version != "prospective-design-set/1":
            raise ValueError("unsupported prospective design-set schema")
        if self.banner != PROSPECTIVE_BANNER or self.prediction_status != PROSPECTIVE_STATUS:
            raise ValueError("design set must carry the exact prospective banner/status")
        if self.claim_boundary != CLAIM_BOUNDARY:
            raise ValueError("design set claim boundary changed")
        if self.target_id != "SMN2_exon7_inclusion" or self.locus != "SMN2":
            raise ValueError("prospective design set is fail-closed outside committed SMN2")
        if self.coordinate_convention != "zero-based half-open within committed context":
            raise ValueError("unsupported candidate coordinate convention")
        try:
            parsed_time = datetime.fromisoformat(self.generated_at)
        except (TypeError, ValueError) as exc:
            raise ValueError("design-set generated_at must be canonical UTC ISO-8601") from exc
        if (parsed_time.tzinfo is None
                or parsed_time.utcoffset() != timezone.utc.utcoffset(parsed_time)
                or parsed_time.isoformat() != self.generated_at):
            raise ValueError("design-set generated_at must be canonical UTC ISO-8601")
        if self.provenance.get("prediction_status") != PROSPECTIVE_STATUS:
            raise ValueError("design-set provenance must carry prospective status")
        corpus_size = self.provenance.get("off_target_corpus_size")
        evaluated = self.provenance.get("off_target_burden_evaluated")
        if (not isinstance(corpus_size, int) or isinstance(corpus_size, bool)
                or corpus_size < 0 or evaluated is not (corpus_size > 0)):
            raise ValueError("off-target evaluation provenance is inconsistent")
        if corpus_size == 0 and (
                self.provenance.get("off_target_burden_discriminating") is not False
                or self.provenance.get("off_target_disclosure")
                    != EMPTY_OFF_TARGET_DISCLOSURE):
            raise ValueError("empty off-target corpus requires the full non-evaluation disclosure")
        identifiers = tuple(row.candidate_id for row in self.candidates)
        if len(identifiers) != len(set(identifiers)):
            raise ValueError("design-set candidate IDs must be unique")
        references = tuple(row.candidate_id for row in self.candidates if row.is_reference)
        expected = () if self.reference.candidate_id is None else (self.reference.candidate_id,)
        if references != expected:
            raise ValueError("reference placement does not match candidate annotations")
        if self.reference.candidate_id is None and self.superior_candidate_ids:
            raise ValueError("superiority cannot be asserted without a scored reference")
        recomputed: tuple[str, ...] = ()
        if self.reference.candidate_id is not None:
            reference = next(row for row in self.candidates
                             if row.candidate_id == self.reference.candidate_id)
            recomputed = tuple(sorted(
                row.candidate_id for row in self.candidates
                if _prediction_axis_dominates(row, reference)))
        annotated = tuple(sorted(row.candidate_id for row in self.candidates
                                 if row.predicted_superior_to_reference))
        if tuple(sorted(self.superior_candidate_ids)) != recomputed or annotated != recomputed:
            raise ValueError(
                "superior set/annotations do not match recomputed three-axis weak "
                "dominance against reference")

    def to_dict(self) -> dict[str, Any]:
        self.validate_for_rendering()
        return _plain(self)
