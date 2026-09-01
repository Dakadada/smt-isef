"""Typed inputs and outputs for reporter-context compression.

Coordinates in this module are reporter-construct coordinates.  They are not
genomic coordinates and are deliberately carried with a distinct domain label.
"""
from __future__ import annotations

from dataclasses import dataclass, field, fields, is_dataclass
from enum import Enum
import hashlib
import json
from typing import Any, Mapping

from design_engine import TargetDefinition, TherapeuticObjective
from physics_adapter import EventDefinition, Interval, TranscriptContext


class ContextDomain(str, Enum):
    MINIGENE_REPORTER = "minigene_reporter"


CLAIM_BOUNDARY = (
    "Compute-only representation of a supplied minigene reporter context; "
    "not an endogenous-genome representation, biological result, experimental "
    "validation, or therapeutic claim."
)
COORDINATE_CONVENTION = "reporter-construct-zero-based-half-open"


def _plain(value: Any) -> Any:
    if is_dataclass(value):
        return {item.name: _plain(getattr(value, item.name)) for item in fields(value)}
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(item) for item in value]
    if isinstance(value, Enum):
        return value.value
    return value


def _hash(value: Any) -> str:
    payload = json.dumps(_plain(value), sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _valid_sha256(value: str) -> bool:
    return len(value) == 64 and all(character in "0123456789abcdef" for character in value)


def _event_intervals(event: EventDefinition) -> tuple[Interval, ...]:
    return event.u1_contacts + event.branchpoint_ppt + event.silencers


@dataclass(frozen=True)
class ReporterSourceContext:
    """Caller-supplied reporter sequence with explicit, non-genomic provenance."""

    reporter_id: str
    locus: str
    sequence: str
    sequence_sha256: str
    sequence_provenance: str
    objective: TherapeuticObjective
    event: EventDefinition
    design_interval: Interval
    domain_label: ContextDomain = ContextDomain.MINIGENE_REPORTER
    coordinate_convention: str = COORDINATE_CONVENTION

    def __post_init__(self) -> None:
        if self.domain_label is not ContextDomain.MINIGENE_REPORTER:
            raise ValueError("reporter source domain must remain minigene_reporter")
        if self.coordinate_convention != COORDINATE_CONVENTION:
            raise ValueError("unsupported reporter coordinate convention")
        if not self.reporter_id.strip() or not self.locus.strip() or not self.sequence_provenance.strip():
            raise ValueError("reporter identity, locus, and sequence provenance are required")
        if not self.sequence or self.sequence != self.sequence.upper() or any(
                base not in "ACGTN" for base in self.sequence):
            raise ValueError("reporter sequence must be uppercase DNA alphabet A/C/G/T/N")
        if not _valid_sha256(self.sequence_sha256) or _hash_sequence(self.sequence) != self.sequence_sha256:
            raise ValueError("reporter source sequence digest mismatch")
        if not isinstance(self.objective, TherapeuticObjective):
            raise TypeError("reporter objective must be typed")
        if self.objective.event_id != self.event.event_id or self.objective.locus != self.locus:
            raise ValueError("reporter event, objective, and locus identities must match")
        intervals = _event_intervals(self.event) + (self.design_interval,)
        if any(interval.end > len(self.sequence) for interval in intervals):
            raise ValueError("reporter interval lies outside the supplied source sequence")


@dataclass(frozen=True, order=True)
class RetainedBlock:
    label: str
    source_interval: Interval

    def __post_init__(self) -> None:
        if not self.label.strip():
            raise ValueError("retained block label is required")


@dataclass(frozen=True)
class ReporterTranscriptContext(TranscriptContext):
    """Engine-compatible transcript context that cannot lose its reporter domain."""

    domain_label: ContextDomain = ContextDomain.MINIGENE_REPORTER
    source_coordinate_convention: str = COORDINATE_CONVENTION
    source_sequence_sha256: str = ""
    source_sequence_provenance: str = ""
    retained_blocks: tuple[tuple[str, int, int], ...] = ()
    source_to_context_mapping: tuple[tuple[int, int], ...] = ()
    artificial_junction_offsets: tuple[int, ...] = ()
    compression_manifest_sha256: str = ""
    claim_boundary: str = CLAIM_BOUNDARY

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.domain_label is not ContextDomain.MINIGENE_REPORTER:
            raise ValueError("engine context domain must remain minigene_reporter")
        if self.source_coordinate_convention != COORDINATE_CONVENTION:
            raise ValueError("unsupported reporter coordinate convention")
        if not _valid_sha256(self.source_sequence_sha256):
            raise ValueError("reporter engine context requires the source sequence digest")
        if not self.source_sequence_provenance.strip():
            raise ValueError("reporter engine context requires source sequence provenance")
        if not self.retained_blocks or not self.source_to_context_mapping:
            raise ValueError("reporter engine context requires compression blocks and mapping")
        if not _valid_sha256(self.compression_manifest_sha256):
            raise ValueError("reporter engine context requires a compression manifest digest")
        if self.claim_boundary != CLAIM_BOUNDARY:
            raise ValueError("reporter engine context claim boundary changed")
        if tuple(offset for _, offset in self.source_to_context_mapping) != tuple(range(len(self.sequence))):
            raise ValueError("reporter engine mapping must exactly cover compressed coordinates")
        if any(offset <= 0 or offset >= len(self.sequence) for offset in self.artificial_junction_offsets):
            raise ValueError("artificial junction offset lies outside the compressed context")
        if tuple(sorted(set(self.artificial_junction_offsets))) != self.artificial_junction_offsets:
            raise ValueError("artificial junction offsets must be unique and ordered")


@dataclass(frozen=True)
class CompressedReporterContext:
    """Validated bridge from reporter source coordinates to an engine target."""

    schema_version: str
    domain_label: ContextDomain
    claim_boundary: str
    source: ReporterSourceContext
    retained_blocks: tuple[RetainedBlock, ...]
    source_to_context_mapping: tuple[tuple[int, int], ...]
    target: TargetDefinition
    input_manifest_sha256: str

    def __post_init__(self) -> None:
        self.validate_for_use()

    def validate_for_use(self) -> None:
        from .compression import _derive

        if self.schema_version != "minigene-context-compression/1":
            raise ValueError("unsupported compressed reporter schema")
        if self.domain_label is not ContextDomain.MINIGENE_REPORTER:
            raise ValueError("compressed context domain must remain minigene_reporter")
        if self.claim_boundary != CLAIM_BOUNDARY:
            raise ValueError("reporter claim boundary changed")
        expected = _derive(self.source, self.retained_blocks)
        if self.source_to_context_mapping != expected.mapping:
            raise ValueError("reporter coordinate mapping is not the deterministic compression output")
        if self.target != expected.target:
            raise ValueError("engine target is not the deterministic compression output")
        if self.input_manifest_sha256 != expected.manifest_sha256:
            raise ValueError("reporter compression input manifest digest mismatch")
        if not isinstance(self.target.context, ReporterTranscriptContext):
            raise TypeError("compressed target requires a reporter-specific engine context")
        if self.target.context.domain_label is not self.domain_label:
            raise ValueError("reporter domain label was not retained into the engine context")

    def to_dict(self) -> dict[str, Any]:
        self.validate_for_use()
        return _plain(self)


def _hash_sequence(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("ascii")).hexdigest()
