"""Typed, serialisable inputs and outputs for effector physics."""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Mapping, Sequence
import hashlib


@dataclass(frozen=True, order=True)
class Interval:
    """Zero-based, half-open transcript interval."""

    start: int
    end: int

    def __post_init__(self) -> None:
        if self.start < 0 or self.end <= self.start:
            raise ValueError("interval must be non-empty and non-negative")

    def overlaps(self, other: "Interval") -> bool:
        return self.start < other.end and other.start < self.end


class CoordinateConvention(str, Enum):
    ZERO_BASED_HALF_OPEN = "zero_based_half_open"


class DuplexOrientation(str, Enum):
    ANTIPARALLEL = "antiparallel"


@dataclass(frozen=True, order=True)
class IntermolecularPair:
    transcript_position: int
    aso_position: int

    def __post_init__(self) -> None:
        if self.transcript_position < 0 or self.aso_position < 0:
            raise ValueError("pair positions must be non-negative")


@dataclass(frozen=True)
class IntermolecularComplexSpec:
    """Complete, auditable definition of a transcript-ASO bound ensemble."""

    transcript_coordinates: Interval
    coordinate_convention: CoordinateConvention
    aso_sequence: str
    transcript_strand_id: str
    aso_strand_id: str
    orientation: DuplexOrientation
    register: int
    chemistry: str
    pairing_constraints: tuple[IntermolecularPair, ...]

    def __post_init__(self) -> None:
        sequence = self.aso_sequence.upper()
        if not sequence or any(base not in "ACGUTN" for base in sequence):
            raise ValueError("ASO sequence must contain only A/C/G/U/T/N")
        if not self.transcript_strand_id or not self.aso_strand_id:
            raise ValueError("both strand IDs are required")
        if self.transcript_strand_id == self.aso_strand_id:
            raise ValueError("transcript and ASO strand IDs must differ")
        if not self.chemistry.strip():
            raise ValueError("ASO chemistry is required")
        if len(self.pairing_constraints) != len(sequence):
            raise ValueError("every ASO base must have an explicit intermolecular pair")
        transcript_positions = {p.transcript_position for p in self.pairing_constraints}
        aso_positions = {p.aso_position for p in self.pairing_constraints}
        if transcript_positions != set(range(self.transcript_coordinates.start,
                                             self.transcript_coordinates.end)):
            raise ValueError("pairing constraints must cover the target exactly")
        if aso_positions != set(range(len(sequence))):
            raise ValueError("pairing constraints must cover the ASO exactly")
        expected = tuple(range(len(sequence) - 1, -1, -1))
        ordered = tuple(p.aso_position for p in sorted(
            self.pairing_constraints, key=lambda pair: pair.transcript_position))
        if self.orientation is DuplexOrientation.ANTIPARALLEL and ordered != expected:
            raise ValueError("pairing constraints do not match antiparallel orientation/register")
        if self.register != self.transcript_coordinates.start:
            raise ValueError("register must identify the first transcript coordinate")

    @property
    def pairing_hash(self) -> str:
        text = ";".join(f"{p.transcript_position}:{p.aso_position}"
                        for p in self.pairing_constraints)
        return hashlib.sha256(text.encode()).hexdigest()


@dataclass(frozen=True)
class TranscriptContext:
    transcript_id: str
    sequence: str
    probing: Mapping[int, float] = field(default_factory=dict)
    eclip_repressive_occupancy: Mapping[str, float] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.sequence or any(x not in "ACGUTNacgutn" for x in self.sequence):
            raise ValueError("sequence must contain only A/C/G/U/T/N")


@dataclass(frozen=True)
class EventDefinition:
    event_id: str
    u1_contacts: tuple[Interval, ...]
    branchpoint_ppt: tuple[Interval, ...] = ()
    silencers: tuple[Interval, ...] = ()

    @property
    def accessible_elements(self) -> tuple[Interval, ...]:
        return self.u1_contacts + self.branchpoint_ppt


@dataclass(frozen=True)
class EffectorSpec:
    effector_id: str
    target: Interval
    bound_complex: IntermolecularComplexSpec
    bulky_rnp: bool = False
    shadow: tuple[Interval, ...] = ()

    def __post_init__(self) -> None:
        if self.bound_complex.transcript_coordinates != self.target:
            raise ValueError("bound complex coordinates must equal the effector target")
        if self.bound_complex.aso_strand_id != self.effector_id:
            raise ValueError("ASO strand ID must equal effector ID")
        if not self.bulky_rnp and self.shadow:
            raise ValueError("only a bulky RNP may declare an excluded-volume shadow")


class ModelClass(str, Enum):
    EQUILIBRIUM_FULL_TRANSCRIPT = "equilibrium_full_transcript"
    COTRANSCRIPTIONAL = "cotranscriptional"
    ECLIP_PROBING_ABLATION = "eclip_probing_ablation"


@dataclass(frozen=True)
class Probability:
    estimate: float
    lower: float
    upper: float

    def __post_init__(self) -> None:
        if not (0 <= self.lower <= self.estimate <= self.upper <= 1):
            raise ValueError("invalid probability interval")


@dataclass(frozen=True)
class ProductiveState:
    name: str
    required_accessible: tuple[Interval, ...]
    required_silencer_free: tuple[Interval, ...]
    target_paired_to_aso: tuple[Interval, ...]
    unavailable_to_factors: tuple[Interval, ...]


@dataclass(frozen=True)
class ClassResult:
    model_class: ModelClass
    role: str
    unbound_state: ProductiveState
    bound_state: ProductiveState
    p_u: Probability
    p_b: Probability
    theta_probabilities: Mapping[float, Probability]
    theta_deltas: Mapping[float, float]
    theta_directions: Mapping[float, int]
    robust_direction: int | None


@dataclass(frozen=True)
class Prediction:
    classes: Mapping[ModelClass, ClassResult]
    primary_direction: int | None
    directional_concordance: bool
    concordance_directions: Mapping[ModelClass, int | None]
    magnitude_claim_allowed: bool
    theta_grid: tuple[float, ...]
    warnings: tuple[str, ...]


class OccupancySource(str, Enum):
    MATCHED_EXPERIMENTAL = "matched_experimental"


@dataclass(frozen=True)
class MagnitudeEvidence:
    provenance: str
    matching_criterion: str
    occupancy_source: OccupancySource
    input_hash: str
    protocol_hash: str

    MATCHING_CRITERION = "same-transcript-effector-event-and-frozen-protocol"

    def __post_init__(self) -> None:
        if not self.provenance.strip() or not self.matching_criterion.strip():
            raise ValueError("evidence provenance and matching criterion are required")
        if self.matching_criterion != self.MATCHING_CRITERION:
            raise ValueError("evidence does not use the frozen matching criterion")
        hashes = (self.input_hash, self.protocol_hash)
        if any(len(value) != 64 or any(char not in "0123456789abcdef" for char in value)
               for value in hashes):
            raise ValueError("evidence hashes must be SHA-256 hex digests")
