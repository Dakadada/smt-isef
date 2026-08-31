"""Backend boundary. Production folders can implement this without adapter coupling."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import hashlib
from typing import Mapping, Protocol

from .types import (CoordinateConvention, DuplexOrientation,
                    IntermolecularComplexSpec, Interval, ModelClass,
                    Probability, TranscriptContext)


class AppliedConstraintType(str, Enum):
    UNBOUND = "unbound"
    INTERMOLECULAR_DUPLEX = "intermolecular_duplex"
    FORCED_UNPAIRED = "forced_unpaired"


@dataclass(frozen=True)
class ConstraintReceipt:
    applied_constraint_type: AppliedConstraintType
    coordinate_convention: CoordinateConvention
    transcript_sequence_hash: str
    aso_sequence_hash: str | None
    strand_hash: str | None
    pairing_hash: str | None
    applied_orientation: DuplexOrientation | None
    applied_register: int | None
    chemistry_hash: str | None
    backend_version: str
    intermolecular_duplex_enforced: bool


@dataclass(frozen=True)
class FoldingResult:
    probability: Probability
    receipt: ConstraintReceipt


@dataclass(frozen=True)
class EnsembleRequest:
    context: TranscriptContext
    model_class: ModelClass
    required_accessible: tuple[Interval, ...]
    required_silencer_free: tuple[Interval, ...]
    bound_complex: IntermolecularComplexSpec | None
    unavailable_to_factors: tuple[Interval, ...]
    max_bp_span: int | None = None
    elongation_rates: tuple[tuple[float, float], ...] = ()
    feature_variant: bool = False


class Folder(Protocol):
    def productive_probability(self, request: EnsembleRequest) -> FoldingResult:
        """Return probability plus proof of the constraints actually applied."""


def _sha(value: str) -> str:
    return hashlib.sha256(value.encode()).hexdigest()


def expected_receipt(request: EnsembleRequest, backend_version: str) -> ConstraintReceipt:
    complex_spec = request.bound_complex
    return ConstraintReceipt(
        AppliedConstraintType.INTERMOLECULAR_DUPLEX if complex_spec else AppliedConstraintType.UNBOUND,
        CoordinateConvention.ZERO_BASED_HALF_OPEN,
        _sha(request.context.sequence.upper()),
        _sha(complex_spec.aso_sequence.upper()) if complex_spec else None,
        _sha(f"{complex_spec.transcript_strand_id}|{complex_spec.aso_strand_id}") if complex_spec else None,
        complex_spec.pairing_hash if complex_spec else None,
        complex_spec.orientation if complex_spec else None,
        complex_spec.register if complex_spec else None,
        _sha(complex_spec.chemistry) if complex_spec else None,
        backend_version,
        complex_spec is not None,
    )


def verified_probability(folder: Folder, request: EnsembleRequest) -> Probability:
    result = folder.productive_probability(request)
    if not isinstance(result, FoldingResult):
        raise ValueError("folder returned no structured constraint receipt")
    if not result.receipt.backend_version.strip():
        raise ValueError("folder did not report a backend version")
    expected = expected_receipt(request, result.receipt.backend_version)
    if result.receipt != expected:
        raise ValueError("folder constraint receipt does not match request; refusing probability")
    return result.probability


class TableFolder:
    """Deterministic tiny/offline backend useful for fixtures and integration tests.

    Keys are ``(model_class.value, bound)``. Requests are retained so callers can
    audit constraint semantics; this is not a thermodynamic approximation.
    """

    def __init__(self, values: Mapping[tuple[str, bool], Probability], backend_version: str = "table-folder/1"):
        self.values = dict(values)
        self.backend_version = backend_version
        self.requests: list[EnsembleRequest] = []

    def productive_probability(self, request: EnsembleRequest) -> FoldingResult:
        self.requests.append(request)
        bound = request.bound_complex is not None
        return FoldingResult(self.values[(request.model_class.value, bound)],
                             expected_receipt(request, self.backend_version))
