"""Blinded-positive recall evaluation; precision is always separate."""

from dataclasses import dataclass
from typing import Iterable, Optional, Sequence, Tuple

from .models import Candidate


@dataclass(frozen=True)
class PositiveSite:
    transcript_id: str
    start: int
    end: int
    label: str = "blinded"
    required_edits: Optional[int] = None
    required_indel_bases: Optional[int] = None

    def __post_init__(self) -> None:
        if self.start < 0 or self.end <= self.start:
            raise ValueError("positive coordinates must be 0-based half-open")


@dataclass(frozen=True)
class PositiveSetMetadata:
    provenance: str
    complete: bool
    blinded: bool


@dataclass(frozen=True)
class DeclaredBounds:
    max_edits: int
    max_indel_bases: int
    min_alignment_length: int
    max_alignment_length: int

    def __post_init__(self) -> None:
        if (self.max_edits < 0 or self.max_indel_bases < 0 or
                self.min_alignment_length < 1 or
                self.min_alignment_length > self.max_alignment_length):
            raise ValueError("invalid declared enumeration bounds")


@dataclass(frozen=True)
class BenchmarkResult:
    recall: Optional[float]
    precision: Optional[float]
    positives: int
    recovered: int
    predictions: int
    target_recall: float
    target_met: Optional[bool]
    missed: Tuple[PositiveSite, ...]
    out_of_bound_misses: Tuple[PositiveSite, ...]
    recall_status: str


def _matches(candidate: Candidate, positive: PositiveSite, tolerance: int) -> bool:
    return (candidate.transcript_id == positive.transcript_id
            and abs(candidate.start - positive.start) <= tolerance
            and abs(candidate.end - positive.end) <= tolerance)


def _out_of_bounds(site: PositiveSite, bounds: Optional[DeclaredBounds]) -> bool:
    if bounds is None:
        return False
    length = site.end - site.start
    return (length < bounds.min_alignment_length or
            length > bounds.max_alignment_length or
            (site.required_edits is not None and site.required_edits > bounds.max_edits) or
            (site.required_indel_bases is not None and
             site.required_indel_bases > bounds.max_indel_bases))


def benchmark(candidates: Sequence[Candidate], positives: Iterable[PositiveSite],
              target_recall: float = 0.90, coordinate_tolerance: int = 0,
              precision_labels_complete: bool = False,
              positive_set_metadata: Optional[PositiveSetMetadata] = None,
              declared_bounds: Optional[DeclaredBounds] = None) -> BenchmarkResult:
    if coordinate_tolerance < 0:
        raise ValueError("coordinate_tolerance must be nonnegative")
    positives = tuple(positives)
    recovered = tuple(p for p in positives if any(
        _matches(c, p, coordinate_tolerance) for c in candidates))
    recall = len(recovered) / len(positives) if positives else None
    # Precision is meaningful only when the entire prediction universe is labeled.
    true_predictions = sum(any(_matches(c, p, coordinate_tolerance) for p in positives)
                           for c in candidates)
    precision = ((true_predictions / len(candidates)) if candidates else 1.0
                 ) if precision_labels_complete else None
    valid_metadata = bool(positive_set_metadata and
                          positive_set_metadata.provenance.strip() and
                          positive_set_metadata.complete and
                          positive_set_metadata.blinded)
    status = ("undefined-empty-positive-set" if not positives else
              "validated" if valid_metadata else "provisional-incomplete-metadata")
    target_met = (recall >= target_recall) if valid_metadata and recall is not None else None
    missed = tuple(p for p in positives if p not in recovered)
    return BenchmarkResult(recall, precision, len(positives), len(recovered),
                           len(candidates), target_recall, target_met, missed,
                           tuple(p for p in missed if _out_of_bounds(p, declared_bounds)),
                           status)
