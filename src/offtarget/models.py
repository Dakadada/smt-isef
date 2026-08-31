"""Typed inputs and outputs for the enumerator."""

from dataclasses import dataclass, field
from typing import Optional, Tuple


@dataclass(frozen=True)
class Transcript:
    """An expressed transcript's complete pre-mRNA, including introns and exons."""

    transcript_id: str
    gene_id: str
    sequence: str
    expression_tpm: float
    nuclear_fraction: float


@dataclass(frozen=True)
class Evidence:
    transcript_id: str
    start: int
    end: int
    source: str = "eCLIP"
    rbp: str = "unknown"
    occupancy: float = 1.0
    effect_sign: Optional[str] = None
    perturbation_supported: bool = False


@dataclass(frozen=True)
class EnumeratorConfig:
    chemistry: str = "2MOE"
    max_edits: int = 8
    max_indel_bases: int = 2
    min_alignment_length: Optional[int] = None
    max_alignment_length: Optional[int] = None
    min_expression_tpm: float = 0.0
    min_score: float = 0.0
    expression_half_saturation_tpm: float = 10.0
    localization_floor: float = 0.1
    rbp_weight: float = 0.12
    unknown_chemistry_penalty: float = 0.75
    allow_unknown_chemistry_heuristic: bool = False

    def __post_init__(self) -> None:
        if self.max_edits < 0 or self.max_indel_bases < 0:
            raise ValueError("enumeration bounds must be nonnegative")
        if ((self.min_alignment_length is None) !=
                (self.max_alignment_length is None)):
            raise ValueError("both alignment length bounds must be declared together")
        if (self.min_alignment_length is not None and
                (self.min_alignment_length < 1 or
                 self.min_alignment_length > self.max_alignment_length)):
            raise ValueError("invalid alignment length window")


@dataclass(frozen=True)
class Candidate:
    rank: int
    transcript_id: str
    gene_id: str
    start: int
    end: int
    target_sequence: str
    mismatches: int
    insertions: int
    deletions: int
    edit_distance: int
    cigar: str
    chemistry: str
    chemistry_parameters_known: bool
    estimated_delta_g: Optional[float]
    energy_estimate_kind: str
    expression_tpm: float
    nuclear_fraction: float
    expression_localization_weight: float
    rbp_evidence: Tuple[Evidence, ...] = field(default_factory=tuple)
    effect_sign: str = "unknown"
    recall_score: float = 0.0
    confidence: float = 0.0
    passes_expression_annotation: bool = True
    passes_score_annotation: bool = True
