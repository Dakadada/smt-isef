"""Locked ASO definitions and sequence-derived target-window placement."""

from __future__ import annotations

from dataclasses import dataclass


ISS_N1_CONTEXT = "GUAAGUCUGCCAGCAUUAUGAAAGUGAAUCUU"


def to_rna(sequence: str) -> str:
    normalized = sequence.upper().replace("T", "U")
    invalid = set(normalized).difference("ACGU")
    if invalid:
        raise ValueError(f"invalid RNA bases: {sorted(invalid)}")
    return normalized


def reverse_complement_rna(sequence: str) -> str:
    return to_rna(sequence).translate(str.maketrans("ACGU", "UGCA"))[::-1]


@dataclass(frozen=True, slots=True)
class ASO:
    aso_id: str
    dna_sequence: str
    target_label: str

    @property
    def rna_sequence(self) -> str:
        return to_rna(self.dna_sequence)

    @property
    def length(self) -> int:
        return len(self.dna_sequence)


@dataclass(frozen=True, slots=True)
class TargetWindow:
    sequence: str
    start_1based: int
    end_1based: int
    requires_human_confirmation: bool


ASOS = {
    "F18": ASO("F18", "TCACTTTCATAATGCTGG", "ISS-N1 core"),
    "F14": ASO("F14", "TTTCATAATGCTGG", "ISS-N1 3-prime anchored"),
    "F10": ASO("F10", "ATAATGCTGG", "ISS-N1 3-prime anchored"),
    "3UP8": ASO("3UP8", "GCTGGCAG", "three-nt-upstream candidate window"),
}


def derive_target(aso_rna: str, context_rna: str = ISS_N1_CONTEXT) -> TargetWindow:
    """Find the unique exact reverse-complement target in the supplied context."""
    context = to_rna(context_rna)
    target = reverse_complement_rna(aso_rna)
    starts = []
    offset = context.find(target)
    while offset >= 0:
        starts.append(offset)
        offset = context.find(target, offset + 1)
    if len(starts) != 1:
        raise ValueError(
            f"expected one exact target for {to_rna(aso_rna)}, found {len(starts)}"
        )
    start = starts[0] + 1
    return TargetWindow(
        sequence=target,
        start_1based=start,
        end_1based=start + len(target) - 1,
        requires_human_confirmation=len(target) == 8,
    )

