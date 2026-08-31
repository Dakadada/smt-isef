"""Checked coordinate conversion for the off-target coordinate contract."""

from dataclasses import dataclass
from typing import Sequence, Tuple


@dataclass(frozen=True)
class GenomicExon:
    """A 0-based, half-open genomic exon interval."""

    start: int
    end: int


def genomic_to_transcript_interval(start: int, end: int, strand: str,
                                   exons: Sequence[GenomicExon]) -> Tuple[int, int]:
    """Map a genomic interval to strand-oriented, spliced transcript indices.

    The interval must be wholly exonic and contiguous in transcript coordinates.
    Exons may be supplied in any order.
    """
    if start < 0 or end <= start or strand not in {"+", "-"}:
        raise ValueError("invalid genomic interval or strand")
    ordered = sorted(exons, key=lambda exon: exon.start,
                     reverse=(strand == "-"))
    if not ordered or any(exon.start < 0 or exon.end <= exon.start for exon in ordered):
        raise ValueError("invalid exon map")
    offsets = []
    cursor = 0
    for exon in ordered:
        offsets.append((exon, cursor))
        cursor += exon.end - exon.start

    covered = 0
    mapped = []
    for exon, offset in offsets:
        left, right = max(start, exon.start), min(end, exon.end)
        if left >= right:
            continue
        covered += right - left
        if strand == "+":
            mapped.append((offset + left - exon.start, offset + right - exon.start))
        else:
            mapped.append((offset + exon.end - right, offset + exon.end - left))
    if covered != end - start:
        raise ValueError("interval contains intronic or unmapped bases")
    tx_start, tx_end = min(x[0] for x in mapped), max(x[1] for x in mapped)
    if tx_end - tx_start != covered:
        raise ValueError("genomic interval is not contiguous in this transcript")
    return tx_start, tx_end
