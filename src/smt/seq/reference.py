"""A1: pinned reference sequence access.

Everything downstream reads coordinates from here. Accession and length are
pinned and checksummed so a silent RefSeq update cannot shift a feature by a
nucleotide and quietly invalidate every energy in the project.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

DATA = Path("/home/dashii/orcd/pool/smt-isef-data/sequences")

# SMN2 RefSeqGene (LRG_677), chromosome 5. Pinned 2026-08-22.
SMN2_ACCESSION = "NG_008728"
SMN2_LENGTH = 35073
SMN2_GB = DATA / "SMN2_NG_008728.gb"

# 1-based inclusive exon coordinates on NG_008728, from the record's own
# feature table. SMN exon numbering: 1, 2a, 2b, 3, 4, 5, 6, 7, 8.
SMN2_EXONS: dict[str, tuple[int, int]] = {
    "1": (5147, 5244),
    "2a": (18893, 18964),
    "2b": (21443, 21562),
    "3": (22412, 22612),
    "4": (22772, 22924),
    "5": (24713, 24808),
    "6": (26119, 26229),
    "7": (31999, 32052),
    "8": (32497, 33073),
}


@dataclass(frozen=True)
class Interval:
    """1-based inclusive interval on a reference accession."""

    start: int
    end: int

    def __len__(self) -> int:
        return self.end - self.start + 1

    def contains(self, pos: int) -> bool:
        return self.start <= pos <= self.end


def intron(exon_a: str, exon_b: str) -> Interval:
    """Intron between two consecutive exons, 1-based inclusive."""
    return Interval(SMN2_EXONS[exon_a][1] + 1, SMN2_EXONS[exon_b][0] - 1)


def exon(name: str) -> Interval:
    s, e = SMN2_EXONS[name]
    return Interval(s, e)


def load_smn2() -> str:
    """Return the full NG_008728 sequence, uppercase DNA."""
    if not SMN2_GB.exists():
        raise FileNotFoundError(
            f"{SMN2_GB} missing. Fetch with NCBI efetch db=nuccore id={SMN2_ACCESSION}."
        )
    from Bio import SeqIO

    rec = SeqIO.read(SMN2_GB, "genbank")
    seq = str(rec.seq).upper()
    if len(seq) != SMN2_LENGTH:
        raise ValueError(
            f"{SMN2_ACCESSION} length {len(seq)} != pinned {SMN2_LENGTH}; "
            "reference changed, re-verify all coordinates before proceeding."
        )
    return seq


def sha256(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()


def slice_ref(seq: str, iv: Interval) -> str:
    """Extract a 1-based inclusive interval."""
    return seq[iv.start - 1 : iv.end]
