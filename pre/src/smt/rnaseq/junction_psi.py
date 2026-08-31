"""Read-level cassette-exon PSI straight from the STAR BAMs, including
multimappers -- the measurement rMATS structurally cannot make for SMN2.

WHY THIS EXISTS. rMATS reported ZERO SMN events (``docs/HANDOFF.md``, GEO
stage 2). That is not a pipeline defect: SMN1 and SMN2 are ~99.9% identical,
89% of SMN2 reads and 99% of SMN1 reads multimap, and rMATS counts UNIQUE
junction reads. So the one exon this whole project is about -- SMN2 exon 7 --
was the one exon the validation could not see. This module counts the same
junctions rMATS counts, with the same normalization, but assigns reads by
what base they carry rather than by whether they mapped uniquely.

WHAT THE PARALOGS ACTUALLY OFFER. Derived by diffing the two loci's exons in
GRCh38 directly (not from memory, not from a citation):

    exon 6 (111 nt)  0 differences   -- carries NO paralog information
    exon 7 ( 54 nt)  1 difference at position 6   c.840C>T   SMN1 C / SMN2 T
    exon 8 (577 nt)  1 difference at position 236 c.*239G>A  SMN1 G / SMN2 A

The consequence is sharp and worth stating plainly: an exon 6->8 SKIPPING
junction read contains no discriminating base at all. Exon 6 is identical
between the paralogs and a 100 nt read cannot reach 236 nt into exon 8. So
the exclusion isoform is unassignable per-read *in principle*, and the
aligner's tie-break decides which locus it lands on. In this dataset every
skipping read landed at SMN1 while every inclusion read landed at SMN2 --
which would fabricate a perfect result if the two loci were counted
separately. THEY MUST BE POOLED.

WHY POOLING IS LEGITIMATE HERE. GSE307431 is GM03813, SMA patient
fibroblasts. Pooling is only valid if SMN1 is absent, so that is MEASURED,
never assumed: ``c840_census`` counts the base at c.840 across both loci in
every sample. See ``docs/HANDOFF.md`` for the observed 0 C / 1072 T.

If a future dataset carries SMN1, this pooling is INVALID and the exclusion
arm becomes genuinely unmeasurable by short reads. ``smn2_psi`` refuses
rather than guessing -- see ``SMN1PresentError``.
"""

from __future__ import annotations

import re
import subprocess
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

from smt.seq.reference import Interval

#: samtools lives in the RNA-seq env (STAR/rMATS), not the thermodynamics env.
SAMTOOLS = Path("/home/dashii/orcd/pool/rna-env/bin/samtools")

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")

#: SAM flags for alignments that are not the read's own primary placement.
#: Counting these would count one read many times, once per multimapping locus.
_NOT_PRIMARY = 0x100 | 0x800


class SMN1PresentError(RuntimeError):
    """Raised when SMN1 reads exist, so the two loci may not be pooled."""


@dataclass(frozen=True)
class Cassette:
    """One skipped-exon event, in rMATS's own terms.

    Coordinates are 1-based inclusive genomic, matching ``Interval``. rMATS
    reports ``exonStart_0base``; convert with ``+1`` before constructing one
    of these so the two agree exactly.
    """

    gene: str
    chrom: str
    upstream: Interval
    cassette: Interval
    downstream: Interval

    @property
    def inc_up_gap(self) -> tuple[int, int]:
        """The intron skipped by an upstream-junction (inclusion) read."""
        return (self.upstream.end + 1, self.cassette.start - 1)

    @property
    def inc_down_gap(self) -> tuple[int, int]:
        """The intron skipped by a downstream-junction (inclusion) read."""
        return (self.cassette.end + 1, self.downstream.start - 1)

    @property
    def skip_gap(self) -> tuple[int, int]:
        """The intron skipped by an exclusion read: straight over the cassette."""
        return (self.upstream.end + 1, self.downstream.start - 1)

    @property
    def region(self) -> str:
        return f"{self.chrom}:{self.upstream.start}-{self.downstream.end}"


@dataclass(frozen=True)
class JunctionCounts:
    """Junction-spanning read counts for one cassette in one sample."""

    inc_up: int
    inc_down: int
    skip: int
    #: Every N gap observed in the queried region, for the audit in
    #: ``all_gaps``-style checks. An exact-coordinate matcher silently
    #: returns zero if the aligner placed a junction one base off, and a
    #: silent zero here reads as "no skipping" -- i.e. as a perfect result.
    gaps: Counter

    @property
    def ijc(self) -> int:
        """rMATS IJC: inclusion junction reads, both junctions summed."""
        return self.inc_up + self.inc_down

    @property
    def sjc(self) -> int:
        return self.skip

    def __add__(self, other: "JunctionCounts") -> "JunctionCounts":
        return JunctionCounts(
            self.inc_up + other.inc_up,
            self.inc_down + other.inc_down,
            self.skip + other.skip,
            self.gaps + other.gaps,
        )


def _parse(pos: int, cigar: str) -> tuple[list[tuple[int, int, int]], list[tuple[int, int]]]:
    """Split one CIGAR into aligned blocks and N gaps.

    Blocks are ``(ref_start, ref_end, query_offset)``, all 1-based inclusive
    on the reference; gaps are ``(ref_start, ref_end)``.
    """
    ref, qpos = pos, 0
    blocks: list[tuple[int, int, int]] = []
    gaps: list[tuple[int, int]] = []
    for n_str, op in _CIGAR.findall(cigar):
        n = int(n_str)
        if op in "M=X":
            blocks.append((ref, ref + n - 1, qpos))
            ref += n
            qpos += n
        elif op == "N":
            gaps.append((ref, ref + n - 1))
            ref += n
        elif op == "D":
            ref += n
        elif op in "IS":
            qpos += n
    return blocks, gaps


def _view(bam: Path, region: str) -> Iterable[list[str]]:
    """Stream primary alignments overlapping ``region`` as split SAM fields."""
    proc = subprocess.Popen(
        [str(SAMTOOLS), "view", str(bam), region],
        stdout=subprocess.PIPE, text=True,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        fields = line.rstrip("\n").split("\t")
        if int(fields[1]) & _NOT_PRIMARY:
            continue
        yield fields
    proc.wait()
    if proc.returncode:
        raise RuntimeError(f"samtools view {bam} {region} exited {proc.returncode}")


def count_junctions(bam: Path, cas: Cassette) -> JunctionCounts:
    """Count inclusion and exclusion junction reads for one cassette.

    Multimappers are INCLUDED -- that is the entire point -- but only at
    their primary placement, so each read contributes at most once.
    """
    want = {cas.inc_up_gap: "inc_up", cas.inc_down_gap: "inc_down",
            cas.skip_gap: "skip"}
    hits: Counter = Counter()
    gaps: Counter = Counter()
    for f in _view(bam, cas.region):
        _, read_gaps = _parse(int(f[3]), f[5])
        for gap in read_gaps:
            gaps[gap] += 1
            key = want.get(gap)
            if key:
                hits[key] += 1
    return JunctionCounts(hits["inc_up"], hits["inc_down"], hits["skip"], gaps)


def read_length(bam: Path, region: str, sample: int = 2000) -> int:
    """Modal SEQ length over the first ``sample`` primary alignments."""
    lengths: Counter = Counter()
    for i, f in enumerate(_view(bam, region)):
        if i >= sample:
            break
        lengths[len(f[9])] += 1
    if not lengths:
        raise ValueError(f"no primary alignments in {region} of {bam}")
    return lengths.most_common(1)[0][0]


def form_lengths(cas: Cassette, read_len: int) -> tuple[int, int]:
    """rMATS ``IncFormLen``/``SkipFormLen`` for a skipped-exon event.

    Reproduced from the rMATS output itself rather than from the paper: for
    the POLR2H event at chr3:184,364,966-184,365,049 (84 nt) with 100 nt
    reads rMATS reports 183/99, and for a 693 nt cassette it reports 198/99.
    Both follow ``min(2*(L-1), exon_len + L - 1)`` and ``L - 1``.
    """
    skip_len = read_len - 1
    return min(2 * skip_len, len(cas.cassette) + skip_len), skip_len


def psi(counts: JunctionCounts, cas: Cassette, read_len: int) -> float | None:
    """rMATS-normalized inclusion level, or None when no junction read exists.

    ``None`` rather than 0.0 or 1.0: a sample with no reads at all is not a
    sample with no inclusion.
    """
    inc_len, skip_len = form_lengths(cas, read_len)
    if counts.ijc == 0 and counts.sjc == 0:
        return None
    inc = counts.ijc / inc_len
    skip = counts.sjc / skip_len
    return inc / (inc + skip)


def pooled_counts(bam: Path, cassettes: Sequence[Cassette]) -> JunctionCounts:
    """Sum counts across paralogous copies of the same event.

    Required for SMN: the exclusion junction carries no paralog-discriminating
    base, so which locus it lands on is an aligner tie-break, not a
    measurement. Only valid when the other paralog is absent -- see
    ``c840_census`` and ``SMN1PresentError``.
    """
    total = JunctionCounts(0, 0, 0, Counter())
    for cas in cassettes:
        total = total + count_junctions(bam, cas)
    return total


def base_census(bam: Path, chrom: str, pos: int, window: Interval) -> Counter:
    """Tally the read base at one genomic position, multimappers included."""
    out: Counter = Counter()
    for f in _view(bam, f"{chrom}:{window.start}-{window.end}"):
        blocks, _ = _parse(int(f[3]), f[5])
        for start, end, qoff in blocks:
            if start <= pos <= end:
                out[f[9][qoff + (pos - start)]] += 1
                break
    return out
