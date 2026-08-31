"""A2: construct definition and reference<->construct coordinate mapping.

A *construct* is a contiguous window of the pinned NG_008728 reference, with
zero or more internal deletion intervals removed. That covers both cases we
care about:

  * the genomic locus (no deletions), and
  * a truncated-intron minigene (introns shortened by deleting an internal
    block), which is a *different molecule* thermodynamically and must be
    folded as such.

Coordinate conventions
----------------------
* Reference space: 1-based inclusive on NG_008728 (same as ``smt.seq.reference``).
* Construct space: 1-based inclusive on the construct's own sequence, which is
  the concatenation of the retained reference segments in 5'->3' order.
* Every ``Interval`` stored in ``Construct.features`` is in REFERENCE space.
  Use :meth:`Construct.feature` to get it in construct space.

Sequences are stored as uppercase DNA (T). Use :meth:`Construct.rna` /
:meth:`Construct.rna_sub` when handing sequence to ViennaRNA or IntaRNA.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Optional

from smt.seq.reference import Interval, exon, intron, load_smn2, slice_ref

# ---------------------------------------------------------------------------
# Literature-anchored features, all verified against NG_008728 in this repo's
# tests. Intron 7 position "+N" == reference coordinate 32052 + N, because
# exon 7 ends at 32052 (smt.seq.reference.SMN2_EXONS["7"]).
# ---------------------------------------------------------------------------

_I7_START = 32052  # exon 7 last base; intron 7 +N == _I7_START + N

# ISS-N1: intron 7 +10..+24, 15 nt. Singh, Androphy & Singh 2006 (Mol Cell Biol
# 26:1333). Sequence check: CCAGCATTATGAAAG.
ISS_N1 = Interval(_I7_START + 10, _I7_START + 24)

# GC-rich core targeted by the 8-mer ASO "3UP8": intron 7 +7..+14, 8 nt.
# Sequence check: CTGCCAGC.
GC_CORE = Interval(_I7_START + 7, _I7_START + 14)

# Exon 7 5' splice site region: exon 7 last 3 nt + intron 7 first 6 nt = 9 nt.
# Sequence check: GGA|GTAAGT (canonical U1 donor).
SS5_EXON7 = Interval(32050, _I7_START + 6)

# SMN1/SMN2 discriminating transition c.840C>T: exon 7 position 6.
# exon 7 starts at 31999, so position 6 == 32004. SMN2 carries T.
C840T_POS = 32004
C840T_SMN2_BASE = "T"

# Default compact-window margins (nt). Chosen to match the span of the classic
# pCI exon 6-8 minigenes: enough intron 6 to carry the exon 7 3' splice site
# context, enough exon 8 to carry the exon 7 3'-side competitor.
DEFAULT_INTRON6_TAIL = 200
DEFAULT_EXON8_HEAD = 100

# Global partition functions are O(n^3) in time and O(n^2) in memory. Above
# this length prefer the compact window or a local (RNAplfold) method.
GLOBAL_PF_SOFT_LIMIT = 2000


@dataclass(frozen=True)
class Construct:
    """A named sequence window on the reference, minus optional deletions.

    Attributes
    ----------
    name:
        Human-readable identifier, used in filenames and IntaRNA sequence ids.
    window:
        Reference interval spanned by the construct (1-based inclusive).
    deletions:
        Reference intervals removed from ``window`` (e.g. a truncated intron).
        Must be disjoint, sorted, and contained in ``window``.
    sequence:
        Construct-space uppercase DNA.
    features:
        Named intervals in REFERENCE space.
    """

    name: str
    window: Interval
    deletions: tuple[Interval, ...]
    sequence: str
    features: dict[str, Interval] = field(default_factory=dict)

    # -- basic accessors ---------------------------------------------------

    def __len__(self) -> int:
        return len(self.sequence)

    def rna(self) -> str:
        """Whole construct as RNA (U for T)."""
        return self.sequence.replace("T", "U")

    def sub(self, iv: Interval) -> str:
        """Construct-space 1-based inclusive slice, DNA alphabet."""
        return self.sequence[iv.start - 1 : iv.end]

    def rna_sub(self, iv: Interval) -> str:
        """Construct-space 1-based inclusive slice, RNA alphabet."""
        return self.sub(iv).replace("T", "U")

    # -- coordinate mapping ------------------------------------------------

    def _segments(self) -> list[tuple[Interval, int]]:
        """Retained reference segments paired with their construct start offset.

        Returns a list of ``(reference_interval, construct_start)`` with
        construct_start 1-based.
        """
        segs: list[tuple[Interval, int]] = []
        cursor = self.window.start
        offset = 1
        for d in self.deletions:
            if d.start > cursor:
                iv = Interval(cursor, d.start - 1)
                segs.append((iv, offset))
                offset += len(iv)
            cursor = d.end + 1
        if cursor <= self.window.end:
            iv = Interval(cursor, self.window.end)
            segs.append((iv, offset))
        return segs

    def to_construct(self, ref_pos: int) -> Optional[int]:
        """Reference 1-based -> construct 1-based, or None if deleted/outside."""
        for iv, off in self._segments():
            if iv.contains(ref_pos):
                return off + (ref_pos - iv.start)
        return None

    def to_reference(self, con_pos: int) -> int:
        """Construct 1-based -> reference 1-based."""
        for iv, off in self._segments():
            if off <= con_pos <= off + len(iv) - 1:
                return iv.start + (con_pos - off)
        raise ValueError(f"construct position {con_pos} outside {self.name} (len {len(self)})")

    def map_interval(self, iv: Interval) -> Interval:
        """Reference interval -> construct interval.

        Raises if the interval is not fully retained and contiguous in the
        construct; a feature split by a deletion is not a feature any more and
        silently returning the fragments would corrupt every energy downstream.
        """
        s = self.to_construct(iv.start)
        e = self.to_construct(iv.end)
        if s is None or e is None:
            raise ValueError(
                f"{iv} is not fully contained in construct {self.name} "
                f"(window {self.window}, deletions {self.deletions})"
            )
        if e - s != len(iv) - 1:
            raise ValueError(f"{iv} is split by a deletion in construct {self.name}")
        return Interval(s, e)

    def feature(self, name: str) -> Interval:
        """Named feature in CONSTRUCT space."""
        if name not in self.features:
            raise KeyError(f"{name!r} not in {self.name}; have {sorted(self.features)}")
        return self.map_interval(self.features[name])

    def ref_feature(self, name: str) -> Interval:
        """Named feature in REFERENCE space."""
        return self.features[name]

    def with_features(self, **extra: Interval) -> "Construct":
        """Copy carrying additional reference-space features."""
        merged = dict(self.features)
        merged.update(extra)
        return Construct(self.name, self.window, self.deletions, self.sequence, merged)


# ---------------------------------------------------------------------------
# builders
# ---------------------------------------------------------------------------


def _build(
    name: str,
    window: Interval,
    deletions: Iterable[Interval] = (),
    seq: Optional[str] = None,
) -> Construct:
    ref = seq if seq is not None else load_smn2()
    dels = tuple(sorted(deletions, key=lambda d: d.start))
    for d in dels:
        if d.start < window.start or d.end > window.end:
            raise ValueError(f"deletion {d} outside window {window}")
    for a, b in zip(dels, dels[1:]):
        if a.end >= b.start:
            raise ValueError(f"overlapping deletions {a} and {b}")

    pieces: list[str] = []
    cursor = window.start
    for d in dels:
        if d.start > cursor:
            pieces.append(slice_ref(ref, Interval(cursor, d.start - 1)))
        cursor = d.end + 1
    if cursor <= window.end:
        pieces.append(slice_ref(ref, Interval(cursor, window.end)))

    features = {
        "ISS_N1": ISS_N1,
        "GC_CORE": GC_CORE,
        "SS5_EXON7": SS5_EXON7,
        "EXON7": exon("7"),
        "INTRON7": intron("7", "8"),
    }
    c = Construct(name, window, dels, "".join(pieces), features)
    # Drop features the deletions destroyed rather than carrying a lie.
    keep = {}
    for k, v in features.items():
        try:
            c.map_interval(v)
        except ValueError:
            continue
        keep[k] = v
    return Construct(name, window, dels, c.sequence, keep)


def genomic_window(seq: Optional[str] = None) -> Construct:
    """Exon 6 start through exon 8 end, full introns 6 and 7 (~7.0 kb).

    Intron 6 alone is 5769 nt, so this window is far above
    ``GLOBAL_PF_SOFT_LIMIT``: use it for local (RNAplfold) accessibility only,
    not for the global constrained partition function.
    """
    return _build("SMN2_ex6_ex8_genomic", Interval(exon("6").start, exon("8").end), (), seq)


def compact_window(
    intron6_tail: int = DEFAULT_INTRON6_TAIL,
    exon8_head: int = DEFAULT_EXON8_HEAD,
    seq: Optional[str] = None,
) -> Construct:
    """Standard minigene-like span: intron 6 3' tail -> exon 8 5' head.

    Default 200 nt of intron 6 + exon 7 (54) + intron 7 (444) + 100 nt of
    exon 8 = 798 nt, which a global partition function handles comfortably.
    This is the default analysis window for :mod:`smt.thermo.constrained`.
    """
    i6 = intron("6", "7")
    e8 = exon("8")
    win = Interval(i6.end - intron6_tail + 1, e8.start + exon8_head - 1)
    return _build(f"SMN2_compact_i6t{intron6_tail}_e8h{exon8_head}", win, (), seq)


def truncated_intron_minigene(
    deletions: Iterable[Interval],
    intron6_tail: int = DEFAULT_INTRON6_TAIL,
    exon8_head: int = DEFAULT_EXON8_HEAD,
    name: str = "SMN2_minigene",
    seq: Optional[str] = None,
) -> Construct:
    """Compact window with explicit internal deletions (intron truncation).

    Deletions are REFERENCE intervals. Use this once the wet lab confirms the
    actual plasmid junctions (see docs/lab-request.md Ask 1); do not guess them.
    """
    i6 = intron("6", "7")
    e8 = exon("8")
    win = Interval(i6.end - intron6_tail + 1, e8.start + exon8_head - 1)
    return _build(name, win, deletions, seq)


# ---------------------------------------------------------------------------
# derived structural features
#
# TODO(SHAPE): TSL2 and ISTL1 coordinates below are DERIVED FROM FOLDING, not
# from literature. They must be cross-checked against Singh lab SHAPE /
# in-cell probing data for SMN2 intron 7 before any of them is treated as a
# real structural element. Until then they are candidates, and are named
# *_candidate in every report.
# ---------------------------------------------------------------------------


def pair_table(structure: str) -> list[int]:
    """Dot-bracket -> 1-based pair table; ``t[i] == j`` (0 if unpaired)."""
    t = [0] * (len(structure) + 1)
    stack: list[int] = []
    for i, ch in enumerate(structure, start=1):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if not stack:
                raise ValueError("unbalanced dot-bracket")
            j = stack.pop()
            t[i] = j
            t[j] = i
    if stack:
        raise ValueError("unbalanced dot-bracket")
    return t


def helices(structure: str, min_len: int = 2) -> list[tuple[Interval, Interval]]:
    """Contiguous stacked helices as (5' arm, 3' arm), construct-local coords.

    A helix is a maximal run of pairs (i,j), (i+1,j-1), ... Only helices of at
    least ``min_len`` base pairs are returned; lonely pairs are noise.
    """
    t = pair_table(structure)
    out: list[tuple[Interval, Interval]] = []
    i = 1
    n = len(structure)
    while i <= n:
        j = t[i]
        if j > i:
            k = 0
            while i + k <= n and t[i + k] == j - k and t[i + k] > i + k:
                k += 1
            if k >= min_len:
                out.append((Interval(i, i + k - 1), Interval(j - k + 1, j)))
            i += k
        else:
            i += 1
    return out


def derive_tsl2_candidate(construct: Construct, scan_nt: int = 60) -> Interval:
    """Candidate terminal stem-loop 2 (TSL2) immediately 3' of the exon 7 5'ss.

    Folds intron 7 +1..+``scan_nt`` in isolation and returns the span of the
    outermost helix, i.e. the hairpin that would sequester the donor site.
    Returned in REFERENCE space.

    TODO(SHAPE): candidate only -- cross-check against Singh lab probing data.
    """
    import RNA

    win = Interval(_I7_START + 1, _I7_START + scan_nt)
    seq = construct.rna_sub(construct.map_interval(win))
    structure, _ = RNA.fold(seq)
    hx = helices(structure)
    if not hx:
        raise ValueError("no helix found in TSL2 scan window")
    five, three = min(hx, key=lambda h: h[0].start)  # outermost = leftmost 5' arm
    return Interval(win.start + five.start - 1, win.start + three.end - 1)


def derive_istl1_candidate(construct: Construct) -> tuple[Interval, Interval]:
    """Candidate ISTL1 long-distance stem: two REFERENCE-space arms.

    Folds the whole of intron 7 and returns the arms of the helix with the
    largest sequence separation that also involves the intron 7 5' region
    (first 40 nt, where the internal stem partner of the terminal stem-loop is
    reported to sit). This reproduces the *shape* of the published ISTL1
    interaction without importing coordinates we cannot verify.

    TODO(SHAPE): candidate only -- cross-check against Singh lab probing data.
    """
    import RNA

    i7 = intron("7", "8")
    local = construct.map_interval(i7)
    seq = construct.rna_sub(local)
    structure, _ = RNA.fold(seq)
    hx = helices(structure, min_len=3)
    # "5' region" = intron 7 +1..+40; span = distance between the two arms.
    cand = [h for h in hx if h[0].start <= 40]
    if not cand:
        cand = hx
    if not cand:
        raise ValueError("no helix found in intron 7 fold")
    five, three = max(cand, key=lambda h: h[1].start - h[0].end)
    return (
        Interval(i7.start + five.start - 1, i7.start + five.end - 1),
        Interval(i7.start + three.start - 1, i7.start + three.end - 1),
    )


def with_derived_features(construct: Construct) -> Construct:
    """Attach TSL2 / ISTL1 candidate features (folding-derived, unverified)."""
    tsl2 = derive_tsl2_candidate(construct)
    arm1, arm2 = derive_istl1_candidate(construct)
    return construct.with_features(TSL2=tsl2, ISTL1_ARM1=arm1, ISTL1_ARM2=arm2)
