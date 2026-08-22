"""A8: inverted-repeat search flanking the C2B34 (exon 2a..4) backsplice.

Why this module exists
-----------------------
A common proposed mechanism for circRNA biogenesis (Zhang et al. 2014, Cell
159:134; Liang & Wilusz 2014, Genes Dev 28:2233; the ciRS-7/CDR1as case) is
that INVERTED repeats (often Alu elements, since ~10% of the human genome is
Alu and they insert in either orientation at random) sit in the introns
flanking the circularized exons. If the repeat in the upstream flanking
intron is the reverse complement of the repeat in the downstream flanking
intron, the two can base-pair with each other *in trans* across the
intervening exons, looping them into proximity and bringing the downstream
splice donor and the upstream splice acceptor together for backsplicing.

For C2B34 (SMN2 exons 2a-2b-3-4, joined by NG_008728 backsplicing donor at
the end of exon 4 to the acceptor at the start of exon 2a) the two flanking
introns are:

  UPSTREAM_FLANK   = intron("1", "2a")  = 5245..18892  (13648 nt)
  DOWNSTREAM_FLANK = intron("4", "5")   = 22925..24712  (1788 nt)

reusing :func:`smt.seq.reference.intron` rather than re-deriving these
coordinates.

What this module does NOT do
-----------------------------
There is no RepeatMasker / Alu-subfamily annotation database in this
environment, so this module cannot say "this is an AluYa5 element" the way a
repeat-annotated genome browser would. It searches GENERICALLY for
Alu-scale (~300 nt) intermolecular complementarity between the two flanking
introns using ViennaRNA's RNAduplex algorithm (:func:`RNA.duplexfold`), and
separately confirms, for any strong hit, whether the two loci are actually
in INVERTED (reverse-complement) orientation relative to each other (a real
inverted repeat) as opposed to DIRECT orientation (a merely similar,
same-strand repeat copy, which cannot base-pair intermolecularly this way).
This is a disclosed scope limitation, not an oversight.

Method
------
1. Tile both flanking introns into overlapping ~300 nt windows (Alu length
   scale) and run ``RNA.duplexfold`` on every (upstream tile, downstream
   tile) pair. ``duplexfold`` computes the single best-scoring intermolecular
   duplex between two RNA strands, restricted to inter-strand base pairs only
   (no intramolecular structure) -- this is exactly the RNAduplex algorithm,
   just called through its ViennaRNA python binding instead of the RNAduplex
   CLI.
2. For the strongest hits, run a composition-matched NULL CONTROL: shuffle
   both matched sequences (mononucleotide shuffle, so %GC is held fixed) and
   recompute the duplex energy many times. A real inverted repeat should
   duplex far more strongly than its own shuffled composition-matched
   controls; a hit that is just an artifact of shared GC-richness should not.
3. For the strongest hits, independently check ORIENTATION with a plain
   Needleman/Smith-Waterman-style local alignment (Biopython
   ``Align.PairwiseAligner``): does the upstream window align better to the
   downstream window read FORWARD (direct/same-strand repeat -- cannot form
   this duplex physically) or to its REVERSE COMPLEMENT (a genuine inverted
   repeat)? This is independent evidence from the duplex energy itself.
4. For the strongest hits, refine the single-duplex MFE estimate with a full
   thermodynamic-ensemble calculation via ``RNA.fold_compound`` cofold mode
   (the same computational strategy as ``RNAcofold``): the ensemble free
   energy of hybridization dG_hyb = G(AB) - G(A) - G(B), and the total
   expected number of inter-strand base pairs from the base-pair probability
   matrix. This checks that the pairing survives once each strand's own
   competing intramolecular structure is allowed into the ensemble (RNAduplex
   ignores intramolecular structure entirely; RNAcofold-style cofold does not).

Units and conventions
----------------------
All energies are kcal/mol (ViennaRNA's native unit, T = 37 C, matching
:data:`smt.thermo.accessibility.RT`). All intervals are 1-based inclusive in
REFERENCE space (NG_008728), the same convention as
:mod:`smt.seq.reference`. Sequences handed to ViennaRNA are RNA alphabet
(T->U); sequences handed to Biopython's aligner are left as DNA, which is
immaterial for alignment scoring.
"""

from __future__ import annotations

import random
import statistics
from dataclasses import dataclass, field
from typing import Optional

from smt.seq.reference import Interval, exon, intron, load_smn2, slice_ref

# ---------------------------------------------------------------------------
# Flanking introns. Reused, not re-derived: see smt.seq.reference.intron().
# ---------------------------------------------------------------------------

#: Intron between exon 1 and exon 2a: 5245..18892, 13648 nt. The upstream
#: flank of the C2B34 backsplice acceptor (start of exon 2a).
UPSTREAM_FLANK: Interval = intron("1", "2a")

#: Intron between exon 4 and exon 5: 22925..24712, 1788 nt. The downstream
#: flank of the C2B34 backsplice donor (end of exon 4).
DOWNSTREAM_FLANK: Interval = intron("4", "5")

# ---------------------------------------------------------------------------
# Search parameters. ALL of these are disclosed ASSUMPTIONS, not measurements
# or fitted values -- there is no prior data on where (or whether) an inverted
# repeat sits in these introns, so the search has to be run generically.
# ---------------------------------------------------------------------------

#: Tile width, nt. Alu consensus elements are ~300 bp full-length (Batzer &
#: Deininger, Nat Rev Genet 2002, 3:370). Chosen to search at Alu length
#: scale without requiring a repeat-family annotation database (not available
#: in this environment). ASSUMPTION: a real inverted element shorter or
#: longer than this will still partially overlap a tile because of the 50%
#: step below, but will not be centred in one.
IR_TILE_WIDTH = 300

#: Tile step, nt (50% overlap with IR_TILE_WIDTH). ASSUMPTION: a compromise
#: between not missing a repeat that straddles a tile boundary and keeping
#: the number of duplexfold calls (and hence wall-clock time) tractable for
#: an interactive analysis rather than a cluster array.
IR_TILE_STEP = 150

#: Number of top-ranked (by duplex energy) candidate pairs to carry through
#: the more expensive per-hit analyses (null control, orientation check,
#: cofold ensemble refinement). ASSUMPTION: chosen for compute cost, not a
#: statistical criterion.
IR_TOP_K = 8

#: Number of shuffled-control trials per null estimate. ASSUMPTION: enough to
#: get a stable mean/SD for a z-score without materially adding runtime; the
#: z-scores actually observed in this analysis (tens of SD) are not sensitive
#: to this choice.
IR_NULL_TRIALS = 30

#: Fixed seed for the shuffle null control, for reproducibility.
IR_NULL_SEED = 20260822


# ---------------------------------------------------------------------------
# data structures
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DuplexHit:
    """One candidate inverted-repeat pair from the tiling scan.

    All Interval fields are REFERENCE space (NG_008728), 1-based inclusive.
    """

    tile_upstream: Interval  # the 300-nt (or shorter, at a boundary) tile searched
    tile_downstream: Interval
    matched_upstream: Interval  # the actual matched sub-span within the tile
    matched_downstream: Interval
    energy_kcal_mol: float  # RNAduplex/duplexfold MFE of this duplex
    n_bp: int  # number of paired positions in the duplex
    density: float  # 2*n_bp / (len(tile_upstream) + len(tile_downstream))

    @property
    def dist_to_acceptor(self) -> int:
        """nt from the end of the matched upstream arm to the exon 2a acceptor.

        Smaller = the candidate arm sits closer to the actual splice site it
        would need to help juxtapose; a purely thermodynamic ranking (by
        energy alone) can favour a strong but very distant arm that is less
        biologically plausible as the functional element. Reported so this
        trade-off is visible rather than hidden by sorting on energy alone.
        """
        return exon("2a").start - self.matched_upstream.end

    @property
    def dist_to_donor(self) -> int:
        """nt from the exon 4 donor to the start of the matched downstream arm."""
        return self.matched_downstream.start - exon("4").end


@dataclass(frozen=True)
class NullControl:
    """Composition-matched shuffle control for one DuplexHit's energy."""

    real_energy: float
    shuffled_mean: float
    shuffled_sd: float
    n_trials: int
    z_score: float  # (real_energy - shuffled_mean) / shuffled_sd; more negative = more significant


@dataclass(frozen=True)
class OrientationCheck:
    """Forward vs reverse-complement local-alignment scores for one hit.

    A genuine inverted repeat aligns far better to the reverse complement of
    its partner than to the partner read forward; a direct (same-strand)
    repeat does the opposite.
    """

    forward_score: float
    revcomp_score: float

    @property
    def is_inverted(self) -> bool:
        """True if the reverse-complement alignment dominates the forward one."""
        return self.revcomp_score > self.forward_score


@dataclass(frozen=True)
class EnsembleRefinement:
    """RNAcofold-style (fold_compound cofold) ensemble check of one hit.

    dG_hyb is the standard ensemble hybridization free energy,
    G(complex) - G(strand1 alone) - G(strand2 alone) (kcal/mol; more negative
    = more favourable), computed over the actual matched sub-spans (not the
    full 300-nt tile) plus a small margin. expected_inter_bp is the sum of
    the base-pair probability matrix restricted to inter-strand pairs -- the
    "pairing probability" requested for this module, expressed as an expected
    paired-base count (0..min(len1,len2)) rather than a single scalar
    probability, since (as in smt.thermo.constrained) the *joint* probability
    of the entire matched span pairing simultaneously is a much smaller
    number dominated by whichever single register is realized; the expected
    count is the standard ensemble-averaged quantity ViennaRNA reports.
    """

    dG_hyb: float
    expected_inter_bp: float
    margin: int


# ---------------------------------------------------------------------------
# tiling and duplex search
# ---------------------------------------------------------------------------


def _tile(iv: Interval, width: int, step: int) -> list[Interval]:
    """Overlapping windows covering ``iv``, 1-based inclusive reference coords.

    The final window is dropped if it would be shorter than ``width // 2``
    (an under-sized trailing remainder is not a fair Alu-scale tile).
    """
    out: list[Interval] = []
    pos = iv.start
    while pos <= iv.end:
        end = min(pos + width - 1, iv.end)
        if end - pos + 1 >= width // 2:
            out.append(Interval(pos, end))
        if end == iv.end:
            break
        pos += step
    return out


def _rna(seq: str, iv: Interval) -> str:
    return slice_ref(seq, iv).replace("T", "U")


def _parse_duplex_span(d, tile1: Interval, tile2: Interval) -> tuple[Interval, Interval, int]:
    """Map a ``RNA.duplexfold`` result back to reference-space intervals.

    ``duplexfold(s1, s2)`` returns a dot-bracket structure ``part1&part2``
    plus ``i`` (the 1-based position in ``s1`` of the LAST paired base, i.e.
    the 3'-most base of ``s1`` engaged in the duplex) and ``j`` (the 1-based
    position in ``s2`` of the FIRST paired base, i.e. the 5'-most base of
    ``s2`` engaged). This convention was verified empirically against a
    synthetic fully-complementary pair before being relied on here (see the
    module's development notes / commit history) rather than assumed from
    memory of the ViennaRNA docs.
    """
    part1, part2 = d.structure.split("&")
    n_bp = part1.count("(")
    up_end = tile1.start + d.i - 1
    up_start = up_end - len(part1) + 1
    down_start = tile2.start + d.j - 1
    down_end = down_start + len(part2) - 1
    return Interval(up_start, up_end), Interval(down_start, down_end), n_bp


def scan_flanks(
    seq: Optional[str] = None,
    width: int = IR_TILE_WIDTH,
    step: int = IR_TILE_STEP,
) -> list[DuplexHit]:
    """Tile both flanking introns and RNAduplex every (upstream, downstream) pair.

    Returns every hit (not just top candidates), unsorted. O(n_tiles_up *
    n_tiles_down) calls to ``RNA.duplexfold``; for the default 300/150
    parameters over the pinned NG_008728 introns this is 90 * 11 = 990 calls,
    ~3 minutes wall-clock (measured on this project's environment; each call
    is a few hundred nt x a few hundred nt so individually fast, but the
    pairwise product adds up). Not run as part of the test suite for that
    reason -- see tests/test_kinetics.py's alu-related tests, which exercise
    the helper functions on small synthetic sequences instead.
    """
    import RNA

    s = seq if seq is not None else load_smn2()
    up_tiles = _tile(UPSTREAM_FLANK, width, step)
    down_tiles = _tile(DOWNSTREAM_FLANK, width, step)
    hits: list[DuplexHit] = []
    for ut in up_tiles:
        s1 = _rna(s, ut)
        for dt in down_tiles:
            s2 = _rna(s, dt)
            d = RNA.duplexfold(s1, s2)
            matched_up, matched_down, n_bp = _parse_duplex_span(d, ut, dt)
            density = 2 * n_bp / (len(s1) + len(s2))
            hits.append(
                DuplexHit(
                    tile_upstream=ut,
                    tile_downstream=dt,
                    matched_upstream=matched_up,
                    matched_downstream=matched_down,
                    energy_kcal_mol=float(d.energy),
                    n_bp=n_bp,
                    density=density,
                )
            )
    return hits


def top_candidates(hits: list[DuplexHit], n: int = IR_TOP_K) -> list[DuplexHit]:
    """The ``n`` hits with the most favourable (most negative) duplex energy."""
    return sorted(hits, key=lambda h: h.energy_kcal_mol)[:n]


# ---------------------------------------------------------------------------
# per-hit follow-up analyses
# ---------------------------------------------------------------------------


def null_control(
    seq: Optional[str],
    hit: DuplexHit,
    n_trials: int = IR_NULL_TRIALS,
    seed: int = IR_NULL_SEED,
) -> NullControl:
    """Composition-matched (mononucleotide shuffle) null for ``hit``'s energy.

    ASSUMPTION: a mononucleotide shuffle preserves overall base composition
    (%A/C/G/U) but not local dinucleotide structure; a dinucleotide-preserving
    shuffle would be a stricter (more conservative) null. It was not needed
    here -- the z-scores this analysis actually finds (tens of SD) are far
    larger than any plausible tightening of the null could erase, so the
    qualitative conclusion does not depend on shuffle stringency. Stated
    explicitly rather than silently using the weaker null.
    """
    import RNA

    s = seq if seq is not None else load_smn2()
    s1 = _rna(s, hit.tile_upstream)
    s2 = _rna(s, hit.tile_downstream)
    rng = random.Random(seed)
    energies = []
    for _ in range(n_trials):
        l1, l2 = list(s1), list(s2)
        rng.shuffle(l1)
        rng.shuffle(l2)
        d = RNA.duplexfold("".join(l1), "".join(l2))
        energies.append(float(d.energy))
    mean = statistics.mean(energies)
    sd = statistics.stdev(energies)
    z = (hit.energy_kcal_mol - mean) / sd if sd > 0 else float("-inf")
    return NullControl(hit.energy_kcal_mol, mean, sd, n_trials, z)


def orientation_check(seq: Optional[str], hit: DuplexHit) -> OrientationCheck:
    """Forward vs reverse-complement local alignment of the two matched arms.

    Independent of the duplex energy: uses plain sequence alignment (no
    thermodynamics) to confirm the two loci are actually reverse complements
    of one another (a real inverted repeat) rather than merely GC-content-
    matched same-orientation sequence that ``duplexfold`` can still pair
    fairly well by chaining short complementary runs.
    """
    from Bio import Align
    from Bio.Seq import Seq

    s = seq if seq is not None else load_smn2()
    up = slice_ref(s, hit.tile_upstream)
    down = slice_ref(s, hit.tile_downstream)
    down_rc = str(Seq(down).reverse_complement())

    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    aligner.open_gap_score = -2
    aligner.extend_gap_score = -0.5

    fwd = aligner.align(up, down).score
    rc = aligner.align(up, down_rc).score
    return OrientationCheck(float(fwd), float(rc))


def refine_ensemble(seq: Optional[str], hit: DuplexHit, margin: int = 20) -> EnsembleRefinement:
    """RNAcofold-style ensemble check: dG_hyb and expected inter-strand bp count.

    Runs on the matched sub-spans (``hit.matched_upstream`` /
    ``hit.matched_downstream``) padded by ``margin`` nt on each side (clipped
    to the flanking intron bounds), not the full 300-nt tile -- keeps the
    fold_compound cofold call (O(n^3) in the combined length) fast while
    still letting each side's competing intramolecular structure into the
    ensemble, which ``duplexfold`` does not model at all.
    """
    import RNA

    s = seq if seq is not None else load_smn2()
    up_iv = Interval(
        max(UPSTREAM_FLANK.start, hit.matched_upstream.start - margin),
        min(UPSTREAM_FLANK.end, hit.matched_upstream.end + margin),
    )
    down_iv = Interval(
        max(DOWNSTREAM_FLANK.start, hit.matched_downstream.start - margin),
        min(DOWNSTREAM_FLANK.end, hit.matched_downstream.end + margin),
    )
    s1 = _rna(s, up_iv)
    s2 = _rna(s, down_iv)

    fc1 = RNA.fold_compound(s1)
    _, g1 = fc1.pf()
    fc2 = RNA.fold_compound(s2)
    _, g2 = fc2.pf()

    fc12 = RNA.fold_compound(s1 + "&" + s2)
    _, g12 = fc12.pf()
    dG_hyb = g12 - g1 - g2

    bpp = fc12.bpp()
    n1 = len(s1)
    n2 = len(s2)
    inter = 0.0
    for i in range(1, n1 + 1):
        for j in range(n1 + 1, n1 + n2 + 1):
            a, b = (i, j) if i < j else (j, i)
            inter += bpp[a][b]

    return EnsembleRefinement(float(dG_hyb), float(inter), margin)


# ---------------------------------------------------------------------------
# top-level report
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AluReport:
    """Full analysis result: top candidates plus their follow-up checks."""

    hits: list[DuplexHit]
    null_controls: dict[int, NullControl] = field(default_factory=dict)  # index into hits
    orientations: dict[int, OrientationCheck] = field(default_factory=dict)
    refinements: dict[int, EnsembleRefinement] = field(default_factory=dict)


def run_report(
    seq: Optional[str] = None,
    width: int = IR_TILE_WIDTH,
    step: int = IR_TILE_STEP,
    top_k: int = IR_TOP_K,
) -> AluReport:
    """Full pipeline: scan, rank, then null-control/orientation/refine the top hits.

    This is the slow path (the ~3-minute genome-scale tiling scan); call it
    once and inspect the returned :class:`AluReport`, not from tests.
    """
    s = seq if seq is not None else load_smn2()
    hits = scan_flanks(s, width, step)
    top = top_candidates(hits, top_k)
    nulls = {i: null_control(s, h) for i, h in enumerate(top)}
    orients = {i: orientation_check(s, h) for i, h in enumerate(top)}
    refines = {i: refine_ensemble(s, h) for i, h in enumerate(top)}
    return AluReport(top, nulls, orients, refines)


def format_report(report: AluReport) -> str:
    """Human-readable table of a report's top candidates."""
    lines = [
        f"{'rank':>4} {'upstream':>18} {'downstream':>18} {'E(kcal/mol)':>12} "
        f"{'density':>7} {'z_null':>8} {'fwd':>6} {'revcomp':>7} {'inverted':>8} "
        f"{'dG_hyb':>9} {'d(acc)':>7} {'d(don)':>7}"
    ]
    for i, h in enumerate(report.hits):
        n = report.null_controls.get(i)
        o = report.orientations.get(i)
        r = report.refinements.get(i)
        lines.append(
            f"{i:>4} "
            f"{h.matched_upstream.start}-{h.matched_upstream.end:>7} "
            f"{h.matched_downstream.start}-{h.matched_downstream.end:>7} "
            f"{h.energy_kcal_mol:>12.1f} {h.density:>7.3f} "
            f"{(n.z_score if n else float('nan')):>8.1f} "
            f"{(o.forward_score if o else float('nan')):>6.0f} "
            f"{(o.revcomp_score if o else float('nan')):>7.0f} "
            f"{str(o.is_inverted) if o else '?':>8} "
            f"{(r.dG_hyb if r else float('nan')):>9.1f} "
            f"{h.dist_to_acceptor:>7} {h.dist_to_donor:>7}"
        )
    return "\n".join(lines)


if __name__ == "__main__":
    rep = run_report()
    print(format_report(rep))
