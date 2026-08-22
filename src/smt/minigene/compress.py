"""A6: algorithmic minigene compression.

Automates the design decision the Singh lab made by hand for SMN2Sup3445
(Ottesen et al. 2024, Nucleic Acids Research 52:3547, DOI 10.1093/nar/gkad1259):
given a locus, decide how much of each intron can be excised while preserving
everything a splicing assay depends on.

RETAIN, UNCONDITIONALLY (never subject to the structural test below)
----------------------------------------------------------------------
* the promoter / upstream flank present in the pinned reference record
* every annotated exon (and every literature-pinned cryptic exon, when the
  caller supplies one via ``extra_keep`` -- STMN2/UNC13A use this)
* a minimal splice-donor window (exon tail + intron head)
* a minimal splice-acceptor window (branch point + polypyrimidine tract + AG
  + exon head), located by scanning for the two consensus elements directly
  (see ``_find_branch_point`` / ``_find_ppt_end``) rather than by a fixed
  offset -- this is the "on principled grounds rather than by hand" part of
  the acceptor side.

EVERYTHING ELSE IN AN INTRON is a CANDIDATE for excision. The candidate
deletion is accepted only if removing it does not shift the PREDICTED LOCAL
STRUCTURE (P_unpaired, via ``smt.thermo.accessibility`` -- no thermodynamics
is reimplemented here) at either of the two regulatory elements bounding it,
by more than an explicit, documented threshold (``CompressionParams.
log2fc_threshold``). This is a MODELLING CHOICE, not a measured constant --
see the parameter's own docstring. If the shift exceeds the threshold at
every candidate deletion size tried, the algorithm falls back, all the way
to full retention if nothing smaller passes.

The growth search
------------------
For each gap between two retained anchors, candidate deletions are tried at
a fixed, documented ladder of "fraction of the gap retained" checkpoints
(``CompressionParams.growth_checkpoints``), from the most aggressive
(checkpoint 0.0: delete the whole gap) to the least (checkpoint 1.0: delete
nothing, i.e. full retention, which always trivially passes since nothing
changed). The first checkpoint (in that order) whose worst-case |log2 fold
change| at either bounding anchor is within threshold is accepted. This is a
monotonicity ASSUMPTION (less deletion is assumed not to make things worse),
not a proof; it is stated here rather than silently relied on.

Exact vs local folding
------------------------
Each gap's test span (bounding-anchor-readout to bounding-anchor-readout) is
folded with the EXACT ensemble ratio (``Method.GLOBAL``, same estimator
``smt.thermo.constrained``/``smt.thermo.occlusion`` already use) when that
span is at or under ``CompressionParams.global_pf_limit`` nt, and with
RNAplfold-style LOCAL folding (``Method.LOCAL``) otherwise. LOCAL folding
cannot see structure formed by base pairs farther apart than
``local_max_bp_span`` nt -- for the multi-kb introns this is unavoidable
(exact O(N^3) folding of a 25 kb intron is not tractable here), and it means
a genuinely long-range interaction beyond that span is invisible to this
test. That is a disclosed limitation, not a silent one; see
``docs/HANDOFF.md``'s repeated span-dependence findings for the SMN2 case,
which is the same class of limitation restated here for compression.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional, Sequence

from smt.seq.construct import Construct
from smt.seq.reference import Interval, slice_ref
from smt.thermo.accessibility import Method, accessibility

# ---------------------------------------------------------------------------
# modelling parameters -- ALL of these are explicit, documented choices, not
# measured constants. See CompressionParams docstring for the threshold in
# particular.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CompressionParams:
    """Every knob the compression algorithm exposes, with its justification.

    None of these numbers come from fitting to a known answer (the
    introns-3-and-7 result this module is validated against). Changing them
    is a legitimate sensitivity analysis; ``smt.minigene.report`` sweeps
    ``log2fc_threshold`` explicitly rather than reporting one hand-picked
    value as if it were derived.
    """

    #: nt of the upstream exon included in the donor readout window. 3+6 nt
    #: (exon+intron) mirrors the classical splice-donor consensus window
    #: (Yeo & Burge 2004, J Comput Biol 11:377 -- MaxEntScan's donor model)
    #: and matches this project's own existing SS5_EXON7 convention
    #: (smt.seq.construct.SS5_EXON7). Not re-derived per gene.
    donor_exon_flank: int = 3
    #: nt of intron retained at the donor (GT + immediate context).
    donor_intron_flank: int = 6

    #: Branch-point search window, nt upstream of the intron's terminal AG.
    #: 18-40 nt is the commonly cited mammalian range (e.g. Corvelo et al.
    #: 2010, PLoS Comput Biol 6:e1001016); widened to 50 here to be
    #: permissive for introns with a weaker/more distal branch. ASSUMPTION,
    #: not independently re-derived from SMN2/STMN2/UNC13A data.
    branch_search_min: int = 18
    branch_search_max: int = 50
    #: Degenerate mammalian branch-point consensus, DNA alphabet (T stands
    #: for the branch-point transcript's U). Textbook consensus, sometimes
    #: written "yUnAy" (Zhang & Chasin 2004, Genes Dev 18:1241); Y=pyrimidine,
    #: N=any, A/T literal. ASSUMPTION: generic consensus, not locus-specific.
    branch_consensus: str = "YTNAY"

    #: Sliding-window size (nt) and minimum pyrimidine fraction for calling a
    #: polypyrimidine tract. Reporting/annotation only (see module docstring
    #: -- the acceptor KEEP boundary is set by the branch point, not the PPT
    #: scan, since the PPT sits between the two and is already covered).
    #: ASSUMPTION: round numbers, not fit to any locus.
    ppt_window: int = 5
    ppt_min_fraction: float = 0.6

    #: nt of the downstream exon included in the acceptor readout window.
    acceptor_exon_flank: int = 3

    #: THE CORE MODELLING CHOICE. |log2(P_after / P_before)| at a bounding
    #: regulatory element must be <= this for a candidate excision to be
    #: accepted. This is NOT a measured constant -- there is no independent
    #: experimental calibration for "how much local structural perturbation
    #: is tolerable before a splice site's behaviour changes" in this
    #: project (per the no-placeholder-calibration rule). 1.0 (i.e. at most a
    #: 2x shift in P_unpaired) is a round, explicitly-flagged choice.
    #: smt.minigene.report sweeps this value; treat any single number here as
    #: illustrative, not derived.
    log2fc_threshold: float = 1.0

    #: Fraction of a candidate gap's length RETAINED (split between both
    #: sides), tried in this fixed order from most-aggressive to
    #: least-aggressive. 0.0 = delete the whole gap; 1.0 = delete nothing
    #: (always trivially passes, guaranteeing the search terminates).
    #: ASSUMPTION: a discretisation choice, not a claim about a smooth
    #: relationship between excision size and structural shift.
    growth_checkpoints: tuple[float, ...] = (0.0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.7, 0.85, 1.0)

    #: A gap's fold span (bounding readout to bounding readout) is folded
    #: EXACTLY (Method.GLOBAL) at or under this length; above it, LOCAL
    #: (RNAplfold) folding is used instead. Set well below
    #: smt.seq.construct.GLOBAL_PF_SOFT_LIMIT (2000 nt): this project's own
    #: A3 sweep smoke run (docs/HANDOFF.md, "97 s / 3 footprints (~32 s
    #: each)") measured ~16 s for a SINGLE exact partition function on a
    #: ~800 nt sequence, and this module needs up to
    #: 2 + 2*len(growth_checkpoints) such calls PER GAP -- at 1800 nt
    #: (~O(N^3) scaling) that is minutes per gap, intractable for the number
    #: of introns/checkpoints this validation runs interactively. 800 nt
    #: keeps every exact call in the tens-of-seconds range. This also means
    #: SMN2 intron 3 (159 nt) and intron 7 (444 nt) -- the two introns this
    #: module is validated against -- are folded EXACTLY over their full
    #: length, which is the highest-fidelity case available; only the large
    #: introns (1, 2a, 6, and others above 800 nt) fall back to LOCAL.
    global_pf_limit: int = 800

    #: RNAplfold window/max_bp_span for the LOCAL fallback. Larger than
    #: smt.thermo.accessibility's defaults (80/40) because the readout
    #: windows here (acceptor windows especially) can exceed 40 nt; NOT a
    #: claim that 300 nt captures every real long-range interaction -- see
    #: module docstring.
    local_window_size: int = 400
    local_max_bp_span: int = 300

    #: Probability floor to avoid log(0) on numerically-zero P_unpaired.
    p_floor: float = 1e-12


DEFAULT_PARAMS = CompressionParams()


# ---------------------------------------------------------------------------
# branch point / polypyrimidine tract detection
# ---------------------------------------------------------------------------


def _branch_score(window: str, consensus: str) -> int:
    score = 0
    for c, w in zip(consensus, window):
        if c == "N":
            score += 1
        elif c == "Y":
            score += 1 if w in "CT" else 0
        else:
            score += 1 if w == c else 0
    return score


def find_branch_point(intron_seq: str, params: CompressionParams = DEFAULT_PARAMS) -> int:
    """0-based offset (within ``intron_seq``) of the best branch-point match.

    Searches the last ``branch_search_max``..``branch_search_min`` nt of the
    intron for the best match to ``params.branch_consensus``. Ties are broken
    toward the 3'-most (closest to the AG) candidate, since the branch point
    is typically the last strong match before the polypyrimidine tract.
    Returns -1 if the intron is too short for the window to fit at all (the
    caller then falls back to full retention -- see ``detect_acceptor``).
    """
    n = len(intron_seq)
    k = len(params.branch_consensus)
    lo = max(0, n - params.branch_search_max)
    hi = n - params.branch_search_min - k
    if hi < lo:
        return -1
    best_score = -1
    best_pos = -1
    for p in range(lo, hi + 1):
        window = intron_seq[p : p + k]
        score = _branch_score(window, params.branch_consensus)
        if score >= best_score:  # >= : ties go to the later (more 3') position
            best_score = score
            best_pos = p
    return best_pos


def find_ppt_end(intron_seq: str, start: int, params: CompressionParams = DEFAULT_PARAMS) -> int:
    """0-based offset where the contiguous pyrimidine-rich run starting at
    ``start`` ends (first position the sliding window drops below
    ``ppt_min_fraction``, or the end of the sequence). Reporting/annotation
    only; see module docstring."""
    n = len(intron_seq)
    pos = start
    w = params.ppt_window
    while pos + w <= n:
        window = intron_seq[pos : pos + w]
        frac = sum(1 for c in window if c in "CT") / w
        if frac < params.ppt_min_fraction:
            break
        pos += 1
    return min(pos, n)


@dataclass(frozen=True)
class AcceptorSite:
    """Result of scanning one intron's 3' end for its acceptor elements.

    All intervals are REFERENCE space. ``keep`` is the minimal always-retained
    block (branch point through the AG); ``branch_point``/``ppt`` are
    reported for the construct spec but do not independently constrain
    ``keep`` (the branch point already anchors it upstream of the PPT).
    """

    intron: Interval
    branch_point: Optional[Interval]  # single consensus-window hit, or None if not found
    ppt: Optional[Interval]
    keep: Interval
    found: bool  # False -> branch search window did not fit; keep == whole intron


def detect_acceptor(seq: str, intron: Interval, params: CompressionParams = DEFAULT_PARAMS) -> AcceptorSite:
    intron_seq = slice_ref(seq, intron)
    branch_off = find_branch_point(intron_seq, params)
    if branch_off < 0:
        # Intron too short for the branch window to fit at all -- nothing to
        # excise beyond what the donor side already claims; keep it whole.
        return AcceptorSite(intron, None, None, intron, found=False)
    k = len(params.branch_consensus)
    branch_iv = Interval(intron.start + branch_off, intron.start + branch_off + k - 1)
    ppt_end_off = find_ppt_end(intron_seq, branch_off + k, params)
    ppt_iv: Optional[Interval] = None
    if ppt_end_off > branch_off + k:
        ppt_iv = Interval(intron.start + branch_off + k, intron.start + ppt_end_off - 1)
    keep = Interval(branch_iv.start, intron.end)
    return AcceptorSite(intron, branch_iv, ppt_iv, keep, found=True)


def detect_donor(exon_end: int, intron: Interval, params: CompressionParams = DEFAULT_PARAMS) -> Interval:
    """Minimal always-retained donor block, REFERENCE space, intron-only
    portion (the exon side is retained unconditionally as part of the exon)."""
    return Interval(intron.start, min(intron.end, intron.start + params.donor_intron_flank - 1))


# ---------------------------------------------------------------------------
# anchors and gaps
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Anchor:
    """One retained landmark bounding a candidate excision gap(s).

    ``keep``: REFERENCE-space interval that must remain in the sequence.
    An anchor in the MIDDLE of an intron (an ``extra_keep`` entry) bounds TWO
    gaps, one on each side, and the regulatory element worth protecting is
    usually different on each side (e.g. a cryptic exon's own acceptor faces
    the upstream gap; its donor, or just its own terminal sequence, faces the
    downstream gap) -- so there are two readouts:

    ``readout_entry``: tested for the gap UPSTREAM of this anchor (the gap
    where this anchor is the RIGHT/3' bound) -- should sit near ``keep.start``.
    ``readout_exit``: tested for the gap DOWNSTREAM of this anchor (the gap
    where this anchor is the LEFT/5' bound) -- should sit near ``keep.end``.

    A donor anchor is only ever a LEFT bound, so only ``readout_exit`` is
    ever read; an acceptor anchor is only ever a RIGHT bound, so only
    ``readout_entry`` is ever read. Both are still populated (with the same
    value) for symmetry/inspection.
    """

    name: str
    keep: Interval
    readout_entry: Interval
    readout_exit: Interval


@dataclass(frozen=True)
class GapMeasurement:
    """One growth-search checkpoint's measurement, kept for full transparency
    and so a threshold sweep can re-decide without recomputing any folds."""

    checkpoint: float
    deletion: Optional[Interval]  # None == nothing deleted at this checkpoint
    p_left_before: float
    p_left_after: float
    p_right_before: float
    p_right_after: float
    log2fc_left: float
    log2fc_right: float
    method: str

    def worst_log2fc(self) -> float:
        return max(abs(self.log2fc_left), abs(self.log2fc_right))


@dataclass(frozen=True)
class GapResult:
    left: Anchor
    right: Anchor
    gap: Optional[Interval]  # None if the anchors already abut (nothing to test)
    measurements: tuple[GapMeasurement, ...]
    chosen_checkpoint: float
    deletion: Optional[Interval]  # REFERENCE space, or None (fully retained)


def _log2fc(p_before: float, p_after: float, floor: float) -> float:
    return math.log2(max(p_after, floor) / max(p_before, floor))


def _accessibility_p(sub_seq: str, iv_local: Interval, method: Method, params: CompressionParams) -> float:
    window_size = min(params.local_window_size, len(sub_seq))
    max_bp_span = min(params.local_max_bp_span, max(len(sub_seq) - 1, 1))
    return accessibility(sub_seq, iv_local, method=method, window_size=window_size, max_bp_span=max_bp_span).p_unpaired


def process_gap(seq: str, left: Anchor, right: Anchor, params: CompressionParams = DEFAULT_PARAMS) -> GapResult:
    """Run the growth search for the candidate excision between two anchors.

    ``seq`` is the FULL natural reference/gene sequence (uppercase DNA);
    slicing to the gap's local test span happens inside this function.
    """
    gap = Interval(left.keep.end + 1, right.keep.start - 1)
    if len(gap) <= 0:
        return GapResult(left, right, None, (), 1.0, None)

    left_readout = left.readout_exit  # left anchor's edge FACING this gap
    right_readout = right.readout_entry  # right anchor's edge FACING this gap
    span = Interval(left_readout.start, right_readout.end)
    natural_sub = slice_ref(seq, span)
    method = Method.GLOBAL if len(natural_sub) <= params.global_pf_limit else Method.LOCAL

    left_local_natural = Interval(left_readout.start - span.start + 1, left_readout.end - span.start + 1)
    right_local_natural = Interval(right_readout.start - span.start + 1, right_readout.end - span.start + 1)
    p_left_before = _accessibility_p(natural_sub, left_local_natural, method, params)
    p_right_before = _accessibility_p(natural_sub, right_local_natural, method, params)

    measurements: list[GapMeasurement] = []
    gap_len = len(gap)

    # Always run the FULL checkpoint ladder (no early break) so a threshold
    # sweep (smt.minigene.report) can re-decide the outcome at any other
    # threshold via rescore_gap() without recomputing a single fold. This
    # costs len(growth_checkpoints) fold evaluations per gap regardless of
    # where the current threshold would have stopped -- an explicit
    # compute-for-honesty tradeoff, not an oversight.
    for r in params.growth_checkpoints:
        retained_each_side = int(round(r * gap_len / 2))
        del_start = gap.start + retained_each_side
        del_end = gap.end - retained_each_side
        deletion: Optional[Interval] = Interval(del_start, del_end) if del_start <= del_end else None

        if deletion is None:
            after_sub = natural_sub
            right_local_after = right_local_natural
        else:
            d_local_start = deletion.start - span.start + 1
            d_local_end = deletion.end - span.start + 1
            after_sub = natural_sub[: d_local_start - 1] + natural_sub[d_local_end:]
            removed_len = d_local_end - d_local_start + 1
            right_local_after = Interval(right_local_natural.start - removed_len, right_local_natural.end - removed_len)
        left_local_after = left_local_natural  # untouched: deletion always lies strictly after it

        p_left_after = _accessibility_p(after_sub, left_local_after, method, params)
        p_right_after = _accessibility_p(after_sub, right_local_after, method, params)

        m = GapMeasurement(
            checkpoint=r,
            deletion=deletion,
            p_left_before=p_left_before,
            p_left_after=p_left_after,
            p_right_before=p_right_before,
            p_right_after=p_right_after,
            log2fc_left=_log2fc(p_left_before, p_left_after, params.p_floor),
            log2fc_right=_log2fc(p_right_before, p_right_after, params.p_floor),
            method=method.value,
        )
        measurements.append(m)

    result = GapResult(left, right, gap, tuple(measurements), 1.0, None)
    chosen_checkpoint, chosen_deletion = rescore_gap(result, params.log2fc_threshold)
    return GapResult(left, right, gap, tuple(measurements), chosen_checkpoint, chosen_deletion)


def rescore_gap(gap_result: GapResult, threshold: float) -> tuple[float, Optional[Interval]]:
    """Re-decide a gap's outcome at a DIFFERENT threshold using already-computed
    measurements -- no folding is repeated. Used by the threshold sweep."""
    if gap_result.gap is None:
        return 1.0, None
    for m in gap_result.measurements:
        if m.worst_log2fc() <= threshold:
            return m.checkpoint, m.deletion
    last = gap_result.measurements[-1]
    return last.checkpoint, last.deletion  # checkpoint 1.0 always passes trivially


# ---------------------------------------------------------------------------
# per-intron and per-gene orchestration
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class IntronPlan:
    name: str  # matches this project's exon-name-based intron labelling
    interval: Interval  # full natural intron, reference space
    donor: Anchor
    acceptor: Anchor
    acceptor_site: AcceptorSite
    extra: tuple[Anchor, ...]
    gaps: tuple[GapResult, ...]

    @property
    def deletions(self) -> tuple[Interval, ...]:
        return tuple(g.deletion for g in self.gaps if g.deletion is not None)

    @property
    def retained_length(self) -> int:
        return len(self.interval) - sum(len(d) for d in self.deletions)

    @property
    def retained_fraction(self) -> float:
        return self.retained_length / len(self.interval)


def compress_intron(
    seq: str,
    name: str,
    exon_a_end: int,
    exon_b_start: int,
    params: CompressionParams = DEFAULT_PARAMS,
    extra_keep: Sequence[tuple] = (),
) -> IntronPlan:
    """Decide retained/excised blocks for one intron.

    ``extra_keep``: literature-pinned ``(keep, entry_readout)`` or
    ``(keep, entry_readout, exit_readout)`` tuples for regions INSIDE this
    intron that must survive wholesale (e.g. a cryptic exon and its own
    splice sites) -- see ``Anchor`` for what entry/exit mean. Empty for a
    blind run (SMN2's validation, item 2); non-empty for STMN2/UNC13A
    (item 3), where the cryptic exon coordinates are already pinned in
    ``smt.seq.targets`` and are not something the algorithm is meant to
    "discover".
    """
    intron = Interval(exon_a_end + 1, exon_b_start - 1)
    donor_iv = detect_donor(exon_a_end, intron, params)
    donor_readout = Interval(exon_a_end - params.donor_exon_flank + 1, donor_iv.end)
    donor = Anchor(f"intron{name}_donor", donor_iv, donor_readout, donor_readout)

    acc_site = detect_acceptor(seq, intron, params)
    acceptor_readout = Interval(acc_site.keep.start, exon_b_start + params.acceptor_exon_flank - 1)
    acceptor = Anchor(f"intron{name}_acceptor", acc_site.keep, acceptor_readout, acceptor_readout)

    extras: list[Anchor] = []
    for i, spec in enumerate(extra_keep):
        keep, entry_readout = spec[0], spec[1]
        # exit_readout defaults to the last few nt of keep (mirrors the donor
        # convention of a short window at the edge facing the downstream
        # gap) when the caller has no pinned feature for that side -- e.g.
        # STMN2's cryptic exon has no annotated 3' regulatory element (it
        # terminates by intrinsic polyadenylation, see smt.seq.targets).
        # ASSUMPTION, flagged: this fallback is a generic "protect this
        # anchor's own terminal sequence" stand-in, not a literature feature.
        if len(spec) > 2 and spec[2] is not None:
            exit_readout = spec[2]
        else:
            w = min(9, len(keep))
            exit_readout = Interval(keep.end - w + 1, keep.end)
        extras.append(Anchor(f"intron{name}_extra{i}", keep, entry_readout, exit_readout))

    anchors = sorted([donor, *extras, acceptor], key=lambda a: a.keep.start)
    # Merge anchors that already overlap/abut (defensive; real pinned data
    # should not need this, but a short intron can push donor and acceptor
    # blocks into contact).
    merged: list[Anchor] = []
    for a in anchors:
        if merged and a.keep.start <= merged[-1].keep.end + 1:
            prev = merged[-1]
            new_keep = Interval(prev.keep.start, max(prev.keep.end, a.keep.end))
            merged[-1] = Anchor(prev.name, new_keep, prev.readout_entry, a.readout_exit)
        else:
            merged.append(a)

    gaps = tuple(process_gap(seq, left, right, params) for left, right in zip(merged, merged[1:]))
    return IntronPlan(name, intron, donor, acceptor, acc_site, extras, gaps)


@dataclass(frozen=True)
class GenePlan:
    gene: str
    seq: str
    exon_order: tuple[str, ...]
    exons: dict[str, Interval]  # reference space
    upstream: Interval  # promoter/flank retained wholesale, reference space
    downstream_excised: Optional[Interval]  # trailing sequence beyond the last exon, dropped wholesale
    introns: dict[str, IntronPlan]
    params: CompressionParams

    def all_deletions(self) -> list[Interval]:
        dels: list[Interval] = []
        for ip in self.introns.values():
            dels.extend(ip.deletions)
        if self.downstream_excised is not None:
            dels.append(self.downstream_excised)
        return sorted(dels, key=lambda d: d.start)

    def natural_length(self) -> int:
        return len(self.seq)

    def transcribed_unit_length(self) -> int:
        """upstream through the last exon -- excludes trailing genomic
        flank beyond the transcript, which is dropped wholesale and isn't a
        "compression" decision at all."""
        last = self.exons[self.exon_order[-1]]
        return last.end - self.upstream.start + 1


def compress_gene(
    gene: str,
    seq: str,
    exon_order: Sequence[str],
    exons: dict[str, Interval],
    params: CompressionParams = DEFAULT_PARAMS,
    extra_keep_by_intron: Optional[dict[str, Sequence[tuple[Interval, Interval]]]] = None,
    keep_upstream_from: int = 1,
) -> GenePlan:
    """Compress every intron of a gene given its (pinned) exon coordinates.

    ``exon_order``: exon names in 5'->3' order, matching keys of ``exons``.
    ``keep_upstream_from``: reference position the retained promoter/flank
    starts at (1 == the whole pinned record upstream of the first exon).
    Nothing after ``keep_upstream_from`` and before the first exon is ever
    excised -- promoter compression is out of scope for this module (see
    module docstring: this is an INTRON compression algorithm).
    """
    extra_keep_by_intron = extra_keep_by_intron or {}
    upstream = Interval(keep_upstream_from, exons[exon_order[0]].start - 1)
    introns: dict[str, IntronPlan] = {}
    for a, b in zip(exon_order, exon_order[1:]):
        ea, eb = exons[a], exons[b]
        plan = compress_intron(seq, a, ea.end, eb.start, params, extra_keep_by_intron.get(a, ()))
        introns[a] = plan

    last = exons[exon_order[-1]]
    downstream = Interval(last.end + 1, len(seq)) if last.end < len(seq) else None

    return GenePlan(gene, seq, tuple(exon_order), dict(exons), upstream, downstream, introns, params)


# ---------------------------------------------------------------------------
# assembling the compressed Construct
# ---------------------------------------------------------------------------


def build_construct(plan: GenePlan, name: Optional[str] = None, extra_features: Optional[dict[str, Interval]] = None) -> Construct:
    """Assemble the compressed sequence as a generic ``smt.seq.construct.
    Construct`` (window = upstream start .. last exon end, deletions = every
    accepted intron excision). Mirrors exactly how ``smt.seq.targets`` builds
    a ``Construct`` by hand for STMN2/UNC13A -- no dependency on
    ``smt.seq.construct._build`` (which is hard-wired to SMN2).
    """
    window = Interval(plan.upstream.start, plan.exons[plan.exon_order[-1]].end)
    dels = tuple(d for d in plan.all_deletions() if d != plan.downstream_excised)
    pieces: list[str] = []
    cursor = window.start
    for d in dels:
        if d.start > cursor:
            pieces.append(slice_ref(plan.seq, Interval(cursor, d.start - 1)))
        cursor = d.end + 1
    if cursor <= window.end:
        pieces.append(slice_ref(plan.seq, Interval(cursor, window.end)))
    sequence = "".join(pieces)

    features: dict[str, Interval] = dict(plan.exons)
    for iname, ip in plan.introns.items():
        features[f"intron{iname}_donor_keep"] = ip.donor.keep
        features[f"intron{iname}_acceptor_keep"] = ip.acceptor.keep
        for extra in ip.extra:
            features[extra.name] = extra.keep
    if extra_features:
        features.update(extra_features)

    c = Construct(name or f"{plan.gene}_compressed", window, dels, sequence, {})
    keep: dict[str, Interval] = {}
    for k, v in features.items():
        try:
            c.map_interval(v)
        except ValueError:
            continue
        keep[k] = v
    return Construct(name or f"{plan.gene}_compressed", window, dels, sequence, keep)
