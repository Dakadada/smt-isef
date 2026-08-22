"""Tests for the A7 deliverable: STMN2/UNC13A cryptic-exon occlusion.

Covers:
* pinned STMN2/UNC13A coordinates and lengths assert correctly, and the
  splice-site dinucleotides (AG acceptors, GT donor) implied by the pinned
  boundaries are actually present in the pinned sequence at those exact
  positions -- an independent structural check, not just a re-statement of
  the numbers in smt.seq.targets;
* the signed figure of merit (smt.thermo.occlusion) is positive for a
  footprint that does the therapeutically desired thing and negative for the
  opposite direction on the SAME underlying physical effect -- i.e. the sign
  convention is exactly a flip, not two different computations;
* the SMN2 FREE regime still reproduces the known positive GC_CORE result
  (regression guard against docs/HANDOFF.md's "RESOLVED 00:55" finding:
  blocking GC_CORE increases A_5ss by ~6-7x at span=100);
* a footprint overlapping the readout window itself is refused, not silently
  scored, in the OCCLUDE regime exactly as it already is in the FREE regime
  (smt.thermo.constrained.a_5ss's guard is shared, not regime-specific);
* smt.validate.cryptic_sweep's shard partition is exact (union == full
  sweep, no overlap) across all four (gene, splice-site) targets;
* the variant-scoring path (smt.thermo.occlusion.variant_effect) runs
  end-to-end on a real pinned SNP (rs12973192) and returns finite numbers.
"""

from __future__ import annotations

import pytest

from smt.seq import targets as T
from smt.seq.construct import GC_CORE, compact_window
from smt.seq.reference import Interval
from smt.thermo.constrained import occlusion_delta
from smt.thermo.occlusion import Direction, SignedOcclusionResult, overlaps, score, variant_effect
from smt.validate import cryptic_sweep as cs

# ---------------------------------------------------------------------------
# pinned coordinates / lengths / structural sanity
# ---------------------------------------------------------------------------


def test_stmn2_pinned_length_and_checksum():
    seq = T.load_stmn2()
    assert len(seq) == T.STMN2_LENGTH
    T.verify_checksums()  # re-derives and asserts both genes' cryptic-exon hashes


def test_unc13a_pinned_length_and_checksum():
    seq = T.load_unc13a()
    assert len(seq) == T.UNC13A_LENGTH


def test_stmn2_exon_table_matches_pinned_tuples():
    for name, (s, e) in T.STMN2_EXONS.items():
        iv = Interval(s, e)
        assert len(iv) == e - s + 1
    # exon 1 and the cryptic exon must both fall inside intron-1's span, i.e.
    # exon1 end < cryptic exon start < exon2 start (the literature's "cryptic
    # exon embedded within intron 1" claim, checked against our own table).
    e1_end = T.STMN2_EXONS["1"][1]
    e2_start = T.STMN2_EXONS["2"][0]
    assert e1_end < T.STMN2_CRYPTIC_EXON.start
    assert T.STMN2_CRYPTIC_EXON.end < e2_start


def test_stmn2_cryptic_exon_length_and_splice_acceptor():
    seq = T.load_stmn2()
    assert len(T.STMN2_CRYPTIC_EXON) == 227  # Melamed et al. 2019 stated length
    # canonical 3' splice site acceptor dinucleotide immediately 5' of the
    # cryptic exon: intron ends "AG" at local 5705-5706.
    assert seq[5704:5706] == "AG"
    # the readout window's own last 2 nt must be that same AG (internal
    # consistency between STMN2_CRYPTIC_3SS and STMN2_CRYPTIC_EXON).
    ss = T.STMN2_CRYPTIC_3SS
    assert seq[ss.end - 5 : ss.end - 3] == "AG"


def test_unc13a_cryptic_exons_between_canonical_exons_20_and_21():
    e20_end = T.UNC13A_EXONS["20"][1]
    e21_start = T.UNC13A_EXONS["21"][0]
    assert e20_end < T.UNC13A_CRYPTIC_EXON_178.start
    assert e20_end < T.UNC13A_CRYPTIC_EXON_128.start
    assert T.UNC13A_CRYPTIC_EXON_128.end < e21_start
    assert T.UNC13A_CRYPTIC_EXON_178.end < e21_start
    # both isoforms share the same downstream (donor-side) boundary
    assert T.UNC13A_CRYPTIC_EXON_128.end == T.UNC13A_CRYPTIC_EXON_178.end == 45940
    # published lengths, recomputed from the pinned intervals independently
    # of the constants that name them
    assert len(T.UNC13A_CRYPTIC_EXON_128) == 128
    assert len(T.UNC13A_CRYPTIC_EXON_178) == 178


def test_unc13a_splice_site_dinucleotides():
    seq = T.load_unc13a()
    # shared donor: GT immediately after the cryptic exon's 3' boundary
    assert seq[45940:45942] == "GT"
    # 128bp-form acceptor: AG immediately before its 5' boundary
    assert seq[45810:45812] == "AG"
    # 178bp-form acceptor: AG immediately before its (further upstream) 5' boundary
    assert seq[45760:45762] == "AG"


def test_unc13a_snp_reference_alleles_match_pinned_sequence():
    seq = T.load_unc13a()
    assert seq[T.RS12608932.local_pos - 1] == T.RS12608932.ref_local
    assert seq[T.RS12973192.local_pos - 1] == T.RS12973192.ref_local
    # rs12973192 is reported to sit inside the cryptic exon itself
    assert T.UNC13A_CRYPTIC_EXON_128.contains(T.RS12973192.local_pos)
    assert T.UNC13A_CRYPTIC_EXON_178.contains(T.RS12973192.local_pos)
    # rs12608932 is reported to sit in the same intron, but NOT inside either
    # cryptic exon isoform
    assert not T.UNC13A_CRYPTIC_EXON_128.contains(T.RS12608932.local_pos)
    assert not T.UNC13A_CRYPTIC_EXON_178.contains(T.RS12608932.local_pos)


# ---------------------------------------------------------------------------
# signed figure of merit: sign convention + SMN2 regression guard
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def gc_core_delta_span100():
    """One real constrained-folding computation, reused by both tests below
    so the sign-convention check and the regression check do not each pay
    for their own partition-function pass."""
    c = compact_window()
    return occlusion_delta(GC_CORE, c, max_bp_span=100)


def test_smn2_gc_core_regression_guard(gc_core_delta_span100):
    """docs/HANDOFF.md, 'RESOLVED 00:55': blocking GC_CORE increases A_5ss at
    span=100 by ~7.1x, delta_A ~= +2.577e-04. This is the FREE-regime
    baseline the OCCLUDE regime is scored against; if this regresses, the
    sign-convention test below is not testing anything meaningful."""
    d = gc_core_delta_span100
    assert d.delta_A > 0
    assert d.delta_A == pytest.approx(2.577e-04, rel=0.05)
    assert d.fold_change == pytest.approx(7.1, rel=0.1)


def test_direction_sign_convention_is_a_pure_flip(gc_core_delta_span100):
    """Same underlying OcclusionDelta, both directions applied: FOM must be
    an exact sign flip of the same delta_A, and the therapeutically-desired
    labelling (`.desired`) must disagree between the two directions for a
    result that is unambiguously nonzero."""
    d = gc_core_delta_span100
    free = SignedOcclusionResult(Direction.FREE, d)
    occlude = SignedOcclusionResult(Direction.OCCLUDE, d)

    assert free.fom == d.delta_A
    assert occlude.fom == -d.delta_A
    assert free.fom == -occlude.fom

    # GC_CORE genuinely frees the SMN2 5'ss (known positive result): scored
    # as desired in the FREE regime, undesired in the OCCLUDE regime.
    assert free.fom > 0
    assert free.desired is True
    assert occlude.fom < 0
    assert occlude.desired is False


def test_direct_overlap_is_refused_in_both_directions():
    """A footprint that swallows the readout window itself must be refused
    (ValueError from smt.thermo.constrained.a_5ss's guard), not silently
    scored -- this holds for OCCLUDE just as much as FREE. See
    smt.thermo.occlusion's module docstring for why routing around this
    guard would produce a physically meaningless number rather than a
    conservative one."""
    c = T.stmn2_cryptic_window()
    site = T.STMN2_CRYPTIC_3SS
    assert overlaps(site, site)  # sanity: a window always overlaps itself
    with pytest.raises(ValueError):
        score(site, Direction.OCCLUDE, "STMN2_CRYPTIC_3SS", c, max_bp_span=100)
    with pytest.raises(ValueError):
        score(site, Direction.FREE, "STMN2_CRYPTIC_3SS", c, max_bp_span=100)


def test_overlaps_helper_matches_a5ss_guard_semantics():
    site = Interval(100, 108)
    assert overlaps(Interval(90, 100), site)  # touches at the left edge
    assert overlaps(Interval(108, 120), site)  # touches at the right edge
    assert overlaps(Interval(103, 105), site)  # fully inside
    assert not overlaps(Interval(80, 99), site)  # strictly left
    assert not overlaps(Interval(109, 130), site)  # strictly right


# ---------------------------------------------------------------------------
# sharding: exact partition, all four targets
# ---------------------------------------------------------------------------


def test_cryptic_sweep_shard_partition_exact():
    items = cs.enumerate_footprints()
    assert len(items) > 0
    nshards = 37  # deliberately not a divisor of len(items), like test_sieve's audit
    parts = [cs.shard_of(items, s, nshards) for s in range(nshards)]

    total = sum(len(p) for p in parts)
    assert total == len(items)  # no duplication, nothing dropped

    seen = set()
    for p in parts:
        keys = {fp.key for fp in p}
        assert not (keys & seen)  # no overlap between shards
        seen |= keys
    assert seen == {fp.key for fp in items}  # union == full sweep

    sizes = [len(p) for p in parts]
    assert max(sizes) - min(sizes) <= 1  # stride partition is balanced

    with pytest.raises(ValueError):
        cs.shard_of(items, nshards, nshards)  # shard index out of range


def test_cryptic_sweep_targets_cover_both_genes():
    genes = {t.gene for t in cs.TARGETS.values()}
    assert genes == {"STMN2", "UNC13A"}


# ---------------------------------------------------------------------------
# SNP variant effect: end-to-end smoke test on rs12973192
# ---------------------------------------------------------------------------


def test_variant_effect_runs_and_returns_finite_numbers():
    c = T.unc13a_cryptic_window()
    ve = variant_effect(
        c,
        ref_pos=T.RS12973192.local_pos,
        alt_base=T.RS12973192.alt_local,
        site_feature="UNC13A_CRYPTIC_DONOR",
        max_bp_span=100,
    )
    assert 0.0 <= ve.ref.A_5ss <= 1.0
    assert 0.0 <= ve.alt.A_5ss <= 1.0
    assert ve.delta_A == pytest.approx(ve.alt.A_5ss - ve.ref.A_5ss)
    import math

    assert math.isfinite(ve.delta_A)
