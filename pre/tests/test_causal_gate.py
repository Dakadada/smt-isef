"""Tests for C4, the pre-declared gate that killed arms C and K.

A gate that cannot fail is not a gate, and a REFUTATION has to be as well
tested as a positive result -- otherwise "we checked and it didn't work"
is unfalsifiable. Covers:

* ``WIDENING_ORDER``'s load-bearing semantic -- that ``max_bp_span = 0`` means
  UNLIMITED and therefore sorts LAST -- re-derived from ViennaRNA itself by
  folding a hairpin whose closing pair spans further than the cap, NOT by
  restating the module's constant. Reading 0 as "smallest" inverts Test 1's
  entire conclusion, so this is the single most dangerous assumption in the
  module;
* ``partner_profile`` reports direction correctly, checked on synthetic
  hairpins built so the answer is known by CONSTRUCTION -- one where the
  readout can only pair downstream, one where it can only pair upstream;
* ``flip_transitions`` counts sign changes at the right step, on hand-built
  ladders including the no-flip and every-step-flips extremes;
* the A3 sweep on disk still yields the 795 full-ladder footprints the report
  was computed from, re-derived from the CSVs rather than restated;
* REGRESSION GUARD ON THE NEGATIVE ITSELF: the exon 7 5'ss remains
  >= 95% upstream-paired at every span. If this ever stops being true the
  refutation must be revisited, because it is the whole basis for the claim
  that a prefix cut has nothing to remove.
"""

from __future__ import annotations

import pytest

from smt.seq.reference import Interval
from smt.validate.causal_gate import (
    DEFAULT_SWEEP_DIR,
    WIDENING_ORDER,
    flip_transitions,
    load_sweep,
    partner_profile,
)

# A hairpin whose closing pair spans ~70 nt: a G-run, a long A loop that cannot
# pair with anything, and the complementary C-run.
LOOP = "A" * 60
LONG_RANGE = "AAAAA" + "GGGGGGG" + LOOP + "CCCCCCC" + "AAAAA"
_G_RUN = Interval(6, 12)
_C_RUN = Interval(len(LONG_RANGE) - 11, len(LONG_RANGE) - 5)


def _pair_span() -> int:
    return _C_RUN.end - _G_RUN.start


def test_max_bp_span_zero_means_unlimited_not_smallest():
    """Re-derived from ViennaRNA: the long-range pair forms at 0, not at 50."""
    span = _pair_span()
    assert span > 50, "fixture must exceed the cap it is testing"

    unlimited = partner_profile(LONG_RANGE, _G_RUN, max_bp_span=0)
    capped = partner_profile(LONG_RANGE, _G_RUN, max_bp_span=50)

    assert unlimited["downstream"] > 0.5, "unlimited span must reach the C-run"
    assert capped["downstream"] < 0.05, "a 50 nt cap must forbid a ~70 nt pair"
    # Hence 0 is the WIDEST setting and must sort last.
    assert WIDENING_ORDER[-1] == 0
    assert WIDENING_ORDER[:-1] == tuple(sorted(WIDENING_ORDER[:-1]))


def test_partner_profile_direction_is_by_construction():
    """The G-run pairs downstream; the C-run pairs upstream. Same molecule."""
    downstream = partner_profile(LONG_RANGE, _G_RUN, max_bp_span=0)
    upstream = partner_profile(LONG_RANGE, _C_RUN, max_bp_span=0)

    assert downstream["downstream"] > 0.9 * downstream["total"]
    assert upstream["upstream"] > 0.9 * upstream["total"]
    # Both readouts describe the SAME pairs, so their totals must agree.
    assert downstream["total"] == pytest.approx(upstream["total"], rel=0.05)


def test_partner_profile_bins_sum_to_downstream_mass():
    prof = partner_profile(LONG_RANGE, _G_RUN, max_bp_span=0)
    assert sum(prof["bins"].values()) == pytest.approx(prof["downstream"], rel=1e-9)
    assert prof["upstream"] + prof["downstream"] == pytest.approx(prof["total"], rel=1e-9)


def test_flip_transitions_locates_the_step():
    order = (50, 100, 150, 300, 0)
    sweep = {
        ("stable", 0): {50: 1.0, 100: 2.0, 150: 3.0, 300: 4.0, 0: 5.0},
        ("flips_early", 0): {50: -1.0, 100: 2.0, 150: 3.0, 300: 4.0, 0: 5.0},
        ("flips_late", 0): {50: 1.0, 100: 2.0, 150: 3.0, 300: 4.0, 0: -5.0},
        ("alternating", 0): {50: 1.0, 100: -1.0, 150: 1.0, 300: -1.0, 0: 1.0},
    }
    trans, n_full, flipping = flip_transitions(sweep, order)

    assert n_full == 4
    assert flipping == 3, "the stable footprint must not count as flipping"
    assert trans[(50, 100)] == 2      # flips_early + alternating
    assert trans[(300, 0)] == 2       # flips_late + alternating
    assert trans[(100, 150)] == 1     # alternating only
    assert trans[(150, 300)] == 1


def test_flip_transitions_ignores_partial_ladders():
    sweep = {("partial", 0): {50: 1.0, 100: -1.0}}
    trans, n_full, flipping = flip_transitions(sweep, WIDENING_ORDER)
    assert (n_full, flipping, sum(trans.values())) == (0, 0, 0)


@pytest.mark.skipif(not DEFAULT_SWEEP_DIR.exists(), reason="A3 sweep not on disk")
def test_sweep_on_disk_still_yields_the_reported_footprint_count():
    """Re-derived from the CSVs, not restated from docs/c4-causal-gate.txt."""
    sweep = load_sweep(DEFAULT_SWEEP_DIR)
    _, n_full, flipping = flip_transitions(sweep)
    assert n_full == 795
    assert 0.70 <= flipping / n_full <= 0.71


def test_exon7_5ss_is_upstream_paired_at_every_span():
    """REGRESSION GUARD ON THE REFUTATION ITSELF.

    Arms C and K were killed because a prefix cut has nothing to remove: the
    readout's partners are upstream. If that ever stops holding, the negative
    is no longer supported and must be re-opened rather than inherited.
    """
    from smt.seq.construct import SS5_EXON7, compact_window

    c = compact_window()
    win = c.window
    readout = Interval(SS5_EXON7.start - win.start + 1, SS5_EXON7.end - win.start + 1)

    for span in WIDENING_ORDER:
        prof = partner_profile(c.sequence, readout, span)
        share = prof["upstream"] / prof["total"]
        assert share >= 0.95, f"span={span}: upstream share fell to {share:.3f}"


def test_readout_has_more_sequence_downstream_than_upstream():
    """The refutation is not a positional artifact -- re-derived from the record."""
    from smt.seq.construct import SS5_EXON7, compact_window

    c = compact_window()
    win = c.window
    upstream_nt = SS5_EXON7.start - win.start
    downstream_nt = win.end - SS5_EXON7.end
    assert downstream_nt > upstream_nt
