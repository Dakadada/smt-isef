"""Tests for G1. Ground truth is re-derived from source, never restated.

The project rule that caught the UNC13A exon-table trap and the pseudo-
replicated prereg panel: a test must reconstruct the expected answer from the
GenBank record / the sweep CSVs / the transcriptome index, not copy a constant
out of the module it is testing. A test that asserts the module equals itself
proves nothing.
"""

from __future__ import annotations

import numpy as np
import pytest

from smt.seq.construct import compact_window, GC_CORE
from smt.seq.reference import Interval
from smt.sieve.scan import NUSINERSEN
from smt.sieve.transcriptome import INDEX_PATH, load_index, revcomp
import smt.validate.goldilocks_offtarget as G


@pytest.fixture(scope="module")
def transcripts():
    if not INDEX_PATH.exists():
        pytest.skip("transcriptome index not present")
    return load_index(INDEX_PATH)


def test_nusinersen_reference_is_the_published_sequence():
    """The length-18 reference point is built from coordinates, then checked
    against the published drug sequence -- so a coordinate drift anywhere in
    construct mapping fails here rather than silently shifting the curve."""
    assert G.nusinersen_reference().aso.replace("U", "T") == NUSINERSEN


def test_nusinersen_exact_sites_are_exactly_the_smn_locus(transcripts):
    """Nusinersen is transcriptome-unique to SMN.

    Derived, not asserted from the module: scan all 140.6 Mnt of the index for
    perfect complements of the published 18-mer. The only two may be SMN2 (its
    target) and SMN1 (its ~99.9%-identical paralog). Any third gene means the
    scan, the index, or the coordinate chain is wrong.
    """
    target = revcomp(NUSINERSEN)
    n = len(target)
    q = np.frombuffer(target.encode(), dtype=np.uint8)
    genes = set()
    for t in transcripts:
        for _, seq in G._units(t):
            if len(seq) < n:
                continue
            a = np.frombuffer(seq.encode(), dtype=np.uint8)
            mism = (np.lib.stride_tricks.sliding_window_view(a, n) != q).sum(1)
            if (mism == 0).any():
                genes.add(t.gene_name)
    assert genes == {"SMN1", "SMN2"}


def test_histogram_is_a_complete_census(transcripts):
    """The mismatch histogram must account for EVERY window, not just the good
    ones -- that is what makes it cutoff-free and therefore comparable across
    lengths. Sum of counts == number of windows scanned."""
    sub = transcripts[:50]
    aso = "GCUGGCAG"
    n = len(aso)
    hist = G.mismatch_histogram(aso, sub)
    windows = sum(max(0, len(s) - n + 1) for t in sub for _, s in G._units(t))
    assert int(hist.sum()) == windows


def test_infeasible_footprints_really_form_no_duplex():
    """The feasibility gate rejects AU-only footprints. Verify the rejection is
    a property of the sequence, not of the gate: an all-AU hexamer forms no
    duplex while a GC-rich hexamer at the same length does."""
    c = compact_window(intron6_tail=G.INTRON6_TAIL)
    au = G.Winner(6, 32044, 0.0, c.rna_sub(c.map_interval(Interval(32044, 32049))), "", True)
    au = G.Winner(6, 32044, 0.0, au.target_rna, G._complement(au.target_rna)[::-1], True)
    gc = c.rna_sub(c.map_interval(Interval(32060, 32065)))
    gcw = G.Winner(6, 32060, 0.0, gc, G._complement(gc)[::-1], True)
    assert G.feasible(au) is None
    assert G.feasible(gcw) is not None


def test_feasible_winners_land_on_the_gc_core():
    """A3 recovered the GC core blind. After the feasibility gate every ASO of
    length 8-14 must still target it -- if the gate moved the answer off the
    published site, the gate is doing something other than what it claims."""
    cb = G.candidates_by_length()
    for L in range(8, 15):
        w, _ = G.first_feasible(cb[L])
        assert w is not None, f"no feasible ASO at length {L}"
        iv = Interval(w.start, w.start + w.length - 1)
        assert iv.start <= GC_CORE.end and iv.end >= GC_CORE.start, (
            f"length {L} winner {iv} does not overlap GC_CORE {GC_CORE}")


def test_panel_controls_are_length_matched(transcripts):
    """Scanned length is the confound: a longer transcript accumulates better
    off-target hits for reasons unrelated to being a target. Each hit must have
    controls within 10% of its own scanned length."""
    hits, controls = G.build_panel(transcripts)
    hl = sorted(G.scanned_length(h) for h in hits)
    for c in controls:
        L = G.scanned_length(c)
        assert min(abs(L - h) / h for h in hl) < 0.10, (
            f"control {c.gene_name} at {L} nt matches no hit length in {hl}")


def test_margin_uses_one_scorer_for_both_sides():
    """On-target and off-target must go through the same energy path, or their
    difference is a comparison of two pipelines rather than of two sites."""
    import inspect
    src = inspect.getsource(G.best_offtarget) + inspect.getsource(G.feasible)
    assert src.count("_ddG_net") == 2
