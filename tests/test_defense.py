"""Tests for the three modules written to answer the attacks on this project.

The statistics here are load-bearing in a way the rest of the codebase is
not: a wrong ``fisher_exact`` would not crash, it would silently produce a
p-value, and a p-value is exactly the thing a reviewer will not re-derive.
So it is checked against values computed independently (scipy where
available, brute-force enumeration otherwise) rather than against itself.
"""

from __future__ import annotations

import itertools
import random

import pytest

from smt.sieve.enrichment import fisher_exact
from smt.validate.prereg_panel import (
    MATCH_TOL, MIN_MEMBER_SEPARATION, MIN_PAIR_SEPARATION, complement,
)
from smt.validate.span_robustness import mean_pairwise_distance


# ------------------------------------------------------------------ Fisher


@pytest.mark.parametrize("table,expected", [
    # Fisher's own tea-tasting experiment: the canonical worked example.
    ((3, 1, 1, 3), 0.4857142857),
    ((4, 0, 0, 4), 0.02857142857),
    # A5's actual observed table: exonic_proximal 2/11 vs background 11/610.
    # Value from scipy.stats.fisher_exact, not from this module's own rounded
    # output -- an expectation transcribed from the report would only prove
    # the report matches itself.
    ((2, 9, 11, 599), 0.020020749859865457),
])
def test_fisher_matches_known_values(table, expected):
    assert fisher_exact(*table) == pytest.approx(expected, rel=1e-6)


def test_fisher_is_symmetric_under_transpose():
    """p([[a,b],[c,d]]) == p([[a,c],[b,d]]) for the two-sided exact test."""
    for a, b, c, d in [(2, 9, 11, 599), (3, 1, 1, 3), (5, 5, 5, 5)]:
        assert fisher_exact(a, b, c, d) == pytest.approx(fisher_exact(a, c, b, d))


def test_fisher_returns_one_for_no_association():
    assert fisher_exact(5, 5, 5, 5) == pytest.approx(1.0)


def test_fisher_never_exceeds_one():
    """Floating-point summation over the tail must not drift above 1."""
    rng = random.Random(1)
    for _ in range(200):
        t = [rng.randint(0, 30) for _ in range(4)]
        if sum(t) == 0:
            continue
        assert 0.0 <= fisher_exact(*t) <= 1.0


def test_fisher_agrees_with_scipy_when_available():
    scipy_stats = pytest.importorskip("scipy.stats")
    rng = random.Random(7)
    for _ in range(50):
        a, b, c, d = (rng.randint(0, 25) for _ in range(4))
        if (a + b) == 0 or (c + d) == 0 or (a + c) == 0 or (b + d) == 0:
            continue
        _, expected = scipy_stats.fisher_exact([[a, b], [c, d]])
        assert fisher_exact(a, b, c, d) == pytest.approx(expected, rel=1e-9)


# ------------------------------------------------- mean pairwise distance


def test_mean_pairwise_distance_matches_brute_force():
    """The O(n log n) form must equal the O(n^2) definition exactly."""
    rng = random.Random(3)
    for _ in range(50):
        xs = [rng.randint(0, 1000) for _ in range(rng.randint(2, 40))]
        brute = sum(abs(i - j) for i, j in itertools.combinations(xs, 2))
        brute /= len(xs) * (len(xs) - 1) / 2
        assert mean_pairwise_distance(xs) == pytest.approx(brute)


def test_clustered_points_score_below_spread_points():
    """Lower == more clustered, which is the direction the p-value assumes."""
    assert mean_pairwise_distance([10, 11, 12, 13]) < mean_pairwise_distance([0, 100, 200, 300])


# ------------------------------------------------------------ ASO panel


def test_complement_is_rna_and_involutive_on_rna():
    assert complement("ACGU") == "UGCA"
    assert complement("ACGT") == "UGCA", "DNA input must still yield RNA"
    assert complement(complement("ACGU")) == "ACGU"


def test_reverse_complement_of_a_known_duplex():
    """5'-GAAGGUCUUC-3' is the panel's top ASO; it must pair with its target."""
    target = "GAAGACCUUC"
    assert complement(target)[::-1] == "GAAGGUCUUC"


def test_panel_constraints_are_ordered_sensibly():
    """Pair separation must exceed member separation, or 'independent pairs'
    would be a weaker constraint than 'distinguishable oligos'."""
    assert MIN_PAIR_SEPARATION > MIN_MEMBER_SEPARATION > 0
    assert 0 < MATCH_TOL < 1
