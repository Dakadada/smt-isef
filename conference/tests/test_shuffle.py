from collections import Counter

from conference.aso_pipeline.shuffle import dinucleotide_shuffles


def _dinucleotides(sequence: str) -> Counter[str]:
    return Counter(sequence[index : index + 2] for index in range(len(sequence) - 1))


def test_dinucleotide_shuffles_are_seeded_and_preserve_exact_edge_multiset():
    first = dinucleotide_shuffles("AUAAUGCUGG", n=100, seed=20260926)
    second = dinucleotide_shuffles("AUAAUGCUGG", n=100, seed=20260926)
    assert first == second
    assert len(first) == 100
    assert all(len(value) == 10 for value in first)
    assert all(_dinucleotides(value) == _dinucleotides("AUAAUGCUGG") for value in first)
    assert len(set(first)) > 1


def test_short_null_keeps_duplicate_draws_instead_of_claiming_100_unique():
    draws = dinucleotide_shuffles("GCUGGCAG", n=100, seed=20260926)
    assert len(draws) == 100
    assert len(set(draws)) < 100

