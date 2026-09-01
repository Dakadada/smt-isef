"""Brute-force completeness tests for the clean-room Feng arm scanner."""

from __future__ import annotations

import ast
from itertools import product
from pathlib import Path
import re

import pytest

from feng.scan.seed_extend import HybridizationHit, SeedLengthPolicy, seed_and_extend


CANONICAL = frozenset("ACGT")
EXHAUSTIVE_ALPHABET = "ATUNR"
SCANNER_PATH = Path(__file__).resolve().parents[1] / "src" / "feng" / "scan" / "seed_extend.py"


def _normalize(sequence: str) -> str:
    return sequence.upper().translate(str.maketrans({"U": "T"}))


def _match(left: str, right: str) -> bool:
    return left in CANONICAL and right in CANONICAL and left == right


def _enumerated_alignment_cost(query: str, target: str) -> tuple[int, int]:
    """Score every global alignment path without using a DP recurrence."""
    query = _normalize(query)
    target = _normalize(target)
    costs: list[tuple[int, int]] = []

    def enumerate_paths(
        query_index: int,
        target_index: int,
        edits: int,
        substitutions: int,
    ) -> None:
        if query_index == len(query) and target_index == len(target):
            costs.append((edits, substitutions))
            return
        if query_index < len(query) and target_index < len(target):
            is_match = _match(query[query_index], target[target_index])
            enumerate_paths(
                query_index + 1,
                target_index + 1,
                edits + (not is_match),
                substitutions + (not is_match),
            )
        if query_index < len(query):
            enumerate_paths(query_index + 1, target_index, edits + 1, substitutions)
        if target_index < len(target):
            enumerate_paths(query_index, target_index + 1, edits + 1, substitutions)

    enumerate_paths(0, 0, 0, 0)
    return min(costs)


def _brute_hamming(
    query: str,
    references: dict[str, str],
    budget: int,
) -> set[HybridizationHit]:
    query = _normalize(query)
    hits: set[HybridizationHit] = set()
    for name, raw_sequence in references.items():
        sequence = _normalize(raw_sequence)
        for start in range(max(0, len(sequence) - len(query) + 1)):
            target = sequence[start : start + len(query)]
            mismatches = sum(not _match(left, right) for left, right in zip(query, target))
            if mismatches <= budget:
                hits.add(
                    HybridizationHit(
                        reference_name=name,
                        start=start,
                        end=start + len(query),
                        mismatches=mismatches,
                        edit_distance=mismatches,
                    )
                )
    return hits


def _brute_edit(
    query: str,
    references: dict[str, str],
    edit_budget: int,
    mismatch_budget: int | None,
) -> set[HybridizationHit]:
    query = _normalize(query)
    minimum_length = max(1, len(query) - edit_budget)
    maximum_length = len(query) + edit_budget
    hits: set[HybridizationHit] = set()
    for name, raw_sequence in references.items():
        sequence = _normalize(raw_sequence)
        for start in range(len(sequence)):
            for length in range(minimum_length, maximum_length + 1):
                end = start + length
                if end > len(sequence):
                    continue
                edit_distance, mismatches = _enumerated_alignment_cost(query, sequence[start:end])
                if edit_distance > edit_budget:
                    continue
                if mismatch_budget is not None and mismatches > mismatch_budget:
                    continue
                hits.add(HybridizationHit(name, start, end, mismatches, edit_distance))
    return hits


def _all_sequences(minimum_length: int, maximum_length: int) -> tuple[str, ...]:
    return tuple(
        "".join(bases)
        for length in range(minimum_length, maximum_length + 1)
        for bases in product(EXHAUSTIVE_ALPHABET, repeat=length)
    )


EXHAUSTIVE_QUERIES = _all_sequences(1, 2)
EXHAUSTIVE_REFERENCES = _all_sequences(0, 3)


def test_exhaustive_small_domain_hamming_completeness() -> None:
    for query in EXHAUSTIVE_QUERIES:
        for reference in EXHAUSTIVE_REFERENCES:
            references = {"toy": reference}
            for budget in range(len(query) + 1):
                expected = _brute_hamming(query, references, budget)
                balanced = set(seed_and_extend(query, references, max_mismatch=budget))
                exhaustive = set(
                    seed_and_extend(
                        query,
                        references,
                        max_mismatch=budget,
                        seed_policy=SeedLengthPolicy.EXHAUSTIVE,
                    )
                )
                assert balanced == exhaustive == expected, (query, reference, budget)


def test_exhaustive_small_domain_edit_completeness() -> None:
    for query in EXHAUSTIVE_QUERIES:
        for reference in EXHAUSTIVE_REFERENCES:
            references = {"toy": reference}
            for budget in range(len(query) + 1):
                expected = _brute_edit(query, references, budget, None)
                balanced = set(seed_and_extend(query, references, max_edit=budget))
                exhaustive = set(
                    seed_and_extend(
                        query,
                        references,
                        max_edit=budget,
                        seed_policy=SeedLengthPolicy.EXHAUSTIVE,
                    )
                )
                assert balanced == exhaustive == expected, (query, reference, budget)


@pytest.mark.parametrize(
    ("query", "target", "expected"),
    [
        ("AU", "AT", (0, 0)),
        ("AC", "AG", (1, 1)),
        ("AC", "A", (1, 0)),
        ("AC", "ATC", (1, 0)),
        ("AN", "AT", (1, 1)),
        ("R", "R", (1, 1)),
    ],
)
def test_independent_alignment_oracle_against_handwritten_truth(
    query: str,
    target: str,
    expected: tuple[int, int],
) -> None:
    assert _enumerated_alignment_cost(query, target) == expected


def test_huge_edit_budget_returns_only_real_intervals() -> None:
    expected = (HybridizationHit("reference", 0, 1, 0, 0),)
    assert seed_and_extend("A", "A", max_edit=1_000_000_000) == expected
    assert seed_and_extend(
        "A",
        "A",
        max_edit=1_000_000_000,
        seed_policy=SeedLengthPolicy.EXHAUSTIVE,
    ) == expected


def test_boundaries_overlaps_and_named_references() -> None:
    hits = seed_and_extend("AAA", {"left": "AAAAA", "right": "TAAA"}, max_mismatch=0)
    assert hits == (
        HybridizationHit("left", 0, 3, 0, 0),
        HybridizationHit("left", 1, 4, 0, 0),
        HybridizationHit("left", 2, 5, 0, 0),
        HybridizationHit("right", 1, 4, 0, 0),
    )


def test_query_longer_than_reference_and_empty_results() -> None:
    assert seed_and_extend("ACGT", "ACG", max_mismatch=4) == ()
    assert seed_and_extend("AAA", "CCC", max_mismatch=0) == ()
    assert seed_and_extend("AAA", "", max_edit=2) == ()
    assert seed_and_extend("ACGT", "ACG", max_edit=1) == (
        HybridizationHit("reference", 0, 3, 0, 1),
    )


def test_ambiguity_is_a_mismatch_not_a_wildcard_and_u_equals_t() -> None:
    assert seed_and_extend("AU", "AT", max_mismatch=0) == (
        HybridizationHit("reference", 0, 2, 0, 0),
    )
    assert seed_and_extend("AN", "ANAT", max_mismatch=1) == (
        HybridizationHit("reference", 0, 2, 1, 1),
        HybridizationHit("reference", 2, 4, 1, 1),
    )
    assert seed_and_extend("NN", "NN", max_mismatch=1) == ()


@pytest.mark.parametrize(
    ("arguments", "error_type"),
    [
        (("", "AC"), ValueError),
        (("AC-", "AC"), ValueError),
        (("AC", "AC"), ValueError),
        (("AC", "AC"), TypeError),
    ],
)
def test_invalid_inputs_are_rejected(arguments, error_type) -> None:
    query, reference = arguments
    if query == "AC" and error_type is ValueError:
        with pytest.raises(error_type):
            seed_and_extend(query, reference, max_mismatch=-1)
    elif query == "AC" and error_type is TypeError:
        with pytest.raises(error_type):
            seed_and_extend(query, reference, max_edit=True)
    else:
        with pytest.raises(error_type):
            seed_and_extend(query, reference)


def test_scanner_contains_no_locus_or_sequence_special_cases() -> None:
    source = SCANNER_PATH.read_text(encoding="utf-8")
    lowered = source.casefold()
    forbidden_loci = ("polr2h", "smn2", "nusinersen")
    assert not any(name in lowered for name in forbidden_loci)

    tree = ast.parse(source)
    strings = (
        node.value
        for node in ast.walk(tree)
        if isinstance(node, ast.Constant) and isinstance(node.value, str)
    )
    hard_coded_sequences = [
        value for value in strings if re.fullmatch(r"[ACGTUN]{8,}", value.upper())
    ]
    assert hard_coded_sequences == []
