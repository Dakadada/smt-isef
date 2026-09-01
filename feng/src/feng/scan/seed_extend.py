"""Complete seed-and-extend candidate scanning for the clean-room Feng arm.

The scanner reports sequence-similarity candidates only. It does not estimate
hybridization thermodynamics or splice activity. Coordinates are zero-based,
half-open offsets into each supplied reference sequence, and the caller must supply
the query in the orientation expected to match the reference.

The default seed policy partitions a query into ``budget + 1`` nonempty, balanced
blocks. An alignment with at most ``budget`` substitutions (Hamming mode) or edit
operations (edit mode) must leave at least one block exactly matched by the
pigeonhole principle, so extending every occurrence of every block cannot miss an
eligible interval. When a nonempty exact seed cannot be guaranteed, the scanner
falls back to exhaustive extension rather than sacrificing completeness.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum


_CANONICAL_BASES = frozenset("ACGT")
_IUPAC_BASES = frozenset("ACGTRYSWKMBDHVN")


class SeedLengthPolicy(str, Enum):
    """Completeness-preserving policies for selecting candidate intervals."""

    BALANCED = "balanced"
    EXHAUSTIVE = "exhaustive"


@dataclass(frozen=True, order=True, slots=True)
class HybridizationHit:
    """One sequence candidate in a named reference.

    ``mismatches`` is the Hamming count in fixed-window mode. In edit mode it is
    the substitution count in the minimum-edit alignment, with fewer substitutions
    used as the deterministic tie-break between equally short alignments.
    """

    reference_name: str
    start: int
    end: int
    mismatches: int
    edit_distance: int


def seed_and_extend(
    query: str,
    reference: str | Mapping[str, str],
    *,
    max_mismatch: int | None = None,
    max_edit: int | None = None,
    seed_policy: SeedLengthPolicy | str = SeedLengthPolicy.BALANCED,
    reference_name: str = "reference",
) -> tuple[HybridizationHit, ...]:
    """Return every reference interval within the declared sequence budget.

    Args:
        query: Nonempty oligo sequence in reference-matching orientation. ``U`` and
            ``T`` compare as equivalent. IUPAC ambiguity symbols are valid input but
            count as mismatches rather than as wildcards.
        reference: One sequence or a mapping of reference names to sequences.
        max_mismatch: Maximum substitutions. When ``max_edit`` is absent this selects
            fixed-length Hamming mode and defaults to zero. In edit mode it is an
            optional cap on substitutions in the optimal alignment.
        max_edit: When provided, select variable-length edit mode and require
            Levenshtein distance at most this value.
        seed_policy: ``balanced`` uses the pigeonhole-complete partition described in
            the module docstring; ``exhaustive`` disables seeding for audit use.
        reference_name: Name assigned when ``reference`` is a single sequence.

    Returns:
        A deterministically sorted tuple of unique :class:`HybridizationHit` values.

    Raises:
        TypeError: If sequences, names, or budgets have invalid types.
        ValueError: If the query is empty, a sequence contains a non-IUPAC symbol, a
            budget is negative, a reference name is empty, or the policy is unknown.

    Hamming mode evaluates query-length windows. Edit mode evaluates every nonempty
    window whose length differs from the query by at most ``max_edit``. Seeding only
    prunes which of those windows require extension; it never changes eligibility.
    """
    normalized_query = _normalize_sequence(query, label="query", allow_empty=False)
    references = _normalize_references(reference, reference_name)
    mismatch_budget = _validate_budget(max_mismatch, "max_mismatch")
    edit_budget = _validate_budget(max_edit, "max_edit")
    try:
        policy = SeedLengthPolicy(seed_policy)
    except (TypeError, ValueError) as error:
        choices = ", ".join(item.value for item in SeedLengthPolicy)
        raise ValueError(f"seed_policy must be one of: {choices}") from error

    if edit_budget is None:
        effective_mismatch_budget = 0 if mismatch_budget is None else mismatch_budget
        hits = _scan_hamming(
            normalized_query,
            references,
            effective_mismatch_budget,
            policy,
        )
    else:
        hits = _scan_edit(
            normalized_query,
            references,
            mismatch_budget,
            edit_budget,
            policy,
        )
    return tuple(sorted(hits))


def _normalize_sequence(sequence: str, *, label: str, allow_empty: bool) -> str:
    if not isinstance(sequence, str):
        raise TypeError(f"{label} must be a string")
    normalized = sequence.upper().replace("U", "T")
    if not normalized and not allow_empty:
        raise ValueError(f"{label} must not be empty")
    invalid = sorted(set(normalized) - _IUPAC_BASES)
    if invalid:
        raise ValueError(f"{label} contains non-IUPAC symbols: {''.join(invalid)}")
    return normalized


def _normalize_references(
    reference: str | Mapping[str, str],
    reference_name: str,
) -> tuple[tuple[str, str], ...]:
    if isinstance(reference, str):
        items = ((reference_name, reference),)
    elif isinstance(reference, Mapping):
        items = tuple(reference.items())
    else:
        raise TypeError("reference must be a sequence string or a name-to-sequence mapping")

    normalized: list[tuple[str, str]] = []
    for name, sequence in items:
        if not isinstance(name, str):
            raise TypeError("reference names must be strings")
        if not name:
            raise ValueError("reference names must not be empty")
        normalized.append(
            (name, _normalize_sequence(sequence, label=f"reference {name!r}", allow_empty=True))
        )
    return tuple(normalized)


def _validate_budget(value: int | None, label: str) -> int | None:
    if value is None:
        return None
    if type(value) is not int:
        raise TypeError(f"{label} must be an integer")
    if value < 0:
        raise ValueError(f"{label} must be nonnegative")
    return value


def _bases_match(left: str, right: str) -> bool:
    return left in _CANONICAL_BASES and right in _CANONICAL_BASES and left == right


def _hamming_distance(query: str, target: str) -> int:
    return sum(not _bases_match(left, right) for left, right in zip(query, target, strict=True))


def _alignment_cost(query: str, target: str) -> tuple[int, int]:
    """Return minimum edit distance and its minimum substitution tie-break."""
    previous = [(column, 0) for column in range(len(target) + 1)]
    for row, query_base in enumerate(query, start=1):
        current = [(row, 0)]
        for column, target_base in enumerate(target, start=1):
            deletion = (previous[column][0] + 1, previous[column][1])
            insertion = (current[column - 1][0] + 1, current[column - 1][1])
            if _bases_match(query_base, target_base):
                diagonal = previous[column - 1]
            else:
                diagonal = (previous[column - 1][0] + 1, previous[column - 1][1] + 1)
            current.append(min(deletion, insertion, diagonal))
        previous = current
    return previous[-1]


def _balanced_seeds(query: str, budget: int) -> tuple[tuple[int, str], ...]:
    if budget >= len(query):
        return ()
    block_count = budget + 1
    short_length, longer_blocks = divmod(len(query), block_count)
    offset = 0
    seeds: list[tuple[int, str]] = []
    for index in range(block_count):
        length = short_length + (index < longer_blocks)
        seeds.append((offset, query[offset : offset + length]))
        offset += length
    return tuple(seeds)


def _exact_occurrences(reference: str, seed: str) -> tuple[int, ...]:
    if any(base not in _CANONICAL_BASES for base in seed):
        return ()
    positions: list[int] = []
    position = reference.find(seed)
    while position != -1:
        positions.append(position)
        position = reference.find(seed, position + 1)
    return tuple(positions)


def _hamming_candidate_starts(
    query: str,
    reference: str,
    budget: int,
    policy: SeedLengthPolicy,
) -> set[int]:
    final_start = len(reference) - len(query)
    if final_start < 0:
        return set()
    seeds = () if policy is SeedLengthPolicy.EXHAUSTIVE else _balanced_seeds(query, budget)
    if not seeds:
        return set(range(final_start + 1))
    starts: set[int] = set()
    for query_offset, seed in seeds:
        for occurrence in _exact_occurrences(reference, seed):
            start = occurrence - query_offset
            if 0 <= start <= final_start:
                starts.add(start)
    return starts


def _scan_hamming(
    query: str,
    references: tuple[tuple[str, str], ...],
    budget: int,
    policy: SeedLengthPolicy,
) -> set[HybridizationHit]:
    hits: set[HybridizationHit] = set()
    for name, sequence in references:
        for start in _hamming_candidate_starts(query, sequence, budget, policy):
            end = start + len(query)
            mismatches = _hamming_distance(query, sequence[start:end])
            if mismatches <= budget:
                hits.add(HybridizationHit(name, start, end, mismatches, mismatches))
    return hits


def _edit_candidate_intervals(
    query: str,
    reference: str,
    budget: int,
    policy: SeedLengthPolicy,
) -> set[tuple[int, int]]:
    minimum_length = max(1, len(query) - budget)
    maximum_length = len(query) + budget
    seeds = () if policy is SeedLengthPolicy.EXHAUSTIVE else _balanced_seeds(query, budget)

    if not seeds:
        return {
            (start, start + length)
            for start in range(len(reference))
            for length in range(
                minimum_length,
                min(maximum_length, len(reference) - start) + 1,
            )
        }

    intervals: set[tuple[int, int]] = set()
    for query_offset, seed in seeds:
        for occurrence in _exact_occurrences(reference, seed):
            for start in range(
                occurrence - query_offset - budget,
                occurrence - query_offset + budget + 1,
            ):
                if start < 0 or start >= len(reference):
                    continue
                for length in range(minimum_length, maximum_length + 1):
                    end = start + length
                    if end <= len(reference):
                        intervals.add((start, end))
    return intervals


def _scan_edit(
    query: str,
    references: tuple[tuple[str, str], ...],
    mismatch_budget: int | None,
    edit_budget: int,
    policy: SeedLengthPolicy,
) -> set[HybridizationHit]:
    hits: set[HybridizationHit] = set()
    for name, sequence in references:
        intervals = _edit_candidate_intervals(query, sequence, edit_budget, policy)
        for start, end in intervals:
            edit_distance, mismatches = _alignment_cost(query, sequence[start:end])
            if edit_distance > edit_budget:
                continue
            if mismatch_budget is not None and mismatches > mismatch_budget:
                continue
            hits.add(HybridizationHit(name, start, end, mismatches, edit_distance))
    return hits
