"""Exact dinucleotide-preserving sequence shuffles via randomized Euler trails."""

from __future__ import annotations

from collections import defaultdict
import random

from .sequences import to_rna


def _one_euler_shuffle(sequence: str, rng: random.Random) -> str:
    outgoing: dict[str, list[str]] = defaultdict(list)
    for left, right in zip(sequence, sequence[1:]):
        outgoing[left].append(right)
    for edges in outgoing.values():
        rng.shuffle(edges)

    stack = [sequence[0]]
    trail: list[str] = []
    while stack:
        node = stack[-1]
        if outgoing[node]:
            stack.append(outgoing[node].pop())
        else:
            trail.append(stack.pop())
    shuffled = "".join(reversed(trail))
    if len(shuffled) != len(sequence):
        raise RuntimeError("Euler shuffle failed to consume every dinucleotide")
    return shuffled


def dinucleotide_shuffles(
    sequence: str, *, n: int = 100, seed: int = 20260926
) -> list[str]:
    """Return ``n`` seeded draws with the exact directed dinucleotide multiset.

    Draws are with replacement. That is essential for short ASOs whose number
    of possible exact dinucleotide shuffles is less than 100; duplicates retain
    their null-sampling multiplicity and are deduplicated only for computation.
    """
    normalized = to_rna(sequence)
    if len(normalized) < 2:
        raise ValueError("dinucleotide shuffling requires at least two bases")
    if n < 1:
        raise ValueError("n must be positive")
    rng = random.Random(seed)
    return [_one_euler_shuffle(normalized, rng) for _ in range(n)]

