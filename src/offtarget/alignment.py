"""Deterministic edit alignment used by the exhaustive scanner."""

from dataclasses import dataclass
from typing import List, Tuple


_COMP = str.maketrans("ACGTUN", "TGCAAN")


def normalize(sequence: str) -> str:
    sequence = sequence.upper().replace("U", "T")
    if not sequence or any(base not in "ACGTN" for base in sequence):
        raise ValueError("sequences must contain only A/C/G/T/U/N")
    return sequence


def reverse_complement(sequence: str) -> str:
    return normalize(sequence).translate(_COMP)[::-1]


@dataclass(frozen=True)
class Alignment:
    distance: int
    mismatches: int
    insertions: int
    deletions: int
    cigar: str


def global_edit_alignment(query: str, target: str) -> Alignment:
    """Levenshtein alignment with deterministic mismatch/indel traceback."""
    query, target = normalize(query), normalize(target)
    rows, cols = len(query) + 1, len(target) + 1
    dp: List[List[int]] = [[0] * cols for _ in range(rows)]
    back: List[List[Tuple[int, int, str]]] = [[(0, 0, "")] * cols for _ in range(rows)]
    for i in range(1, rows):
        dp[i][0], back[i][0] = i, (i - 1, 0, "D")
    for j in range(1, cols):
        dp[0][j], back[0][j] = j, (0, j - 1, "I")
    for i in range(1, rows):
        for j in range(1, cols):
            op = "M" if query[i - 1] == target[j - 1] else "X"
            choices = [
                (dp[i - 1][j - 1] + (op == "X"), 0, i - 1, j - 1, op),
                (dp[i - 1][j] + 1, 1, i - 1, j, "D"),
                (dp[i][j - 1] + 1, 2, i, j - 1, "I"),
            ]
            cost, _, pi, pj, chosen = min(choices)
            dp[i][j], back[i][j] = cost, (pi, pj, chosen)
    ops, i, j = [], len(query), len(target)
    while i or j:
        i, j, op = *back[i][j][:2], back[i][j][2]
        ops.append(op)
    ops.reverse()
    runs = []
    for op in ops:
        if runs and runs[-1][1] == op:
            runs[-1] = (runs[-1][0] + 1, op)
        else:
            runs.append((1, op))
    return Alignment(
        dp[-1][-1], ops.count("X"), ops.count("I"), ops.count("D"),
        "".join("{}{}".format(n, op) for n, op in runs),
    )
