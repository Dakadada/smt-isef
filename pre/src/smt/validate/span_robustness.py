"""Is the max_bp_span sign-flip a defect of the observable, or of near-zero noise?

THE ATTACK THIS ANSWERS. 87.0% of A7 footprints flip the sign of their
predicted effect across the ``max_bp_span`` ladder. If the sign of a
prediction depends on a computational parameter, the observable is not
well-defined and no single-span number means anything. That is close to a
falsification, not a caveat, and it deserves a test rather than a disclaimer.

TWO HYPOTHESES, both tested here, and the first one LOSES:

  H-noise  Flipping is concentrated among near-zero effects -- footprints
           that barely change anything, whose sign is numerical noise about
           zero. If true, the 87% is a statement about noise and the strong
           predictions are safe.
           >>> REFUTED. The flip rate is 85.7% even among footprints with
           |fom| >= 1e-2, versus 86.9% at 1e-3 and 96.7% at 1e-4. Magnitude
           buys essentially nothing. The problem is real.

  H-struct The span-ROBUST minority is not a random leftover but marks
           specific positions whose control of splice-site accessibility
           survives every folding regime. If true, the framework still makes
           well-defined predictions -- just far fewer of them.
           >>> SUPPORTED. Robust footprints are significantly clustered in
           start position at all four targets (permutation test below).

WHAT THIS BUYS, stated honestly: it does NOT rescue the 87%. It narrows the
framework's scope to the span-robust subset (4-18% of footprints depending
on target) and shows that subset is structured rather than arbitrary. Any
A7 claim must be restricted to it, and any single-span number remains
meaningless.
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

from smt.validate.cryptic_report import DEFAULT_DIR, SPANS, load

#: Declared before running. 2000 gives a permutation-p floor of 1/2001.
N_PERMUTATIONS = 2000
SEED = 0

#: Magnitude bins for H-noise, declared up front and all reported.
BINS = ((0, 1e-4), (1e-4, 1e-3), (1e-3, 1e-2), (1e-2, float("inf")))


def mean_pairwise_distance(xs: list[int]) -> float:
    """Mean |xi - xj| over all pairs, in O(n log n). Lower == more clustered."""
    xs = sorted(xs)
    total, prefix = 0, 0
    for i, x in enumerate(xs):
        total += i * x - prefix
        prefix += x
    n = len(xs)
    return total / (n * (n - 1) / 2)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DEFAULT_DIR)
    ap.add_argument("--permutations", type=int, default=N_PERMUTATIONS)
    args = ap.parse_args()

    prints, _ = load(args.out)
    full = [p for p in prints if len(p.fom) == len(SPANS)]
    print(f"## span robustness, n = {len(full)} full-ladder footprints\n")

    print("### H-noise: does flipping concentrate among weak effects?")
    print(f"  {'max |fom| bin':>22}{'n':>7}{'flip rate':>11}")
    for lo, hi in BINS:
        group = [p for p in full if lo <= max(abs(v) for v in p.fom.values()) < hi]
        if not group:
            continue
        flip = sum(1 for p in group if not p.consistent)
        print(f"  {f'[{lo:.0e}, {hi:.0e})':>22}{len(group):>7}{100*flip/len(group):>10.1f}%")
    print("  H-noise is REFUTED if the rate is flat across bins. Read the column.")

    print(f"\n### H-struct: are span-robust footprints positionally clustered?")
    print(f"  Permutation test, {args.permutations} draws, seed {SEED}: the robust")
    print("  set's mean pairwise start distance vs. equal-sized samples drawn")
    print("  from the SAME target's scanned footprints. Lower == more clustered.")
    rng = random.Random(SEED)
    print(f"  {'target':<22}{'robust':>8}{'of':>7}{'observed':>10}{'null med':>10}{'p':>9}")
    for target in sorted({p.target for p in full}):
        group = [p for p in full if p.target == target]
        robust = [p for p in group if p.desired_everywhere]
        if len(robust) < 10:
            print(f"  {target:<22}{len(robust):>8}{len(group):>7}"
                  f"{'too few to test':>29}")
            continue
        starts = [p.start for p in robust]
        pool = [p.start for p in group]
        observed = mean_pairwise_distance(starts)
        null = sorted(mean_pairwise_distance(rng.sample(pool, len(starts)))
                      for _ in range(args.permutations))
        pval = (sum(1 for v in null if v <= observed) + 1) / (args.permutations + 1)
        print(f"  {target:<22}{len(robust):>8}{len(group):>7}{observed:>10.1f}"
              f"{null[len(null)//2]:>10.1f}{pval:>9.4f}")
    print("\n  Clustering does NOT rescue the 87%. It shows the robust minority")
    print("  marks specific positions rather than being a random leftover, so")
    print("  A7 claims are restricted to that minority and no further.")


if __name__ == "__main__":
    main()
