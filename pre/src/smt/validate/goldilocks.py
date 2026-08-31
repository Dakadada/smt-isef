"""G: the Goldilocks length -- most on-target effect, least off-target affinity.

Answers the spec bullet "determine the precise oligonucleotide length that
unravels TSL2 and ISTL1 without inducing off-target affinity", and needs no
new on-target compute: the A3 sweep already scored every (start, length, span)
footprint across ISS-N1 and its output is on disk.

DIRECTION. SMN2 exon 7 wants MORE 5'ss accessibility, so this is the FREE
regime and a footprint is desirable when delta_A > 0. (The A3 sweep predates
the signed `fom` field, so the sign convention is applied here explicitly
rather than read from a column that does not exist.)

RANKED ON THE WORST SPAN, never the best -- the standing project rule. A
footprint that only helps at one max_bp_span is not a candidate, and 87% of
footprints flip sign across the ladder (see A7).

HONEST EXPECTATION, stated before the numbers: shorter ASOs may win on
off-target simply because shorter seeds match less often. If the trade-off
curve is monotonic in length with no interior optimum, THERE IS NO GOLDILOCKS
POINT and that is the result -- it does not get bins tuned until one appears.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

SWEEP_DIR = Path("/home/dashii/orcd/pool/smt-isef-data/sweep")
SPANS = (0, 50, 100, 150, 300)

#: The published answer A3 recovered blind. Used ONLY to verify this loader
#: reproduces it -- never to select or rank anything.
KNOWN_8MER = (32059, 8)

#: Scan positions at either boundary to exclude as edge artifacts.
#:
#: MEASURED, NOT ASSUMED. Every length 6-10 had its best span-robust footprint
#: at exactly 32040 -- the first scanned position -- and dropping the first
#: three positions cuts those winners by 2-3 orders of magnitude (7.9e-03 ->
#: 9.1e-06 at length 8). A footprint at the boundary sits against the end of
#: the folding window, so its accessibility change is partly an artifact of
#: where the window stops rather than of the molecule. Both views are reported.
EDGE_MARGIN = 3


def load(sweep_dir: Path, nshards: int = 200) -> dict[tuple[int, int], dict[int, float]]:
    """(start, length) -> {span: delta_A}, from the A3 shards."""
    out: dict[tuple[int, int], dict[int, float]] = defaultdict(dict)
    missing = []
    for shard in range(nshards):
        p = sweep_dir / f"sweep_{shard:04d}_of_{nshards:04d}.csv"
        if not p.exists():
            missing.append(shard)
            continue
        with p.open(newline="") as fh:
            for row in csv.DictReader(fh):
                if not row.get("delta_A"):
                    continue
                out[(int(row["start"]), int(row["length"]))][int(row["max_bp_span"])] = \
                    float(row["delta_A"])
    if missing:
        raise FileNotFoundError(f"{len(missing)} A3 shards missing: {missing[:5]}")
    return out


def worst_span(ladder: dict[int, float]) -> float:
    return min(ladder.values())


def robust(ladder: dict[int, float]) -> bool:
    """Desirable (FREE regime: delta_A > 0) at EVERY span in the ladder."""
    return len(ladder) == len(SPANS) and all(v > 0 for v in ladder.values())


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", type=Path, default=SWEEP_DIR)
    ap.add_argument("--top", type=int, default=3)
    args = ap.parse_args()

    data = load(args.sweep)
    print(f"## G -- on-target axis from the A3 sweep. {len(data)} footprints.\n")

    by_len: dict[int, list] = defaultdict(list)
    for (start, length), ladder in data.items():
        if robust(ladder):
            by_len[length].append((worst_span(ladder), start))

    lo = min(s for s, _ in data)
    print(f"  {'len':>4}{'robust':>9}{'best (all)':>16}{'at':>8}"
          f"{'best (no edge)':>18}{'at':>8}")
    for length in sorted(by_len):
        cands = sorted(by_len[length], reverse=True)
        inner = [c for c in cands if c[1] > lo + EDGE_MARGIN - 1]
        edge_s = f"{inner[0][0]:.3e}" if inner else "NONE"
        edge_p = f"{inner[0][1]}" if inner else "-"
        print(f"  {length:>4}{len(cands):>9}{cands[0][0]:>16.3e}{cands[0][1]:>8}"
              f"{edge_s:>18}{edge_p:>8}")
    zero = [l for l in range(6, 21) if l not in by_len]
    if zero:
        print(f"\n  LENGTHS WITH ZERO SPAN-ROBUST FOOTPRINTS: {zero}")
        print("  Nusinersen is an 18-mer. If 18 is in that list, the clinically")
        print("  effective length has NO span-robust footprint in this scan --")
        print("  consistent with long ASOs acting by DIRECT occupancy, which this")
        print("  framework explicitly refuses to model. Report it that way, not")
        print("  as evidence that nusinersen should be shorter.")

    print("\n### VERIFICATION -- does this loader reproduce A3's blind answer?")
    ks, kl = KNOWN_8MER
    ranked8 = [s for _, s in sorted(by_len.get(kl, []), reverse=True)]
    if ks in ranked8:
        r = ranked8.index(ks) + 1
        print(f"  the published 8-mer at {ks} ranks {r}/{len(ranked8)} among "
              f"span-robust 8-mers.")
        print("  PASS" if r <= 5 else "  WARNING: expected it at or near the top.")
    else:
        print(f"  FAIL: {ks}+{kl} is not span-robust here. The loader or the")
        print("  sign convention is wrong; do not trust anything above.")

    print("\n  Lengths are ranked on their WORST span. Off-target scoring (G1,")
    print("  targeted panel) attaches to these winners; until it runs, this is")
    print("  the on-target axis ALONE and implies nothing about a trade-off.")


if __name__ == "__main__":
    main()
