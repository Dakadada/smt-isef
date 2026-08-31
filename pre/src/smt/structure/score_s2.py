"""S2: is the accessibility ceiling visible to a model that has never heard of it?

THE CLAIM UNDER TEST. A7's headline finding is that three of four cryptic
targets sit AT their accessibility ceiling -- the site cannot be occluded much
further because it is already largely paired -- and that the 6.5e5-fold spread
in achievable effect across targets is baseline accessibility, not footprint
quality. That claim rests entirely on `A_5ss`, a ViennaRNA quantity.

Boltz-2 is a learned potential with no `max_bp_span` and no nearest-neighbour
table. If the four targets' UNBOUND 3D structures rank in the same order as
their `A_5ss` baselines, the ceiling is a property of the molecules. If they do
not, that is a falsification signal and it goes in the record either way.

NO NEW GPU TIME. The four unbound predictions (`unbound_{target}`) were already
produced as controls in the S1 run, five diffusion samples each. This module
only scores them.

TWO CORRECTIONS TO THE PLAN, both measured rather than argued
--------------------------------------------------------------
1. The plan's S2 text says the expected order is "STMN2 most exposed,
   UNC13A_DONOR most buried". That is WRONG. By the project's own standing
   worst-span rule the order is

       STMN2_3SS 7.109e-03 > ACCEPTOR_128 9.669e-05
                            > DONOR 1.899e-05 > ACCEPTOR_178 3.446e-06

   so ACCEPTOR_178 is the most buried target and DONOR is third. The
   expectation is taken from the data, not from the prose.

2. The baseline rank order is ITSELF span-dependent -- at span 50 DONOR is the
   MOST accessible target, at span 150 ACCEPTOR_178 is. There is no single
   "the" A_5ss ordering to compare against. This module therefore pins the
   reference to the WORST span, the same rule that governs every other ranking
   in the project, and reports the per-span orderings alongside so the reader
   can see how much the reference itself moves. A test whose reference standard
   wobbles is a weaker test, and saying so is not optional.

POWER, STATED BEFORE THE RESULT. There are four targets. Under an exact
permutation test over all 4! = 24 orderings, a PERFECT rank match gives a
one-sided p of 1/24 = 0.042 and nothing can ever do better. S2 is therefore a
direction check with essentially no power to detect partial agreement, and it
is reported as such -- a corroboration or a falsification, never a p-value to
lean on.
"""

from __future__ import annotations

import argparse
import csv
import glob
import statistics
from itertools import permutations
from pathlib import Path

from smt.seq.reference import Interval
from smt.structure.gate import base_pairs
from smt.structure.score_s1 import S_DIR, target_windows
from smt.validate.cryptic_sweep import TARGETS

SWEEP_DIR = Path("/home/dashii/orcd/pool/smt-isef-data/cryptic_sweep")
PANEL_CSV = Path("/home/dashii/smt-isef/docs/prereg-panel.csv")


def baselines(sweep_dir: Path = SWEEP_DIR) -> dict[str, dict[int, float]]:
    """target -> {span: A_baseline}, re-derived from the A7 shards.

    The baseline depends only on (construct, estimator, feature, span) and never
    on the footprint, so every row for one (target, span) must carry the same
    value. That invariant is ASSERTED here rather than assumed: a non-constant
    baseline means baseline hoisting broke and every delta in the sweep is
    suspect.
    """
    seen: dict[str, dict[int, set]] = {}
    for p in sorted(glob.glob(str(sweep_dir / "*.csv"))):
        with open(p, newline="") as fh:
            for r in csv.DictReader(fh):
                if not r.get("A_baseline"):
                    continue
                t, s = r["target"], int(r["max_bp_span"])
                seen.setdefault(t, {}).setdefault(s, set()).add(
                    round(float(r["A_baseline"]), 12))
    out: dict[str, dict[int, float]] = {}
    for t, spans in seen.items():
        out[t] = {}
        for s, vals in spans.items():
            if len(vals) != 1:
                raise ValueError(
                    f"{t} span {s} has {len(vals)} distinct baselines {sorted(vals)[:4]}; "
                    "baseline hoisting is broken, do not trust any delta")
            out[t][s] = next(iter(vals))
    return out


def worst_span_baseline(b: dict[str, dict[int, float]]) -> dict[str, float]:
    return {t: min(v.values()) for t, v in b.items()}


def unbound_scores(out_dir: Path, panel_csv: Path = PANEL_CSV) -> dict[str, list[float]]:
    """target -> per-seed fraction of splice-site nucleotides base-paired.

    Same geometric pair detector and same window/site reconstruction as S1, so
    this number is the one whose accuracy was measured at F1 0.971 on the gate.
    """
    wins = target_windows(panel_csv)
    out: dict[str, list[float]] = {}
    for pdir in sorted(out_dir.rglob("predictions/unbound_*")):
        if not pdir.is_dir():
            continue
        tk = pdir.name[len("unbound_"):]
        if tk not in TARGETS:
            raise KeyError(f"{pdir.name} does not name a known target")
        window, site = wins[tk]
        lo = site.start - window.start
        vals = []
        for m in sorted(pdir.glob(f"{pdir.name}_model_*.cif")):
            partners = {i for p in base_pairs(m) for i in p}
            vals.append(sum(1 for i in range(lo, lo + len(site)) if i in partners)
                        / len(site))
        if vals:
            out[tk] = vals
    return out


def spearman(x: list[float], y: list[float]) -> float:
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        for pos, i in enumerate(order):
            r[i] = pos
        return r
    rx, ry = rank(x), rank(y)
    n = len(x)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = sum((a - mx) ** 2 for a in rx) ** 0.5
    dy = sum((b - my) ** 2 for b in ry) ** 0.5
    return num / (dx * dy) if dx and dy else 0.0


def exact_p(x: list[float], y: list[float]) -> tuple[float, float]:
    """(rho, one-sided exact permutation p) for rho <= observed.

    Direction is pre-specified: MORE baseline accessibility should mean LESS
    base-pairing at the site in 3D, i.e. a NEGATIVE correlation. Exhaustive
    over all n! relabelings -- n is 4, so there is no reason to approximate.
    """
    obs = spearman(x, y)
    perms = list(permutations(range(len(y))))
    hits = sum(1 for p in perms if spearman(x, [y[i] for i in p]) <= obs)
    return obs, hits / len(perms)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=S_DIR / "s1_out")
    ap.add_argument("--sweep", type=Path, default=SWEEP_DIR)
    ap.add_argument("--panel", type=Path, default=PANEL_CSV)
    a = ap.parse_args()

    b = baselines(a.sweep)
    worst = worst_span_baseline(b)
    scores = unbound_scores(a.out, a.panel)

    print("## S2 -- accessibility ceiling, checked against an orthogonal 3D model\n")

    print("### The reference standard moves with span (reported, not hidden)")
    spans = sorted(next(iter(b.values())).keys())
    print(f"  {'span':>6}  most accessible -> least accessible")
    for s in spans:
        order = sorted(b, key=lambda t: b[t][s], reverse=True)
        print(f"  {s:>6}  {' > '.join(order)}")
    print(f"  {'WORST':>6}  {' > '.join(sorted(worst, key=worst.get, reverse=True))}"
          "   <-- the reference used, per the standing worst-span rule\n")

    missing = [t for t in worst if t not in scores]
    if missing:
        print(f"  MISSING unbound predictions for {missing}; S2 cannot be graded.")
        return

    print("### Measured")
    print(f"  {'target':<24}{'A_5ss (worst span)':>20}{'3D paired frac':>16}"
          f"{'seed spread':>13}{'n':>4}")
    ts = sorted(worst, key=worst.get, reverse=True)
    for t in ts:
        v = scores[t]
        sd = statistics.stdev(v) if len(v) > 1 else 0.0
        print(f"  {t:<24}{worst[t]:>20.3e}{statistics.mean(v):>16.3f}"
              f"{sd:>13.3f}{len(v):>4}")

    x = [worst[t] for t in ts]
    y = [statistics.mean(scores[t]) for t in ts]
    rho, p = exact_p(x, y)
    print(f"\n### Result\n  Spearman rho = {rho:+.3f}, one-sided exact p = {p:.4f} "
          f"(best attainable {1/24:.4f})")

    pred = ts
    obs = sorted(scores, key=lambda t: statistics.mean(scores[t]))
    print(f"  predicted order (most exposed first): {' > '.join(pred)}")
    print(f"  observed 3D order  (least paired first): {' > '.join(obs)}")

    if rho <= -0.9:
        print("\n  CORROBORATED: the 3D model reproduces the A_5ss baseline ordering.")
        print("  The ceiling is a property of the molecules, not of ViennaRNA.")
    elif rho >= 0.5:
        print("\n  FALSIFIED: the orthogonal model orders these targets the OTHER way.")
        print("  A7's ceiling interpretation does not survive an independent instrument.")
    else:
        print("\n  INCONCLUSIVE: no clear agreement or disagreement at n = 4.")
    print("\n  n = 4 targets. This is a direction check, not a powered test; a")
    print("  perfect match could not have produced p < 0.042. Do not report the")
    print("  p-value without the n.")


if __name__ == "__main__":
    main()
