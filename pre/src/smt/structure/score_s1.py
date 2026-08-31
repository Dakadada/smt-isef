"""Grade S1: does an orthogonal 3D model separate span-robust from span-flipping?

THE TEST, as pre-registered in docs/preregistration.md (P1) and restated here
because S1 is the same hypothesis on a different instrument:

    Pairs are MATCHED on best-span |fom| and on length, so a single-span
    analysis ranks the two members equally and predicts no difference BY
    CONSTRUCTION. The framework predicts the span-robust member occludes the
    splice site more. Boltz-2 has no concept of max_bp_span at all, so if it
    separates them, span-robustness is a property of the molecule and not of
    ViennaRNA's parameterization.

THE MEASURE is the fraction of splice-site nucleotides that are base-paired in
the predicted structure -- the 3D analogue of 1 - A_5ss. Pairs are detected by
the SAME geometric code that graded the gate (smt.structure.gate.base_pairs),
so the number being compared here is the one whose accuracy was measured at
F1 0.971 against deposited structures.

SEED SPREAD GATES EVERY PAIR. Five diffusion samples per input give a spread;
a pair whose two members differ by less than their own pooled spread is
recorded as UNRESOLVED and excluded from the test, because that difference is
not distinguishable from sampling noise. Excluding them is declared here,
before the numbers, and the count of exclusions is always reported.
"""

from __future__ import annotations

import argparse
import csv
import statistics
from pathlib import Path

from smt.seq.reference import Interval
from smt.structure.gate import base_pairs
from smt.structure.inputs import WINDOW_PAD
from smt.validate.cryptic_sweep import TARGETS

S_DIR = Path("/home/dashii/orcd/pool/smt-isef-data/boltz_s")


def target_windows(panel_csv: Path) -> dict[str, tuple[Interval, Interval]]:
    """(window, site) in CONSTRUCT space, rebuilt exactly as inputs.py did.

    Recomputed from the same panel and the same WINDOW_PAD rather than stored,
    so the two can never silently disagree about where the splice site sits.
    """
    rows = list(csv.DictReader(panel_csv.open()))
    out: dict[str, tuple[Interval, Interval]] = {}
    for tk in {r["target"] for r in rows}:
        c, t = TARGETS[tk].construct_fn(), TARGETS[tk]
        site = c.map_interval(t.site)
        fps = [c.map_interval(Interval(int(r["start"]), int(r["end"])))
               for r in rows if r["target"] == tk]
        lo = max(1, min([site.start] + [f.start for f in fps]) - WINDOW_PAD)
        hi = min(len(c.sequence), max([site.end] + [f.end for f in fps]) + WINDOW_PAD)
        out[tk] = (Interval(lo, hi), site)
    return out


def site_paired_fraction(model: Path, window: Interval, site: Interval) -> float:
    """Fraction of splice-site nucleotides that are base-paired, chain A only."""
    pairs = base_pairs(model)
    partners = {i for p in pairs for i in p}
    lo = site.start - window.start          # 0-based index of the site in chain A
    idx = range(lo, lo + len(site))
    return sum(1 for i in idx if i in partners) / len(site)


def collect(out_dir: Path, panel_csv: Path) -> dict[str, list[float]]:
    wins = target_windows(panel_csv)
    scores: dict[str, list[float]] = {}
    for part in sorted(out_dir.glob("part*")):
        for pdir in sorted(part.rglob("predictions/*")):
            if not pdir.is_dir():
                continue
            name = pdir.name
            tk = next((t for t in TARGETS if t in name), None)
            if tk is None:
                continue
            window, site = wins[tk]
            vals = [site_paired_fraction(m, window, site)
                    for m in sorted(pdir.glob(f"{name}_model_*.cif"))]
            if vals:
                scores[name] = vals
    return scores


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=S_DIR / "s1_out")
    ap.add_argument("--panel", type=Path, default=Path("docs/prereg-panel.csv"))
    args = ap.parse_args()

    scores = collect(args.out, args.panel)
    rows = list(csv.DictReader(args.panel.open()))
    print(f"## S1 -- splice-site occlusion in 3D. {len(scores)} predictions scored.\n")

    by_pair: dict[tuple[str, str], dict[str, str]] = {}
    for r in rows:
        by_pair.setdefault((r["target"], r["pair"]), {})[r["role"]] = (
            f"s1_{r['target']}_p{r['pair']}_{r['role']}_{r['start']}_{r['length']}")

    print(f"  {'target':<22}{'pair':>5}{'robust':>9}{'flip':>9}"
          f"{'diff':>9}{'spread':>9}  verdict")
    diffs, unresolved = [], 0
    for (tk, pair), roles in sorted(by_pair.items()):
        rn, fn = roles.get("ROBUST"), roles.get("FLIP")
        if rn not in scores or fn not in scores:
            continue
        rv, fv = scores[rn], scores[fn]
        mr, mf = statistics.mean(rv), statistics.mean(fv)
        spread = max(statistics.pstdev(rv), statistics.pstdev(fv))
        diff = mr - mf
        if abs(diff) < spread:
            unresolved += 1
            verdict = "UNRESOLVED (< seed spread)"
        else:
            diffs.append(diff)
            verdict = "robust more occluded" if diff > 0 else "flip more occluded"
        print(f"  {tk:<22}{pair:>5}{mr:>9.3f}{mf:>9.3f}{diff:>+9.3f}"
              f"{spread:>9.3f}  {verdict}")

    print(f"\n  usable pairs {len(diffs)}, unresolved {unresolved}")
    if len(diffs) < 6:
        print("  FEWER THAN 6 USABLE PAIRS -- underpowered, per the decision rule")
        print("  fixed in docs/preregistration.md. Report the count, claim nothing.")
        return
    wins_ = sum(1 for d in diffs if d > 0)
    print(f"  robust more occluded in {wins_}/{len(diffs)} pairs; "
          f"mean difference {statistics.mean(diffs):+.4f}")
    print("  Primary readout is the graded paired test over these differences,")
    print("  one-sided in the direction P1 registered. A null result here puts")
    print("  P1's core premise in trouble and MUST be reported as such.")


if __name__ == "__main__":
    main()
