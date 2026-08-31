"""G: the assembled Goldilocks trade-off -- on-target effect vs off-target risk.

Joins the on-target axis (`smt.validate.goldilocks`, from the A3 sweep) to the
off-target axis (`smt.validate.goldilocks_offtarget`) and states whether an
interior optimum exists. See `goldilocks_offtarget` for why a fixed absolute
mismatch budget is NOT a length-comparable criterion and what replaces it.

The headline quantity is the SPECIFICITY MARGIN

    margin = ddG_net(best panel off-target) - ddG_net(on-target)      kcal/mol

positive = the intended duplex is stronger than the best unintended one. A raw
off-target energy is not a specificity measure, because a longer ASO binds
everything more tightly including its real target.

WHAT COUNTS AS A GOLDILOCKS POINT, declared before looking: an INTERIOR maximum
in margin -- some length strictly better than both its neighbours and better
than the endpoints. A curve that rises or falls monotonically has no optimum,
and that is a reportable result, not a reason to re-bin.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from smt.validate.goldilocks_offtarget import OUT_DIR, PANEL_HITS


def load(out_dir: Path) -> list[dict]:
    rows = [json.loads(p.read_text()) for p in sorted(out_dir.glob("len_*.json"))]
    if not rows:
        raise FileNotFoundError(f"no G1 outputs in {out_dir}")
    return sorted(rows, key=lambda r: r["length"])


def summarize(r: dict) -> dict:
    panel = r.get("panel", {})
    hits = {k: v for k, v in panel.items() if v["group"] == "hit"}
    ctrl = {k: v for k, v in panel.items() if v["group"] == "control"}
    best_all = min(panel.values(), key=lambda v: v["ddG_net"]) if panel else None
    best_hit = min(hits.values(), key=lambda v: v["ddG_net"]) if hits else None
    best_ctl = min(ctrl.values(), key=lambda v: v["ddG_net"]) if ctrl else None
    ont = r["on_target_ddG_net"]
    return {
        "length": r["length"],
        "start": r["start"],
        "aso": r["aso"],
        "gc": r["gc_fraction"],
        "on_target": ont,
        "best_off": best_all["ddG_net"] if best_all else None,
        "best_off_gene": best_all["gene_name"] if best_all else None,
        "margin": (best_all["ddG_net"] - ont) if best_all else None,
        "best_hit": best_hit["ddG_net"] if best_hit else None,
        "best_hit_gene": best_hit["gene_name"] if best_hit else None,
        "best_ctl": best_ctl["ddG_net"] if best_ctl else None,
        "polr2h": panel.get("POLR2H", {}).get("ddG_net"),
        "polr2h_mm": panel.get("POLR2H", {}).get("mismatches"),
        "exact": r["mismatch_histogram"][0],
        "mm1": r["mismatch_histogram"][1],
        "worst_fom": r["worst_fom"],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=OUT_DIR)
    a = ap.parse_args()
    rows = load(a.out)

    print("## G -- the Goldilocks trade-off, both axes\n")

    infeas = [r for r in rows if not r.get("feasible")]
    skipped = [(r["length"], r["infeasible_skipped"]) for r in rows
               if r.get("infeasible_skipped")]
    print("### Feasibility pre-gate (applied before any off-target number was seen)")
    print("  A_5ss scores the CONSEQUENCE of occluding a window, never whether an")
    print("  oligo can occupy it. These top-ranked A3 footprints form NO duplex at")
    print("  all against their own perfect complement, so their predicted")
    print("  accessibility gain is unreachable by any real oligo:")
    if skipped:
        for L, xs in skipped:
            for x in xs:
                gc = sum(c in "GC" for c in x["aso"]) / len(x["aso"])
                print(f"    len {L:>2} @ {x['start']}  aso {x['aso']:<12} GC {gc:>4.0%}  REJECTED")
    else:
        print("    (none)")
    if infeas:
        print(f"  Lengths with NO feasible span-robust footprint: "
              f"{[r['length'] for r in infeas]}")
    print()

    s = [summarize(r) for r in rows if r.get("feasible")]

    print("### Axis 1 -- transcriptome-wide uniqueness (exhaustive, no cutoff)")
    print(f"  {'len':>4}{'ASO':>18}{'GC':>6}{'exact sites':>14}{'1 mismatch':>13}")
    for r in s:
        print(f"  {r['length']:>4}{r['aso']:>18}{r['gc']:>6.0%}"
              f"{r['exact']:>14,}{r['mm1']:>13,}")

    print("\n### Axis 2 -- thermodynamic affinity on the validated panel (kcal/mol)")
    print(f"  {'len':>4}{'on-target':>11}{'best off':>10}{'gene':>10}"
          f"{'MARGIN':>9}{'POLR2H':>9}{'mm':>4}")
    for r in s:
        p = f"{r['polr2h']:>9.2f}" if r["polr2h"] is not None else f"{'n/a':>9}"
        pm = f"{r['polr2h_mm']:>4}" if r["polr2h_mm"] is not None else f"{'-':>4}"
        print(f"  {r['length']:>4}{r['on_target']:>11.2f}{r['best_off']:>10.2f}"
              f"{r['best_off_gene']:>10}{r['margin']:>9.2f}{p}{pm}")
    print("  MARGIN = best off-target minus on-target; positive = intended duplex")
    print("  is the stronger one. This is the specificity measure.")

    print("\n### Is there a Goldilocks point?")
    m = [(r["length"], r["margin"]) for r in s]
    best = max(m, key=lambda t: t[1])
    interior = [t for t in m if t[0] != m[0][0] and t[0] != m[-1][0]]
    is_interior = best[0] not in (m[0][0], m[-1][0])
    diffs = [m[i + 1][1] - m[i][1] for i in range(len(m) - 1)]
    monotonic = all(d >= 0 for d in diffs) or all(d <= 0 for d in diffs)
    print(f"  best margin: {best[1]:+.2f} kcal/mol at length {best[0]}")
    print(f"  margin monotonic in length: {monotonic}")
    if is_interior and not monotonic:
        print(f"  INTERIOR OPTIMUM at length {best[0]} -- a Goldilocks point exists.")
    elif monotonic:
        print("  NO GOLDILOCKS POINT: the margin curve is monotonic in length.")
        print("  Reported as such. Bins were not tuned until an optimum appeared.")
    else:
        print(f"  Best margin is at an ENDPOINT (length {best[0]}), so the scanned")
        print("  range does not contain an interior optimum. The true optimum, if")
        print("  any, lies outside 6-18 and this scan cannot locate it.")

    print("\n### Friendly fire on POLR2H")
    print("  POLR2H is the field's canonical nusinersen off-target and the A5")
    print("  blind prediction (rMATS dPSI -0.920, FDR 0).")
    for r in s:
        if r["polr2h"] is None:
            continue
        verdict = "SAFE" if r["polr2h"] - r["on_target"] > 0 else "AT RISK"
        print(f"    len {r['length']:>2}: POLR2H {r['polr2h']:>7.2f} vs on-target "
              f"{r['on_target']:>7.2f}  ->  {verdict}")

    print("\n### Caveats that travel with these numbers")
    print("  - The off-target panel is 4 validated hits + 20 length-matched")
    print("    controls, not the transcriptome. Axis 1 is transcriptome-wide;")
    print("    axis 2 is not. A panel margin is not a transcriptome guarantee.")
    print("  - Length 18 is nusinersen's published site, NOT a scan winner: no")
    print("    18-mer has a span-robust footprint at all (see goldilocks). It is")
    print("    a reference point, and a length-18 row is not evidence about what")
    print("    an optimised 18-mer could do.")
    print("  - Every feasible ASO from 6 to 14 nt lands on the GC core at")
    print("    32059-32062, which is the site A3 recovered blind. The length axis")
    print("    and the position axis are therefore NOT independent here.")


if __name__ == "__main__":
    main()
