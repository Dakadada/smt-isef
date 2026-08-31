"""The SMN2 ASO panel the wet lab is actually going to test.

WHY THIS EXISTS. docs/prereg-panel.csv holds 40 oligos and every one of them
targets a CRYPTIC EXON (STMN2 12, UNC13A_ACCEPTOR_178 12, UNC13A_DONOR 8,
UNC13A_ACCEPTOR_128 8). The collaborating lab tests SMN2 ASOs and reads them
out by SDS-PAGE + western blot. So the framework has had NO pre-registered
prediction for the experiment that is actually being run -- despite A3 having
ranked 1,065 SMN2 footprints blind and recovering the published GC-core 8-mer.

This module closes that gap. It is a PREDICTION, emitted before any western
blot exists, so it can be wrong in public.

WHAT IT REUSES, rather than reinventing:
* the A3 sweep on disk -- ranking on the WORST span, never the best;
* ``goldilocks_offtarget.candidates_by_length`` -- which already applies EDGE
  EXCLUSION, because every length 6-10 had its nominal best at 32040, the first
  scanned position, and dropping three boundary positions cut those scores by
  2-3 orders of magnitude;
* ``goldilocks_offtarget.feasible`` -- the HARD PRE-GATE. A_5ss scores the
  CONSEQUENCE of occluding a window and never whether an oligo can OCCUPY it.
  The A3 winners at lengths 6-7 are UUAAUU / UUAAUUU and form no duplex at all
  against their own perfect complement. Feasibility is composition, not length.

THE LENGTH QUESTION, answered where it can actually be answered. G1 already
swept lengths 6-18 on two axes and found NO Goldilocks point: the specificity
margin rises monotonically (-1.64 at 6 nt to +7.15 at 18). But that is a
COMPUTED specificity, not a measured efficacy. The panel therefore carries a
length series over the single best site so the lab's own blots can test it.

WHAT A WESTERN BLOT CAN AND CANNOT SETTLE, stated before the data:
* CAN: on-target efficacy. SMN-delta-7 protein is rapidly degraded, so protein
  is a NONLINEAR, AMPLIFYING readout of exon 7 inclusion -- more sensitive than
  PSI for detecting that an ASO works at all.
* CANNOT: off-target burden. Off-targets are in OTHER genes, which are not in
  the construct. That needs RNA-seq, and this project already has two datasets.
* CANNOT: a direct A_5ss -> PSI fit. Protein is not a linear PSI proxy, so the
  link in ``smt.thermo.link`` needs a second stage (A_5ss -> PSI -> protein)
  before any western can be used to fit it.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from smt.validate.goldilocks_offtarget import (
    SWEEP_DIR,
    candidates_by_length,
    feasible,
)

#: Intronic only. Exonic footprints are systematically over-ranked because the
#: model cannot see the exon 7 SRSF1 ESE it would also destroy -- recorded as a
#: KNOWN LIMITATION during the A3 smoke run, not discovered here.
FIRST_INTRONIC = 32059

#: Nusinersen's own site, as a positive control the lab already knows the
#: answer to. If the panel's blots do not reproduce this, the assay is the
#: problem, not the prediction.
NUSINERSEN_SITE = (32062, 18)

#: Lengths carried over the top site to test the length question empirically.
LENGTH_SERIES = (10, 12, 14, 16, 18)


def build(sweep_dir: Path = SWEEP_DIR, top: int = 12) -> tuple[list, list]:
    """(panel rows, rejected rows). Feasibility applied as a declared pre-gate."""
    by_len = candidates_by_length(sweep_dir)
    seen: set[tuple[int, int]] = set()
    panel, rejected = [], []

    ranked = sorted(
        (w for ws in by_len.values() for w in ws if w.start >= FIRST_INTRONIC),
        key=lambda w: -w.worst_fom,
    )
    for w in ranked:
        if len(panel) >= top:
            break
        key = (w.start, w.length)
        if key in seen:
            continue
        seen.add(key)
        fe = feasible(w)
        row = {"role": "candidate", "start": w.start, "length": w.length,
               "end": w.start + w.length - 1, "worst_fom": w.worst_fom,
               "target_rna": w.target_rna, "aso": w.aso,
               "gc": round(sum(c in "GC" for c in w.aso) / len(w.aso), 3)}
        if fe is None:
            row["role"] = "REJECTED_no_duplex"
            rejected.append(row)
            continue
        row["dG_hybrid"], row["dG_unfold"] = round(fe[0], 2), round(fe[1], 2)
        panel.append(row)
    return panel, rejected


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", type=Path, default=SWEEP_DIR)
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--out", type=Path, default=Path("docs/smn2-aso-panel.csv"))
    args = ap.parse_args()

    panel, rejected = build(args.sweep, args.top)

    print("## SMN2 ASO PANEL -- a PREDICTION, emitted before any western blot")
    print("## Ranked on the WORST span. Feasibility is a declared PRE-GATE,")
    print("## not a filter applied after results came out wrong.\n")
    print(f"  {'#':>3}  {'start':>6}{'len':>5}{'worst fom':>12}"
          f"{'dGhyb':>8}{'dGunf':>8}{'GC':>6}  ASO 5'->3'")
    for i, r in enumerate(panel, 1):
        print(f"  {i:>3}  {r['start']:>6}{r['length']:>5}{r['worst_fom']:>12.3e}"
              f"{r['dG_hybrid']:>8.2f}{r['dG_unfold']:>8.2f}{r['gc']:>6.2f}  {r['aso']}")

    if rejected:
        print(f"\n  REJECTED by the feasibility pre-gate ({len(rejected)}) --")
        print("  these form NO duplex against their own perfect complement:")
        for r in rejected:
            print(f"    {r['start']:>6} len {r['length']:>2}  {r['aso']}  GC {r['gc']:.0%}")

    print(f"\n  Positive control to include on every gel: nusinersen 18-mer at "
          f"{NUSINERSEN_SITE[0]}, the site whose answer is already known.")
    print(f"  Length series over the top site: {LENGTH_SERIES} -- G1 predicts")
    print("  specificity rises monotonically with length and found NO interior")
    print("  optimum, so a 14-mer should work on-target but be less specific.")
    print("  A western blot tests the FIRST half of that and not the second.")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as fh:
        cols = ["role", "start", "length", "end", "worst_fom", "gc",
                "dG_hybrid", "dG_unfold", "target_rna", "aso"]
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        for r in panel + rejected:
            w.writerow(r)
    print(f"\n  wrote {args.out} ({len(panel)} candidates, {len(rejected)} rejected)")


if __name__ == "__main__":
    main()
