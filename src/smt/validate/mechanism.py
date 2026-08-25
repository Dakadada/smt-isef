"""Arm M: WHY does occluding the GC core free the exon 7 5' splice site?

A3 recovered the GC-core 8-mer blind and the mechanism was never asked for.
C4 made the gap sharp: the 5'ss is held by a LOCAL UPSTREAM hairpin (dominant
partners ref 32030-32042, p ~ 0.99 at -10..-12, inside exon 7), while GC_CORE
is 8 nt DOWNSTREAM and does not touch it -- yet occluding GC_CORE raises A_5ss
6-7x. Something indirect connects them and nobody has said what.

Pre-registered in docs/preregistration.md ("Arm M") BEFORE any number here.
M1: GC_CORE occlusion LOWERS P(5'ss : 32030-32042) -- partner competition.
M2: ISS_N1 occlusion does NOT, because its partners are downstream.
M3: across 780 footprints, partner overlap with 32030-32042 predicts rank.

If M1's null holds, A3 is an empirical hit with no mechanism and the standing
stopping rule applies: Arms U and E do not run.
"""

from __future__ import annotations

import argparse

from smt.seq.reference import Interval

SPANS = (50, 100, 150, 300, 0)
PARTNER_BLOCK = Interval(32030, 32042)   # from C4, measured before this arm


def _to_construct(iv: Interval, win: Interval) -> Interval:
    """Reference -> construct space. _blocked_fc indexes the SEQUENCE, so a
    reference interval passed straight through is out of bounds and segfaults
    ViennaRNA rather than raising."""
    return Interval(iv.start - win.start + 1, iv.end - win.start + 1)


def block_pairing(seq, win, a: Interval, b: Interval, blocked, span: int) -> float:
    """Summed pair probability between two reference-space blocks."""
    from smt.thermo.constrained import _blocked_fc

    fc = _blocked_fc(seq, [_to_construct(iv, win) for iv in blocked], span)
    fc.pf()
    bpp = fc.bpp()
    lo1, hi1 = a.start - win.start + 1, a.end - win.start + 1
    lo2, hi2 = b.start - win.start + 1, b.end - win.start + 1
    total = 0.0
    for i in range(lo1, hi1 + 1):
        for j in range(lo2, hi2 + 1):
            x, y = (i, j) if i < j else (j, i)
            if x != y:
                total += bpp[x][y]
    return total


def main() -> None:
    from smt.seq.construct import GC_CORE, ISS_N1, SS5_EXON7, compact_window

    ap = argparse.ArgumentParser()
    ap.parse_args()
    c = compact_window()
    seq, win = c.sequence, c.window

    print("## ARM M -- the mechanism of A3, pre-registered before these numbers")
    print(f"## readout {SS5_EXON7.start}-{SS5_EXON7.end}, partner block "
          f"{PARTNER_BLOCK.start}-{PARTNER_BLOCK.end} (from C4)\n")
    print("### M1/M2 -- P(5'ss : partner block) under each occlusion")
    print(f"  {'span':>6}{'baseline':>11}{'GC_CORE':>11}{'dP':>9}"
          f"{'ISS_N1':>11}{'dP':>9}")
    rows = []
    for span in SPANS:
        base = block_pairing(seq, win, SS5_EXON7, PARTNER_BLOCK, (), span)
        gc = block_pairing(seq, win, SS5_EXON7, PARTNER_BLOCK, (GC_CORE,), span)
        iss = block_pairing(seq, win, SS5_EXON7, PARTNER_BLOCK, (ISS_N1,), span)
        rows.append((span, base, gc, iss))
        print(f"  {span if span else 'inf':>6}{base:>11.4f}{gc:>11.4f}"
              f"{gc-base:>+9.4f}{iss:>11.4f}{iss-base:>+9.4f}")

    worst_gc = max(r[2] - r[1] for r in rows)      # worst = least favourable
    worst_iss = max(r[3] - r[1] for r in rows)
    print(f"\n  worst-span dP, GC_CORE: {worst_gc:+.4f}   ISS_N1: {worst_iss:+.4f}")
    print("  M1 predicted GC_CORE dP < 0 at the worst span.")
    print(f"  M1: {'SUPPORTED' if worst_gc < -0.05 else 'NULL -- see the stopping rule'}")
    print(f"  M2 predicted ISS_N1 does NOT show the same decrease.")
    print(f"  M2: {'SUPPORTED' if worst_iss >= -0.05 else 'NOT SUPPORTED'}")


if __name__ == "__main__":
    main()
