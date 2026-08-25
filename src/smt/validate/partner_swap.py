"""Arm N: does occluding ISS-N1 FREE the 5' splice site, or just move its jailer?

POST-HOC, AND SAID SO. This hypothesis was formed after Arm M's predictions
failed, and the project's stopping rule would have blocked it. It runs by
explicit user override, recorded in docs/preregistration.md ("Arm N"). The
post-hoc origin is a permanent caveat on the result.

THE PARADOX. C4 showed the exon 7 5'ss is held almost entirely by one upstream
block (32030-32042): 6.17 of 6.19 total pairing at span 50. Arm M then showed
occluding ISS_N1 nearly abolishes that pairing (6.12 -> 0.29), while GC_CORE --
the footprint that actually works -- barely touches it. If ISS_N1 destroyed the
jailer, A_5ss should rocket. Instead ISS_N1 is the footprint whose sign INVERTS
across the span ladder, unexplained since the project's first week.

N1: ISS_N1 occlusion REDISTRIBUTES the 5'ss's pairing upstream -> downstream
    rather than removing it. Total stays high; the partner changes.
N2: GC_CORE occlusion LOWERS total pairing -- genuine freeing. This would be
    the mechanism A3 has never had.
N3: the span where the downstream share moves is the span where the sign flips.
"""

from __future__ import annotations

import argparse

SPANS = (50, 100, 150, 300, 0)


def pairing_profile(seq, win, readout, blocked, span):
    """Total / upstream / downstream pairing mass for ``readout``."""
    from smt.thermo.constrained import _blocked_fc
    from smt.seq.reference import Interval

    conv = lambda iv: Interval(iv.start - win.start + 1, iv.end - win.start + 1)
    fc = _blocked_fc(seq, [conv(iv) for iv in blocked], span)
    fc.pf()
    bpp = fc.bpp()
    lo, hi = conv(readout).start, conv(readout).end
    up = dn = 0.0
    for i in range(lo, hi + 1):
        for j in range(1, len(seq) + 1):
            if i == j:
                continue
            x, y = (i, j) if i < j else (j, i)
            p = bpp[x][y]
            if p <= 0:
                continue
            if j < i:
                up += p
            else:
                dn += p
    return {"total": up + dn, "up": up, "dn": dn}


def main() -> None:
    from smt.seq.construct import GC_CORE, ISS_N1, SS5_EXON7, compact_window

    argparse.ArgumentParser().parse_args()
    c = compact_window()
    seq, win = c.sequence, c.window

    print("## ARM N -- does ISS-N1 occlusion free the 5'ss, or swap its jailer?")
    print("## POST-HOC hypothesis, run by explicit user override of the")
    print("## stopping rule. Pre-registered before these numbers.\n")

    cases = (("baseline", ()), ("GC_CORE", (GC_CORE,)), ("ISS_N1", (ISS_N1,)))
    prof = {}
    print(f"  {'span':>6}  " + "".join(
        f"{n + ' total':>16}{'dn%':>7}" for n, _ in cases))
    for span in SPANS:
        row = f"  {span if span else 'inf':>6}  "
        for name, blocked in cases:
            p = pairing_profile(seq, win, SS5_EXON7, blocked, span)
            prof[(span, name)] = p
            share = 100 * p["dn"] / p["total"] if p["total"] else 0.0
            row += f"{p['total']:>16.3f}{share:>6.1f}%"
        print(row)

    print("\n### N1 -- ISS_N1: redistribution or freeing?")
    print(f"  {'span':>6}{'total retained':>16}{'dn share shift':>17}")
    worst_ret, worst_shift = 1e9, 1e9
    for span in SPANS:
        b, i = prof[(span, "baseline")], prof[(span, "ISS_N1")]
        ret = i["total"] / b["total"] if b["total"] else 0.0
        sb = 100 * b["dn"] / b["total"] if b["total"] else 0.0
        si = 100 * i["dn"] / i["total"] if i["total"] else 0.0
        worst_ret, worst_shift = min(worst_ret, ret), min(worst_shift, si - sb)
        print(f"  {span if span else 'inf':>6}{100*ret:>15.1f}%{si-sb:>+16.1f}pp")
    print(f"\n  worst span: total retained {100*worst_ret:.1f}% "
          f"(N1 needs >= 50%), downstream shift {worst_shift:+.1f}pp "
          f"(N1 needs >= +20pp)")
    print(f"  N1: {'SUPPORTED' if worst_ret >= 0.5 and worst_shift >= 20 else 'NULL'}")

    print("\n### N2 -- GC_CORE: does total pairing actually DROP?")
    worst_gc = -1e9
    for span in SPANS:
        b, g = prof[(span, "baseline")], prof[(span, "GC_CORE")]
        d = g["total"] - b["total"]
        worst_gc = max(worst_gc, d)
        print(f"  {span if span else 'inf':>6}{d:>+12.3f}")
    print(f"\n  worst-span change {worst_gc:+.3f} (N2 needs < 0)")
    print(f"  N2: {'SUPPORTED -- this is A3 mechanism' if worst_gc < 0 else 'NULL'}")


if __name__ == "__main__":
    main()
