"""N4: does 'the 5'ss's partner moved downstream' predict the sign of delta_A?

Arm N found, on ISS_N1 alone, that occlusion drives the exon 7 5'ss's pairing
from ~2% downstream to 95-97% downstream at spans 100/300/inf, and that the
recorded sign of delta_A matches a 90%-downstream rule 5/5 across the ladder.
That is p ~ 0.1 on one footprint -- suggestive, not a result.

This sweeps the rule over every A3 footprint. If it holds, the 87% span
sign-flip -- the project's largest open wound -- has a structural explanation:
the sign is set by WHERE the 5'ss's partner goes when you occlude, and the span
controls whether the downstream partner is reachable at all.

Pre-registered as N4 in docs/preregistration.md, threshold fixed at 90% from
N3's observed separation and NOT re-fit here.

Sharded, skip-if-done, atomic .tmp rename -- the standing preemption pattern.
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

from smt.seq.reference import Interval

SPANS = (50, 100, 150, 300, 0)
DEFAULT_OUT = Path("/home/dashii/orcd/pool/smt-isef-data/swap_sweep")

#: From N3's observed separation (95-97% vs 15-50%). NOT re-fit on this sweep.
DOWNSTREAM_RULE = 90.0


def footprints() -> list[tuple[int, int]]:
    """The same grid A3 swept: starts 32040..32110, lengths 6..20."""
    return [(s, L) for s in range(32040, 32111) for L in range(6, 21)]


def profile(seq, win, readout, footprint, span):
    from smt.thermo.constrained import _blocked_fc

    conv = lambda iv: Interval(iv.start - win.start + 1, iv.end - win.start + 1)
    fc = _blocked_fc(seq, [conv(footprint)], span)
    fc.pf()
    bpp = fc.bpp()
    r = conv(readout)
    up = dn = 0.0
    for i in range(r.start, r.end + 1):
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
    return up, dn


def main() -> None:
    from smt.seq.construct import SS5_EXON7, compact_window

    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, required=True)
    ap.add_argument("--nshards", type=int, default=200)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    dest = args.out / f"swap_{args.shard:04d}_of_{args.nshards:04d}.csv"
    if dest.exists():
        print(f"{dest} exists, skipping")
        return

    c = compact_window()
    seq, win = c.sequence, c.window
    mine = [f for k, f in enumerate(footprints()) if k % args.nshards == args.shard]

    tmp = dest.with_suffix(".tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["start", "length", "end", "max_bp_span",
                    "up", "dn", "total", "dn_share", "predicted_sign"])
        for start, length in mine:
            fp = Interval(start, start + length - 1)
            # A footprint overlapping the readout is refused upstream by a_5ss;
            # carry it with empty scores rather than dropping it silently.
            overlap = not (fp.end < SS5_EXON7.start or fp.start > SS5_EXON7.end)
            for span in SPANS:
                if overlap:
                    w.writerow([start, length, fp.end, span, "", "", "", "", ""])
                    continue
                up, dn = profile(seq, win, SS5_EXON7, fp, span)
                tot = up + dn
                share = 100 * dn / tot if tot else 0.0
                w.writerow([start, length, fp.end, span, f"{up:.6f}",
                            f"{dn:.6f}", f"{tot:.6f}", f"{share:.3f}",
                            "-" if share >= DOWNSTREAM_RULE else "+"])
    os.replace(tmp, dest)
    print(f"shard {args.shard}: {len(mine)} footprints -> {dest}")


if __name__ == "__main__":
    main()
