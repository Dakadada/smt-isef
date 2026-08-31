"""Draw what A_5ss actually measures: a splice site opening or closing
because an ASO bound somewhere ELSE on the molecule.

This is the mechanism figure. The framework's whole claim is that an
oligo can bind 300 nt away from a splice site and change whether the
spliceosome can reach it, through nothing but a shift in how the RNA
folds. That is a secondary-structure claim, so this draws secondary
structure -- not a docked pose. THE PROJECT NEVER DOCKED ANYTHING and a
3D figure implying otherwise would be a fabrication.

Both panels are folded by the same code path the sweep uses
(``smt.thermo``), at a max_bp_span drawn from the same ladder, so the
picture cannot disagree with the table.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import RNA

from smt.seq.reference import Interval
from smt.validate.cryptic_sweep import TARGETS

#: Drawn from the sweep's own ladder, never a value picked for looks.
DEFAULT_SPAN = 150


def fold(seq: str, span: int, blocked: Interval | None):
    """MFE structure + per-base unpaired probability, optionally with a
    footprint forced unpaired -- the same hard constraint the sweep applies."""
    md = RNA.md()
    md.max_bp_span = span
    fc = RNA.fold_compound(seq, md)
    if blocked is not None:
        for i in range(blocked.start, blocked.end + 1):
            fc.hc_add_up(i, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
    structure, _ = fc.mfe()
    fc.exp_params_rescale(_)
    fc.pf()
    bpp = np.array(fc.bpp())
    n = len(seq)
    paired = np.zeros(n + 1)
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            p = bpp[i][j]
            paired[i] += p
            paired[j] += p
    return structure, np.clip(1.0 - paired[1:], 0.0, 1.0)


def draw_structure(ax, seq, structure, site, footprint, title):
    # ViennaRNA returns a tuple of COORDINATE proxies, not an array, and it
    # is one longer than the sequence (a trailing sentinel).
    raw = RNA.simple_xy_coordinates(structure)
    x = np.array([p.X for p in raw[: len(seq)]])
    y = np.array([p.Y for p in raw[: len(seq)]])
    pt = RNA.ptable(structure)
    for i in range(1, len(seq) + 1):
        j = pt[i]
        if j > i:
            ax.plot([x[i - 1], x[j - 1]], [y[i - 1], y[j - 1]],
                    color="#C9CBD1", lw=0.5, zorder=1)
    ax.plot(x, y, color="#6E7480", lw=0.8, zorder=2)
    def mask(iv):
        m = np.zeros(len(seq), dtype=bool)
        m[iv.start - 1: iv.end] = True
        return m
    ax.scatter(x, y, s=3, color="#B8BCC4", zorder=3, linewidths=0)
    m = mask(site)
    ax.scatter(x[m], y[m], s=34, color="#B24726", zorder=5, linewidths=0,
               label="cryptic splice site (the readout)")
    if footprint is not None:
        m = mask(footprint)
        ax.scatter(x[m], y[m], s=26, color="#0F8A72", zorder=4, linewidths=0,
                   label="ASO footprint (forced unpaired)")
    ax.set_title(title, fontsize=10, pad=8)
    ax.set_aspect("equal")
    ax.axis("off")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="STMN2_3SS")
    ap.add_argument("--start", type=int, default=5716, help="footprint start, reference coords")
    ap.add_argument("--length", type=int, default=10)
    ap.add_argument("--span", type=int, default=DEFAULT_SPAN)
    ap.add_argument("--out", type=Path, default=Path("docs/figures"))
    args = ap.parse_args()

    t = TARGETS[args.target]
    c = t.construct_fn()
    seq = c.sequence.replace("T", "U")
    site = c.map_interval(t.site)
    fp_ref = Interval(args.start, args.start + args.length - 1)
    fp = c.map_interval(fp_ref)

    s_base, u_base = fold(seq, args.span, None)
    s_block, u_block = fold(seq, args.span, fp)

    a_base = float(np.prod(u_base[site.start - 1: site.end]))
    a_block = float(np.prod(u_block[site.start - 1: site.end]))

    fig = plt.figure(figsize=(13, 8.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.1, 1.0], hspace=0.28, wspace=0.06)
    ax1, ax2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    draw_structure(ax1, seq, s_base, site, None, "unoccluded")
    draw_structure(ax2, seq, s_block, site, fp,
                   f"ASO bound at {fp_ref.start}, {args.length} nt")
    ax2.legend(loc="lower right", fontsize=8, frameon=False)

    ax3 = fig.add_subplot(gs[1, :])
    idx = np.arange(1, len(seq) + 1)
    lo, hi = max(1, site.start - 220), min(len(seq), site.end + 220)
    sl = slice(lo - 1, hi)
    ax3.plot(idx[sl], u_base[sl], color="#6E7480", lw=1.4, label="unoccluded")
    ax3.plot(idx[sl], u_block[sl], color="#0F8A72", lw=1.4, label="ASO bound")
    ax3.axvspan(site.start, site.end, color="#B24726", alpha=0.16, lw=0)
    ax3.axvspan(fp.start, fp.end, color="#0F8A72", alpha=0.13, lw=0)
    # Stagger the two annotations: at close separations they collide.
    ax3.annotate("cryptic splice site", xy=((site.start + site.end) / 2, 1.06),
                 fontsize=8, color="#B24726", ha="center", annotation_clip=False)
    ax3.annotate("ASO", xy=((fp.start + fp.end) / 2, 0.995 if abs(fp.start - site.end) < 60 else 1.06),
                 fontsize=8, color="#0F8A72", ha="center", annotation_clip=False)
    ax3.set_ylabel("P(base unpaired)", fontsize=9)
    ax3.set_xlabel("position in construct (nt)", fontsize=9)
    ax3.set_ylim(0, 1.12)
    ax3.legend(fontsize=8, frameon=False, loc="lower left")
    ax3.spines[["top", "right"]].set_visible(False)

    sep = abs(fp.start - site.end) if fp.start > site.end else abs(site.start - fp.end)
    fig.suptitle(
        f"{args.target}   max_bp_span = {args.span}\n"
        f"A_5ss  {a_base:.3e} -> {a_block:.3e}   "
        f"({a_block/a_base:.2f}x)   ASO sits {sep} nt from the splice site "
        f"and never touches it",
        fontsize=11, y=0.98)
    args.out.mkdir(parents=True, exist_ok=True)
    path = args.out / f"occlusion_{args.target}_{args.start}_{args.length}_span{args.span}.png"
    fig.savefig(path, dpi=170, bbox_inches="tight", facecolor="white")
    print(f"A_5ss unoccluded {a_base:.6e}  occluded {a_block:.6e}  "
          f"ratio {a_block/a_base:.4f}")
    print(f"ASO-to-splice-site separation: {sep} nt")
    print(f"written: {path}")


if __name__ == "__main__":
    main()
