"""A7 discovery sweep: tile candidate effector footprints around the STMN2
and UNC13A cryptic splice sites and score each with the generalized, signed
OCCLUDE readout (``smt.thermo.occlusion``), over a ``max_bp_span`` LADDER --
never a single span. ``docs/HANDOFF.md`` records that ISS-N1's sign flips
with span for SMN2; the same span-dependence reproduces here (see
``docs/HANDOFF.md``'s A7 entry), so a single-span result would be exactly as
misleading here as it was there.

Modeled directly on ``smt.validate.length_sweep`` (the SMN2 / FREE-regime
sweep): same stride-partition sharding, same idempotent skip-if-exists +
atomic ``.tmp`` rename per shard, same "annotate but never use in scoring"
discipline. The only structural difference is that there are FOUR
(gene, splice-site) targets here instead of one (STMN2's single 3' acceptor;
UNC13A's shared donor and its two alternative 3' acceptors), so a "footprint"
in the sweep is qualified by which target it is being scored against.

Footprints that overlap the target's own readout window are not filtered out
(unlike ``length_sweep``) -- they are INCLUDED with ``direct_overlap=1`` and
no score, rather than silently absent from the table. See
``smt.thermo.occlusion`` module docstring for exactly why this framework
cannot score that case (a splice site occluded by direct ASO occupancy is a
different physical mechanism -- competitive hybridization, not structural
propagation -- than what the hard-constraint partition-function ratio
models), and why faking a number there would be worse than flagging it.
"""

from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from smt.seq import targets as T
from smt.seq.construct import Construct
from smt.seq.reference import Interval
from smt.thermo.constrained import A5ssResult, a_5ss
from smt.thermo.occlusion import Direction, overlaps, score

#: Effector lengths in nt. Same range as length_sweep.py's SMN2 sweep (6 is
#: below any published active ASO, 20 is above nusinersen's 18).
LENGTHS = tuple(range(6, 21))

#: max_bp_span ladder. 0 == unrestricted global folding.
SPANS = (0, 50, 100, 150, 300)

#: How far (nt) each target's scan range extends beyond its own splice-site
#: readout window, on both sides. Matches the scale of length_sweep.py's
#: ISS-N1-centered scan (~70 nt total range there); here applied per target
#: so each of the four splice sites gets its own local flanking scan rather
#: than one scan across a whole shared window.
SCAN_MARGIN = 150


@dataclass(frozen=True)
class Target:
    """One (gene, splice-site) pair to sweep footprints against."""

    gene: str
    site_feature: str
    construct_fn: Callable[[], Construct]
    site: Interval  # reference space


TARGETS: dict[str, Target] = {
    "STMN2_3SS": Target("STMN2", "STMN2_CRYPTIC_3SS", T.stmn2_cryptic_window, T.STMN2_CRYPTIC_3SS),
    "UNC13A_DONOR": Target(
        "UNC13A", "UNC13A_CRYPTIC_DONOR", T.unc13a_cryptic_window, T.UNC13A_CRYPTIC_DONOR
    ),
    "UNC13A_ACCEPTOR_128": Target(
        "UNC13A", "UNC13A_CRYPTIC_ACCEPTOR_128", T.unc13a_cryptic_window, T.UNC13A_CRYPTIC_ACCEPTOR_128
    ),
    "UNC13A_ACCEPTOR_178": Target(
        "UNC13A", "UNC13A_CRYPTIC_ACCEPTOR_178", T.unc13a_cryptic_window, T.UNC13A_CRYPTIC_ACCEPTOR_178
    ),
}


@dataclass(frozen=True)
class Footprint:
    target_key: str
    start: int  # reference coord, 1-based inclusive
    length: int

    @property
    def interval(self) -> Interval:
        return Interval(self.start, self.start + self.length - 1)

    @property
    def key(self) -> str:
        return f"{self.target_key}_{self.start}_{self.length}"


def _scan_range(target: Target) -> tuple[int, int]:
    """(lo, hi) footprint-start bounds for ``target``, clamped to its
    construct's own window so every enumerated footprint is fully containable."""
    c = target.construct_fn()
    lo = max(c.window.start, target.site.start - SCAN_MARGIN)
    hi = min(c.window.end, target.site.end + SCAN_MARGIN)
    return lo, hi


def enumerate_footprints() -> list[Footprint]:
    """Full sweep across every target, in a fixed deterministic order."""
    out: list[Footprint] = []
    for key, target in TARGETS.items():
        lo, hi = _scan_range(target)
        for s in range(lo, hi + 1):
            for length in LENGTHS:
                if s + length - 1 > hi:
                    continue
                out.append(Footprint(key, s, length))
    return out


def shard_of(items: list[Footprint], shard: int, nshards: int) -> list[Footprint]:
    """Stride partition. Union over all shards == items, with no overlap."""
    if not 0 <= shard < nshards:
        raise ValueError(f"shard {shard} out of range for nshards {nshards}")
    return items[shard::nshards]


FIELDS = [
    "target", "gene", "site_feature", "start", "length", "end", "max_bp_span",
    "direct_overlap", "A_baseline", "A_blocked", "delta_A", "fom", "desired",
    "fold_change",
]


def build_baselines(
    constructs: dict[str, Construct], spans=SPANS
) -> dict[tuple[str, int], A5ssResult]:
    """Precompute the unoccluded A_5ss for every (target, span) pair.

    The baseline depends on (construct, estimator, site feature, span) and NOT
    on the footprint, so across a sweep it is the same number every time. It is
    also not cheap: measured on these constructs it is 50.0-50.5% of each
    baseline+blocked pair at every span. Computing it once per (target, span)
    -- 4 x 5 = 20 partition functions instead of 17820 x 5 = 89100 -- halves
    the whole sweep. ``occlusion_delta`` validates what it is handed, so a
    mismatched cache entry raises rather than silently scoring against the
    wrong reference.
    """
    return {
        (key, span): a_5ss(
            constructs[key], None, ss5_feature=TARGETS[key].site_feature, max_bp_span=span
        )
        for key in TARGETS
        for span in spans
    }


def score_footprint(
    fp: Footprint,
    constructs: dict[str, Construct],
    spans=SPANS,
    baselines: dict[tuple[str, int], A5ssResult] | None = None,
) -> list[dict]:
    """Score one footprint across the span ladder against its own target.

    ``fom`` (see ``smt.thermo.occlusion``) is positive iff the footprint does
    the therapeutically desired thing: OCCLUDE direction throughout, since
    every target here is a cryptic splice site whose accessibility we want
    to DECREASE. Rows where the footprint directly overlaps the readout
    window carry ``direct_overlap=1`` and empty score fields -- see module
    docstring.
    """
    target = TARGETS[fp.target_key]
    construct = constructs[fp.target_key]
    direct = overlaps(fp.interval, target.site)
    rows = []
    for span in spans:
        if direct:
            rows.append(
                {
                    "target": fp.target_key, "gene": target.gene, "site_feature": target.site_feature,
                    "start": fp.start, "length": fp.length, "end": fp.interval.end,
                    "max_bp_span": span, "direct_overlap": 1,
                    "A_baseline": "", "A_blocked": "", "delta_A": "", "fom": "", "desired": "",
                    "fold_change": "",
                }
            )
            continue
        r = score(
            fp.interval, Direction.OCCLUDE, target.site_feature, construct,
            max_bp_span=span,
            baseline=None if baselines is None else baselines[(fp.target_key, span)],
        )
        rows.append(
            {
                "target": fp.target_key, "gene": target.gene, "site_feature": target.site_feature,
                "start": fp.start, "length": fp.length, "end": fp.interval.end,
                "max_bp_span": span, "direct_overlap": 0,
                "A_baseline": f"{r.baseline.A_5ss:.6e}", "A_blocked": f"{r.blocked.A_5ss:.6e}",
                "delta_A": f"{r.delta_A:.6e}", "fom": f"{r.fom:.6e}", "desired": int(r.desired),
                "fold_change": f"{r.fold_change:.4f}" if r.baseline.A_5ss > 0 else "inf",
            }
        )
    return rows


def run_shard(shard: int, nshards: int, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    final = out_dir / f"cryptic_sweep_{shard:04d}_of_{nshards:04d}.csv"
    if final.exists():
        print(f"shard {shard}: output exists, skipping -> {final}")
        return final

    constructs = {key: target.construct_fn() for key, target in TARGETS.items()}
    items = shard_of(enumerate_footprints(), shard, nshards)
    print(f"shard {shard}/{nshards}: {len(items)} footprints x {len(SPANS)} spans")
    # 20 partition functions up front, reused by every footprint in the shard.
    baselines = build_baselines(constructs)
    print(f"  {len(baselines)} baselines precomputed", flush=True)

    # Unique per process, NOT a fixed "<shard>.csv.tmp". The plan allows
    # overflowing this array onto a second partition, and skip-if-done is
    # checked at start, so the same shard genuinely can be in flight twice.
    # With a shared tmp path those two writers would interleave into one file
    # and os.replace would publish the mixture -- a corrupt shard that looks
    # complete. With a unique tmp path each writes its own and the last
    # os.replace wins atomically, which is harmless: the shard is
    # deterministic, so both would have written identical bytes.
    tmp = final.with_suffix(f".csv.tmp.{os.getpid()}")
    with tmp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for n, fp in enumerate(items, 1):
            for row in score_footprint(fp, constructs, baselines=baselines):
                w.writerow(row)
            if n % 10 == 0:
                fh.flush()
                print(f"  {n}/{len(items)}", flush=True)
    # Atomic: a kill mid-write never leaves a file the skip-check trusts.
    os.replace(tmp, final)
    print(f"shard {shard}: wrote {final}")
    return final


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--shard", type=int, required=True)
    p.add_argument("--nshards", type=int, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    run_shard(a.shard, a.nshards, a.out)


if __name__ == "__main__":
    main()
