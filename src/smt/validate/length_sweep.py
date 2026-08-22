"""A3 falsification test: does the model recover the known ISS-N1 answer blind?

Tiles hypothetical effector footprints of length 6-20 nt across and around
ISS-N1, scores each by its effect on exon 7 5'ss accessibility, and writes the
ranked result. The published ground truth -- that an 8-mer over the GC-rich
core (intron 7 +7..+14) potently restores exon 7 inclusion -- is NOT used
anywhere in the scoring. It is only consulted afterwards, when grading.

Reporting rule: A_5ss is sensitive to ``max_bp_span``, and for some footprints
the sign of the effect flips with it. Every footprint is therefore scored over
a span ladder and the ladder is reported. A single-span result is not a result.

Sharding is idempotent: each shard writes its own file, skips it if already
present, and writes atomically via a .tmp rename. A preempted shard costs at
most the work of that one shard.
"""

from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass
from pathlib import Path

from smt.seq.construct import compact_window, GC_CORE, ISS_N1
from smt.seq.reference import Interval
from smt.thermo.constrained import a_5ss

#: Footprint 5' start positions to scan, reference coords. Spans the exon 7
#: 3' end through well past ISS-N1 (intron 7 ends at 32496).
SCAN_START = 32040
SCAN_END = 32110

#: Effector lengths in nt. 6 is below any published active ASO; 20 is above
#: nusinersen's 18. The known answer (8) sits in the interior, not at an edge.
LENGTHS = tuple(range(6, 21))

#: max_bp_span ladder. 0 == unrestricted global folding.
SPANS = (0, 50, 100, 150, 300)

#: Window size. intron 6 tail length; exon 8 head uses the module default.
INTRON6_TAIL = 200


@dataclass(frozen=True)
class Footprint:
    start: int  # reference coord, 1-based inclusive
    length: int

    @property
    def interval(self) -> Interval:
        return Interval(self.start, self.start + self.length - 1)

    @property
    def key(self) -> str:
        return f"{self.start}_{self.length}"


def enumerate_footprints() -> list[Footprint]:
    """Full sweep, in a fixed deterministic order. Sharding partitions this."""
    return [
        Footprint(s, L)
        for s in range(SCAN_START, SCAN_END + 1)
        for L in LENGTHS
    ]


def shard_of(items: list[Footprint], shard: int, nshards: int) -> list[Footprint]:
    """Stride partition. Union over all shards == items, with no overlap."""
    if not 0 <= shard < nshards:
        raise ValueError(f"shard {shard} out of range for nshards {nshards}")
    return items[shard::nshards]


def score(fp: Footprint, construct, spans=SPANS) -> list[dict]:
    """Score one footprint across the span ladder.

    delta_A > 0 means the footprint frees the 5'ss, which is the direction a
    therapeutic steric blocker must produce.
    """
    rows = []
    for span in spans:
        base = a_5ss(construct, None, max_bp_span=span).A_5ss
        blk = a_5ss(construct, fp.interval, max_bp_span=span).A_5ss
        rows.append(
            {
                "start": fp.start,
                "length": fp.length,
                "end": fp.interval.end,
                "max_bp_span": span,
                "A_baseline": f"{base:.6e}",
                "A_blocked": f"{blk:.6e}",
                "delta_A": f"{blk - base:.6e}",
                "fold_change": f"{blk / base:.4f}" if base > 0 else "inf",
                # Annotations for grading only -- never used in scoring.
                "overlaps_gc_core": int(
                    fp.start <= GC_CORE.end and fp.interval.end >= GC_CORE.start
                ),
                "overlaps_iss_n1": int(
                    fp.start <= ISS_N1.end and fp.interval.end >= ISS_N1.start
                ),
            }
        )
    return rows


FIELDS = [
    "start", "length", "end", "max_bp_span", "A_baseline", "A_blocked",
    "delta_A", "fold_change", "overlaps_gc_core", "overlaps_iss_n1",
]


def run_shard(shard: int, nshards: int, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    final = out_dir / f"sweep_{shard:04d}_of_{nshards:04d}.csv"
    if final.exists():
        print(f"shard {shard}: output exists, skipping -> {final}")
        return final

    construct = compact_window(intron6_tail=INTRON6_TAIL)
    items = shard_of(enumerate_footprints(), shard, nshards)
    print(f"shard {shard}/{nshards}: {len(items)} footprints x {len(SPANS)} spans")

    tmp = final.with_suffix(".csv.tmp")
    with tmp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for n, fp in enumerate(items, 1):
            for row in score(fp, construct):
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
