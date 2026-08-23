"""Aggregate the A7 cryptic sweep shards into a ranked candidate table.

THE ONE RULE THIS MODULE ENFORCES: a candidate is reported across the WHOLE
``max_bp_span`` ladder, never at one span. ``docs/HANDOFF.md`` records that
ISS-N1's sign flips with span for SMN2, and the same span-dependence
reproduces at these cryptic loci. Reporting the span that agrees with
expectation would be the single easiest way to manufacture a result here, so
the ranking is built on SIGN CONSISTENCY ACROSS ALL FIVE SPANS first and on
effect size second. Footprints whose sign flips are not discarded -- they are
reported in their own section, because "this locus is span-dependent" is a
finding about the physics, not a nuisance to hide.

``fom`` is already signed for the desired direction (OCCLUDE at every target
here, since these are cryptic splice sites whose accessibility we want to
DECREASE), so fom > 0 means the footprint does the therapeutically desired
thing. See ``smt.thermo.occlusion``.

Rows with ``direct_overlap=1`` carry no scores by construction: a splice site
occluded by direct ASO occupancy is competitive hybridization, not structural
propagation, and this framework does not model it. They are counted and
reported, never silently dropped and never imputed.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

DEFAULT_DIR = Path("/home/dashii/orcd/pool/smt-isef-data/cryptic_sweep")

#: Must match smt.validate.cryptic_sweep.SPANS.
SPANS = (0, 50, 100, 150, 300)


@dataclass(frozen=True)
class Footprint:
    """One (target, start, length) scored across the whole span ladder."""

    target: str
    gene: str
    start: int
    length: int
    end: int
    fom: dict[int, float]
    fold_change: dict[int, float]

    @property
    def signs(self) -> set[int]:
        return {1 if v > 0 else -1 if v < 0 else 0 for v in self.fom.values()}

    @property
    def consistent(self) -> bool:
        """True iff every span agrees on the sign AND all five spans are present."""
        return len(self.fom) == len(SPANS) and len(self.signs) == 1

    @property
    def desired_everywhere(self) -> bool:
        return self.consistent and self.signs == {1}

    @property
    def worst_fom(self) -> float:
        """The weakest span's fom -- the honest headline for a candidate.

        Ranking on the best span would report each footprint at whichever
        span flatters it most, which is cherry-picking one level down.
        """
        return min(self.fom.values())


def load(out_dir: Path, nshards: int = 800) -> tuple[list[Footprint], dict]:
    by_key: dict[tuple, dict] = defaultdict(dict)
    meta: dict[tuple, tuple] = {}
    audit = {"shards_found": 0, "missing": [], "rows": 0, "direct_overlap": 0}
    for shard in range(nshards):
        path = out_dir / f"cryptic_sweep_{shard:04d}_of_{nshards:04d}.csv"
        if not path.exists():
            audit["missing"].append(shard)
            continue
        audit["shards_found"] += 1
        with path.open(newline="") as fh:
            for row in csv.DictReader(fh):
                audit["rows"] += 1
                if row["direct_overlap"] == "1":
                    audit["direct_overlap"] += 1
                    continue
                key = (row["target"], int(row["start"]), int(row["length"]))
                meta[key] = (row["gene"], int(row["end"]))
                by_key[key][int(row["max_bp_span"])] = (
                    float(row["fom"]), float(row["fold_change"]) if row["fold_change"] != "inf" else float("inf")
                )
    prints = [
        Footprint(target=k[0], gene=meta[k][0], start=k[1], length=k[2], end=meta[k][1],
                  fom={s: v[0] for s, v in spans.items()},
                  fold_change={s: v[1] for s, v in spans.items()})
        for k, spans in by_key.items()
    ]
    return prints, audit


def report(prints: list[Footprint], audit: dict, top: int = 10) -> None:
    print("## A7 cryptic sweep -- aggregate")
    print(f"shards {audit['shards_found']}  rows {audit['rows']}  "
          f"direct_overlap {audit['direct_overlap']}  scored footprints {len(prints)}")
    if audit["missing"]:
        m = audit["missing"]
        print(f"!! {len(m)} SHARDS MISSING -- this table is INCOMPLETE: "
              f"{m[:10]}{' ...' if len(m) > 10 else ''}")
    else:
        print("all shards present")

    partial = [p for p in prints if len(p.fom) != len(SPANS)]
    if partial:
        print(f"!! {len(partial)} footprints lack a full span ladder -- excluded from ranking")

    by_target: dict[str, list[Footprint]] = defaultdict(list)
    for p in prints:
        by_target[p.target].append(p)

    for target in sorted(by_target):
        group = [p for p in by_target[target] if len(p.fom) == len(SPANS)]
        good = sorted((p for p in group if p.desired_everywhere),
                      key=lambda p: p.worst_fom, reverse=True)
        flips = [p for p in group if not p.consistent]
        wrong = [p for p in group if p.consistent and p.signs == {-1}]
        print(f"\n### {target}  ({group[0].gene if group else '?'})")
        print(f"  {len(group)} footprints with a full ladder: "
              f"{len(good)} desired at every span, {len(wrong)} wrong-signed at every span, "
              f"{len(flips)} SIGN-FLIPPING ({100*len(flips)/len(group):.1f}%)")
        if not good:
            print("  NO footprint achieves the desired sign across the whole ladder.")
            print("  That is the result. Do not report a single-span winner instead.")
            continue
        header = "  " + f"{'start':>7}{'len':>5}{'end':>7}" + "".join(f"{'fom@'+str(s):>13}" for s in SPANS)
        print(header)
        for p in good[:top]:
            row = f"  {p.start:>7}{p.length:>5}{p.end:>7}"
            row += "".join(f"{p.fom[s]:>13.4e}" for s in SPANS)
            print(row)
        print(f"  ranked by WORST span, not best. Best of these at its weakest span: "
              f"{good[0].worst_fom:.4e}")

    print("\n## cross-target comparison, best ladder-robust candidate per target")
    print("  A sign-consistent footprint whose WEAKEST span is ~0 is not a")
    print("  candidate, it is a rounding artifact. These numbers are printed")
    print("  side by side so that cannot be read past.")
    best: list[tuple[str, Footprint]] = []
    for target in sorted(by_target):
        good = [p for p in by_target[target]
                if len(p.fom) == len(SPANS) and p.desired_everywhere]
        if good:
            best.append((target, max(good, key=lambda p: p.worst_fom)))
    for target, p in sorted(best, key=lambda t: t[1].worst_fom, reverse=True):
        print(f"  {target:<22} {p.start}+{p.length:<3} worst-span fom {p.worst_fom:.4e}")
    if len(best) > 1:
        ratio = best[0][1].worst_fom / max(best[-1][1].worst_fom, 1e-300)
        strongest = max(best, key=lambda t: t[1].worst_fom)
        weakest = min(best, key=lambda t: t[1].worst_fom)
        print(f"  {strongest[0]} is {strongest[1].worst_fom / max(weakest[1].worst_fom, 1e-300):.3g}x "
              f"{weakest[0]} at the weakest span.")

    total = [p for p in prints if len(p.fom) == len(SPANS)]
    if total:
        flip_frac = sum(1 for p in total if not p.consistent) / len(total)
        print(f"\n## span-dependence, all targets: {100*flip_frac:.1f}% of footprints "
              f"flip sign across the ladder.")
        print("Read that as the headline caveat, not a footnote: at these loci the "
              "sign of the predicted effect is span-dependent for a large fraction "
              "of footprints, exactly as recorded for ISS-N1 in docs/HANDOFF.md.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DEFAULT_DIR)
    ap.add_argument("--nshards", type=int, default=800)
    ap.add_argument("--top", type=int, default=10)
    args = ap.parse_args()
    prints, audit = load(args.out, args.nshards)
    report(prints, audit, args.top)


if __name__ == "__main__":
    main()
