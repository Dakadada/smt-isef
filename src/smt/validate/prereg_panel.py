"""Pre-registered matched-pair ASO panel: the experiment that can falsify A7.

WHY THIS EXISTS. The framework's link function ``A_5ss -> PSI`` is
deliberately unfit, because no real PSI exists yet. That is the correct call
-- fitting a link to synthetic data would be worthless -- but it means the
central claim is currently untested rather than supported, and a reviewer is
right to say so.

A pre-registration answers that without any PSI: the prediction is written
down, committed, and timestamped BEFORE the measurement exists, so it cannot
be adjusted afterwards. What is registered here is deliberately the
project's own weakest point.

THE PREDICTION. A7 found that 87% of footprints flip the sign of their
predicted effect across the ``max_bp_span`` ladder, and that the robust
minority is positionally clustered (``smt.validate.span_robustness``). The
framework therefore claims that SPAN ROBUSTNESS IS A REAL DESIGN CRITERION,
not a bookkeeping convenience. That is falsifiable:

    Take pairs of footprints MATCHED on their best-span |fom| -- so a
    single-span analysis would rank them equally -- where one holds its
    sign across the whole ladder and the other does not. The framework
    predicts the span-robust member of each pair suppresses the STMN2
    cryptic exon MORE than its partner.

    Null: the two members are indistinguishable. If the null holds, span
    robustness carries no information, worst-span ranking was a mistake,
    and A7's candidate list should not be used.

This is a strong test because a single-span analysis -- the naive thing, and
what the field would do -- predicts NO difference within a pair by
construction. Only the ladder-based framework predicts one.

DECISION RULE, fixed here in advance:
  * paired one-sided sign test across the pairs, alpha = 0.05;
  * a pair counts only if both members are measurable above the assay's own
    noise floor, declared by the wet lab BEFORE unblinding;
  * fewer than 6 usable pairs -> report as underpowered, claim nothing.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from smt.validate.cryptic_report import DEFAULT_DIR, SPANS, Footprint, load
from smt.seq.reference import Interval
from smt.validate.cryptic_sweep import TARGETS

#: Matched on best-span |fom| to within this relative tolerance.
MATCH_TOL = 0.15
#: Only pairs of equal length -- an ASO length difference is a confound.
DEFAULT_TARGET = "STMN2_3SS"

#: Minimum separation (nt) between the ROBUST members of different pairs.
#:
#: WITHOUT THIS THE TEST IS PSEUDO-REPLICATED. Ranking robust footprints by
#: worst-span fom puts the same site at 7, 8, 9 and 10 nt in the top four
#: slots, so a "six pair" sign test would really be one site measured six
#: times and the n would be a fiction. Pairs must interrogate DISTINCT
#: positions to count as independent observations.
MIN_PAIR_SEPARATION = 25

#: A ROBUST member qualifies for the PRIMARY test only if its weakest span
#: retains at least this fraction of its own best-span magnitude.
#:
#: Sign-consistency alone is not enough. A footprint can hold its sign across
#: the whole ladder at |fom| ~ 1e-7, which is the framework predicting NO
#: effect rather than a robust one; pairing it against a partner and then
#: expecting it to win would be testing a prediction the framework never
#: made. This is declared here, before the panel is generated, and pairs
#: below it are still printed -- as SECONDARY, never dropped.
PRIMARY_RETENTION = 0.01

#: Minimum separation between the two MEMBERS of a pair. Members 1 nt apart
#: are nearly the same oligo and cannot be distinguished by any assay.
MIN_MEMBER_SEPARATION = 10


def complement(seq: str) -> str:
    return seq.translate(str.maketrans("ACGUTacgut", "UGCAAugcaa"))


def aso_for(construct, fp: Footprint) -> str:
    """The ASO is the reverse complement of the footprint's target RNA.

    ``Footprint.start`` is a REFERENCE coordinate (see
    ``smt.validate.cryptic_sweep.Footprint``), so it must be mapped into
    construct space before slicing. ``map_interval`` raises rather than
    returning fragments if the footprint is not contiguous in the construct,
    which is the behaviour we want: a split footprint is not a footprint.
    """
    rna = construct.rna_sub(construct.map_interval(Interval(fp.start, fp.end)))
    return complement(rna)[::-1]


def best_mag(fp: Footprint) -> float:
    return max(abs(v) for v in fp.fom.values())


def build_pairs(prints: list[Footprint], target: str, n_pairs: int) -> list[tuple]:
    group = [p for p in prints if p.target == target and len(p.fom) == len(SPANS)]
    robust = sorted((p for p in group if p.desired_everywhere),
                    key=lambda p: min(p.fom.values()), reverse=True)
    flipping = [p for p in group if not p.consistent]
    used: set[tuple] = set()
    chosen_starts: list[int] = []
    pairs = []
    for r in robust:
        if len(pairs) >= n_pairs:
            break
        if any(abs(r.start - s) < MIN_PAIR_SEPARATION for s in chosen_starts):
            continue
        target_mag = best_mag(r)
        cands = [f for f in flipping
                 if f.length == r.length
                 and (f.start, f.length) not in used
                 and abs(best_mag(f) - target_mag) <= MATCH_TOL * target_mag]
        if not cands:
            continue
        # closest match on magnitude, so the pair is as confounded-free as possible
        cands = [f for f in cands
                 if abs(f.start - r.start) >= MIN_MEMBER_SEPARATION]
        if not cands:
            continue
        partner = min(cands, key=lambda f: abs(best_mag(f) - target_mag))
        used.add((partner.start, partner.length))
        used.add((r.start, r.length))
        chosen_starts.append(r.start)
        pairs.append((r, partner))
    return pairs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DEFAULT_DIR)
    ap.add_argument("--target", default=DEFAULT_TARGET)
    ap.add_argument("--all-targets", action="store_true")
    ap.add_argument("--pairs", type=int, default=8)
    ap.add_argument("--csv", type=Path, default=None)
    args = ap.parse_args()

    prints, _ = load(args.out)
    targets = list(TARGETS) if args.all_targets else [args.target]
    all_rows, primary, secondary = [], 0, 0
    for tname in targets:
        pairs = build_pairs(prints, tname, args.pairs)
        construct = TARGETS[tname].construct_fn()
        primary, secondary = emit(pairs, construct, tname, all_rows, primary, secondary)

    print(f"\n## POWER, computed after the panel and before any measurement")
    print(f"  PRIMARY pairs (worst span retains >= {100*PRIMARY_RETENTION:.0f}% of best): {primary}")
    print(f"  SECONDARY pairs (sign-consistent but predicted effect ~0):    {secondary}")
    if primary < 6:
        print(f"  {primary} primary pairs is BELOW the 6-pair minimum declared in this")
        print(f"  module's docstring. A one-sided sign test on {primary} pairs cannot")
        print(f"  reach p<0.05 even if EVERY pair goes the predicted way "
              f"(min p = {0.5**primary:.4f}).")
        print("  STATE THAT UP FRONT. The panel needs either more targets that clear")
        print("  the retention bar, or a graded readout (measure the SIZE of the")
        print("  difference within each pair, not just its direction) -- a paired")
        print(f"  Wilcoxon over all {primary + secondary} pairs has real power where a "
              f"sign test on {primary}")
        print("  does not. Recommend the graded readout to the wet lab.")
    if args.csv and all_rows:
        with args.csv.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(all_rows[0]))
            w.writeheader()
            w.writerows(all_rows)
        print(f"\n  written: {args.csv}  ({len(all_rows)} oligos)")


def emit(pairs, construct, tname, rows, primary, secondary):
    print(f"\n## Pre-registered matched-pair panel -- {tname}")
    print(f"construct {construct.name}; {len(pairs)} pairs matched on best-span "
          f"|fom| to within {100*MATCH_TOL:.0f}% and on length exactly.")
    print(f"ROBUST members are >= {MIN_PAIR_SEPARATION} nt apart, so the pairs are "
          f"independent positions and not one site measured repeatedly.")
    print("PREDICTION: within every pair, ROBUST suppresses the cryptic exon more")
    print("than FLIP. A single-span analysis predicts no difference by construction.\n")
    for i, (r, f) in enumerate(pairs, 1):
        retention = min(r.fom.values()) / best_mag(r)
        tier = "PRIMARY" if retention >= PRIMARY_RETENTION else "secondary"
        if tier == "PRIMARY":
            primary += 1
        else:
            secondary += 1
        print(f"  pair {i:>2} [{tier}]  robust worst/best = {100*retention:.2f}%")
        for role, p in (("ROBUST", r), ("FLIP", f)):
            signs = "".join("+" if p.fom[s] > 0 else "-" for s in SPANS)
            aso = aso_for(construct, p)
            print(f"     {role:<7} start {p.start} len {p.length:>2}  "
                  f"best|fom| {best_mag(p):.3e}  worst {min(p.fom.values()):+.3e}  "
                  f"signs {signs}  ASO 5'-{aso}-3'")
            rows.append({"target": tname, "pair": i, "tier": tier, "role": role,
                         "start": p.start, "length": p.length, "end": p.end,
                         "best_abs_fom": f"{best_mag(p):.6e}",
                         "worst_fom": f"{min(p.fom.values()):.6e}",
                         "retention": f"{retention:.6f}",
                         "sign_ladder": signs, "aso_5to3": aso})
        print()
    return primary, secondary


if __name__ == "__main__":
    main()
