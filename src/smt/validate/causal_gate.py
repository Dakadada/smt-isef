"""C4: is the max_bp_span sign flip caused by NOT-YET-TRANSCRIBED downstream RNA?

THE HYPOTHESIS THIS KILLS. Splice-site selection is co-transcriptional
(Bentley 2014, Nat Rev Genet 15:163). ``max_bp_span`` is an ISOTROPIC band
around the diagonal -- it limits reach 5' and 3' equally -- whereas a real
nascent transcript is CAUSAL and one-sided: when Pol II sits at position L,
sequence 3' of L does not exist and cannot pair with anything.

So the proposal was: the span ladder's 70-87% sign flipping is an artifact of
letting downstream sequence that has not been transcribed yet reach back and
sequester the exon 7 5' splice site. Replace the span cap with a prefix cut at
L, and the observable becomes well-defined and physically parameterised (L in
nucleotides transcribed, convertible to seconds at Pol II's measured rate).

``smt.thermo.constrained``'s own docstring motivates the span ladder exactly
this way: global folding "over-predicts long-range helices relative to what a
co-transcriptionally folding transcript can reach". The ladder was adopted as a
STAND-IN for co-transcriptional folding and the process itself was never
modelled.

>>> REFUTED, on data that was already on disk plus five partition functions.
    Both tests below fail, and neither is close.

TEST 1 -- where in the widening order do the flips sit?
    Reachability is monotone in span, so a partner at distance d first becomes
    available at the step that crosses d. SMN2's ISTL1 long-range interaction
    reaches ref ~32496, +438 nt from the readout, which ONLY the 300 -> inf
    step can reach. That step accounts for 1.3% of flips. 62.6% happen at
    50 -> 100. The proposed mechanism predicts the flips in the wrong place.

TEST 2 -- where are the readout's pairing partners, actually?
    98.0-99.7% of the exon 7 5'ss's pairing probability is with UPSTREAM
    partners at every span. The dominant partners sit 10-20 nt upstream, inside
    exon 7 (ref 32030-32042, p ~ 0.99 at -10..-12) -- a local intra-exonic
    hairpin, fully formed as soon as exon 7 is transcribed. The largest single
    downstream partner contributes 0.2%.

WHY THAT IS FATAL, not a tuning problem. A prefix cut at any L past the readout
removes at most ~2% of the readout's pairing. There is nothing downstream to
remove. The structure occluding the exon 7 5'ss is upstream and already exists
when the site appears, so co-transcriptional folding cannot change it.

NOT A POSITIONAL ARTIFACT: the readout sits at construct position 252 of 798,
with MORE sequence downstream (538 nt) than upstream (251 nt).

RECORDED BUT NOT PURSUED HERE. At span=inf, ISS_N1 is 77.0% downstream-paired
while GC_CORE is 76.3% upstream-paired. That asymmetry is a plausible lead on
why full-length ISS-N1 occlusion inverts while the GC-core 8-mer does not
(docs/HANDOFF.md, "OPEN: span-dependence of the ISS-N1 result is not
understood"). It is a NEW hypothesis, discovered after this gate failed, and it
is deliberately NOT tested here -- swapping a fresh hypothesis in to rescue a
refuted arm is how a gate stops being a gate. It gets its own pre-registration
or it does not get run.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

from smt.seq.reference import Interval

DEFAULT_SWEEP_DIR = Path("/home/dashii/orcd/pool/smt-isef-data/sweep")

#: The span ladder in order of INCREASING reach. ``max_bp_span = 0`` disables
#: the cap in ViennaRNA (``smt.thermo.constrained._blocked_fc`` only sets
#: ``md.max_bp_span`` when the value is truthy), so 0 is UNLIMITED and sorts
#: LAST, not first. Reading it as "smallest" inverts the whole analysis.
WIDENING_ORDER = (50, 100, 150, 300, 0)

#: Partner-distance buckets, declared before the measurement.
DISTANCE_BINS = ((0, 50), (50, 100), (100, 150), (150, 300), (300, None))


def load_sweep(sweep_dir: Path, nshards: int = 200) -> dict[tuple[int, int], dict[int, float]]:
    """(start, length) -> {max_bp_span: delta_A} from the A3 sweep CSVs."""
    by: dict[tuple[int, int], dict[int, float]] = defaultdict(dict)
    for shard in range(nshards):
        path = sweep_dir / f"sweep_{shard:04d}_of_{nshards:04d}.csv"
        if not path.exists():
            continue
        with path.open(newline="") as fh:
            for row in csv.DictReader(fh):
                key = (int(row["start"]), int(row["length"]))
                by[key][int(row["max_bp_span"])] = float(row["delta_A"])
    return dict(by)


def flip_transitions(
    sweep: dict[tuple[int, int], dict[int, float]],
    order: tuple[int, ...] = WIDENING_ORDER,
) -> tuple[Counter, int, int]:
    """Count sign changes at each step of the widening order.

    Returns ``(transitions, n_full_ladder, n_flipping)``. A footprint counts as
    flipping if its sign is not constant across the whole ladder.
    """
    full = {k: v for k, v in sweep.items() if all(s in v for s in order)}
    transitions: Counter = Counter()
    flipping = 0
    for values in full.values():
        signs = [values[s] > 0 for s in order]
        if len(set(signs)) > 1:
            flipping += 1
        for i in range(len(order) - 1):
            if signs[i] != signs[i + 1]:
                transitions[(order[i], order[i + 1])] += 1
    return transitions, len(full), flipping


def partner_profile(
    seq: str,
    readout: Interval,
    max_bp_span: int,
    bins: tuple = DISTANCE_BINS,
) -> dict:
    """Where does ``readout``'s pairing probability go, by partner distance?

    ``readout`` is in CONSTRUCT space (1-based inclusive on ``seq``). Returns
    ``{"total", "upstream", "downstream", "bins": {label: mass}}`` where the
    bins cover DOWNSTREAM distances only; upstream mass is reported whole,
    because the hypothesis under test is about downstream reach.
    """
    import RNA

    md = RNA.md()
    if max_bp_span:
        md.max_bp_span = max_bp_span
    fc = RNA.fold_compound(seq, md)
    _, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()
    bpp = fc.bpp()

    out = {"total": 0.0, "upstream": 0.0, "downstream": 0.0,
           "bins": {_bin_label(lo, hi): 0.0 for lo, hi in bins}}
    for i in range(readout.start, readout.end + 1):
        for j in range(1, len(seq) + 1):
            if i == j:
                continue
            a, b = (i, j) if i < j else (j, i)
            p = bpp[a][b]
            if p <= 0:
                continue
            out["total"] += p
            if j < i:
                out["upstream"] += p
                continue
            out["downstream"] += p
            d = j - i
            for lo, hi in bins:
                if d >= lo and (hi is None or d < hi):
                    out["bins"][_bin_label(lo, hi)] += p
                    break
    return out


def _bin_label(lo: int, hi) -> str:
    return f"<{hi}" if lo == 0 else (f">{lo}" if hi is None else f"{lo}-{hi}")


def dominant_partners(seq: str, readout: Interval, max_bp_span: int,
                      top: int = 12, floor: float = 1e-4) -> list[tuple[int, float]]:
    """The ``top`` construct positions carrying most of ``readout``'s pairing."""
    import RNA

    md = RNA.md()
    if max_bp_span:
        md.max_bp_span = max_bp_span
    fc = RNA.fold_compound(seq, md)
    _, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()
    bpp = fc.bpp()
    mass: Counter = Counter()
    for i in range(readout.start, readout.end + 1):
        for j in range(1, len(seq) + 1):
            if i == j:
                continue
            a, b = (i, j) if i < j else (j, i)
            p = bpp[a][b]
            if p > floor:
                mass[j] += p
    return mass.most_common(top)


def main() -> None:
    from smt.seq.construct import GC_CORE, ISS_N1, SS5_EXON7, compact_window

    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep", type=Path, default=DEFAULT_SWEEP_DIR)
    ap.add_argument("--nshards", type=int, default=200)
    args = ap.parse_args()

    print("## C4 -- is span-driven sign flipping caused by downstream reach?")
    print("## PRE-DECLARED GATE. Both tests below had to pass for arms C and K")
    print("## to proceed. Read the verdict at the bottom.\n")

    sweep = load_sweep(args.sweep, args.nshards)
    transitions, n_full, flipping = flip_transitions(sweep)
    print(f"### TEST 1 -- flips by widening step (n = {n_full} full-ladder footprints)")
    print(f"  overall flip rate {100 * flipping / n_full:.1f}%")
    print("  max_bp_span = 0 is UNLIMITED and therefore sorts LAST.\n")
    print(f"  {'step':>14}{'flips':>8}{'share':>9}")
    for i in range(len(WIDENING_ORDER) - 1):
        a, b = WIDENING_ORDER[i], WIDENING_ORDER[i + 1]
        n = transitions[(a, b)]
        label = f"{a} -> {b if b else 'inf'}"
        print(f"  {label:>14}{n:>8}{100 * n / n_full:>8.1f}%")
    print("\n  ISTL1's distal arm is +438 nt from the readout and is reachable")
    print("  ONLY at the 300 -> inf step. Read that row.\n")

    c = compact_window()
    win = c.window
    to_construct = lambda iv: Interval(iv.start - win.start + 1, iv.end - win.start + 1)
    ro = to_construct(SS5_EXON7)
    print(f"### TEST 2 -- readout partner distances (construct {ro.start}-{ro.end}"
          f" of {len(c.sequence)})")
    print(f"  {'span':>6}{'total':>9}", end="")
    for lo, hi in DISTANCE_BINS:
        print(f"{_bin_label(lo, hi):>9}", end="")
    print(f"{'upstream':>10}")
    for span in WIDENING_ORDER:
        prof = partner_profile(c.sequence, ro, span)
        tot = prof["total"] or 1.0
        print(f"  {span if span else 'inf':>6}{prof['total']:>9.3f}", end="")
        for lo, hi in DISTANCE_BINS:
            print(f"{100 * prof['bins'][_bin_label(lo, hi)] / tot:>8.1f}%", end="")
        print(f"{100 * prof['upstream'] / tot:>9.1f}%")

    print("\n  dominant partners at span = inf:")
    for j, p in dominant_partners(c.sequence, ro, 0):
        print(f"    construct {j:>3}  ref {j + win.start - 1}  {j - ro.start:+5d} nt  p={p:.4f}")

    print("\n### Recorded, NOT pursued -- footprint-site pairing polarity")
    for name, iv in (("GC_CORE", GC_CORE), ("ISS_N1", ISS_N1)):
        prof = partner_profile(c.sequence, to_construct(iv), 0)
        tot = prof["total"] or 1.0
        print(f"  {name:<8} ref {iv.start}-{iv.end}  total {prof['total']:>7.3f}"
              f"  upstream {100 * prof['upstream'] / tot:>5.1f}%"
              f"  downstream {100 * prof['downstream'] / tot:>5.1f}%")
    print("  A NEW hypothesis found after this gate failed. Not tested here.")

    print("\n### VERDICT: C4 FAILS. Arms C and K do not proceed.")
    print("  Flips are not where downstream reachability predicts, and the")
    print("  readout has almost no downstream partners to remove. A prefix cut")
    print("  cannot fix an observable held by UPSTREAM structure.")


if __name__ == "__main__":
    main()
