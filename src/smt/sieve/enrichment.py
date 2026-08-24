"""Is the sieve's hit rate ENRICHMENT, or is it the base rate?

THE ATTACK THIS ANSWERS. A5 produced four striking validated predictions
(POLR2H -0.920, MAP1S -0.879, ANKRD28 -0.619, CCSER2 -0.236, all FDR<1e-6)
against a transcriptome-wide correlation of rho = 0.0025. Those two facts
together are exactly what selection from a large candidate pool looks like:
36,613 SE events were measured, so a handful of large-effect "hits" is
expected by chance and proves nothing on its own.

``docs/HANDOFF.md`` already names the confound: nearly every significant
event in this dataset is NEGATIVE, so "predict skipping" is right by base
rate. The only way to tell a real candidate-generator from that base rate is
to compare the predicted set against a BACKGROUND MATCHED ON THE SAME
MEASUREMENT, and to state the null explicitly.

WHAT IS BEING TESTED, declared before the numbers are computed so it cannot
be adjusted afterwards:

  H0  Genes carrying a sieve prediction of a given category are no more
      likely to show a large-magnitude measured dPSI than genes measured by
      the same rMATS run without such a prediction.
  H1  ``exonic_proximal`` predictions specifically ARE enriched. This is not
      a post-hoc slice: all four validated hits are exonic_proximal and all
      four failed inclusion-direction predictions are intronic, which is a
      pattern visible in the handoff before this module existed.

WHAT WAS CHECKED AND REJECTED, recorded so "did you try X?" has an answer:
rMATS also emits JCEC counts (junction + exon-body reads, more power per
event). Switching to them raises the significant-gene count from 602 to 616,
about 2%, which does not move the binding constraint (the overlap between
predicted genes and measurable events) and would inflate the Bonferroni
denominator for nothing. It was measured, not assumed, and not used.

Every category is reported whether or not it supports H1, and the primary
statistic is a two-sided Fisher exact test against the measured background,
not a p-value hunted across thresholds. Thresholds are declared in
``MAGNITUDES`` up front and ALL of them are reported.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

SIEVE = Path("/home/dashii/orcd/pool/smt-isef-data/sieve/annotated.csv")
RMATS_DIR = Path("/home/dashii/orcd/scratch/smt-rna/rmats/out")

#: ALL five rMATS event classes, not just SE.
#:
#: This is not a second bite at the apple after SE alone came out
#: underpowered -- it is the analysis that matches what the sieve actually
#: predicts. The sieve scores SPLICE-SITE ACCESSIBILITY. A change in
#: accessibility shows up as a retained intron or a shifted 3'/5' site just
#: as readily as a skipped cassette exon, so restricting to SE discards
#: 10,306 measured events for no stated reason. Both analyses are reported
#: side by side and every p-value is Bonferroni-corrected over every test
#: printed, so the expansion cannot buy significance by itself.
EVENT_TYPES = ("SE", "A3SS", "A5SS", "MXE", "RI")

#: Declared up front, all reported. No threshold is chosen after seeing p.
MAGNITUDES = (0.2, 0.5, 0.8)
FDR_MAX = 0.05


def load_measured(rmats_dir: Path, events: tuple[str, ...]) -> dict[str, float]:
    """Strongest significant |dPSI| per gene, over the given event classes."""
    best: dict[str, float] = {}
    for event in events:
        path = rmats_dir / f"{event}.MATS.JC.txt"
        if not path.exists():
            raise FileNotFoundError(path)
        _accumulate(path, best)
    return best


def _accumulate(path: Path, best: dict[str, float]) -> None:
    with path.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            gene = row["geneSymbol"].strip('"')
            try:
                fdr = float(row["FDR"])
                dpsi = float(row["IncLevelDifference"])
            except (ValueError, KeyError):
                continue
            if fdr > FDR_MAX:
                continue
            if gene not in best or abs(dpsi) > abs(best[gene]):
                best[gene] = dpsi


def load_predictions(path: Path) -> dict[str, dict[str, float]]:
    """Strongest |signed_effect| per (gene, category)."""
    out: dict[str, dict[str, float]] = defaultdict(dict)
    with path.open() as fh:
        for row in csv.DictReader(fh):
            gene, cat = row["gene_name"], row["category"]
            try:
                eff = float(row["signed_effect"])
            except ValueError:
                continue
            if cat not in out[gene] or abs(eff) > abs(out[gene][cat]):
                out[gene][cat] = eff
    return out


def fisher_exact(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher exact p for [[a,b],[c,d]], by exact enumeration."""
    from math import comb
    n = a + b + c + d
    row1, col1 = a + b, a + c
    def prob(x: int) -> float:
        return comb(row1, x) * comb(n - row1, col1 - x) / comb(n, col1)
    observed = prob(a)
    lo = max(0, col1 - (n - row1))
    hi = min(row1, col1)
    # 1e-9 relative slack: exact-arithmetic ties must count as "as extreme".
    return min(1.0, sum(p for x in range(lo, hi + 1)
                        if (p := prob(x)) <= observed * (1 + 1e-9)))


def analyse(measured: dict[str, float], predicted: dict[str, dict[str, float]],
            label: str, collect: list) -> None:
    print(f"\n## {label}")
    print(f"genes with a significant event (FDR<{FDR_MAX}): {len(measured)}")
    neg = sum(1 for v in measured.values() if v < 0)
    print(f"BASE RATE CONFOUND, stated first: {neg}/{len(measured)} "
          f"({100*neg/len(measured):.1f}%) of significant events are NEGATIVE. "
          f"'Predict skipping' is right this often by chance alone.")

    categories = sorted({c for cats in predicted.values() for c in cats})
    for mag in MAGNITUDES:
        big = {g for g, v in measured.items() if abs(v) >= mag}
        print(f"\n### |dPSI| >= {mag}   ({len(big)}/{len(measured)} genes, "
              f"base rate {100*len(big)/len(measured):.1f}%)")
        print(f"  {'category':<20}{'n meas':>8}{'n hit':>7}{'rate':>9}"
              f"{'base':>9}{'fold':>8}{'raw p':>11}")
        for cat in categories:
            genes = {g for g, cats in predicted.items() if cat in cats} & set(measured)
            if len(genes) < 5:
                continue
            hit = len(genes & big)
            bg = set(measured) - genes
            bg_hit = len(bg & big)
            rate, base = hit / len(genes), bg_hit / len(bg)
            pval = fisher_exact(hit, len(genes) - hit, bg_hit, len(bg) - bg_hit)
            collect.append((label, mag, cat, len(genes), hit, rate, base, pval))
            fold = rate / base if base else float("inf")
            print(f"  {cat:<20}{len(genes):>8}{hit:>7}{100*rate:>8.1f}%"
                  f"{100*base:>8.1f}%{fold:>8.2f}{pval:>11.4g}")

    genes = [(g, cats["exonic_proximal"], measured[g])
             for g, cats in predicted.items()
             if "exonic_proximal" in cats and g in measured]
    print(f"\n### Within exonic_proximal: does the RANK carry information? (n={len(genes)})")
    if len(genes) >= 8:
        genes.sort(key=lambda t: t[1])
        half = len(genes) // 2
        strong = sorted(t[2] for t in genes[:half])
        weak = sorted(t[2] for t in genes[-half:])
        print(f"  median measured dPSI, most-negative-predicted half: "
              f"{strong[len(strong)//2]:+.3f}")
        print(f"  median measured dPSI, least-negative-predicted half: "
              f"{weak[len(weak)//2]:+.3f}")
        print("  If those are the same, only the CATEGORY LABEL is doing work and")
        print("  the signed magnitude orders nothing. Say that plainly.")
    else:
        print("  too few genes to rank. Report the n, not a p-value.")


def power_note(collect: list, ntests: int) -> None:
    """Turn the negative into a specification: how many genes would it take?

    A null result is only useful if it says what would settle the question.
    This reports the smallest predicted-gene count n at which the OBSERVED
    effect size would clear the corrected threshold, by exact Fisher on the
    same background rate -- i.e. assuming the observed rate is the true rate.
    That is an optimistic floor, and is labelled as one.
    """
    best = min(collect, key=lambda c: c[7])
    label, mag, cat, ngen, hit, rate, base, pval = best
    alpha = 0.05 / ntests
    print(f"## WHAT WOULD SETTLE IT -- power for the strongest signal seen")
    print(f"  {cat} at |dPSI|>={mag}: {hit}/{ngen} = {100*rate:.1f}% vs "
          f"base {100*base:.1f}%  (fold {rate/base:.1f})")
    bg_n, bg_hit = 10000, round(10000 * base)
    for n in range(ngen, 2001):
        k = round(n * rate)
        if fisher_exact(k, n - k, bg_hit, bg_n - bg_hit) < alpha:
            print(f"  If {100*rate:.1f}% is the TRUE rate, this clears "
                  f"alpha={alpha:.2g} at n = {n} predicted genes with a")
            print(f"  measurable event -- {n/ngen:.0f}x the {ngen} available here. "
                  f"That is the experiment to ask for,")
            print(f"  and it is an OPTIMISTIC floor: it assumes the observed "
                  f"rate is exactly right.")
            return
    print(f"  No n <= 2000 clears alpha={alpha:.2g} at this effect size.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sieve", type=Path, default=SIEVE)
    ap.add_argument("--rmats", type=Path, default=RMATS_DIR)
    ap.add_argument("--prior-tests", type=int, default=0,
                    help="tests already run against the SAME samples in another "
                         "grouping; added to the Bonferroni denominator. Two "
                         "groupings of one dataset are NOT independent "
                         "replication, and pretending otherwise is how a "
                         "re-analysis manufactures significance.")
    args = ap.parse_args()

    predicted = load_predictions(args.sieve)
    print(f"genes carrying any sieve prediction: {len(predicted)}")
    collect: list = []
    analyse(load_measured(args.rmats, ("SE",)), predicted,
            "SE ONLY -- the original analysis", collect)
    analyse(load_measured(args.rmats, EVENT_TYPES), predicted,
            "ALL FIVE EVENT CLASSES -- what the sieve actually predicts", collect)

    n = len(collect) + args.prior_tests
    if args.prior_tests:
        print(f"\n## MULTIPLE TESTING -- {len(collect)} tests printed above, PLUS "
              f"{args.prior_tests} already")
        print("spent on the same samples under a different grouping. Two groupings")
        print("of ONE dataset are not independent replication, so the denominator")
        print(f"is {n}, not {len(collect)}.")
    else:
        print(f"\n## MULTIPLE TESTING -- {n} tests were printed above, so every one "
              f"of them is corrected")
    print("Bonferroni over ALL tests in this report, including the ones that")
    print("did not support the hypothesis. Nothing is excluded from the count.")
    survivors = [c for c in collect if c[7] * n < 0.05]
    print(f"  {'label':<46}{'|dPSI|':>7}{'category':<20}{'raw p':>11}{'corrected':>12}")
    for label, mag, cat, ngen, hit, rate, base, pval in sorted(collect, key=lambda c: c[7])[:6]:
        print(f"  {label[:44]:<46}{mag:>7}{cat:<20}{pval:>11.4g}{min(1.0, pval*n):>12.4g}")
    print()
    power_note(collect, n)
    if survivors:
        print(f"\n  {len(survivors)} test(s) SURVIVE correction.")
    else:
        print(f"\n  NOTHING SURVIVES Bonferroni correction over {n} tests.")
        print("  That is the result. The four validated hits remain individually")
        print("  real and individually striking, but they are NOT demonstrated")
        print("  enrichment, and the sieve may not be described as validated.")


if __name__ == "__main__":
    main()
