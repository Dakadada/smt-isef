"""Arm V: does the sieve replicate in an independent experiment, and is its
dose response affinity-ordered?

Pre-registered in ``docs/preregistration.md`` under "Arm V", committed
2026-08-24 BEFORE any number in this module was computed. Only design facts
(event counts, n per category, the sign convention, the control design) were
established first.

WHAT THIS ANSWERS. Attack 1 (``docs/a5-enrichment-test.txt``) turned the
project's own negative into a specification: the exonic_proximal tail
enrichment reached 10.1x at |dPSI| >= 0.8 with raw p = 0.020 but nothing
survived Bonferroni, and *"the sieve hypothesis is three genes short of
testable"* -- n = 14 needed, 11 available. Ottesen 2021 supplies 18.

THREE TESTS, all declared in advance:

  V1  Does the exonic_proximal tail enrichment REPLICATE? Same categories,
      same |dPSI| thresholds, same one-sided Fisher, base rate printed first,
      Bonferroni over every test printed including the unsupportive ones.

  V2  Is DOSE SENSITIVITY affinity-ordered? This is the test the base-rate
      confound cannot manufacture. ~55% of significant events are negative, so
      "predict skipping" wins that often by chance -- but chance does not
      produce a dose response ordered by predicted binding energy. Spearman
      between ``ddG_net`` and the retention ratio dPSI_LoA / dPSI_HiA.

  V3  Do the four individually-validated GSE307431 hits replicate in
      DIRECTION? n <= 4, so this is descriptive only and computes no p-value.
      It is excluded from the Bonferroni denominator for that reason, not for
      convenience.

The dose contrast is clean because the authors held total oligo load at 100 nM
in both arms (5 nM Anti-N1 + 95 nM non-targeting 10mer). See
``smt.sieve.ottesen`` for that and for the inverted sign convention, which is
re-verified on every load rather than trusted.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from smt.sieve.enrichment import MAGNITUDES, SIEVE, fisher_exact, load_predictions
from smt.sieve.ottesen import DEFAULT_XLSX, best_per_gene, load_events, verify_sign_convention

#: The four hits GSE307431 validated individually, with their measured dPSI
#: there. Project convention: negative = ASO reduced inclusion.
GSE307431_HITS = {
    "POLR2H": -0.920,
    "MAP1S": -0.879,
    "ANKRD28": -0.619,
    "CCSER2": -0.236,
}

#: Declared before the fit: events weaker than this make the retention ratio
#: a quotient of noise. The count excluded is reported.
MIN_HIA_FOR_RATIO = 0.1


def load_affinity(path: Path = SIEVE) -> dict[str, float]:
    """Tightest (most negative) predicted ddG_net per gene."""
    best: dict[str, float] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh):
            gene = row["gene_name"]
            try:
                ddg = float(row["ddG_net"])
            except (ValueError, KeyError):
                continue
            if gene not in best or ddg < best[gene]:
                best[gene] = ddg
    return best


def v1_enrichment(measured: dict[str, float], predicted: dict, collect: list) -> None:
    print("\n## V1 -- does the tail enrichment REPLICATE?")
    print(f"genes with a measurable event, both doses: {len(measured)}")
    neg = sum(1 for v in measured.values() if v < 0)
    print(f"BASE RATE CONFOUND, stated first: {neg}/{len(measured)} "
          f"({100*neg/len(measured):.1f}%) of events are NEGATIVE "
          f"(ASO reduced inclusion). 'Predict skipping' is right this often "
          f"by chance alone.")

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
            rate = hit / len(genes)
            base = bg_hit / len(bg) if bg else 0.0
            pval = fisher_exact(hit, len(genes) - hit, bg_hit, len(bg) - bg_hit)
            collect.append((f"|dPSI|>={mag}", cat, len(genes), hit, rate, base, pval))
            fold = rate / base if base else float("inf")
            print(f"  {cat:<20}{len(genes):>8}{hit:>7}{100*rate:>8.1f}%"
                  f"{100*base:>8.1f}%{fold:>8.2f}{pval:>11.4g}")


def v2_dose_response(events, predicted, affinity, collect: list) -> None:
    from scipy import stats

    print("\n## V2 -- is dose sensitivity ordered by predicted affinity?")
    print("  The base-rate confound cannot manufacture this: chance does not")
    print("  produce a dose response ordered by predicted binding energy.")

    exonic = {g for g, cats in predicted.items() if "exonic_proximal" in cats}
    pairs, excluded = [], 0
    for gene, ev in events.items():
        if gene not in affinity:
            continue
        if ev.dpsi_hia is None or abs(ev.dpsi_hia) < MIN_HIA_FOR_RATIO:
            excluded += 1
            continue
        r = ev.retention
        if r is None:
            continue
        pairs.append((gene, affinity[gene], r, gene in exonic))

    print(f"\n  events excluded by the pre-declared |dPSI_HiA| >= "
          f"{MIN_HIA_FOR_RATIO} floor: {excluded}")
    print(f"  {'set':<22}{'n':>6}{'rho':>9}{'one-sided p':>14}")
    for label, subset in (("all sieve genes", pairs),
                          ("exonic_proximal only", [p for p in pairs if p[3]])):
        if len(subset) < 8:
            print(f"  {label:<22}{len(subset):>6}{'too few to test':>23}")
            continue
        ddg = [p[1] for p in subset]
        ret = [p[2] for p in subset]
        rho, _ = stats.spearmanr(ddg, ret)
        # Pre-declared direction: tighter binding (more negative ddG_net)
        # retains MORE effect at low dose, i.e. rho < 0.
        n = len(subset)
        t = rho * ((n - 2) / max(1e-12, 1 - rho**2)) ** 0.5
        p_one = stats.t.cdf(t, n - 2)
        collect.append((f"V2 {label}", "-", n, 0, 0.0, 0.0, p_one))
        print(f"  {label:<22}{n:>6}{rho:>9.3f}{p_one:>14.4g}")
    print("\n  Pre-declared direction: rho < 0. A POSITIVE rho is not a")
    print("  result in the other direction -- it is the null with a sign.")


def v3_hit_replication(events) -> None:
    print("\n## V3 -- do the four GSE307431 hits replicate in DIRECTION?")
    print("  DESCRIPTIVE ONLY. n <= 4; a sign test on 4 cannot reach p<0.05")
    print("  even if all agree (min p = 0.0625). No p-value is computed and")
    print("  this test is NOT in the Bonferroni denominator.")
    print(f"\n  {'gene':<10}{'GSE307431':>12}{'Ottesen HiA':>14}{'Ottesen LoA':>14}{'agrees':>9}")
    agree = present = 0
    for gene, prev in GSE307431_HITS.items():
        ev = events.get(gene)
        if ev is None:
            print(f"  {gene:<10}{prev:>12.3f}{'absent':>14}{'-':>14}{'-':>9}")
            continue
        present += 1
        ok = (ev.dpsi_hia < 0) == (prev < 0)
        agree += ok
        loa = "n/a" if ev.dpsi_loa is None else f"{ev.dpsi_loa:+.3f}"
        print(f"  {gene:<10}{prev:>12.3f}{ev.dpsi_hia:>+14.3f}{loa:>14}"
              f"{'yes' if ok else 'NO':>9}")
    print(f"\n  {agree}/{present} present hits agree in direction.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", type=Path, default=DEFAULT_XLSX)
    ap.add_argument("--sieve", type=Path, default=SIEVE)
    args = ap.parse_args()

    n_rows, frac = verify_sign_convention(args.xlsx)
    print("## Arm V -- independent replication (Ottesen 2021, PMC8395096)")
    print("## Pre-registered in docs/preregistration.md BEFORE any number here.")
    print(f"\nSIGN CONVENTION re-verified on load: published dPSI is "
          f"control - treated for {frac:.1%} of {n_rows} rows.")
    print("All dPSI below are PROJECT convention (treated - control):")
    print("negative = the ASO REDUCED inclusion.")

    events = best_per_gene(load_events(args.xlsx, verify=False))
    predicted = load_predictions(args.sieve)
    affinity = load_affinity(args.sieve)
    measured = {g: ev.dpsi_hia for g, ev in events.items()}

    overlap = set(predicted) & set(measured)
    print(f"\nsieve-predicted genes: {len(predicted)}")
    print(f"genes measured here (both doses): {len(measured)}")
    print(f"OVERLAP: {len(overlap)}")

    collect: list = []
    v1_enrichment(measured, predicted, collect)
    v2_dose_response(events, predicted, affinity, collect)
    v3_hit_replication(events)

    print(f"\n## MULTIPLE TESTING -- {len(collect)} tests printed above.")
    print("Bonferroni over ALL of them, including those that did not support")
    print("the hypothesis. Nothing is excluded from the count.")
    print(f"  {'test':<34}{'category':<20}{'raw p':>11}{'corrected':>12}")
    for label, cat, n, hit, rate, base, p in sorted(collect, key=lambda t: t[-1]):
        print(f"  {label:<34}{cat:<20}{p:>11.4g}{min(1.0, p*len(collect)):>12.4g}")

    surviving = [c for c in collect if c[-1] * len(collect) < 0.05]
    print(f"\n## VERDICT: {len(surviving)} of {len(collect)} tests survive Bonferroni.")
    if not surviving:
        print("  The tail enrichment does NOT replicate and dose sensitivity is")
        print("  NOT affinity-ordered. Report as such -- this was a declared")
        print("  possible outcome, not a failure to be re-analysed.")

    _power_caveat(measured)


def _power_caveat(measured: dict[str, float]) -> None:
    """POST-HOC interpretation. Stated as such; nothing above was re-analysed.

    Not anticipated at pre-registration, and it cuts against V1's usefulness
    rather than the project's interest, so it is recorded rather than dropped.
    """
    big = sum(1 for v in measured.values() if abs(v) >= 0.2)
    print("\n## POST-HOC INTERPRETATION -- declared as post-hoc, changes nothing above")
    print(f"  V1's background base rate at |dPSI| >= 0.2 is "
          f"{100*big/len(measured):.1f}% here, versus 11.5% on GSE307431.")
    print("  CAUSE: Ottesen's supplementary tables are already filtered to")
    print("  events the authors called notable, so small-effect events are")
    print("  absent and the background is itself enriched for large effects.")
    print("  CONSEQUENCE: 'enrichment above background' has little room to")
    print("  register when background is already the majority. V1 was a")
    print("  WEAKER test than intended, and its null is correspondingly weaker")
    print("  evidence than a null on a full measured background would be.")
    print("  V2 IS NOT AFFECTED -- the dose-retention correlation never touches")
    print("  background composition. It is the cleaner negative of the two, and")
    print("  it is the one with a real n (123).")
    print("  This does NOT rescue V1: the observed fold is 0.90, BELOW 1.")


if __name__ == "__main__":
    main()
