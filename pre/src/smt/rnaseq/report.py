"""R3: measure SMN2 exon 7 inclusion from the GSE307431 BAMs.

Run: ``python -m smt.rnaseq.report`` (needs the rna-env interpreter, which is
where samtools lives).

ORDER MATTERS, and it is the point of this module. Three things happen, and
the first two are permitted to stop the third:

  1. c.840 CENSUS over all 20 samples, both loci. Pooling SMN1 with SMN2 is
     only legitimate if SMN1 is absent. That is a premise, so it is measured
     first and refused loudly if it fails.
  2. POSITIVE CONTROL. POLR2H must reproduce rMATS's dPSI of -0.920 through
     this exact counting code. If it does not, the counter is wrong and the
     SMN2 number means nothing -- so the SMN2 number is not printed.
  3. SMN2 exon 7 PSI, pooled across both loci, F18MOE vs ScrMOE.

A PSI of exactly 1.000 is a BOUNDARY, not a measurement: it means zero
skipping junction reads were seen, which is not the same as zero skipping.
Those samples are reported with a binomial upper bound on the skipping
fraction so the depth is visible next to the number.
"""

from __future__ import annotations

import math
from collections import Counter
from pathlib import Path

from smt.rnaseq import targets as T
from smt.rnaseq.junction_psi import (
    SMN1PresentError, base_census, count_junctions, pooled_counts, psi,
    read_length,
)
from smt.seq.reference import Interval

#: Sample accessions, in the order slurm/sra_runs.txt lists them.
ALL_SAMPLES = tuple(f"SRR{35298699 + i}" for i in range(20))

#: Tolerance on the POLR2H control. This counter includes multimappers and
#: applies no anchor-length filter, while rMATS uses unique reads only, so
#: the counts differ by ~20% and the PSIs by a few points. 0.05 accepts that
#: while still failing on any counter that has the event wrong.
CONTROL_TOL = 0.05


def wilson_upper(k: int, n: int, z: float = 1.96) -> float:
    """Upper 95% bound on a binomial proportion (Wilson). ``n == 0`` -> 1.0."""
    if n == 0:
        return 1.0
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return min(1.0, centre + half)


def c840_census() -> Counter:
    """Base at c.840 across every sample and BOTH loci, multimappers included."""
    total: Counter = Counter()
    print("c.840 census -- C means SMN1 is present, T means SMN2")
    for acc in ALL_SAMPLES:
        bam = Path(T.bam_path(acc))
        per: Counter = Counter()
        for cas, pos in ((T.SMN2_EXON7, T.C840_SMN2), (T.SMN1_EXON7, T.C840_SMN1)):
            per += base_census(bam, cas.chrom, pos, cas.cassette)
        total += per
        print(f"  {acc}  C={per['C']:>4}  T={per['T']:>4}  A={per['A']:>3}  G={per['G']:>3}")
    print(f"  POOLED  C={total['C']}  T={total['T']}  A={total['A']}  G={total['G']}")
    return total


def group_psi(cassettes, accessions) -> tuple[list[float], list[tuple[str, int, int, float]]]:
    rows, values = [], []
    for acc in accessions:
        bam = Path(T.bam_path(acc))
        counts = pooled_counts(bam, cassettes) if len(cassettes) > 1 \
            else count_junctions(bam, cassettes[0])
        length = read_length(bam, cassettes[0].region)
        value = psi(counts, cassettes[0], length)
        if value is None:
            raise ValueError(f"{acc}: no junction reads at all -- cannot form a PSI")
        rows.append((acc, counts.ijc, counts.sjc, value))
        values.append(value)
    return values, rows


def positive_control() -> float:
    print("\npositive control -- POLR2H, which rMATS calls at dPSI -0.920")
    f18, rows_a = group_psi((T.POLR2H_EXON3,), T.F18MOE)
    scr, rows_b = group_psi((T.POLR2H_EXON3,), T.SCRMOE)
    for label, rows in (("F18MOE", rows_a), ("ScrMOE", rows_b)):
        for acc, ijc, sjc, value in rows:
            print(f"  {label}  {acc}  IJC={ijc:>4}  SJC={sjc:>4}  PSI={value:.3f}")
    dpsi = sum(f18) / len(f18) - sum(scr) / len(scr)
    delta = abs(dpsi - T.POLR2H_RMATS_DPSI)
    print(f"  dPSI = {dpsi:+.3f}   rMATS = {T.POLR2H_RMATS_DPSI:+.3f}   |diff| = {delta:.3f}")
    if delta > CONTROL_TOL:
        raise RuntimeError(
            f"POLR2H control failed: {dpsi:+.3f} vs rMATS {T.POLR2H_RMATS_DPSI:+.3f}. "
            "The counter is wrong; any SMN2 number it produces is meaningless."
        )
    print("  CONTROL PASSES")
    return dpsi


def smn2_exon7() -> float:
    print("\nSMN2 exon 7, pooled across both SMN loci")
    f18, rows_a = group_psi(T.SMN_EXON7_LOCI, T.F18MOE)
    scr, rows_b = group_psi(T.SMN_EXON7_LOCI, T.SCRMOE)
    for label, rows in (("F18MOE", rows_a), ("ScrMOE", rows_b)):
        for acc, ijc, sjc, value in rows:
            bound = wilson_upper(sjc, ijc + sjc)
            note = f"  (skipping <= {bound:.3f} at 95%)" if sjc == 0 else ""
            print(f"  {label}  {acc}  IJC={ijc:>4}  SJC={sjc:>4}  PSI={value:.3f}{note}")
    dpsi = sum(f18) / len(f18) - sum(scr) / len(scr)
    print(f"  mean F18MOE {sum(f18)/len(f18):.3f}   mean ScrMOE {sum(scr)/len(scr):.3f}")
    print(f"  dPSI (F18MOE - ScrMOE) = {dpsi:+.3f}")
    return dpsi


def main() -> None:
    census = c840_census()
    if census["C"]:
        raise SMN1PresentError(
            f"{census['C']} reads carry C at c.840, so SMN1 is expressed here. "
            "The two loci may NOT be pooled, and the exon 6->8 skipping junction "
            "carries no paralog-discriminating base, so the exclusion arm is not "
            "measurable from these reads. Report that, do not force a number."
        )
    print("  SMN1 absent -- pooling the two loci is legitimate for this dataset")
    positive_control()
    smn2_exon7()


if __name__ == "__main__":
    main()
