# Post-mortem: why was each null null?

Every null in this project has been recorded. None has been asked the question
that decides what it means: **was the hypothesis false, or was the test unable
to see it?** Those have opposite consequences. A decisive null is knowledge and
closes a question. An underpowered null closes nothing and must not be cited as
evidence of absence.

Power computed at alpha = 0.05, 80% power, one-sided where the prediction was
directional.

## Classification

| arm | n | detectable effect | observed | verdict |
|---|---|---|---|---|
| **N4** sign rule | 3975 | 1.76 pp lift over base | 60 SD BELOW base | **DECISIVE** |
| **A5** global sieve | 426 | \|rho\| >= 0.120 | 0.0025 | **DECISIVE** |
| **V2** dose, all genes | 123 | \|rho\| >= 0.223 | -0.007 | PARTIAL |
| **V1** tail enrichment | 18 | fold >= 1.49 | 0.90 | **INVALID TEST** |
| **V2** dose, exonic_prox | 18 | \|rho\| >= 0.566 | +0.020 | **UNINFORMATIVE** |
| **G1** off-target panel | 24 | \|rho\| >= 0.495 | n/a | UNDERPOWERED |
| **S2** ceiling in 3D | 4 | \|rho\| >= 0.986 | -0.200 | **UNINFORMATIVE** |
| **C4 / C / K** | n/a | structural measurement | 98% upstream | **DECISIVE** |
| **S5** gate | 4 | n/a -- tool capability | 0/2 drug-like | TOOL LIMIT |
| **M1 / M2** | n/a | structural measurement | inverted | **INFORMATIVE, not null** |

## The two findings that change what we may say

### V1 could not have replicated the original effect IN PRINCIPLE
V1 tried to replicate GSE307431's **5.97x** exonic_proximal tail enrichment.
Ottesen's supplementary tables are pre-filtered to events the authors called
notable, so V1's background base rate is **55.8%** rather than GSE307431's
11.5%. The maximum achievable fold is therefore **1/0.558 = 1.79x**. A 5.97x
effect is arithmetically unreachable.

**V1 is not a failed replication. It is an invalid replication design**, and it
must not be cited as evidence against the enrichment. The honest statement is
that the enrichment has still been tested exactly once.

The fix, if the question is worth re-asking: a background of ALL measured genes
rather than only the authors' notable ones -- which requires the raw counts,
not the supplementary tables.

### N4 is decisive and closes its question
n = 3975 cells, powered to detect a 1.76 percentage-point lift. Observed
accuracy 29.6% against a 72.3% base rate -- 60 SD below. No plausible effect
was missed. The 90%-downstream rule does not predict the sign of delta_A.

OPEN, and to be resolved at the next session: the rule's accuracy (29.6%) is
close to the MINORITY class rate (27.7%), which suggests it is behaving as a
near-constant "+" predictor rather than being anti-correlated. **The confusion
matrix must be computed before any claim about anti-correlation is made.**

## Completed 2026-08-25 — the three remaining nulls, plus two audit findings

### S1 — genuinely null, adequately powered
Sign test 10/16, p = 0.2272; it would have needed 12/16 (75%) to reach p<0.05.
Powered for a large effect and did not see one. **DECISIVE for a large effect,
silent about a small one.**

### P — NOT a power failure. A strict self-imposed threshold.
    hit 1 upstream  z = +2.57   one-sided p = 0.0051
    hit 1 partner   z = +2.41   one-sided p = 0.0080
    pre-declared threshold z >= 2.58  =  p <= 0.0049
At a conventional alpha = 0.05 (z >= 1.645) **both arms clear it easily**. The
"NOT ESTABLISHED" verdict is an artifact of a threshold set roughly ten times
stricter than convention. The threshold was declared in advance and is NOT
being moved — but any citation of this null must say that the bar was 0.005,
not 0.05, or it misrepresents the evidence.

### A6 — a real confound found, and the hit survives it
The intron-length confound was never quantified. It is real:
**Spearman(intron length, retention) = -0.833, p = 0.0102** — shorter introns
are retained more.

But intron 7 beats the trend decisively. Intron 3 is the SHORTEST (159 nt) and
retains 0.421; intron 7 is **2.8x longer** (444 nt) and retains 1.000. A
length-only model fit on the other seven introns predicts intron 7 = 0.274;
observed 1.000, **residual +0.726**. The blind hit stands, now quantified
rather than asserted. Intron 3 remains a weak test for the reason already in
the record.

### AUDIT FINDING — a bounded scan reported as universal
"The nusinersen 18-mer has exactly two exact-match sites in 140.6 Mnt" is
**arithmetically correct** (72,863,420 nt mRNA + 67,709,953 nt intron flanks =
140,573,373). But the index covers only **3.8% of annotated intronic space** —
~1.71 Gnt of intron interior was never searched, and `transcriptome.py`'s own
docstring declares that assumption. An ASO acts on pre-mRNA, so the uniqueness
claim is bounded by the index and is not a pre-mRNA-wide statement.

It is also a NEAR-MATCH claim: axis 1 counted exact and 1-mismatch sites, while
POLR2H — the validated off-target at dPSI -0.920 — sits at **four** mismatches.

### AUDIT FINDING — a published figure that was an estimate
The run-log artifact carried "~120k partition functions". The record supports
~105k (A3 10,650 + A7 89,100 rows, 6,300 of which carried empty scores, +
N4 5,325 + smaller arms). Corrected in the artifact.

## Still to do

- Confusion matrix for N4 (above).
- V2's PARTIAL verdict needs stating precisely: it rules out |rho| >= 0.223 but
  a real dose effect of |rho| ~ 0.15 would have been missed. The dose test is
  therefore not yet decisive and should not be described as one.
- S1 (p = 0.4545, 16/20) and P (z = +2.57/+2.41 vs threshold 2.58) are not yet
  in this table and need the same treatment. P in particular sits JUST under
  its pre-declared threshold, which is a power question, not a verdict.
- A6 intron 3 (0.421) -- already noted in the record as a weak test by
  construction, since a 159 nt intron is nearly covered by the published
  design's own flanking rule. Confirm and add.

## The rule this establishes

An arm reported as a null must carry its **detectable effect size** alongside
its p-value. "No effect" and "no power" are different claims and the project
has been writing both the same way.
