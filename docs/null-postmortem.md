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
