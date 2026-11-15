# Arm E StericRuler length-titration report

> Claim boundary: predictive, pre-registration-eligible, and hypothesis-driven. Every value in this package is a sequence-derived model output tagged `candidate`. Nothing here is a measurement or a wet-lab outcome.

## What Arm E claims

- [candidate] Over a symmetric truncation/extension ladder of the true nusinersen 18-mer (even lengths, centred on the ISS-N1 target `chr5:70076584-70076601`, GRCh38 v50), Arm E emits one pre-declared **therapeutic window** value per length: predicted on-target potency divided by predicted off-target burden.
- [candidate] The predicted window-vs-length curve and the directional hypothesis are captured in a preregistration payload whose SHA-256 is recorded here, to be sealed and compared against a tandem wet-lab length series.
- [candidate] Directional hypothesis under test: a shorter mer improves the therapeutic window: the predicted potency-over-off-target-burden ratio favours a length below the 18-nt nusinersen anchor.

## Methods

- [candidate] Ladder frame: each even length `L` takes the `L`-nt window centred on the 18-nt target, trimming (`L<18`) or extending (`L>18`) one nucleotide per end; the ASO is the reverse complement of that genomic window. The `L=18` rung reproduces nusinersen exactly (fail-closed in-code assertion).
- [candidate] On-target axes: IntaRNA duplex ΔG and RNAup accessibility-corrected ΔG form the potency backbone. IntaRNA's `no_duplex` (0.0) and sub-seed `None` readings are never averaged as real binding. RNAduplex/RNAcofold (a perfect reverse-complement duplex by construction) and RNAplfold accessibility (already folded into RNAup) are carried as secondary/diagnostic columns and are excluded from the backbone to avoid double-counting.
- [candidate] Off-target axis: near-complement hit count from the seed-and-extend scan. A missing scan is reported as `None`, never as zero burden.
- [candidate] Endpoint definitions (all candidate, pre-registration-eligible):

- [candidate] **potency**: negative mean of the usable on-target binding energies (IntaRNA on-target dg and RNAup accessibility-corrected dg); IntaRNA no_duplex/None excluded; RNAduplex/RNAcofold/RNAplfold accessibility excluded from the backbone (secondary/diagnostic)
- [candidate] **off_target_burden**: count of near-complement off-target hits; None (not zero) when the off-target scan was not run
- [candidate] **therapeutic_window**: potency / (off_target_burden + OFF_TARGET_PSEUDOCOUNT); OFF_TARGET_PSEUDOCOUNT=1

- [candidate] Coverage: 15 rungs scored; 0 have a computable therapeutic window (the rest await the off-target scan and/or GPU potency axes).
- [candidate] Preregistration formula version `armE-window-v2`; payload SHA-256 `b75de7356db19ad1a58b0e9c03e6443db6614a9f25e11612d4840eca7e2ba978`. The lock is an unsealed DRAFT.
- [candidate] Software: Python 3.11.16, NumPy 1.26.4, pandas 3.0.5, matplotlib 3.10.9.
- [candidate] File-access audit: this generator consumes only the in-memory candidate window rows. It reads no sealed material and reuses no earlier-arm code.

## Per-mer summary [candidate]

| length | aso_5to3 | intarna_dg | intarna_no_duplex | rnaup_dg | potency | off_target_burden | therapeutic_window | window_incomplete |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 2 | AT | nan | None | nan | nan | None | None | True |
| 4 | CATA | nan | None | 0.1 | nan | None | None | True |
| 6 | TCATAA | nan | None | -2.8 | nan | None | None | True |
| 8 | TTCATAAT | -4.8 | False | -4.8 | 4.8 | None | None | True |
| 10 | TTTCATAATG | -8.29 | False | -8.28 | 8.285 | None | None | True |
| 12 | CTTTCATAATGC | -14.18 | False | -14.18 | 14.18 | None | None | True |
| 14 | ACTTTCATAATGCT | -17.46 | False | -17.45 | 17.455 | None | None | True |
| 16 | CACTTTCATAATGCTG | -22.4 | False | -22.39 | 22.395 | None | None | True |
| 18 | TCACTTTCATAATGCTGG | -27.0 | False | -26.99 | 26.994999999999997 | None | None | True |
| 20 | TTCACTTTCATAATGCTGGC | -31.27 | False | -31.27 | 31.27 | None | None | True |
| 22 | ATTCACTTTCATAATGCTGGCA | -33.84 | False | -33.83 | 33.835 | None | None | True |
| 24 | GATTCACTTTCATAATGCTGGCAG | -38.94 | False | -38.94 | 38.94 | None | None | True |
| 26 | AGATTCACTTTCATAATGCTGGCAGA | -42.19 | False | -41.79 | 41.989999999999995 | None | None | True |
| 28 | AAGATTCACTTTCATAATGCTGGCAGAC | -45.67 | False | -42.36 | 44.015 | None | None | True |
| 30 | TAAGATTCACTTTCATAATGCTGGCAGACT | -48.47 | False | -42.27 | 45.370000000000005 | None | None | True |

## Figure caption

- [candidate] `arm_e_length_axis`: on-target potency, therapeutic window, and off-target burden over the symmetric nusinersen length ladder. The dotted line marks the 18-nt drug anchor. Panels whose axis has not been run are annotated as pending; values are model outputs, not measurements.

## Limitations

- [candidate] CPU-only potency (duplex energy + accessibility) rises with length by construction; the directional test lives in the full therapeutic window (÷ off-target burden) and the GPU potency axes, not in raw binding energy.
- [candidate] RNAduplex/RNAcofold score a perfect reverse-complement duplex and carry little signal beyond composition; they are diagnostic only.
- [candidate] The off-target burden counts near-complement sequence sites and does not model expression, chromatin, or delivery; it is a sequence-specificity proxy.

## What Arm E does not claim

- Arm E does not claim that any predicted therapeutic-window value is a measured potency, a measured specificity, or a wet-lab outcome.
- Arm E does not claim that the sequence-derived CPU scores (duplex energy, accessibility, self-structure, off-target count) establish mechanism, causality, or efficacy.
- Arm E does not claim a clinical response, a dosing recommendation, or that any rung is safe or effective in cells, animals, or patients.
- Arm E does not claim that a shorter or longer oligo is superior to nusinersen; the ladder is scored on its own pre-declared axes only.
- Arm E does not claim the preregistration is sealed: the reported lock is an unsealed DRAFT until the GPU potency axes are included and the seal is cut with explicit sign-off.
- Arm E does not claim the directional prediction is confirmed; confirmation requires the tandem wet-lab length series, read only after the seal.
