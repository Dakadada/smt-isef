# Discriminating the C2B34 mechanism — proposed experiment

## Why this experiment exists

Nusinersen treatment is reported to raise the C2B34 circRNA up to ~1000-fold in
a subset of SMA patients, correlating with better motor outcomes. The obvious
explanation does not work: C2B34 spans exons 2a–2b–3–4, which are transcribed
*before* intron 7. An ASO bound at ISS-N1 cannot retroactively change a splicing
decision that already happened upstream.

So rather than assert a mechanism, the kinetic model encodes three and asks
which experiment tells them apart.

| | Mechanism | Defining signature |
|---|---|---|
| **H1** | Restored SMN feeds back on snRNP biogenesis, shifting spliceosome availability globally | **Lag.** Near-zero response immediately after dosing; depends on accumulated SMN protein, not ASO occupancy |
| **H2** | At high concentration the ASO engages upstream sites and alters Alu duplex formation directly | **No lag.** Dose-response tracks occupancy immediately |
| **H3** | Backsplicing rate is unchanged; circles accumulate because linear decay changes | **Flat synthesis flux.** `k_back × N` does not move with dose |

## What the model says about separability

Power to correctly identify the true mechanism, from simulated experiments at
6 doses × 8 timepoints, 200 replicates:

| σ (relative measurement error) | H1 | H2 | H3 |
|---|---|---|---|
| 2%  | 1.000 | 1.000 | 0.995 |
| 5%  | 1.000 | 0.935 | 0.935 |
| 10% | 1.000 | 0.800 | 0.780 |
| 20% | 1.000 | 0.620 | 0.655 |
| 40% | 0.970 | 0.620 | 0.580 |

**H1 is easy.** Its lag is unmistakable — 97% power even at 40% measurement
error. A coarse time course settles whether the effect is feedback-driven.

**H2 vs H3 is the hard pair.** Separating them needs σ ≲ 5% for >90% power.
That is out of reach for gel densitometry and within reach of ddPCR or qPCR.

Rate constants in the model are symbolic; the numbers above come from an
illustrative parameter set explicitly labelled as *not* a measurement. They
indicate the shape of the requirement, not its exact value.

## Proposed protocol

**Time course** (this alone resolves H1). Fixed dose, harvest at roughly
0.5, 1, 2, 5, 10, 50, 200, 1000 in model time units — mapped to real hours by
the transfection/treatment kinetics of the system in use. A near-flat early
response with a later rise indicates H1.

**Dose-response** (needed for H2 vs H3). At least 6 doses spanning zero to
saturating, harvested at a late timepoint where the system has settled.

**Measure separately, in the same RNA prep:**
- C2B34 circular species, via divergent primers across the backsplice junction
- Linear SMN2 transcript, exon-7-included and exon-7-skipped
- SMN protein (H1 predicts the circle response tracks protein, not dose)

**Critical prep constraint:** circles are **not polyadenylated**. Oligo-dT
priming will miss them entirely. Use random hexamers or gene-specific primers.
This is the same constraint already flagged in `docs/lab-request.md`, Ask 4.

**To test H3 specifically:** measure circle synthesis flux and decay separately
— a pulse-chase, or an actinomycin D chase to isolate decay. H3 predicts
synthesis flux is unchanged by dose while decay shifts.

## What each result would mean

- Lag present, response tracks SMN protein → **H1**
- No lag, response tracks dose immediately, synthesis flux rises → **H2**
- No lag, synthesis flux flat, decay changes → **H3**
- None of these → the model is missing a mechanism, which is itself worth knowing

Reproduce the table above with
`smt.kinetics.hypotheses.separability_power(...)`.
