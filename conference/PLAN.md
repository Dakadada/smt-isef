# Conference program: computational prediction of splice-switching ASO behavior

> Plan of record for two manuscripts (IEEE BIBE 2026, BIOS 2026). Scope = **ASO only**
> (dCas13d / SpliceRUSH / fusion-protein / cross-modality all EXCLUDED). New work in `conference/`.

## Context & locked framing

The research contribution **is the computational model**. Wet-lab experiments are an **independent
held-out validation layer**. ASO size/length is **established background, not a novelty claim**
(Singh 2009 etc.) — it is used as a *controlled perturbation axis* to stress-test the model.

Paper narrative order (both papers): **computational problem → model development → in-silico
candidate evaluation → predictions → validation on independent wet-lab data.** Never "experiments
showed X, so we explained it"; always "the framework generated predictions, which were then tested
against super-minigene and endogenous experiments."

**Central objective:** *Can an integrated computational framework predict the efficacy and
specificity of splice-switching ASOs across sequence, target position, ASO length, and biological
context?* — and does integrating features beat sequence-complementarity / binding-energy-only?

**The headline modeling problem = context transfer.** The size series measured in **both
super-minigene and endogenous** SMN is the key validation matrix: does a relationship learned in the
simplified reporter survive the move to the endogenous transcript? Report ASOs as:
**context-stable / context-amplified / context-dependent** — far stronger than separate PSI values.

- **BIBE (primary — methods):** *Structure-Aware Computational Prediction of Splice-Switching ASO
  Efficacy and Specificity Across Molecular and Biological Contexts.* Ablation ladder: seq-only →
  +thermo → +structure/accessibility → +transcriptome specificity → +positional → full; then test
  on experiment; then test super-minigene→endogenous transfer.
- **BIOS (biology):** *Which sequence/structural/positional properties determine whether an
  ASO-mediated splice effect is stable across biological context?* Landscape = target position ×
  ASO seq/length × RNA structure × context → measured outcome. The positional (E3) series
  integrated with the context comparison; biology = why some regions are robustly modulatory and
  others highly context-dependent.

Feature set to integrate (all defensible): ASO sequence, length, target position, local RNA
sequence, predicted target accessibility, RNA secondary structure, hybridization ΔG, GC%, mismatch
architecture, transcriptome-wide partial complementarity, splice-site proximity, local regulatory
context, and super-minigene-vs-endogenous context.

---

## Data reality & dependencies

Excluding the dCas13d screen removes the only *dense* perturbation set → this is a **mechanistic-
feature + public-data framework externally validated on a small wet-lab series**, NOT a data-hungry
deep net. Foundation-model embeddings (RiNALMo/LAMAR) → **optional stretch only** (overfit risk;
absent from repo, Explore-confirmed).

| Dependency | Role | Status |
|---|---|---|
| **Wet-lab densitometry (Singh Lab)** — 10 CSVs at repo root | **Held-out validation** | **IN HAND** but raw. Endogenous size: `...7-28 ASO size comparison...endogenous...csv`, `...7-27 truncated ASOs endogenous...csv`. Super-minigene size: `...7-7-26 ASO 5nM super minigene.csv`, `...7-9-26 ASO 5nm super minigene.csv`. Baseline: `...8-1011-26...no ASO treatments.csv`. **Excluded:** dCas13d (`7-12`/`7-17`) + fusion (`8-4`). Needs reconstruction (below). |
| **Lane→sample→ASO map + amplicon sizes** | Turns densitometry into PSI | **MISSING — need from user.** CSV `Label`/`MW` cols are all `-`; no sample map inside. Blocks band-identity locking. |
| **Public RNA-seq reanalysis PSI** (`rawread_reanalysis`; GSE307431 etc.) | In-repo public train/eval PSI | PRESENT (`artifacts/judge-demo/.../evidence.json`). |
| **Holgersen 2021 S1/S4–S6** (81 PS-MOE 16–20 nt) | Off-target benchmark | Not in repo. Default = obtain if accessible, else public rMATS. |
| **Literature walks** — Hua 2007 exon-7, Singh 2006/2009 length | Positional/length known-answers | Partial (`harmonized_aso_positions.csv`, qualitative); PDFs need numeric extraction. |

Context transfer is now fully covered by data on disk: **super-minigene (7-7, 7-9) ↔ endogenous
(7-27, 7-28)**, both ASO size series. Replicate/rerun status of the two dates per context is unknown
until the sample map resolves it.

**Densitometry reality (from reading the flagship CSVs):** columns = `Lane#,Band#,Label,Peak Value,
Peak Rf,MW,Raw Volume,Vol Ratio,Cal. Volume`. `Label` and `MW` are unset (`-`); lane 1 carries
~17–19 bands (likely ladder) vs 5–10 in sample lanes; recurring products ~Rf 0.56/0.40 (7-28) and
~0.57/0.41 (7-27) with extra peaks in some lanes. So **band identity must be locked from expected
amplicon size + ladder calibration, never top-2-by-Rf.** 7-27 and 7-28 share "samples 1-9" on
consecutive days → treat replicate/rerun status as **unknown until the sample map resolves it; do
not count as independent n.**

**Prior negative result = an asset:** A5 sieve validated on real rMATS gives **POLR2H −0.920,
MAP1S −0.879 blind**, but **fails genome-wide (ρ=0.0025)** → specificity head is framed as a
**strong-exonic-off-target candidate generator, not a genome-wide oracle.** And **SMN2 on-target
ΔPSI is not short-read-rMATS measurable** (SMN1/SMN2 ~99.9% identical) → on-target ground truth =
**minigene + literature**, which is exactly why the minigene validation context matters.

---

## Reuse map (do NOT rebuild — pinned by Explore). Real science lives in `pre/` (`smt`) + `feng/`

**BIOS atlas / occlusion engine (mostly reuse):**
- `pre/src/smt/thermo/constrained.py` — `a_5ss()`, `occlusion_delta()`, `scan_footprints()` (5′ss
  accessibility under occlusion). Extend the scan across the whole region = the atlas.
- `pre/src/smt/thermo/occlusion.py` — `score()`, `variant_effect()` (signed FOM).
- `pre/src/smt/thermo/link.py` — `psi_from_a5ss()`, `delta_psi()` (**fittable a_5ss→ΔΨ logistic
  link; refuses if unfit**) — the physics-based efficacy signal to calibrate.

**Efficacy signal (feng):** `feng/src/feng/titration/gpu_scores.py` — `pangolin_exon7_delta()`
(on-target ΔΨ), `pangolin_off_target_burden()`.

**Feature/thermo/structure:** `feng/.../titration/evaluate.py` (`rnaduplex_dg/rnacofold_dg/rnaup_dg/
intarna_on_target/rnafold_mfe/rnaplfold_accessibility`); `pre/.../thermo/{hybridization,
accessibility}.py`. **No Tm scorer anywhere → new (trivial, nearest-neighbor).**

**Length ladder:** `feng/.../titration/{ladder.py:build_ladder, evaluate.py:evaluate_ladder,
report.py:build_arm_e_report}` (schema `arm_e_per_mer_table`); `pre/.../validate/length_sweep.py` +
`goldilocks.py` (span-robust read).

**Specificity / sieve:** `pre/src/smt/sieve/{transcriptome,scan,polarity,enrichment}.py`
(ddG_net, category, signed_effect; `enrichment.py` already Fishers predictions vs measured rMATS);
feng `scan/seed_extend.py`. Indexes live in `orcd/pool` / scratch — **do not scan.**

**Coordinates (two frames to unify):** pre `smt/seq/{reference,construct}.py` (NG_008728 space;
ISS_N1, GC_CORE, SS5_EXON7, c.840); feng genomic GRCh38 anchors in `gpu_scores.py`/`ladder.py`.
`feng/src/feng/atlas/arm_a.py:window_alignment()` already harmonizes ASO windows.

**Minigene:** `pre/src/smt/minigene/{compress.py:compress_gene, spec.py, report.py}`.

**The learned-model home is an empty stub:** `feng/src/feng/gate/{features,model,power}.py` and
`feng/src/feng/scan/intarna.py` are 2-line placeholders → natural home for the ΔPSI head. sklearn is
**NOT installed** (per `AGENTS.md`); a learned head needs it added, else numpy/physics-link only.

**Caution:** repo-root `src/` (`design_engine` Pareto, `benchmark/harness` baselines,
`rawread_reanalysis`) is a separate effector-design/benchmark layer (includes excluded modalities) —
reuse `rawread_reanalysis` (public PSI) and the baseline-contract idea selectively; do not treat it
as the core.

## Genuinely new work
Unified SMN2 coordinate frame; nt-resolution occlusion **atlas** (extend `scan_footprints`); learned
**ΔPSI efficacy head** (or calibrated physics link) in the feng `gate` stub; **specificity head**
(reframed sieve); **efficacy×specificity Pareto**; **counterfactual ±nt/±length scan**;
**context-transfer classifier** (stable/amplified/dependent); Tm feature; (stretch) uncertainty,
generative/Bayesian-opt design, FM embeddings.

---

## Phasing (MVP → ceiling; submittable at every gate)

- **P0 — Unified coordinate frame.** Reconcile pre/feng coordinate systems; place all ASOs (size
  series, E3 walk, literature) + regulatory elements on one axis. *Now, no new data.*
- **P1 — Occlusion atlas (BIOS core).** Extend `scan_footprints`/`a_5ss` to signed susceptibility
  S(x) at nt resolution; overlay structure/accessibility/RBP to explain discontinuities;
  regression-guard vs A3 GC-core. *Now.*
- **P2 — Feature matrix + efficacy head (BIBE core).** Assemble the full feature set; fit ΔPSI head
  on **public reanalysis + literature, blind to wet-lab**; ablation ladder vs baselines.
- **P3 — Specificity head + off-target landscape** (reframed sieve; seq-only → +thermo → +structure
  → +splice-context ablation).
- **P4 — Pareto + counterfactual** (efficacy×specificity; size series sampled on it; ±nt/±length
  neighborhoods explaining why each tested ASO sits where it does).
- **P5 (stretch)** — uncertainty; generative/BO design; FM embeddings.
- **Pexp — Experimental reconstruction (parallel track; gated on the sample map).** Parse the 10
  densitometry CSVs → apply the user's lane→sample→ASO map → **lock inclusion/skipping band identity
  by expected amplicon size + ladder calibration** (never Rf/top-2) → compute **replicate-level**
  `PSI = incl/(incl+skip)×100` and ΔPSI vs the matched control → keep replicates separate, mark
  reruns. Produces the MASTER band-identity + PSI tables (raw and derived kept separate; formula
  recorded per value). Nothing here feeds model fitting — it is the sealed test set.
- **Validation gate (held-out).** Evaluate **frozen** models on the reconstructed PSI; true external
  test if n permits, else **leave-one-ASO-out / leave-one-condition-out**. Headline =
  super-minigene→endogenous transfer + the stable/amplified/dependent classification.

## Validation rigor (non-negotiable)
Never fit to all experimental points then claim they validate. Freeze features + splits + the
a_5ss→ΔΨ fit **before** reading wet-lab PSI (reuse `pre/` sealed-prediction discipline). LOACO/LOCO
when n is small. Per-file pytest inline (full feng suite OOMs under the 10G slice). Regression-guard
the atlas against A3 GC-core recovery.

## Decisions locked
- Context = **both minigene + endogenous** (transfer is the headline); both now on disk.
- Validation sets: **primary** = 7-28 size comparison + 7-27 truncated (endogenous) + 7-7/7-9 super
  minigene; **secondary** = E3 positional walk; no-treatment = baseline only; dCas13d + fusion
  **excluded**.
- Off-target set (default) = **both**: public rMATS reanalysis now, Holgersen if accessible.
- Sequencing (default) = **shared backbone (P0–P3), then split** BIBE (methods) / BIOS (biology).

## What I still need from you to execute Pexp / the transfer headline
1. **Lane→sample→ASO map** for 7-27, 7-28, 7-7, 7-9 (which lane = which ASO/size/control; which lane
   is the ladder) — the CSVs carry none.
2. **Expected RT-PCR product sizes** (inclusion vs skipping amplicon lengths) for the endogenous and
   super-minigene assays, to lock band identity against the ladder.
3. Replicate structure: within each context, are the two dates (7-7 vs 7-9; 7-27 vs 7-28) independent
   biological replicates, technical replicates, or reruns?
