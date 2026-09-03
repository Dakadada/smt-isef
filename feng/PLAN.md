# feng/ — the new prospective arm

## Context

`docs/HANDOFF.md` holds the old plan: a compute-native cis-vs-trans steric-splicing
benchmark plus a generative design engine, honestly capped at "SMN2 feasibility +
timestamped prospective predictions." It is committed through task #11 at 182/182 tests,
with #16's adversarial gate still deferred.

You asked for a new arm in `feng/` — a Feng-Zhang-shaped move that turns one result
(nusinersen) into a platform — and proposed four directions: a pharmacogenomic responder
map, an exhaustive ASO design landscape, a transcriptome-wide off-target sieve, and a
transferable grammar of splice-silencer targeting.

Two Codex `gpt-5.6-sol` passes (an adversarial kill-check and a literature scout) found
that **all four, as stated, are already occupied.** But the scout found something better
in the same space, and it is more Feng-Zhang-shaped than the original four. This plan
keeps your four-arm platform structure and replaces each arm with a version that
survives a hostile judge.

### Decisions taken
- Timeline: **~12 months** (Sept 2026 → Aug 2027), with a showable Fall 2026 checkpoint.
- `feng/` is **clean-room** — no imports from `src/` or `pre/src/`. `pre/` may be read
  and its published logic re-derived, never copied.
- Build budget: **months 1–2 infrastructure, months 3–12 science.** Real ViennaRNA/
  IntaRNA stack, real sequences. No fixture-backend results.
- Data: public-first. A controlled-access direction may appear only behind an explicit
  early gate on whether your institutional route actually materializes.
- Claim posture: hunt for the genuinely open niche, hold the honesty line structurally.

---

## What the adversarial pass killed (do not re-litigate)

| Original arm | Verdict | Killing prior art |
|---|---|---|
| 1. Responder map | **CROWDED, data-contingent** | SMN1/SMN2 are >99.9% identical, copy-variable, gene-converting. gnomAD/1kGP short-read calls cannot be assigned to a copy. Only Paraphase-class HiFi data counts. Signoria 2024 already reported heterogeneous patient-specific nusinersen response across 45 fibroblast lines. |
| 2. Exhaustive landscape | **DEAD** | Hua 2007 exon walk + microwalks; Hua 2008 intron tiling (found the drug's silencer); Singh 2010 ISS-N1 microwalk incl. 8–15-mers; Recinos 2024 SpliceRUSH; **Celaj 2023 BigRNA already ranked all 437 local 18-mers, nusinersen top-3.** Tools: ASOG, PFRED, eSkip-Finder, Sfold, OligoWalk. |
| 3. Off-target sieve | **DEAD as a "first"** | Ottesen 2026 (PMC12805893 / GSE307431) ran transcriptome-wide RNA-seq on F18MOE — sequence- and chemistry-identical to nusinersen — and orthogonally validated 8 exon-skipping events incl. POLR2H exon 2. Holgersen 2021: naive edit-distance-3 = **14% sensitivity at 99% FDR**. |
| 4. Transferable grammar | **DEAD** | n=1 positive. Plus ESEfinder/RESCUE-ESE/PESX/HSF/SpliceAid2, RBPmap/ATtRACT/oRNAment/eCLIP, SpliceAI/Pangolin/MMSplice/SPLAM, Vex-seq/MaPSy/MFASS, FAS-ISS, SpliceRUSH, BigRNA, DANGO. |

---

## The new program

**The move:** stop mapping *where* to put an ASO. Start inferring **the operator an
effector imposes**. A splicing element has a latent regulatory field; an 18-nt duplex, a
25-nt dCas13d guide plus a bulky RNP, and a persistent U7 snRNP each apply a *different
steric shadow* to that same field. Nobody has derived that transfer law — and the
1,937-guide SpliceRUSH screen plus the historical SMN2 ASO walks make it derivable now.

That is the platform claim: **a cheap programmable RNP tiling screen becomes a
chemically actionable ASO design map.** One molecule → a general tool, with a hard
falsifiable core: *predict one effector from another.*

### Arm A — StericScope (HEADLINE)
*Replaces your arm 2. Novelty 5 / Feasibility 5 / Legibility 5.*

**Claim.** A latent regulatory field convolved with an effector-specific steric-shadow
kernel predicts splice modulation by the *other* modality.

**Why it is open.** SpliceRUSH asserts a shared steric mechanism and notes guide
orientation/footprint effects, but never infers footprint width, asymmetry, or excluded
volume, and never predicts ASO activity from dCas13d activity. DANGO's 21-nt-deletion
proxy and BigRNA's virtual ASOs do not compare *physical* effectors.

**Primary endpoint.** Harmonize each historical ASO study to within-study standardized
effect. Fit the latent field + dCas13d kernel using **no ASO labels**. Then, under
contiguous-region blocked validation, test Spearman ρ between predicted and measured ASO
effects with a one-sided circular-shift permutation that preserves local autocorrelation.

**Pre-declared success (all four):** ρ ≥ 0.50; permutation p < 0.01; ρ improvement ≥ 0.15
over guide-window averaging, accessibility-only, and DANGO-style deletion baselines;
nusinersen's interval in the top 5% of eligible windows.

**Kill criterion.** Fewer than 25 independent recoverable ASO windows; or leave-block-out
ρ < 0.30; or the fitted kernel is unstable across bootstrap resamples. → downgrade to a
descriptive cross-effector comparison.

**Terminal artifact.** `StericScope` — effector kernels with uncertainty, a reproducible
SMN2 benchmark, an interactive locus map, and an API turning a dCas13d tiling result into
ranked ASO windows.

**Data (all public):** SpliceRUSH 1,937-guide source data + PRJNA1052769; Hua 2007/2008
and Singh 2010 walk tables from published supplements; GSE307431/PRJNA1321660. BigRNA as
a baseline, never as training labels.

**Compute:** 200–2,000 CPU core-hours; <20 GPU-hours for optional neural kernels.

**GPU tier** *(design + hooks now; build deferred to a follow-up plan).* A neural
effector-kernel / DL latent-field fit augmenting the analytic kernel, wired to the shared
`feng/src/feng/gpu/` infra (`mit_preemptable`, checkpoint + requeue).

### Arm B — RBPShadow (mechanism)
*Replaces your arm 4 — the real "grammar" is which RBP gets displaced.*

**Claim.** A steric ASO is a programmable **local knockout** of the RBP occupying its
binding interval, so ASO-induced ΔPSI direction should match that RBP's knockdown
direction.

**Why it is open.** Scharner used *predicted* RBPmap motifs and saw enrichment, but never
connected ASO responses to *measured* eCLIP peaks with matched RBP-KD ΔPSI. ENCODE now
supplies 223 eCLIP datasets and 472 RBP-knockdown RNA-seq profiles — including HepG2,
where 48 of Holgersen's experiments were run.

**Primary endpoint.** Restrict to active HepG2 sites with adequate coverage and exactly
one high-confidence nuclear splicing-RBP eCLIP assignment. Two-sided exact binomial test
on sign concordance, with oligo-cluster bootstrap CIs.
**Success:** concordance ≥70%, lower 95% bound >50%, p < 0.01, replicated in a secondary
multi-RBP model and one non-Holgersen dataset.
**Kill:** <50 uniquely assigned events, concordance <60%, or the effect vanishes against
nearby non-overlapping eCLIP peaks as negative controls.

**Data:** ENCODE eCLIP ENCSR456FVU; HepG2 KD-RNA-seq ENCSR369TWP; K562 ENCSR795JHH;
processed splicing ENCSR413YAF / ENCSR870OLK; Holgersen *processed* tables; GSE307431.
**Compute:** 300–2,000 CPU core-hours on processed ENCODE events; <20 GPU-hours.

**GPU tier** *(design + hooks now; build deferred to a follow-up plan).* A deep RBP-binding
predictor for hnRNPA1 displacement over eCLIP/KD, wired to the shared `feng/src/feng/gpu/`
infra (`mit_preemptable`, checkpoint + requeue).

### Arm C — Bind→Act gate (safety)
*Replaces your arm 3 — the honest, useful version.*

**Claim.** Duplex formation and splice action are **separate gates**. Complementarity
predicts binding; a distinct molecular gate predicts whether binding changes splicing.

**Why it is open.** Scharner found several *inactive* sites had stronger predicted binding
than active off-targets; Holgersen found only ~35% of perfect exonic matches changed PSI
(16% for non-designed oligos). Neither integrated measured RBP occupancy, RBP-knockdown
direction, splice-site commitment, accessibility, and nuclear pre-mRNA availability into a
held-out-ASO actionability model evaluated by leave-one-ASO-out.

**Primary endpoint.** Define active/inactive by **coverage-aware power analysis before
feature generation**. Nested leave-one-ASO-out only. Compare against a preregistered
baseline of edit distance + mismatch position + duplex/accessibility energy; ASO-cluster
bootstrap on the AUCPR difference.
**Success:** absolute AUCPR gain ≥0.10 (or ≥2×), bootstrap p < 0.01, ≥70% sensitivity at
≤20% empirical FDR *within the candidate set*, directional replication on GSE307431 and
Scharner (neither used for tuning).
**Kill:** >30% of negative labels flip under reasonable coverage/power thresholds, or
external AUCPR gain <0.05.

**Compute:** 5,000–20,000 CPU core-hours for IntaRNA/RNAplfold candidate generation.

**GPU tier** *(design + hooks now; build deferred to a follow-up plan).* A DL actionability
model for the bind→act leave-one-ASO-out harness, wired to the shared `feng/src/feng/gpu/`
infra (`mit_preemptable`, checkpoint + requeue).

### Arm D — The locked nusinersen benchmark (the floor)
*The null-safe terminal validation. Guarantees a real result even if A–C all fail.*

A preregistered public benchmark: does RNA-accessibility-aware duplex scoring rank the 8
experimentally supported nusinersen off-target exons above matched detectable exons better
than sequence distance alone?

Per positive exon, build a stratum of ~100 background internal exons matched **using
ScrMOE controls only** on baseline PSI, junction depth, expression, exon length, GC.
Rank under baseline (min 18-nt Hamming distance to nusinersen's reverse complement) vs
test (best IntaRNA interaction energy over ±100-nt windows, accessibility-aware).
ΔR = mean percentile improvement over the 8.
Test: one-sided matched-stratum label permutation, 10⁶ permutations, frozen seed.
**Success:** ΔR ≥ 0.10 **and** p ≤ 0.05. Mandatory sensitivity analysis excluding POLR2H
(already known to you — the evaluation is *locked*, not blind).
**Kill:** any positive unmappable to the locked GENCODE universe; <50 eligible background
exons per positive; clean-room reanalysis fails to recover published direction for ≥6 of 8.

A secondary panel tests whether a fixed sequence score assigns F18MOE and F18OMe the same
risk despite different observed event sets — a **model transport failure**, not a learned
chemistry model.

**GPU tier** *(design + hooks now; build deferred to a follow-up plan).* A GPU alternative
score (MD / DL) wired to the shared `feng/src/feng/gpu/` infra, added **only as a secondary,
preregistered score defined BEFORE the FREEZE tag**. It must **not** alter the locked primary
IntaRNA score or the sealed Phase-9 permutation test; it is a preregistered add-on ranked
alongside, never a replacement for, the locked benchmark.

### Arm E — StericRuler (length titration)
*New fifth arm. Symmetric length-titration of the true nusinersen 18-mer, run in tandem
with a wet-lab length series. Pre-FREEZE; freeze-safe (reads no sealed ΔPSI).*

**Claim.** A symmetric length ladder over the true nusinersen 18-mer (even lengths 2→30,
GRCh38 v50 anchor `chr5:70076584-70076601`, ±1 nt/end from the center) predicts a per-length
**therapeutic window** = predicted on-target potency ÷ off-target burden. Pre-registered
directional prediction — **"shorter is better"** (window favors a sub-18 length) — locked
before the wet-lab readout.

**Why it is open.** No study titrates length symmetrically around the real drug in the
correct genomic frame and scores potency *against* specificity per length. CPU-thermo
potency alone rises monotonically with length, so whether a shorter mer trades potency for
specificity *favorably* is genuinely unknown until the off-target and GPU axes enter the
window.

**Primary endpoint.** The per-mer therapeutic-window curve vs length and a directional
verdict (`supports` / `refutes` shorter-better), preregistration-locked (currently **DRAFT,
unsealed**) before the wet-lab length series is read.

**Success:** window non-monotonic / favoring a sub-18 length **and** rank-concordant with
the wet-lab series.
**Kill:** window monotonically increasing with length (refutes "shorter better"), or wet-lab
discordant.

**Terminal artifact.** `StericRuler` — the locked per-length ladder, the window curve with
its axis definitions, a length-axis figure, and the prereg hash, comparable against the blot.

**Data (all public):** nusinersen anchor + GRCh38 v50 genome. CPU scorers (IntaRNA on-target,
Vienna RNAduplex/RNAcofold/RNAup/RNAplfold/RNAfold) + off-target seed-extend transcriptome
scan. GPU axes: GROMACS MD binding free energy, Pangolin/SpliceAI deep-splice (on-target
potency + off-target splicing), hnRNPA1 displacement. Reads **no sealed ΔPSI** — freeze-safe;
its score/window definitions are preregistered at the eventual Phase-8 FREEZE.

**Compute:** CPU core done (`feng/src/feng/titration/`); GPU tier on `mit_preemptable` via the
shared `feng/src/feng/gpu/` infra (checkpoint + requeue).

### Appendix — variant × ASO decomposition (kill-gated)
*Your arm 1, reduced to what survives.*

The factorial decomposition — `cis effect = f(ALT,−ASO) − f(REF,−ASO)` vs
`interaction = [f(ALT,+ASO) − f(ALT,−ASO)] − [f(REF,+ASO) − f(REF,−ASO)]` — is
**methodologically unpublished** for natural copy-specific SMN2 variants. But with no
variant ± nusinersen measurements it separates two terms *inside your own model*, not a
biological interaction.

**Two-day recon, then a hard gate:** continue only if ≥2 independent public copy-specific
SMN2 haplotypes carry a non-reference allele inside the 18-nt footprint. Otherwise publish
a target-site conservation audit and stop.

---

## Phase plan

| Months | Phase | Deliverable | Gate |
|---|---|---|---|
| 1–2 | **Infrastructure** | `feng/` package, tests, schemas, own `feng/env/environment.yml` (py3.11, ViennaRNA 2.7.2, IntaRNA 3.4.1, biopython, numpy/scipy/pandas, pyarrow, matplotlib, pytest, STAR, samtools, bcftools, bedtools, rMATS-turbo 4.3.0). Data manifest with URLs + MD5/SHA-256. GENCODE v50 + STAR index. Re-derived seed-and-extend scanner with exhaustive toy-reference tests. | All imports + STAR/IntaRNA/rMATS smoke tests pass. Scanner has zero seed-related false negatives. Otherwise stop science and fix infra. |
| 2 | **Variant recon** (2 days) | SMN2 target-haplotype audit from public Paraphase 107-sample data + HPRC r2. | ≥2 independent copy-specific non-reference target alleles, or the appendix dies. |
| 3–4 | **Ground truth** | Clean-room full-transcriptome alignment of the 20 GSE307431 libraries; control-only detectable-exon universe; digitized Hua/Singh walk tables with provenance; SpliceRUSH guide-level effects separated from reporter enrichment. **Treatment contrasts written to a sealed directory, not inspected.** | ≥80% unique mapping; no chemistry/control swap; SMN2 target junction moves in the expected direction. Catastrophic batch confounding kills Arm D. |
| 5–7 | **Arm A build** | Latent field + effector kernel fits; blocked-region CV harness; baselines (guide-window averaging, accessibility-only, DANGO-style deletion). | ≥25 independent ASO windows recoverable, kernel bootstrap-stable. |
| 6–8 | **Arm C build** (parallel) | Transcriptome candidate generation; power-aware active/inactive labels **defined before features**; leave-one-ASO-out harness. | >30% negative-label instability → Arm C dies. |
| 5–8 | **Arm E build** (parallel) | StericRuler symmetric length ladder + CPU scorers + therapeutic-window endpoint + shared GPU tier (`feng/src/feng/gpu/`). Prereg **DRAFT** (unsealed). | CPU potency alone rises with length → the window verdict needs the off-target + GPU axes; score/window defs sealed at FREEZE. |
| 8 | **FREEZE** | Preregistration: event universe, score definitions, matching tolerances, endpoints, tests, thresholds, exclusions, and the null-result text — all written out. Git tag + timestamped OSF/Zenodo deposit of code/config/data hashes. | No terminal contrast or rank result inspected before freeze. |
| 9 | **Terminal execution** | One shot. Unlock sealed contrasts. Run Arm A validation, Arm D permutation test, Arm C external replication. | **No refitting.** A bug found after unblinding invalidates the confirmatory label; the correction is explicitly exploratory. |
| 10 | **Arm B** | RBPShadow concordance test + atlas. | <50 uniquely assigned events → mechanistic claim dropped, panel stays descriptive. |
| 11 | **Robustness** | Specification curve: mismatch ceilings 4–6, junction-coverage thresholds, POLR2H exclusion, isoform/multi-mapping audit. Variant appendix if it passed. | Contradictory robustness is *reported*, never optimized away. |
| 12 | **Package + audit** | Poster, methods diagram, one-command judge demo, reproducibility capsule, adversarial limitations table. Independent rerun from manifests with hash verification. | Every headline number traces to a frozen table and script. Any unexplained manual step removes the "fully reproducible" claim. |

**Fall 2026 checkpoint (end of month 4):** the ground-truth phase is a showable result on
its own — a clean-room reproduction of a published drug's off-target events from raw public
reads, with a locked event universe.

**Compute envelope:** STAR index ~32 CPU-h (40–60 GB RAM); alignment array 240–640 CPU-h;
rMATS 100–300 CPU-h; scanner shards 50–150 CPU-h; IntaRNA 100–400 CPU-h (Arm D) and
5,000–20,000 CPU-h (Arm C, sharded); GPU only for Arm A's optional neural kernels
(<20 h). Scratch 0.3–1 TB. All heavy work via sbatch, checkpointed, preemptable-safe.

---

## `feng/` layout

```
feng/
  env/environment.yml          # own stack, not pre/env's
  README.md                    # claim boundary, clean-room declaration
  PREREGISTRATION.md           # frozen at month 8
  src/feng/
    data/        manifest.py, fetch.py, hashes.py
    align/       star.py, rmats.py, universe.py       # detectable-exon universe
    scan/        seed_extend.py, intarna.py           # re-derived, not copied
    field/       latent.py, kernels.py                # Arm A
    rbp/         encode.py, concordance.py            # Arm B
    gate/        power.py, features.py, model.py      # Arm C
    bench/       strata.py, permute.py, terminal.py   # Arm D
    titration/   ladder.py, evaluate.py, window.py, report.py   # Arm E
    gpu/         env.py, checkpoint.py                 # shared GPU tier (all arms)
    honesty/     claims.py, linter.py, states.py
  tests/
  terminal/                    # SEALED until month 8 freeze
  artifacts/
```

## Verification

- `pytest` in the `feng` env; scanner verified against brute-force toy references with
  zero seed-related false negatives (POLR2H must **not** be special-cased).
- Smoke: STAR alignment of `tests/fixtures/tiny_R1.fastq`; IntaRNA and rMATS one-shot runs.
- Clean-room check: a manifest test that **fails** if anything under `feng/` imports or
  path-references `src/`, `pre/src/`, `design_engine`, `offtarget`, `physics_adapter`, or
  `design_panel`.
- Freeze check: a test that fails if feature/model code can read `feng/terminal/` before
  the freeze tag exists.
- End-to-end: one command regenerates every headline number from the data manifest.
- **Standing gate (CLAUDE.md):** every arm goes build → Codex sol adversarial review → fix,
  before advancing.

## Honesty architecture

Verbatim on the artifact:

> This study benchmarks candidate hybridization-associated splice off-targets in one
> fibroblast RNA-seq experiment. It does not measure clinical safety, CNS exposure,
> patient response, or therapeutic efficacy.

> F18MOE has the sequence and nominal 2′-MOE/phosphorothioate chemistry of nusinersen, but
> 6 µM nucleofection in fibroblasts is not clinical intrathecal dosing.

> Matched background exons are detectable comparison sites, not proven biological
> negatives; precision and false-discovery rate are therefore not estimated.

> IntaRNA and ViennaRNA use RNA thermodynamic parameters. Their scores are RNA-surrogate
> interaction/accessibility values, not calibrated thermodynamics for 2′-MOE, 2′-OMe, or
> PMO oligonucleotides.

> An inferred kernel width is a statistical property of these data, not a measured
> molecular dimension.

Enforced in code: allowed result states limited to `candidate`, `paper_supported`,
`clean_reanalysis_reproduced`. No fields named `clinical_response`, `patient_responder`,
`safe`, `efficacious`, `superior_to_nusinersen`. Report generation refuses to print
"precision" or "FDR" when the negative-label manifest is incomplete. Empty output renders
as "no candidates detected under this specification," never "zero risk." Chemistry is
mandatory metadata; pooling chemistries is prohibited. A claim linter blocks "first",
"genome-wide safety map", "clinical response", "beats nusinersen", "generalizable" from
final figures unless explicitly overridden and logged.

## If the headline fails

Arm A failing is not project failure — it is the finding that **steric effects do not
transfer across effector geometries**, which contradicts SpliceRUSH's stated shared-mechanism
premise and matters to anyone using dCas13d screens for ASO design. Arm D's locked benchmark
delivers a real result either way. A negative result *with* preregistration, raw-read
reproduction, and terminal hashes is a strong entry; the same result without them is
ordinary and weak.

## Traps (explicit do-nots)

- Do not train a chemistry model on F18MOE/F18OMe/F20PMO — one sequence, one site, and
  F20PMO is a 20-mer with a neutral backbone and different delivery.
- Do not compare recall/FDR numerically against Holgersen; the label universes differ.
- Do not use gnomAD or 1kGP target-region calls as copy-specific SMN2 evidence.
- Do not count HPRC assemblies as SMA patients.
- Do not wait on dbGaP/EGA/author-mediated data (EGAS50000001222 is controlled;
  Holgersen raw reads are proprietary; the 31-patient ONT gene-conversion set needs an IRB
  and a data-sharing agreement). "An application was submitted" is not data.
- Do not download whole-genome gnomAD; one indexed chr5 interval suffices.
- Do not stack SpliceAI/Pangolin/SPLAM/IntaRNA as an ensemble because they exist — each
  needs a defined causal role and a locked ablation.
- Do not treat Stage-A targeted junction counts as transcriptome-wide off-target discovery.
- Do not copy-and-rename code from `src/` or `pre/src/`. Published logic may be
  re-derived; code may not be moved.
- Do not repair a failed terminal model post hoc and present it as confirmatory.
