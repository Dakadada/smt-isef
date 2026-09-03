# feng/ — execution taskboard

Created 2026-09-01. Detailed constraints, evidence, and kill criteria live in the plan
(`~/.claude/plans/yeah-so-just-talk-reflective-moon.md`) and in `feng/HANDOFF.md`.

**State: greenfield — nothing started.** Timeline ~12 months (Sept 2026 → Aug 2027).
Every arm goes build → Codex sol adversarial review → fix/log → advance (standing gate).

Status keys: `[ ]` todo · `[~]` in progress · `[x]` done · `[!]` blocked/gated ·
`[-]` killed by its own criterion.

---

## Phase 0 — Pre-flight (advisory, cheap, do first)

- [ ] **P0.1 Arm A overlap probe (~1 day).** From SpliceRUSH source (PRJNA1052769) +
  Hua 2007/2008 & Singh 2010 supplement walk tables: confirm the dCas13d field and the
  ASO walks co-locate on the same SMN2 coordinate frame with **≥25 recoverable ASO
  windows** under guide-covered field. **Gate:** thin overlap → Arm A pre-degraded to a
  descriptive cross-effector comparison before infra is built.

## Phase 1 — Infrastructure (months 1–2)

- [x] **1.1** `feng/` package + `env/environment.yml` (py3.11, ViennaRNA 2.7.2,
  IntaRNA 3.4.1, biopython, numpy/scipy/pandas, pyarrow, matplotlib, pytest, STAR,
  samtools, bcftools, bedtools, rMATS-turbo 4.3.0). *(package skeleton + env FILE done;
  env not yet solved/built. + `pyproject.toml`, CLI entry point.)*
- [x] **1.2** Copy frozen plan into `feng/PLAN.md`; write `README.md` (claim boundary +
  clean-room declaration).
- [x] **1.3** Data manifest with URLs + MD5/SHA-256 (`data/manifest.py`, `fetch.py`,
  `hashes.py`). *(schema + hash-verify path done; no real URLs/downloads yet.)*
- [x] **1.4** GENCODE v50 + STAR index. *Scripts adversarial-reviewed + fixed (4 high +
  1 med). GENCODE v50 fetched + MD5-verified (genome 3.0G, GTF 4.4G in
  `$SCRATCH/feng/gencode/v50`). **STAR index BUILT + verified** (job 21773261, COMPLETED
  32min): 30G at `$SCRATCH/feng/star/gencode-v50-sjdb149`, all 15 required files +
  completion marker (STAR 2.7.11b, sjdb149). Requeue-safe sbatch, atomic publish confirmed.*
- [x] **feng conda env** BUILT at `/home/dashii/orcd/pool/feng-env` via `module load
  miniforge/25.11.0-0` + `mamba env create -p` (pkgs cache on scratch). Tools verified:
  STAR 2.7.11b, IntaRNA 3.4.1, ViennaRNA 2.7.2, samtools 1.23, bcftools 1.23.1,
  bedtools 2.31.1, rmats.py 4.3.0. **91 passed** under feng-env python.
- [x] **1.5** Re-derived seed-and-extend scanner (`scan/seed_extend.py`) with exhaustive
  brute-force toy-reference tests — zero seed-related false negatives, no POLR2H special
  case. *Done + adversarial-reviewed clean (2 medium findings fixed: unbounded fallback
  bounded to remaining ref length; randomized tests replaced with 26,520-config exhaustive
  enumeration + independent Levenshtein/Hamming oracle). Pigeonhole seed = budget+1 blocks.*
- [x] **1.6** Clean-room manifest test (fails on any import/path-ref to `src/`,
  `pre/src/`, `design_engine`, `offtarget`, `physics_adapter`, `design_panel`).
- [x] **1.7** Freeze-seal test (feature/model code fails if it can read
  `feng/terminal/` before the freeze tag exists). Freeze tag recognized as
  `feng-freeze` / `feng-freeze-*`.
- [x] **Honesty layer** (standing): `honesty/states.py` (allowed states only),
  `honesty/claims.py` + `linter.py` (forbidden strings/fields blocked, override+log
  path). Tested.
- [x] **GATE-1 smoke tests** (`tests/test_smoke_tools.py` + tiny fixtures): STAR tiny
  index+align, IntaRNA one-shot duplex, rMATS one-shot — all PASS under feng-env.
  (IntaRNA fixture seqs fixed to a verified complementary pair, E=-34.78.)
- **GATE 1: PASSED.** imports + guard + honesty + scanner + STAR/IntaRNA/rMATS smoke all
  green (**94 passed** under feng-env); GENCODE v50 STAR index built + verified.
  **PHASE 1 COMPLETE.**

## Phase 2 — Variant recon (2 days, month 2) — kill-gated

- [ ] **2.1** SMN2 target-haplotype audit from public Paraphase 107-sample data + HPRC
  r2. **Gate:** ≥2 independent copy-specific non-reference alleles inside the 18-nt
  footprint → appendix lives; else publish a target-site conservation audit and stop
  the appendix.

## Phase 3–4 — Ground truth (months 3–4) — Fall 2026 checkpoint

- [ ] **3.1** Clean-room full-transcriptome alignment of the 20 GSE307431 libraries
  (`align/star.py`, `universe.py`).
- [ ] **3.2** Control-only detectable-exon universe.
- [ ] **3.3** Digitized Hua/Singh walk tables with provenance.
- [ ] **3.4** SpliceRUSH guide-level effects separated from reporter enrichment.
- [ ] **3.5** Treatment contrasts written to a **sealed** dir, not inspected.
- **GATE 3:** ≥80% unique mapping; no chemistry/control swap; SMN2 target junction moves
  in expected direction. Catastrophic batch confounding **kills Arm D**.
- **★ CHECKPOINT (end of month 4):** clean-room reproduction of a published drug's
  off-target events from raw public reads, with a locked event universe — showable
  standalone.

## Phase 5–7 — Arm A build: StericScope (HEADLINE)

- [ ] **A.1** Latent field + effector kernel fit (`field/latent.py`, `kernels.py`),
  dCas13d only, no ASO labels.
- [ ] **A.2** Blocked-region (contiguous) CV harness.
- [ ] **A.3** Baselines: guide-window averaging, accessibility-only, DANGO-style
  deletion proxy.
- **GATE A:** ≥25 independent recoverable ASO windows; kernel bootstrap-stable;
  leave-block-out ρ ≥ 0.30. Else → descriptive cross-effector comparison.

## Phase 6–8 — Arm C build: Bind→Act gate (parallel with A)

- [ ] **C.1** Transcriptome candidate generation (`gate/`, IntaRNA/RNAplfold; sharded,
  5,000–20,000 CPU core-hours — confirm it earns cost).
- [ ] **C.2** Power-aware active/inactive labels **defined before features**
  (`gate/power.py`).
- [ ] **C.3** Nested leave-one-ASO-out harness (`gate/model.py`) vs preregistered
  edit-distance + mismatch-position + duplex/accessibility baseline.
- **GATE C:** >30% negative-label instability under reasonable power thresholds → Arm C
  dies.

## Phase 5–8 — Arm E build: StericRuler (length titration, parallel; pre-FREEZE)

- [x] **E.1** Symmetric ladder (`titration/ladder.py`) — even 2→30 around the true
  nusinersen 18-mer, GRCh38 v50 anchor `chr5:70076584-70076601`, fail-closed L=18 ==
  nusinersen. *Done + committed `b08fae7` (12 tests). Opus-fork gate: ADVANCE.*
- [x] **E.2** Per-mer CPU scorers (`titration/evaluate.py`) — IntaRNA on-target (3 honest
  states), Vienna RNAduplex/RNAcofold/RNAup/RNAplfold(genomic-context)/RNAfold, off-target
  seed-extend. *Done + committed `b08fae7` (19 tests). MED-1 accessibility fix applied.*
- [x] **E.3** Therapeutic-window endpoint (`titration/window.py`) — potency ÷ off-target
  burden + preregistration-lock mechanism (**DRAFT, no seal cut**). *Done + committed
  `b08fae7` (32 tests). Opus-fork gate: ADVANCE; MED-1 verdict-plateau fix applied.*
- [x] **E.4** Artifact builder (`titration/report.py`) — README + per-mer table + length-
  axis figure + provenance, load-bearing honesty lint. *Done + committed `b08fae7`
  (11 tests).*
- [x] **E.5** Shared GPU infra (`gpu/env.py`, `gpu/checkpoint.py`,
  `scripts/templates/gpu_preemptable.sbatch`, documented-not-run `scripts/setup_gpu_env.sh`).
  *Done (24 tests); NOT yet committed. Reusable by all arms.*
- [ ] **E.6** Off-target transcriptome scan → populate off-target burden (the decisive
  CPU axis; runnable today; makes `therapeutic_window` real + gives the directional verdict).
- [ ] **E.7** GPU env install (`feng-gpu-env`: cu124 torch + Pangolin/SpliceAI + OpenMM;
  quota-checked, manually run) + wire G/H scorers (GROMACS MD, deep-splice, hnRNPA1).
- [!] **E.8** Prereg **seal** (`locked=true`) — gated: needs user sign-off + all axes
  (off-target + GPU) present. Currently DRAFT.
- **GATE E:** window non-monotonic / favoring sub-18 **and** wet-lab-concordant → SUCCESS;
  monotonic-increasing or wet-lab-discordant → KILL. **HONEST DRAFT SIGNAL:** CPU-thermo
  potency alone `refutes shorter-better` (rises with length) — verdict pends off-target + GPU.
  Codex sol durable gate deferred to Sep 7; Opus-fork fallback gates used on E.1–E.4.

## Phase 8 — FREEZE

- [ ] **8.1** `PREREGISTRATION.md`: event universe, score defs, matching tolerances,
  endpoints, tests, thresholds, exclusions, and the null-result text — all written out.
- [ ] **8.2** Git tag + timestamped OSF/Zenodo deposit of code/config/data hashes.
- **GATE 8:** no terminal contrast or rank result inspected before freeze.

## Phase 9 — Terminal execution (one shot)

- [ ] **9.1** Unlock sealed contrasts. Run Arm A validation, Arm D permutation test,
  Arm C external replication. **No refitting** — a post-unblind bug invalidates the
  confirmatory label; corrections are explicitly exploratory.
  - Arm A success (all four): ρ ≥ 0.50; perm p < 0.01; ρ gain ≥ 0.15 over all three
    baselines; nusinersen interval in top 5% of eligible windows.
  - Arm D success: ΔR ≥ 0.10 **and** p ≤ 0.05 (10⁶ matched-stratum permutations, frozen
    seed); mandatory POLR2H-excluded sensitivity.
  - Arm C success: AUCPR gain ≥ 0.10 (or ≥2×), bootstrap p < 0.01, ≥70% sensitivity at
    ≤20% empirical FDR within the candidate set, directional replication on GSE307431 +
    Scharner.

## Phase 10 — Arm B: RBPShadow (mechanism)

- [ ] **B.1** ENCODE eCLIP + KD-RNA-seq ingest (`rbp/encode.py`); restrict to active
  HepG2 sites with one high-confidence nuclear splicing-RBP assignment.
- [ ] **B.2** Sign-concordance test (`rbp/concordance.py`, exact binomial, oligo-cluster
  bootstrap CIs) + atlas.
- **GATE B:** <50 uniquely assigned events → mechanistic claim dropped, panel stays
  descriptive. Success: concordance ≥70%, lower 95% bound >50%, p < 0.01, replicated.

## Phase 11 — Robustness

- [ ] **11.1** Specification curve: mismatch ceilings 4–6, junction-coverage thresholds,
  POLR2H exclusion, isoform/multi-mapping audit. Contradictory robustness is *reported*,
  never optimized away.
- [ ] **11.2** Variant appendix decomposition — only if Phase 2 gate passed.

## Phase 12 — Package + audit

- [ ] **12.1** Poster, methods diagram, one-command judge demo, reproducibility capsule,
  adversarial limitations table.
- [ ] **12.2** Independent rerun from manifests with hash verification.
- **GATE 12:** every headline number traces to a frozen table + script; any unexplained
  manual step removes the "fully reproducible" claim.

---

## Standing / cross-cutting

- [ ] Codex `gpt-5.6-sol` adversarial review after **each** arm build; fix or log
  findings before advancing.
- [ ] Honesty linter + allowed-result-states enforcement wired before any report is
  generated (`honesty/`).
- [ ] All heavy compute via `sbatch`, checkpointed, preemptable-safe.
- [x] Shared GPU tier (`gpu/`) + `gpu_preemptable.sbatch` template + documented installer.
  *Done (E.5); not yet committed.*
- [ ] **GPU-slot A** (design-only now): neural effector-kernel / DL latent-field fit
  augmenting Arm A's kernel, wired to `gpu/`. Build deferred to a follow-up plan.
- [ ] **GPU-slot B** (design-only now): deep RBP-binding predictor for hnRNPA1 displacement
  over eCLIP/KD, wired to `gpu/`. Build deferred.
- [ ] **GPU-slot C** (design-only now): DL actionability model for the bind→act
  leave-one-ASO-out harness, wired to `gpu/`. Build deferred.
- [!] **GPU-slot D** (design-only now): GPU alternative score (MD/DL) — **secondary,
  preregistered, defined BEFORE the FREEZE tag; must NOT alter the locked primary IntaRNA
  score or the sealed Phase-9 test.** Build deferred.
