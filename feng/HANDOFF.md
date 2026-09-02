# HANDOFF — feng/ (the new prospective arm)

**Created: 2026-09-01. State: greenfield. Nothing built yet. Plan is frozen.**

This file is the entry point for anyone (Claude or Codex) picking up the `feng/`
arm. The old `src/`/`pre/` program is a separate, committed effort — see
`docs/HANDOFF.md` for it. `feng/` does not touch that code.

---

## ★ LIVE CHECKPOINT — 2026-09-02 ~09:25 EDT (three fronts in flight)

**STOP-ALL deadline this session: 10:05 AM EDT** (user-set). Checkpoint + halt then.

Key unblock: **the login node HAS internet.** Every prior "environment-limited"
NO-GO / staging failure was the *Codex sandbox* network block, NOT the login node.
Fetch directly Claude-side (curl) or via sbatch on a networked partition.

Three parallel fronts launched, all via `codex:codex-rescue` -> background Codex Tasks:

1. **Arm A descriptive atlas** — Codex `task-mtk4nwmu-x37mb4`. Implement
   ARM_A_DESCRIPTIVE_SPEC.md; fix OPEN finding #4 (+9nt S4 frame offset) first;
   descriptive only, no inferential pooling. RUNNING.
2. **Phase 3 rMATS** — Codex `task-mtk4opx5-stoslt`. Resolve GSM9224134-9224153 ->
   condition mapping (each ASO chemistry vs control, NO per-chemistry model), write
   `feng/scripts/phase3_rmats.sbatch` on the 20/20 finished STAR BAMs. RUNNING.
3. **Phase 2 HapSMA audit** — Codex `task-mtk4pg5v-4dc6p0`.
   Audit MOESM5 Table S3 for donor-independent SMN2 target-site ALTs in T2T interval
   `chr5:71408792-71408809`. This is the extraction the prior NO-GO couldn't do
   (sandbox had no net); data now STAGED. RUNNING.

**Phase 2 data staged Claude-side** at `$SCRATCH/feng/phase2/hapsma/`:
MOESM3/4/5/6/7.xlsx (HapSMA Genome Medicine Additional files 3-7; MOESM5=Table S3 the
decisive per-haplotype variant table; sha256 recorded in session). Springer URL base:
`static-content.springer.com/esm/art%3A10.1186%2Fs13073-025-01448-2/MediaObjects/`.
**Paraphase `SMN_phased_data` repo probed: all 538 paths are `HG*` = GIAB/HPRC
controls, NOT SMA patients** (confirms README cohort-label discrepancy). So Paraphase
bamlets are control-conservation evidence only; HapSMA MOESM5 is the sole patient
source. Not staging the control bamlets (200+ LFS BAMs; HPC courtesy).

Resume after /clear: check `/codex:status <task-id>` for each; verify outputs on disk
(atlas figures/code under feng/; `feng/scripts/phase3_rmats.sbatch` +
`feng/artifacts/phase3/sample_conditions.tsv`; `feng/artifacts/phase2_variant_recon/
EXTRACTION_HAPSMA.md`). Standing gate: adversarial-review each before advancing.
Nothing committed this session.

---

## What `feng/` is

A Feng-Zhang-shaped second arm: turn the single nusinersen result into a platform.
Full plan lives at:

- `~/.claude/plans/yeah-so-just-talk-reflective-moon.md` (canonical, frozen)

Copy it into `feng/PLAN.md` before starting real work so the repo is self-contained.

**The move.** Stop mapping *where* to place an ASO. Infer **the operator an effector
imposes**: a splicing element has a latent regulatory field; an 18-nt ASO duplex, a
25-nt dCas13d guide + bulky RNP, and a persistent U7 snRNP each cast a *different
steric shadow* on that same field. Derive the transfer law. Platform claim: a cheap
programmable RNP tiling screen becomes a chemically actionable ASO design map, with a
hard falsifiable core — **predict one effector from another.**

## The four arms (+ appendix)

**Framing (2026-09-01, user):** there is NO single headline. All four arms A–D are
**concurrently headline** — four independent falsifiable cores, each standing on its own,
not a hero-arm-plus-fallback. "FLOOR" on D means guaranteed-to-land, not lesser. This
reframe is load-bearing for the Arm-A RED verdict below: Arm A failing its transfer-law
form does not sink the project, and it also means Arm A is worth *contesting* rather than
quietly degrading (see "Arm A extraction" — user chose to contest the independence bar).

| Arm | Name | Role | Falsifiable core |
|---|---|---|---|
| A | **StericScope** | headline (contested — see RED gate below) | Latent field × effector steric-shadow kernel predicts the *other* modality's splice modulation. Fit on dCas13d, no ASO labels; test on ASO. |
| B | RBPShadow | headline (mechanism) | Steric ASO = local RBP knockout; ASO ΔPSI direction matches that RBP's KD direction (ENCODE eCLIP + KD-RNA-seq). |
| C | Bind→Act gate | headline (safety) | Binding and action are separate gates; held-out actionability model beats an edit-distance/energy baseline (leave-one-ASO-out). |
| D | Locked nusinersen benchmark | headline (guaranteed to land) | Accessibility-aware duplex scoring ranks the 8 supported nusinersen off-target exons above matched exons better than sequence distance. Stands alone even if A–C fail. |
| Appx | variant × ASO decomposition | kill-gated | Your old arm 1, reduced to what survives a 2-day recon gate. |

The original four proposed arms are all **dead / crowded** (responder map, exhaustive
landscape, off-target sieve, transferable grammar). Do not re-litigate — killing prior
art is tabulated in the plan.

## Current state (checkpoint 2026-09-01, Phase 1 in progress)

Working style: substantive work delegated to Codex `gpt-5.6-sol` via the
`codex:codex-rescue` subagent (thread `01a05eb9-90cc-7092-8036-719af336f70f`); Claude
orchestrates, verifies on disk, runs the standing adversarial gate. **All work
uncommitted** — user chose "hold until milestone" (commit at end of Phase 1).

DONE + adversarial-reviewed clean:
- **Package skeleton** per PLAN layout: `feng/src/feng/{data,align,scan,field,rbp,gate,
  bench,honesty}/`, `tests/`, `terminal/` (sealed placeholder), `artifacts/`,
  `pyproject.toml`, CLI entry point (`python -m feng --help` prints claim-boundary
  usage, does no pipeline work).
- **Honesty layer** (load-bearing): `honesty/states.py` (allowed states only:
  candidate / paper_supported / clean_reanalysis_reproduced), `honesty/claims.py` +
  `linter.py` (forbidden strings/fields blocked, override+log path). Tested.
- **Guard tests**: `test_cleanroom.py` (fails on any ref to src/, pre/src/,
  design_engine, offtarget, physics_adapter, design_panel), `test_freeze_seal.py`
  (freeze tag `feng-freeze` / `feng-freeze-*`), `test_honesty.py`, `test_data.py`,
  `test_cli.py`.
- **Data manifest scaffolding**: `data/manifest.py` (schema), `hashes.py` (MD5/SHA-256),
  `fetch.py` (hash-verify path real; no real URLs/downloads yet).
- **Seed-and-extend scanner** (`scan/seed_extend.py`) — the correctness keystone.
  Re-derived, pigeonhole seed = budget+1 blocks, zero seed-related false negatives.
  Reviewed clean after fixing 2 medium findings: unbounded exhaustive fallback bounded
  to remaining ref length; randomized tests replaced with **26,520-config exhaustive**
  enumeration + a genuinely independent Levenshtein/Hamming oracle. No POLR2H/locus
  special-casing (AST-guarded).
- **Test status: 84 passed** (via `/home/dashii/orcd/pool/smt-isef-env/bin/python -m
  pytest feng/tests/` — the OLD-arm pool env used only as a pytest runner; feng's own
  env is being built).

DONE since:
- **feng conda env BUILT** at `/home/dashii/orcd/pool/feng-env` (py3.11). Built directly
  (Codex sandbox could not write conda's pkgs dir): `module load miniforge/25.11.0-0` →
  `mamba env create -p ... -f feng/env/environment.yml`, pkgs cache on scratch. Tool
  versions verified: STAR 2.7.11b, IntaRNA 3.4.1, ViennaRNA 2.7.2, samtools 1.23,
  bcftools 1.23.1, bedtools 2.31.1, rmats 4.3.0, pandas 3.0.5 / numpy 1.26.4 / scipy /
  biopython 1.88 / pyarrow / matplotlib.
- **`.gitignore` fixed**: bare `data` rule was hiding the whole `feng/src/feng/data`
  source package (clean-checkout would fail imports) — added `!feng/src/feng/data/`
  negations. (Adversarial-review finding.)
- **1.4 scripts written + under review**: `scripts/fetch_gencode.sh`,
  `scripts/star_index.sbatch` (account=mit_general, partition=mit_preemptable, requeue-
  safe), manifest wired. Adversarial review found 4 high + 1 med; the gitignore one is
  fixed, the other 4 (TOFU digests, decompress self-validation, weak STAR completion
  predicate, personal-path defaults) are being fixed by Codex now. Authoritative GENCODE
  v50 MD5s pinned: genome.fa.gz=da1a11258be075cfa7af718162c894e7,
  annotation.gtf.gz=289b91e5e95e8b0450d223246f10a12e.

DONE since:
- **GATE-1 smoke tests wired + green under feng-env** (`tests/test_smoke_tools.py` +
  tiny fixtures): STAR tiny index+align, IntaRNA one-shot duplex, rMATS one-shot.
  IntaRNA fixture fixed to a verified complementary pair. **94 passed** under feng-env
  (91 + 3 smokes). GATE 1 tests CLOSED.
- **GENCODE v50 fetched + verified** into `/home/dashii/orcd/scratch/feng/gencode/v50/`:
  genome.fa (3.0G, MD5 da1a1125…, SHA-256 b760d18d…) + annotation.gtf (4.4G,
  MD5 289b91e5…, SHA-256 89bbad69…), provenance sidecars written.
- **STAR index sbatch SUBMITTED — job 21773261** (mit_preemptable, requeue-safe,
  16 cpu / 60G / 2h). Output → `$SCRATCH/feng/star/gencode-v50-sjdb149`, logs →
  feng/logs/. Submit env: `SCRATCH=/home/dashii/orcd/scratch
  FENG_ENV_PREFIX=/home/dashii/orcd/pool/feng-env`.

DONE since:
- **Phase 0 Arm A overlap probe — PARTIAL** (`artifacts/phase0_arm_a_overlap/README.md`).
  Data sources verified + pinned (SpliceRUSH PMC11076525 / SRA PRJNA1052769; Hua 2007
  PLoS Biol; Hua 2008 PMC2427210; Singh 2010 RNA j.). Key findings: SpliceRUSH per-guide
  effects are NOT tabulated (must re-derive from raw reads / Suppl Data 4); coordinate
  frames differ 3 ways (construct-relative + hg19 vs walk-specific SMN2 numbering);
  the ≥25 ASO windows come from the historical walks (SpliceRUSH itself has only ~3 ASOs).
  **Verdict AMBER** — promising, but the ≥25 count is UNMEASURED and gated on extraction.
  No number claimed. Execution TODOs listed in the memo.

BLOCKER:
- **Codex usage limit hit (~2026-09-01 evening; resets ~9:25 PM).** Substantive
  delegation paused; Claude doing the work directly until reset. Resumable Codex thread:
  `01a05eb9-90cc-7092-8036-719af336f70f`.

### ★ PHASE 1 COMPLETE (2026-09-01)

STAR index BUILT + verified (job 21773261, COMPLETED 32min): 30G at
`$SCRATCH/feng/star/gencode-v50-sjdb149`, all 15 required files + completion marker
(STAR 2.7.11b, sjdb149). GATE 1 PASSED. Committing the milestone on a feng branch.

### Arm A extraction — IN PROGRESS (user chose this as the next major step)

- **dCas13d field DONE (Claude, 2026-09-01).** Dataset S4 directly tabulates per-guide
  position (`start/end pos`, `Region` E6/I6/E7/I7/…) + precomputed combined `Z-score` —
  no enrichment/raw-SRA reprocessing needed. Extracted 1,937 guides (**1,625 SMN2**) →
  `$SCRATCH/feng/phase0/splicerush_dcas13d_field.csv`. Coverage DENSE: **842 guides in
  intron 7**, 1,032 tiling E7/I7 (minigene pos 296–835) — exactly the ASO-walked region.
  **Verdict upgraded AMBER → GREEN-leaning.** (Details in
  `artifacts/phase0_arm_a_overlap/README.md`.)
- **ASO-walk harmonization DONE + adversarially reviewed (2026-09-01). VERDICT: RED.**
  Digitized Hua 2007 / Hua 2008 / Singh 2010 walks, harmonized to S4 minigene frame, in
  `artifacts/phase0_arm_a_overlap/harmonized_aso_positions.csv` (108 records) + README
  "ASO-walk harmonization" + "Adversarial review" sections. First pass claimed GREEN (98
  "distinct exact" windows); Codex sol adversarial review **retracted it as claim inflation**
  — the defensible count of *independent* (pairwise non-overlapping) recoverable windows is
  **11 < 25 kill threshold**. So **Arm A's overlap kill-gate is RED**: StericScope cannot be
  built as a held-out transfer-law fit on this overlap; it degrades to a descriptive
  cross-effector comparison up front. Arm D (locked nusinersen benchmark) remains the floor.
  - **OPEN review findings (logged in README):** (3) effects only qualitatively "recoverable"
    — no standardized numeric outcome variable yet, so even 11 is coordinate-only; (4) S4
    coordinate frame is **off by +9 nt** vs recorded ASO coords (must fix before any downstream
    modeling; doesn't move the gate — coverage stays 98/98); (5) Hua 2008 intron-6 microwalk
    (17 ASOs) not yet digitized; (6) exact-interval dedup masks Hua-vs-Singh outcome conflicts
    → need a preregistered study-specific outcome/aggregation model.
  - **DECISION (user, 2026-09-01): CONTESTED the RED — contest LOST. RED STANDS.**
    (`artifacts/phase0_arm_a_overlap/CONTEST.md`, Codex sol.) Pairwise-non-overlap IS the wrong
    independence unit (spatial effective dimension of the 115-window walk = 12.1, not 11), so
    the narrow contest claim is right — but immaterial. Against the plan's FROZEN bar (≥25
    effective obs; 80% power for ρ=0.5 at one-sided α=0.01): spatial effective-N of the
    numerically-recoverable subset = **10.9** (sensitivity ell=3–10nt → 7.7–16.9, all <25);
    minimum detectable Spearman ≈ **0.85** at 80% power; power at the frozen ρ=0.5 ≈ **19%**.
    And **0 analysis-ready standardized numeric ASO outcomes on disk** (86/115 windows are
    digitizable from figures but not poolable across 3 studies / assays / chemistries; study is
    spatially confounded with region). So Arm A CANNOT be an adequately-powered held-out transfer
    test on this overlap.
  - **RESOLUTION: Arm A → DESCRIPTIVE** cross-effector comparison (not held-out transfer). Under
    the four-concurrent-headlines framing this does not sink the project; B, C, D stand on their
    own falsifiable cores. **Only path to revive Arm-A-as-transfer:** add genuinely new,
    independent numeric ASO datasets beyond Hua 2007/2008 + Singh 2010 (raise effective-N ≥25).
    Findings #3 (numeric outcome) and #6 (study-specific outcome model) remain OPEN but are moot
    for the transfer claim unless that new data arrives.
  - **Arm A descriptive SPEC written + honesty-verified (2026-09-01):**
    `artifacts/phase0_arm_a_overlap/ARM_A_DESCRIPTIVE_SPEC.md` — a single-frame SMN2 field atlas
    (raw dCas13d Z tracks + per-study ASO lanes, direction-categorical, 11 display blocks,
    study-separated tables), inferential pooling prohibited, self-specs a forbidden-string linter
    scan, revival data/power target set. States limited to `candidate`/`paper_supported`. NOT yet
    implemented (that's real code + figures — next Arm-A build step).

### Phase 2 variant recon — DONE (recon memo), verdict NEEDS-EXTRACTION (2026-09-01)

`artifacts/phase2_variant_recon/README.md` (Codex sol). Provisional verdict **NEEDS-EXTRACTION**.
- **GO:** ≥2 unrelated donors with public, callable, read-backed, **copy-resolved** SMN2
  haplotypes carrying an ALT inside the verified 18-nt footprint.
- **NO-GO:** <2 qualify, coordinates ambiguous, or the threshold needs unphased calls,
  same-donor copies, HPRC-as-patients, or controlled data.
- **Biggest uncertainty / the extraction to run:** whether the public **HapSMA per-haplotype
  table or Paraphase bamlets** yield two donor-independent target-site ALTs after coordinate
  reconciliation + read-level QC.
  - **Extraction ATTEMPTED 2026-09-02 -> UNRESOLVED (environment-limited, NOT a true NO-GO).**
    Codex reconciled the 18-nt footprint (GRCh38 SMN2 `chr5:70076584-70076601`; Paraphase SMN1
    `chr5:70952004-70952021`; HapSMA T2T `chr5:71408792-71408809`; seq `CCAGCATTATGAAAGTGA`) but
    could **not stage** the HapSMA workbook / Paraphase bamlets in its sandbox (network-blocked,
    same limit that hit the PMC helper). It returned "NO-GO -- insufficient public evidence" with
    **0** ALTs, explicitly flagged "not a biological zero-variant claim." **NOT settled** -- data
    unreachable, not audited. **Redo the fetch on a networked node (Claude-side or sbatch), not in
    Codex's sandbox,** before any GO/NO-GO stands.

PENDING:
- Phase 2 extraction REDO on a networked node (HapSMA/Paraphase actually staged + audited) ->
  real GO/NO-GO. Prior NO-GO is environment-limited only.
- Arm A descriptive atlas implementation (per the spec above). Unblocked, no heavy compute.

### Phase 3 -- GSE307431 clean-room alignment: COMPLETE + verified (2026-09-02)

**STAR alignment DONE: 20/20 samples, 0 failures.** Verified on disk: 20 sample dirs +
20 `.done` sentinels at `$SCRATCH/feng/phase3/star/<run>/`, each with a sorted BAM
(~2.2 G) + `ReadsPerGene.out.tab` (~1.9 M); no empty BAMs; 107 GB total footprint. FASTQs
staged + md5-verified at `$SCRATCH/feng/phase3/fastq/` (40 files). Monitor closed.
**NEXT: rMATS ΔPSI** — needs the GSM->condition mapping (3 ASO chemistries F18OMe/F18MOE/F20PMO
vs control; 20 samples from GEO GSM9224134-9224153 / manifest `sample_accession`) to define
the group contrasts. Not yet scripted. Do-not still applies: off-target splicing OK, NO
per-chemistry model. History of the 5 staging submits below is kept for the record.

<details><summary>Phase 3 submission history (resolved)</summary>

Target **GSE307431 / PRJNA1321660** -- off-target transcriptome of 3 ISS-N1 ASO chemistries
(F18OMe/F18MOE/F20PMO) + controls in SMA fibroblasts (GM03813). **20 samples, paired-end**
Illumina (~2x100bp, 25-35M pairs each). Feeds **B/RBPShadow** (ASO off-target splicing).
Do-not: reanalyze off-target splicing, but do NOT build a per-chemistry model (one seq/site).
- Run manifest (URLs + MD5s): `artifacts/phase3/ena_manifest.tsv` (20 runs, ~60 GB gz).
- Compute est ~80-120 core-hr (STAR dominates); peak ~40 G RAM/task, ~150-200 GB SCRATCH
  (700 GB free). Small job, NOT the Arm-C tail.
- Plan: ENA direct download + MD5 verify (sbatch on networked node -- Codex sandbox can't fetch)
  then STAR array (20), `--readFilesCommand zcat --quantMode GeneCounts`, sorted BAM, index
  `gencode-v50-sjdb149`, mit_preemptable + requeue-safe, then rMATS.
- **SUBMITTED 2026-09-02 (Claude-authored, no Codex).** Live handles:
  - Scripts: `feng/scripts/phase3_stage.sbatch` (partition `mit_data_transfer`, has internet),
    `feng/scripts/phase3_star_array.sbatch` (array `1-20%10`, `mit_preemptable`).
  - Jobs (CURRENT): **staging `21796885`**, **STAR array `21796886`** held `afterok:21796885`.
  - Submit env (REQUIRED to resubmit): `SCRATCH=/home/dashii/orcd/scratch`,
    `FENG_ENV_PREFIX=/home/dashii/orcd/pool/feng-env`, `FENG_P3_MAXPAR=2`; submit from repo root
    `/home/dashii/smt-isef`. (37/40 files already `.ok`-verified on scratch; resubmit fetches only
    the rest.)
  - Progress monitor: Monitor task `bgbw9e1af` (persistent; breaks on STAR-array terminal or
    dependency-never-satisfied).
  - Stage script: `-sS`, per-file 12-attempt retry loop w/ `-C -` resume (env `FENG_P3_ATTEMPTS`),
    up to `FENG_P3_MAXPAR` concurrent (default 4; running the tail at 2). Idempotent + resumable.
  - **Four prior submits FAILED + superseded** (staging skips verified `.ok` files each time):
    (1) `21791561`: curl 52 (ENA empty-reply) — added retry loop. (2) `21792534`:
    `--no-progress-meter` unknown to curl 4.4.20 (rc=2) — swapped to `-sS`. (3) `21793149`:
    serial ~1 MB/s would exceed 12 h wall (user can't raise TimeLimit — admin-only) —
    parallelized 4-way. (4) `21793773`: got 37/40 but 3 files hit curl rc=18 (partial transfer,
    ENA truncation worsened by 4-way) past 6 retries — bumped to 12 attempts + MAXPAR=2.
    (5) `21796885`/`21796886`: fetched the last 3 files, STAR ran 20/20 clean — COMPLETE.

</details>
  - Outputs: FASTQs -> `$SCRATCH/feng/phase3/fastq/`; STAR -> `$SCRATCH/feng/phase3/star/<run>/`
    (per-sample `.done` sentinel; both scripts idempotent/requeue-safe -- safe to resubmit).
  - Resume after /clear: if jobs still queued/running, just re-arm a monitor on 21791561/21791562;
    if they died, re-`sbatch` both (with the env above) -- staging skips verified files, STAR
    skips `.done` samples. Next after alignment: rMATS ΔPSI (not yet scripted).

Staged data (scratch, not repo): `$SCRATCH/feng/phase0/` — SpliceRUSH supplements
(moesm4.zip → Dataset S1–S4, moesm6.xlsx source data) + extracted dCas13d field CSV.
`$SCRATCH/feng/gencode/v50/` — GENCODE FASTA+GTF. `$SCRATCH/feng/star/gencode-v50-sjdb149`
— STAR index. `openpyxl` added to env (+ installed) for xlsx parsing.

Repo working tree also has untracked `work/` and Stage-A `.out` logs from the OLD arm —
unrelated to `feng/`; leave them alone.

## Runner / env notes

- feng env prefix (being built): `/home/dashii/orcd/pool/feng-env`.
- Interim pytest runner (pure-Python units only): `/home/dashii/orcd/pool/smt-isef-env/bin/python`.
- Load tooling with `module load miniforge/25.11.0-0` (provides conda+mamba).
- Quota (2026-09-01): HOME 42%, SCRATCH 28%, POOL 14% — room to build.

## Immediate next steps

1. Confirm feng-env built + tool versions; wire + run GATE-1 smoke (tiny fixtures).
2. GENCODE v50 download + STAR index sbatch (1.4). Then GATE 1 is closed → Phase 2.
3. **Pre-flight feasibility probe (advisory).** Arm A lives or dies on *overlap*: do the
   SpliceRUSH 1,937-guide dCas13d field and the historical ASO walks (Hua 2007/2008,
   Singh 2010) co-locate on the same SMN2 frame with ≥25 recoverable ASO windows? ~1 day
   from published supplements + SpliceRUSH source; de-risks months 5–7. Thin overlap →
   Arm A degrades to descriptive up front.

## Working style / ownership

- **Clean-room.** No imports or path-references from `src/`, `pre/src/`,
  `design_engine`, `offtarget`, `physics_adapter`, `design_panel`. Published logic may
  be *re-derived*; code may **not** be moved or renamed in. A manifest test must fail on
  any such reference.
- **Freeze discipline.** `feng/terminal/` is sealed until the month-8 preregistration
  freeze tag. Feature/model code must fail a test if it can read `terminal/` pre-freeze.
  No terminal contrast or rank result is inspected before freeze; month-9 is one shot,
  no refitting.
- **Standing gate (CLAUDE.md).** Every arm goes build → Codex `gpt-5.6-sol` adversarial
  review → fix/log → advance. Do not move on with open unlogged findings.
- **Delegation.** Substantive work goes to Codex via the `codex:codex-rescue` subagent;
  Claude orchestrates and relays. Only trivial fixes done inline.
- **HPC courtesy.** Shared quota-limited slice. All heavy work via `sbatch`,
  checkpointed, preemptable-safe. Arm C (5,000–20,000 CPU core-hours) is the fat tail —
  shard it and confirm it earns its cost before committing.

## Do-nots (see plan "Traps" for the full list)

- No chemistry model on F18MOE/F18OMe/F20PMO (one sequence, one site).
- No gnomAD/1kGP target-region calls as copy-specific SMN2 evidence; no HPRC assemblies
  as patients.
- No waiting on controlled-access data (dbGaP/EGA/author-mediated). "Application
  submitted" is not data.
- No stacking SpliceAI/Pangolin/SPLAM/IntaRNA as an ensemble without a defined causal
  role + locked ablation.
- No post-hoc repair of a failed terminal model presented as confirmatory.

## Verification (what "done" means for `feng/`)

- `pytest` green in the `feng` env; scanner verified against brute-force toy references,
  zero seed-related false negatives (POLR2H must **not** be special-cased).
- STAR/IntaRNA/rMATS one-shot smoke runs pass.
- Clean-room import test passes; freeze-seal test passes.
- One command regenerates every headline number from the data manifest.

## Honesty architecture (enforced in code)

Result states limited to `candidate`, `paper_supported`, `clean_reanalysis_reproduced`.
No fields named `clinical_response`, `patient_responder`, `safe`, `efficacious`,
`superior_to_nusinersen`. Report generation refuses "precision"/"FDR" when the
negative-label manifest is incomplete. Empty output renders "no candidates detected
under this specification," never "zero risk." Chemistry is mandatory metadata; pooling
chemistries is prohibited. Claim linter blocks "first", "genome-wide safety map",
"clinical response", "beats nusinersen", "generalizable" unless overridden and logged.
Full verbatim disclaimers are in the plan's "Honesty architecture" section.
