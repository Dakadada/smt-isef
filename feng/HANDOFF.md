# HANDOFF — feng/ (the new prospective arm)

**Created: 2026-09-01. State: greenfield. Nothing built yet. Plan is frozen.**

This file is the entry point for anyone (Claude or Codex) picking up the `feng/`
arm. The old `src/`/`pre/` program is a separate, committed effort — see
`docs/HANDOFF.md` for it. `feng/` does not touch that code.

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

| Arm | Name | Role | Falsifiable core |
|---|---|---|---|
| A | **StericScope** | HEADLINE | Latent field × effector steric-shadow kernel predicts the *other* modality's splice modulation. Fit on dCas13d, no ASO labels; test on ASO. |
| B | RBPShadow | mechanism | Steric ASO = local RBP knockout; ASO ΔPSI direction matches that RBP's KD direction (ENCODE eCLIP + KD-RNA-seq). |
| C | Bind→Act gate | safety | Binding and action are separate gates; held-out actionability model beats an edit-distance/energy baseline (leave-one-ASO-out). |
| D | Locked nusinersen benchmark | FLOOR | Accessibility-aware duplex scoring ranks the 8 supported nusinersen off-target exons above matched exons better than sequence distance. Guarantees a real result even if A–C fail. |
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

PENDING (next, on Codex reset ~9:25 PM):
- Arm A probe execution TODOs (Dataset S4 enrichment → dCas13d field; Hua 2007/2008 +
  Singh 2010 walk digitization; coordinate harmonization; ≥25-window count → finalize
  GREEN/AMBER/RED). Data staged at `$SCRATCH/feng/phase0/`; plan in
  `artifacts/phase0_arm_a_overlap/README.md`.
- Phase 2 variant recon (SMN2 haplotype audit, kill-gated).
- Phase 3 clean-room alignment of GSE307431 (now unblocked — STAR index ready).

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
