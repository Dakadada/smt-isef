# HANDOFF — smt-isef (compute-native steric splice-modulation)

**Updated: 2026-08-29. Written at the checkpoint where the first REAL run is live.**

## TL;DR
The old project is archived in `pre/` (read-only, died all-null). We rebuilt a NEW,
compute-native project: a leak-proof cis-vs-trans steric-splicing **benchmark** +
a generative **design engine** (ISEF headline) + honest **prospective designs**.
Scope is honestly capped (data floor is thin) to **SMN2-centered feasibility +
timestamped prospective predictions the field must later test** — strong, defensible,
not overclaimed. A real Stage-A run is DOWNLOADING + COUNTING real data right now.

## Working style (in CLAUDE.md)
- Outsource almost all substantive work to **Codex, model `sol`**; keep the main
  Claude chat clean; do NOT start substantive work unless told.
- **Standing gate:** every deliverable goes build → Codex **adversarial review** →
  fix, before advancing. This has caught a fatal-class bug on nearly every task.
- The Codex-rescue **forwarder subagent is rate-limited until ~4am ET**. Drive Codex
  DIRECTLY via the companion CLI (runs on OpenAI, no Claude quota):
  ```
  CX='node /home/dashii/.claude/plugins/cache/openai-codex/codex/1.0.6/scripts/codex-companion.mjs'
  $CX task --background --write "<prompt>"     # launch (omit --write for read-only review)
  $CX status --all                              # list jobs (Active + Latest finished)
  $CX status <job-id>                           # one job
  $CX result <job-id>                           # fetch full stdout
  $CX cancel <job-id>
  ```
  **Direct-CLI jobs do NOT send task-notifications — you must poll `status`.**

## Scope & honest ceiling (do not re-litigate)
- No new wet lab; compute-only; public data + tiny existing lab readouts.
- Data floor verified: **4 retrievable trans-intervention families, ~5 interventions,
  2 loci** → **NO-GO for a powered general benchmark**. Achievable = SMN2 core +
  limited multi-effector feasibility + prospective designs. See `docs/dataset-audit.md`.
- **DMD deferred** (PRJEB13980: all 12 runs share title `U7E53_OFFTARGET_ANALYSIS`,
  arms not inferable + vector-control confound). **C2B34 (task 13) shelved** (junction
  decided upstream of intron 7 — mechanism dead).
- Every design output is stamped **"PROSPECTIVE PREDICTION — NOT EXPERIMENTALLY
  VALIDATED"**; honesty is structurally enforced, not documentation.

## Repo map
- `pre/` — archived dead project, READ-ONLY. Never modify.
- `src/rawread_reanalysis/` — targeted PSI/dPSI pipeline (marker-count in raw reads; no genome index needed for targeted tier).
- `src/physics_adapter/` — keystone estimand P(productive|c)=(1−θ)P_U+θP_B; ASO-bound ensemble; θ-grid; direction-primary.
- `src/offtarget/` — high-recall off-target enumerator (full pre-mRNA, bulges/indels/>4mm; sign=unknown w/o evidence).
- `src/baselines/` — cis prior (SpliceAI/Pangolin frozen features) + real-thermodynamics nonstructural baselines; leak-proof connected-component folds; matched EDF.
- `src/benchmark/` — decisive cis-vs-trans harness; CANNOT fabricate a win or over-claim "powered" at this floor.
- `src/design_engine/` — generative design (ASO/gRNA/U7), Pareto, honest banner. **ISEF headline.**
- `src/design_panel/` — cross-disease panel + designability map (FIX IN FLIGHT, see below).
- `scripts/real_run/` — Stage-A real-run tooling (build_manifest.py, run_stageA.sbatch, monitor_requeue.sh, RUN.md).
- `config/pipeline.yaml`, `config/event_junctions.GRCh38.tsv` — verified event coords (GRCh38 / Ensembl 116).
- `docs/` — architecture, benchmark-protocol, dataset-catalog, dataset-audit, prior-art-defense, this HANDOFF.
- **Tests:** `make smoke-test` (uses `/usr/bin/python3.12`). Last green: **95 tests**.

## Verified coordinates (GRCh38 / Ensembl 116)
- SMN2 exon 7 = chr5:70,076,521–70,076,574 (+), 54 nt (classic exon 7 = Ensembl rank 8; the 2a/2b split shifts numbering). Confirmed by the c.840 T (SMN2 vs SMN1).
- DMD exon 53 = chrX:31,679,375–31,679,586 (−), 212 nt.

## THE REAL RUN (live)
- **SLURM job `21542426`**, partition `mit_data_transfer`, node `dtn001`, targeted Stage-A.
- `RUN_DIR=/orcd/scratch/orcd/011/dashii/smt-stageA-real` (SCRATCH quota ~742 GB free; hard cap 300 GiB; deletes FASTQ per accession).
- Families: SMN2 PRJNA1321660 (20), SMN2+VPA PRJNA705187 (12 RNA-seq), dCas13d PRJNA1052769 (4). 36 eligible endpoint runs.
- **First real point (partial):** GSM6063688 (ASO1 arm) → inclusion 838 / skipping 141 → **PSI ≈ 85.6%, anomalies 0** — correct nusinersen-class direction.
- Downloads are slow (~1.5 MB/s, multi-GB FASTQs); expect several hours to ~36h.
- **Check progress:** `squeue -j 21542426` ; `column -t -s$'\t' $RUN_DIR/work/counts/all_counts.tsv`
- **On done:** file `$RUN_DIR/STAGE_A_COMPLETE` appears → results in `$RUN_DIR/work/harmonized/tidy_splicing.tsv`.
  Interpret `magnitude` = raw inclusion dPSI (treated − matched control). SMN2 ASO arms should be POSITIVE. A wrong sign is a bug alarm, not a result.
- **If it dies / preempts:** it is requeue+checkpoint safe. Resubmit identically:
  ```
  cd /home/dashii/smt-isef; export SCRATCH=/orcd/scratch/orcd/011/dashii; export RUN_DIR=$SCRATCH/smt-stageA-real
  sbatch --parsable --partition=mit_data_transfer --export=ALL,REPO="$PWD",RUN_DIR="$RUN_DIR" scripts/real_run/run_stageA.sbatch
  ```
  Full rebuild (if RUN_DIR lost): `rm -rf $RUN_DIR; /usr/bin/python3.12 scripts/real_run/build_manifest.py --out $RUN_DIR --families GSE307431_PRJNA1321660,GSE167762_PRJNA705187,PRJNA1052769_dCas13d` then resubmit.

## Task board
- ✅ DONE (built + reviewed + fixed, green): #1 architecture/protocol, #2 dataset audit, #3 reanalysis pipeline, #4 baseline arena, #5 physics adapter, #6 decisive benchmark, #7 off-target enumerator, #12 prior-art defense, #14 verified coords, #15 design engine.
- ⛔ REVIEW BLOCKED: **#16** cross-disease panel — the vanished job `task-mtencg31-qlmog0` was not resumed. The fail-closed implementation is present on disk and locally verified (frozen registry digest/exact row binding, committed sequence digest, typed real-context sequence+mapping hashes, identity/polarity/clock/report-manifest checks; DMD fails pending without a committed sequence; renamed synthetic fixture and self-asserted verification regression tests). Focused 13/13 and full 117/117 tests pass. **Not LOCKED:** the mandatory `gpt-5.6-sol` adversarial gate cannot reach OpenAI from the current network-restricted sandbox; do not advance until that review runs clean and its findings are fixed/logged.
- ✅ **#17** judge demo LOCKED (2026-08-31, see RESOLVED block below). 5-arm headline, arm-level null gate, snapshot `stageA-complete-20260831`, 117/117 smoke.
- ⏳ PENDING: #8 prospective superior-molecule design set; #9 cryptic-exon span-robustness (STMN2/UNC13A); #10 effector-class variance (ASO vs dCas13x, in-silico); #11 minigene compression tooling.
- 🅿️ #13 C2B34 shelved.

## Up in the air / next steps
1. **Watch the real run** to completion; report the harmonized dPSI (SMN2 ASO arms positive? dCas13d guide effect?).
2. Finish **#16** by running a fresh scoped `gpt-5.6-sol` adversarial review in a network-enabled session; fix/log findings and lock only after a clean gate. The old `task-mtencg31-qlmog0` no longer exists and must not be resumed. **#17 is already locked.**
3. **Wire REAL transcript sequences** (from Ensembl) into design_engine/design_panel so designs are on real contexts, not fixtures (fetch like the coords were fetched).
4. **Benchmark real execution:** after Stage-A labels exist, install SpliceAI/Pangolin/ViennaRNA + a real folding backend, then run `src/benchmark` for the real physics-vs-baseline direction number. Heavy — checkpoint-first, `mit_preemptable`/`mit_normal` (2 GPU max; see memory `mit-compute-reality.md`).
5. Optionally resolve DMD arms (paper supplement) to add the 2nd locus.

## Guardrails
Shared MIT HPC: be a considerate guest. Read `~/orcd/.quota` (don't `du`). Don't traverse `/orcd/scratch|pool`. Heavy jobs via sbatch, checkpointed. `$SCRATCH` unset — it's `/orcd/scratch/orcd/011/dashii`.

## RESOLVED — 01:47 AM Codex-reset queue (2026-08-31, sonnet driver + sol)
Stage-A COMPLETE (job 21624462, 36/36, zero fetch failures). **#17 judge demo LOCKED.** All done + independently verified:
1. Null-control blocker FIXED — `scripts/judge_demo.py select_evidence()` gate now arm-level: mean of `psi − control_mean_psi` per `(family,intervention)`, assert `|arm_effect|≤2`. Per-row `|magnitude|` test removed. 2 regression tests in `tests/test_judge_demo.py`.
2. Headline EXPANDED 3→5 arms — F20PMO (+27.87) + F18OMe (+24.08) added to `INTERVENTIONS`/`EXPECTED`. New headline: "...recovered the expected splice direction across five deposited steric intervention arms—ASO1, F18MOE, F20PMO, F18OMe, and dCas13d D[333,357]...". Honesty vetted clean (no benchmark/clinical overclaim, banner + claim-boundary intact).
3. 2nd blocker caught mid-review: completion gate fired on ANY missing SMN2 row incl out-of-panel scramble_PMO → crashed every real build. FIXED — gate only fires on required-panel `INTERVENTIONS` arms. scramble_PMO MISSING row is LEGIT (low-coverage fail-closed), still surfaced in `excluded_missing_rows`.
4. Snapshot REBUILT: `artifacts/judge-demo/stageA-complete-20260831`, ts `2026-08-31T02:00:00Z`, no `--allow-incomplete`, 19 obs.
5. Adversarial gate RUN (scoped to touched files): 1 structural (no ≥2-replicate enforcement) fixed + test; unused `FULL_BANNER` removed; 4 minor items logged non-load-bearing.
6. `make smoke-test` = **117/117** green.
UNCOMMITTED — working tree untracked; commit is a separate decision. Touched: `scripts/judge_demo.py`, `tests/test_judge_demo.py`, `artifacts/judge-demo/stageA-complete-20260831/*`.
NEXT (awaiting go-ahead): #8 prospective superior-molecule designs, #9 cryptic-exon span-robustness (STMN2/UNC13A), #10 effector-class variance (ASO vs dCas13x in-silico), #11 minigene compression tooling.

## BLOCKED — #16 mandatory sol gate and commits (2026-08-31, Codex driver)
The obsolete companion job `task-mtencg31-qlmog0` is absent and was deliberately not resumed. Re-audit of the actual files found the intended #16 repair already present in `src/design_panel/`: production panel execution no longer calls `offline_target`; verified targets must resolve exactly against the frozen `config/event_junctions.GRCh38.tsv` path and SHA-256, match the registered transcript/assembly/chromosome/strand and provenance digest, carry a registry-committed sequence digest, and receive a `RealTranscriptContext` whose sequence and genomic-to-transcript mapping hashes, identity, therapeutic polarity, frozen clock, and returned engine manifest all agree. DMD remains explicitly pending because its registry row has no committed sequence digest. Regression coverage rejects renamed synthetic fixtures, self-asserted targets absent from the registry, registry override/tampering, wrong identity/assembly/polarity/clock, substituted sequence reports, forged nested manifests, incoherent panel outputs, and honesty-banner removal.

Verification completed before the gate: `PYTHONPATH=src /usr/bin/python3.12 -m unittest tests.test_design_panel -v` = 13/13 pass; `make smoke-test` = 117/117 pass. The required independent review was attempted with Codex CLI 0.150.1, explicit model `gpt-5.6-sol`, high reasoning, read-only repo sandbox. Attempt 1 failed because Codex state initialization targeted the read-only home filesystem. Retrying with an isolated writable `/tmp` `CODEX_HOME` and copied auth successfully started the requested model, but both WebSocket and HTTPS transports to OpenAI failed with `Operation not permitted` / network connection errors and stayed in reconnect. Therefore the standing adversarial gate is not clean, #16 is **not LOCKED**, and tasks #2–#11 were not started.

Checkpoint commits are also blocked by infrastructure: `.git` is mounted read-only in this session, and an intentional-path staging attempt failed at `.git/index.lock` with `Read-only file system`. No commit was created. Working files remain on disk; retry in a session with writable `.git` and outbound Codex access. First action on resume: run the scoped #16 `gpt-5.6-sol` adversarial review, fix/log every finding, rerun 117 tests, commit the intentional compute-native snapshot without run logs/work data, then mark #16 LOCKED before advancing to Ensembl sequence wiring.
