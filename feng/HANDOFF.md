# HANDOFF — feng/ (the new prospective arm)

Entry point for anyone (Claude or Codex) picking up the `feng/` arm. The old
`src/`/`pre/` program is a separate, committed effort — see `docs/HANDOFF.md`.
`feng/` does not touch that code.

> **History trimmed 2026-09-09.** The blow-by-blow Arm D publish-atomicity review
> saga (R1–R4, H1–H5) and the Arm E hoff bug#1–#4 / site-15 saga are SUPERSEDED
> (both arms' compute is out) and live in git history of this file (working tree at
> `f0b0d48` + uncommitted). Earlier dated checkpoints (09-01…09-04) are in
> `feng/HANDOFF_ARCHIVE.md`. This file = LIVE state + evergreen orientation only.

---
## ⏩ RESUME AFTER /clear — LIVE: Phase-1 compute DONE, review in flight (2026-09-10)
GOAL: get Phase-1 GPU compute out. Compute COMPLETE + validated. Remaining: review → seal → commit.

- **JOB 22343195 — ALL 15 COMPLETED (2026-09-10 ~01:00).** L10 (753k sites, task_5) was
  the long pole, exit 0:0. `export SCRATCH=/orcd/scratch/orcd/011/dashii`.
- **VALIDATION PASS (2026-09-10, no burden VALUE read):** 15/15 `.done` markers (0-byte,
  expected) + 15 non-empty JSON in `$SCRATCH/feng/arm_e/hoff`, all state=candidate +
  terminal=True + manifest_sha256 + input_digest + masking_method='N-ablation sensitivity
  proxy'. L02/04/06 pangolin_off_target_incomplete=True n_sites=0; L08..L30 incomplete=False
  (L10 n_sites=753089). KEY NAMES: `pangolin_off_target_incomplete`, `terminal` (not `incomplete`).
- **REVIEW DONE (task keqt8dxnt, sol high) → SEAL-READY: NO, 6 blockers.** NO Pangolin
  rerun needed (verified): #4 digest gap = provenance only (d=50/context=5000 fixed consts,
  results deterministic → re-stamp digest, seconds); #1 silent-zero never fired (writer always
  co-writes int site_index+burden). Other blockers = code/text/commit-sequence.
- **FIX IN FLIGHT:** MCP Codex gpt-5.6-sol high, task **k6d0g6kgb**, workspace-write,
  freeze-safe (burden byte-identical, no GPU rerun). Fixing all 6: (1) strict int resume,
  (2) exact bool/count validation, (3) prereg text = complete-zero vs uninformative-empty,
  (4) bind d+context+scorer-hash into input_digest + re-stamp 15 payloads (*.json.prebak
  backups), (5) self-contained CODE commit (Arm E launcher+scores+tests together),
  (6) prereg stays DRAFT, seal in a SEPARATE later commit.
- **COMMIT PLAN (revised, resolves #5+#6):** commit-1 = ALL endpoint code + Arm D support +
  tests + re-stamped payloads + prereg-as-DRAFT (self-contained). commit-2 = flip prereg
  DRAFT→SEALED stamping commit-1 SHA + job 22343195 + manifest sha.
- **FIX DONE (task k6d0g6kgb):** all 6 fixed, 556 pass / 5 skip / 0 fail. BUT restamp #4
  did NOT run (sandbox blocked SCRATCH) — 15 payloads still OLD digest, no *.prebak yet.
  Manual restamp cmd (run w/ SCRATCH write): `TMPDIR=/dev/shm PYTHONPATH=feng/src
  /home/dashii/orcd/pool/feng-env/bin/python -m feng.titration.restamp_arm_e_hoff
  --scratch-root /orcd/scratch/orcd/011/dashii` (expect: re-stamped 15 payloads + 12 ckpts,
  input_digest bytes only).
- **PROVENANCE-HONESTY ISSUE (assessed, no rerun):** payload `inputs` dict recorded genome/
  manifest/weight/pangolin-version + scan_burden_sha256, but NO scorer-code hash at launch
  (that omission = blocker #4). Sealing binds CURRENT gpu_scores.py SHA. Claim: burden formula
  max(0,gain)+max(0,-loss) was PRE-launch fix (job ran it); fix pass touched only non-behavioral
  gpu_scores.py (naming, commutative by_id order, finite guards) → binding current SHA HONEST iff
  prereg states scorer-hash reconstructed post-hoc. Terra charged to confirm/refute (Q2 equivalence).
- **TERRA RE-REVIEW DONE (task k48vy0kh2) → ADVANCE-TO-SEAL: NO.** Numbers SOUND (reorder
  commutative on finite results, rename inert — equivalence CONFIRMED for the actual run). But
  Q2=CANNOT-CONFIRM: no launch source snapshot → a code-AUTHENTICATED seal needs recompute from
  an immutable commit; binding current SHA is only a weaker retrospective claim. Q3=OVERCLAIM
  (prereg "fixed-code" wording). +MAJOR restamp not transactional/idempotent. +MINOR prereg
  self-contradiction on DRAFT-in-commit.
- **★ DECISION: PATH A (user, 2026-09-10)** — code-authenticated recompute from an immutable
  sealed commit. Restamp ABANDONED (fresh run writes correct digest). Recompute ~20h (L10 long
  pole) + re-freeze. Gold-standard: SHA genuinely authenticates runtime code.
- **PATH A PREP v1 SUPERSEDED:** the single-file hash/read-only approach below failed whole-tree,
  child-process, launcher, and crash-recovery authentication review; do not launch that version.
- **PREP RE-REVIEW (task ka3mptfiu, terra HIGH) → LAUNCH-READY: NO.** 3 blockers: (1) chmod≠
  immutable + child worker re-imports w/o check = TOCTOU; (2) only gpu_scores.py bound, but burden
  depends on whole import closure (Check#1=NO); (3) supersession not crash-atomic/locked. +MAJOR
  launcher doesn't prove own bytes from SEALED_COMMIT. +MINOR digest=inputs not execution.
- **AUTH-FIX VERIFIED landed (02:04, Codex reset over).** kc1dv9cxc failed on usage limit at
  summary stage but code landed: source_auth.py + submit_arm_e_hoff.sh + test_source_auth.py
  present, compile + bash -n clean, auth tests GREEN (test_source_auth 5/5, test_arm_e_launch
  56/56). NOTE: full suite in ONE process gets OOM-killed (slice MemoryMax=10G) — run per-file,
  never the whole suite at once.
- **AUTH REDESIGN RE-REVIEW DONE (terra round-2):** LAUNCH-READY: NO. Supersession + launcher/
  digest-split now CLEAN (2 old blockers closed). 3 remain: #1 TOCTOU (rehash moment-in-time),
  #2 pre-bootstrap import escape (inherited PYTHONPATH/.pth/sitecustomize), #3 Pangolin runtime
  code unauthenticated (only version+weight hashes bound, pangolin/model.py code executes → swap
  file, digest unchanged).
- **★ DECISION: RECOMMENDATION A (user, 2026-09-11)** — fix #3 (real drift gap, bind Pangolin
  PACKAGE CODE hash into input_digest + worker-verify), and SCOPE #1/#2 in prereg threat-model
  (defend code-drift + at-rest substitution; live-write TOCTOU / hostile interpreter = out of
  scope for single-user prereg). NOT chasing immutable-substrate/in-memory-import (exotic threat,
  heavy eng, zero science gain — result identical, only guarantee-strength differs).
- **#3 FIX DONE (thread 01a09070…, sol high):** 6 files — source_auth.py (stdlib pangolin
  package-hash helper, no import), arm_e_hoff.py (provenance + passes expected hash to workers),
  gpu_scores.py (_worker_main verifies BEFORE numpy/torch/pangolin import), prereg (code-hash
  identity + explicit TOCTOU/hostile-startup out-of-scope, STAYS DRAFT), +3 test files. input_digest
  now binds pangolin_package_sha256 + .py manifest alongside version+weights+feng/src tree+commit.
  Tests per-file GREEN: test_source_auth 8/0/0, test_arm_e_launch 60/0/0, test_titration_gpu_scores
  24/1/0. Burden formula + scored values untouched/unread.
- **#3 STANDING-GATE REVIEW DONE (task k30f3rlww, sol high) → LAUNCH-READY: NO, 2 BLOCKERs.**
  Both at-rest / code-drift (NO live-write TOCTOU) → IN-SCOPE for the A claim, must fix.
  B1: only .py hashed, but CPython can run cached .pyc (install ships model.cpython-311.pyc) or
  higher-precedence model*.so → at-rest swap leaves hash unchanged, runs substituted code.
  B2: expected hash measured fresh from the SAME live install then verified against it = circular,
  no committed trust anchor → substitution before a fresh task / drift between array tasks yields
  a NEW valid digest and passes. MAJOR: prereg "pins package code" currently false until B1+B2 close.
  Confirmed GOOD: worker verify ordered before numpy/torch/pangolin, fail-closed on absence,
  digest-sensitive+deterministic, resume via fresh workers.
- **#3 BLOCKER FIX DONE (thread 01a09070…, sol high):** both closed. PINNED_PANGOLIN_PACKAGE_SHA256
  = e5fc7aa31ea91408632b16b8fb094f6271445186f5ccea50d2083e70e698eba4 in authenticated source_auth.py
  (inside feng/src tree hash). Manifest now covers .py/.pyc/.pyo/.pyd + EXTENSION_SUFFIXES (6 files
  current) + rejects sourceless/extension/non-SourceFileLoader. Parent+worker+resume compare live ==
  pinned before scoring, both bound into input_digest. Pin seeded via stdlib-only scan under feng-gpu-env
  (no import/GPU). Tests per-file GREEN: test_source_auth 10/0/0, test_arm_e_launch 61/0/0,
  test_titration_gpu_scores 23/1/0. Burden untouched/unread.
- **#3 BLOCKER RE-REVIEW DONE (task kwo3noz6u, terra HIGH):** B1 CLOSED. B2 = closes-on-commit
  (source_auth.py/arm_e_hoff.py/launcher still UNTRACKED → pin not yet in committed tree; commit-1
  fixes it — not a code defect). 2 MINORs: M1 raw-.pyc mtime-header → benign recompile = false
  fail-closed; M2 seed manifest/command not recorded. Prereg accurate, TOCTOU retained, DRAFT.
- **#3 MINOR FIX IN FLIGHT:** MCP Codex gpt-5.6-sol HIGH, task **k72bi9ag3**, fix thread (01a09070…),
  workspace-write, freeze-safe. M1: strip 16-byte CPython header on .pyc/.pyo before hashing (benign
  recompile inert, body-tamper still caught), full-byte .py/.so/.pyd, re-seed pin. M2: committed
  feng/src/feng/pangolin_package_seed_provenance.json (manifest+interpreter+command+version, inside
  tree hash). Benign-recompile/body-tamper tests, per-file. Prereg DRAFT.
- **NEXT after MINOR fix clean:** commit-1 (CODE, feng-arm, immutable sealed — file set HANDOFF
  lines 79-88, ADD source_auth.py + pangolin_package_seed_provenance.json) → this MECHANICALLY
  closes B2 → light post-commit verify (committed tree hash == pin authenticated, launcher runs
  from commit) → **PAUSE for user launch go/no-go** (recompute ~20h GPU shared slice + freeze).
  Launch via submit_arm_e_hoff.sh <commit-1> <run-id> → re-arm state monitor (L10 long pole).
- **AUTH-FIX IMPLEMENTED (pending verify + fresh review):** whole `feng/src` is manifested from the selected
  commit's Git object store and re-hashed by a stdlib-only bootstrap in every parent/worker process
  before any `feng` import. Digest-bound identity is `source_tree_sha256` + `executed_git_commit`;
  run/job/task/launcher hash are separate execution provenance. The submit wrapper streams the
  committed launcher blob, whose executing bytes are checked at task start. A full-lifetime per-length
  `flock` plus recovery-aware moves and atomic marker rename make supersession crash-safe.
- **COMMIT FILE SETS (PATH A):** CODE commit = feng/src/feng/titration/{arm_e_scan,
  arm_e_hoff,gpu_scores}.py + feng/src/feng/bench/{scores,strata,match,permute}.py +
  feng/src/feng/{source_auth.py,gpu/checkpoint.py} +
  feng/scripts/{arm_e_reference,arm_e_scan,arm_e_hoff,bench_match}.sbatch +
  feng/scripts/submit_arm_e_hoff.sh + setup_gpu_env.sh + templates/gpu_preemptable.sbatch +
  feng/tests/{test_arm_e_launch,test_titration_gpu_scores,test_bench_scores,test_bench_universe,
  test_bench_match,test_bench_permute,test_gpu_infra,test_source_auth}.py +
  feng/PREREG_ARM_E_AXIS_H_GPU_PHASE1.md as DRAFT + this HANDOFF update. No prior payload belongs
  in CODE. SEAL commit later changes the same prereg DRAFT→SEALED and stamps the CODE commit SHA,
  authenticated recompute job id, source-tree SHA, launcher SHA, manifest SHA, and date.
- **RESUME STEPS (PATH A, after auth-fix review):**
  1. Re-review auth-fix diff — terra HIGH fresh (whole-tree/parent/worker/launcher authentication,
     supersession transaction, provenance split, and prereg framing).
  2. On clean: **commit-1 (CODE, feng-arm)** = immutable sealed producing commit. Note its SHA.
  3. Launch only through `submit_arm_e_hoff.sh <commit-1> <new-run-id>`; it streams the committed
     launcher blob and passes its SHA plus the sealed commit. Confirm old H-off payloads are moved
     into the run's per-length superseded directories so the run starts clean.
  4. Monitor (re-arm state monitor; L10 ~20h long pole). On 15/15: VALIDATE (no burden read) —
     digest binds the complete source-tree identity matching commit-1; execution provenance records
     the run/job/task and committed launcher SHA.
  5. **commit-2** = flip prereg DRAFT→SEALED stamping commit-1 SHA + recompute job id + manifest
     sha. **PHASE-1 OUT.**
  6. Do NOT run Phase-9, do NOT read burden/scored ranking (freeze). Phase-2 MD DEFERRED.
- OLD job 22343195 outputs are SUPERSEDED by the Path-A recompute (numbers were sound, but not
  code-authenticated). Move to superseded/, do not delete, do not read burden.
- NOTE: MCP Codex tasks DO fire a harness notification on completion now (both keqt8dxnt +
  k6d0g6kgb did) — the old "poll jobs json, no notification" rule is stale.

## UNCOMMITTED INVENTORY — two separate commits (after seal)
- **commit-A (Phase-1 blocker fix + prereg):** `gpu_scores.py` + `arm_e_hoff.py` +
  `test_titration_gpu_scores.py` + `test_arm_e_launch.py` + `PREREG_ARM_E_AXIS_H_GPU_PHASE1.md`.
  = the 5-blocker fix (readout last-exonic-base, N-ablation-proxy reframe + masking_method,
  directional burden max(0,gain)+max(0,-loss), strict checkpoint resume, weight+version
  provenance) + site-15 `>` fix.
- **commit-B (Arm D publish-atomicity chain + Arm E env fixes):** `strata.py` / `permute.py` /
  `bench_match.sbatch` / `test_bench_universe.py` / `test_bench_match.py` / `test_bench_permute.py`
  (H1–H5 versioned-bundle + `current` pointer atomic publish) + `setup_gpu_env.sh`
  (pandas+pyarrow+pysam) + `gpu_preemptable.sbatch` (SCRATCH self-heal) + `scores.py` N-pad.
- Also published to scratch (compute OUT, not git): Arm D scored_strata fresh key
  `55eadcfd…` (`current` pointer); rebuilt `detectable_universe_v50.parquet` + `strata_v50`.

## PRE-FREEZE CLEANUP (DONE 2026-09-09, uncommitted — folds into commit-A)
Codex gpt-5.6-sol, suite 536 pass / 5 skip / 0 fail, diff clean:
(1) MED-2 e2e footprint test `test_offtarget_edge_window_preserves_genomic_padding_and_footprint_offset`
    (test_bench_scores.py); (2) non-palindrome tie test `test_offtarget_nonpalindromic_orientation_tie_is_scored`
    + strict `rc>direct` reject (test_titration_gpu_scores.py); (3) relabel unused
    `hybridization_query_5to3` → `..._provenance` + comment in gpu_scores.py (inert rename,
    not in input_digest/burden). Still PREREG the `rc<=direct` eligibility rule at FREEZE.
NOTE: (3) touched gpu_scores.py while 22343195 pending tasks 9-15 outstanding — rename of an
UNUSED field, burden identical whichever version a pending task imports; no provenance change.

---
## ★ ARMS — compute status (2026-09-08/09)
- **Arm C: DONE — NO-GO.** Committed `f0b0d48`. Model AUC 0.838 / AP 0.042 vs best baseline
  duplex_dg 0.797 / edit 0.040 — misses +0.05 bar on both margins. Duplex energy largely
  explains action. Closed.
- **Arm D: COMPUTE ENDPOINT REACHED.** Chain 22336666 (universe) → 22336667 (match) →
  22336668 (IntaRNA array 0-31, all 32 COMPLETED) all-terminal 2026-09-08. scored_strata
  PUBLISHED to fresh cache key `55eadcfd…` via new atomic publish (`current` pointer, per-task
  --validate-output gated). Prior 09-05 key `0209bf8e…` stale, untouched. Ranking NOT read.
  Publish-atomicity fix chain H1–H5 CLOSED (round-4 re-review ADVANCE / SAFE-TO-REBUILD=YES).
- **Arm E hoff: RECOMPUTE IN FLIGHT (job 22343195 above).** Prior 15/15 .done endpoint is
  SUSPECT under structural review findings #3 (burden sign) / #4 (checkpoint silent-zero) →
  recompute MANDATORY. 54 pre-fix artifacts discarded.
- Freeze remains UNCUT, seal UNOPENED, Phase-9 untouched.

## ★ GPU PHASE-1 GATE — where it stands
- Phase-1 off-target burden = the Arm E hoff endpoint (compute scale-up of built H-off, not
  new science). PATH B (recompute at H200 scale, fixed code) chosen after structural review.
- Structural review (thread 01a0840a) returned SEAL-READY=NO with 5 blockers; all 5 FIXED
  (thread 01a08411) + re-review SEAL-READY=YES (thread 01a08424). Recompute launched (22343195).
- PREREG DRAFT written: `feng/PREREG_ARM_E_AXIS_H_GPU_PHASE1.md` (3 candidate defs frozen
  verbatim: steric N-ablation-proxy masking, exon-7 readout at last exonic base chr5:70076574,
  per-mer directional off-target burden = TOTAL not mean, n_sites is a count). STATUS DRAFT →
  flip SEALED on all-terminal validate.
- Phase-2 MD = GO-WITH-HEAVY-CAVEATS, DEFER (not before Phase-1 readout). Full plan
  `feng/ARM_E_GPU_PLAN.md`. GPL-3.0 Pangolin vendoring flag still open.

## ★ KETTERLE ARM — spun out to `../ketterle/`
Standalone (Wolfgang Ketterle): in-silico ASO specificity. Bivalent coincidence/AND-gate =
NO-GO on SMN2 exon-7. Pivoting to a different mechanism. Full checkpoint `ketterle/HANDOFF.md`.

## ★ OPERATIONAL NOTES (evergreen)
- **Codex route:** bare aliases `sol/luna/terra` are REJECTED on this account
  (`'sol' model is not supported`). Use fully-qualified `gpt-5.6-sol` (+ `-terra`/`-luna`).
  MCP codex (`mcp__codex__codex`, sandbox=workspace-write) WRITES fine — the old
  "Codex write-broken" blocker is CLEARED.
- **MCP Codex tasks fire NO harness notification** — poll
  `~/.claude/plugins/data/codex-openai-codex/state/smt-isef-1074aefacd3b422d/jobs/<task>.json`
  `status`. bg-isolation guard disabled (`.claude/settings.json` worktree.bgIsolation=none)
  so HANDOFF stays editable in the main checkout.
- **Codex app-server wedge** ("failed to clean up stale arg0 temp dirs / initialize sqlite"):
  `rm -rf ~/.codex/tmp/arg0/codex-arg0*` (needs classifier override; residual .nfs* stubs harmless).
- **Model tiering** (memory codex-review-model-policy): dev/fixes=sol high; FIRST review/arm=
  sol high; re-reviews=terra (HIGH effort when delicate); escalate to sol only on significant
  structural findings; monitoring=luna. Never cut round count. Claude fallback: sol→opus5 deep,
  terra→opus5 light, luna→haiku4.5; opus5 worker-only, main chat stays opus 4.8.

---
> Superseded dated checkpoints (2026-09-01…09-04) in `feng/HANDOFF_ARCHIVE.md`.
---

## What `feng/` is

A Feng-Zhang-shaped second arm: turn the single nusinersen result into a platform.
Full plan: `~/.claude/plans/yeah-so-just-talk-reflective-moon.md` (canonical, frozen).
Copy into `feng/PLAN.md` before real work so the repo is self-contained.

**The move.** Stop mapping *where* to place an ASO. Infer **the operator an effector
imposes**: a splicing element has a latent regulatory field; an 18-nt ASO duplex, a
25-nt dCas13d guide + bulky RNP, and a persistent U7 snRNP each cast a *different steric
shadow* on that same field. Derive the transfer law. Platform claim: a cheap programmable
RNP tiling screen becomes a chemically actionable ASO design map, with a hard falsifiable
core — **predict one effector from another.**

## The four arms (+ appendix)

**Framing (user):** there is NO single headline. All four arms A–D are concurrently
headline — four independent falsifiable cores, each standing on its own. "FLOOR" on D
means guaranteed-to-land, not lesser.

| Arm | Name | Role | Falsifiable core |
|---|---|---|---|
| A | **StericScope** | headline (contested) | Latent field × effector steric-shadow kernel predicts the *other* modality's splice modulation. Fit on dCas13d, no ASO labels; test on ASO. |
| B | RBPShadow | headline (mechanism) | Steric ASO = local RBP knockout; ASO ΔPSI direction matches that RBP's KD direction (ENCODE eCLIP + KD-RNA-seq). |
| C | Bind→Act gate | headline (safety) | Binding and action are separate gates; held-out actionability model beats an edit-distance/energy baseline (leave-one-ASO-out). **[Arm C closed NO-GO.]** |
| D | Locked nusinersen benchmark | headline (guaranteed to land) | Accessibility-aware duplex scoring ranks the 8 supported nusinersen off-target exons above matched exons better than sequence distance. Stands alone even if A–C fail. |
| Appx | variant × ASO decomposition | kill-gated | Old arm 1, reduced to what survives a 2-day recon gate. |

The original four proposed arms are all dead/crowded (responder map, exhaustive landscape,
off-target sieve, transferable grammar). Do not re-litigate.

## Runner / env notes
- feng env prefix: `/home/dashii/orcd/pool/feng-env`. GPU env: `~/orcd/pool/feng-gpu-env`
  (Pangolin 1.0.2 + 12/12 weights, pandas/pyarrow/pysam installed).
- Interim pure-Python pytest runner: `/home/dashii/orcd/pool/smt-isef-env/bin/python`.
- Load tooling: `module load miniforge/25.11.0-0`.
- Tests need `TMPDIR=/dev/shm` (rnaplfold/sigterm/hoff paths).

## Working style / ownership
- **Clean-room.** No imports/path-refs from `src/`, `pre/src/`, `design_engine`, `offtarget`,
  `physics_adapter`, `design_panel`. Logic may be re-derived; code may NOT be moved/renamed in.
- **Freeze discipline.** `feng/terminal/` sealed until the month-8 preregistration freeze tag.
  No terminal contrast or rank result inspected before freeze; month-9 is one shot, no refitting.
- **Standing gate (CLAUDE.md).** Every arm: build → Codex gpt-5.6-sol adversarial review →
  fix/log → advance. No moving on with open unlogged findings.
- **Delegation.** Substantive work → Codex via `codex:codex-rescue`; Claude orchestrates.
  Only trivial fixes inline. Give workers scoped self-contained prompts, never full conversation.
- **HPC courtesy.** Shared quota-limited slice. Heavy work via `sbatch`, checkpointed,
  preemptable-safe.

## Do-nots (see plan "Traps")
- No chemistry model on F18MOE/F18OMe/F20PMO (one sequence, one site).
- No gnomAD/1kGP target-region calls as copy-specific SMN2 evidence; no HPRC assemblies as patients.
- No waiting on controlled-access data (dbGaP/EGA/author-mediated).
- No stacking SpliceAI/Pangolin/SPLAM/IntaRNA as an ensemble without a defined causal role + locked ablation.
- No post-hoc repair of a failed terminal model presented as confirmatory.

## Verification (what "done" means)
- `pytest` green in feng env; scanner verified against brute-force toy refs, zero seed-related
  false negatives (POLR2H must NOT be special-cased).
- STAR/IntaRNA/rMATS one-shot smoke runs pass.
- Clean-room import test passes; freeze-seal test passes.
- One command regenerates every headline number from the data manifest.

## Honesty architecture (enforced in code)
Result states limited to `candidate`, `paper_supported`, `clean_reanalysis_reproduced`. No
fields named `clinical_response`, `patient_responder`, `safe`, `efficacious`,
`superior_to_nusinersen`. Report generation refuses "precision"/"FDR" when the negative-label
manifest is incomplete. Empty output renders "no candidates detected under this specification,"
never "zero risk." Chemistry is mandatory metadata; pooling chemistries prohibited. Claim linter
blocks "first", "genome-wide safety map", "clinical response", "beats nusinersen",
"generalizable" unless overridden and logged.
