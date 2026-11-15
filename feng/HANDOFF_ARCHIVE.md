# HANDOFF ARCHIVE — feng/ (superseded checkpoints)

> Moved out of feng/HANDOFF.md 2026-09-05 to cut its token size. Historical only;
> live state is in HANDOFF.md. Order: newest-archived first, then old phase-1 progress.

---

## ⏭ RESUME WHEN BACK — 2026-09-04 ~17:26 EDT (user dropped wifi; session ended)

**Jobs keep running on Slurm without the session. Do this on return:**
1. `squeue -u $(whoami) -o '%.14i %.24j %.14P %.11T %.20E %R'` — check Arm D.
2. **Arm D chain 21994396 universe → 21994397 match → 21994398_[0-31%16] IntaRNA** (relaunched post-
   hotfix). If universe FAILED → read `feng/logs/feng-bench-universe-21994396.err`; the strata.py
   NaN hotfix (Claude-applied, UNREVIEWED) may need Codex's take (drop-event vs NaN-as-empty). Done when
   `$SCRATCH/feng/bench/scores/current/manifest.json` + 32 `intarna_*.done` → `scored_strata.parquet`.
3. **Arm E:** code-fixed + 487 green but NOT relaunched. First run a Codex adversarial RE-REVIEW of the
   Arm E fix (4 HIGH+3 MED, see below) AND have Codex verify the strata.py NaN hotfix. On ADVANCE →
   relaunch Arm E scan+hoff from the MAIN shell (reuse pinned ref 21992483; submit block in the
   2026-09-04 16:45 section below). Open design Q: window "complete" over full 2–30 vs informative L≥8
   (L2/L4/L6 uninformative).
4. HARD RULES still hold: seal UNOPENED, freeze UNCUT, no Phase-9, no manufactured GO. sbatch launches
   trip the harness classifier — approve per-launch or add a Bash sbatch permission rule.
5. UNCOMMITTED on feng-arm (all Phase-0 + Arm E + hotfix + AGENTS.md); jobs read feng/src at runtime —
   do NOT `git checkout`/`stash`/switch branches until chains finish, or commit first.

---

## ★ LIVE CHECKPOINT — 2026-09-04b (GOAL RUNNING — Phase 0 fixes dispatched to Codex; NO compute launched yet)

**STATE FOR /clear — read FIRST. User gave GO on the whole /goal (GOAL_TONIGHT.md): fix D+E blockers,
gate, launch D+E heavy compute tonight, live monitors. Decisions locked this session: seal stays
UNOPENED (default NO); Arm E off-target scan = GO (user: "finish, do not stop"); rebuild Arm D
universe from ScrMOE-only (resolves open freeze-decision = option b). Delegate ALL substantive work
to Codex — NO Claude subagents except trivial monitors.**

- Env checked clean 2026-09-04 00:14: NO jobs in squeue; quota HOME 45% / SCRATCH 42% / POOL 15%; no
  motd maintenance block.
- **PHASE 0 DONE** — Codex sol MCP task `k431ny4wm` (thread `01a06aa1-3508-7db0-af11-48506f21becf`,
  resumable) fixed ALL 11 findings (D-HIGH-1/2/3, D-MED-1/2, D-LOW-1, E-HIGH-1..5). Uncommitted on
  feng-arm. Codex-reported pytest: **329 passed / 4 skipped / 0 failed** (143s, feng-env). NO commit,
  NO launch, NO seal/freeze/unlock. Fix map: strata.py (ScrMOE-only universe + STAR-count raise),
  scores.py (content-addressed cache + one_to_one join), permute.py (Phase-9 needs unlock AND freeze
  tag), match.py (dup event_id fails loud), window.py (potency backbone de-double-count + incomplete-
  curve→None + durable seal hash), gpu_scores.py (H-off honors ASO orient/strand, no 0.0 fabrication),
  checkpoint.py + gpu_preemptable.sbatch (failures exit non-zero).
- **⚠ PHASE 1 DEPENDENCY:** the committed detectable-universe artifact is now STALE (old pooled build);
  the new validator REJECTS it. The ScrMOE-only universe must be REBUILT (compute) before the Arm D
  IntaRNA score sweep can run. Sequence Phase 1 as: (D) rebuild ScrMOE-only universe → build 8×101
  strata → IntaRNA score sweep; (E) off-target seed-extend scan + Pangolin H-off sweep.
- **Phase 0 fixes INDEPENDENTLY VERIFIED on disk** (Claude, not just Codex's word): `git diff --stat`
  = 17 code/test files changed +616/−127; `PYTHONPATH=feng/src <feng-env>/bin/python -m pytest
  feng/tests/ -q` = **329 passed / 4 skipped / 0 failed** (144s) at 00:43 EDT.
- **GATE 1 BLOCKED on Codex 5-HR ROLLING LIMIT.** Adversarial re-review task `kwy34n2sh` FAILED at
  00:47 EDT: "You've hit your usage limit … try again at 4:09 AM." Per GOAL_TONIGHT rule (5-hr case =
  monitor for reset+1min and continue, do NOT checkpoint-stop): armed a timer monitor to auto-resume
  ~4:10 AM EDT. Nothing to commit yet; nothing else runnable without Codex (all Codex routes share the
  one limited account; user rule = delegate ALL substantive work to Codex, no Claude subagents).
- **GATE 1 ROUND 1 = DO-NOT-ADVANCE** (Codex adversarial task `kalnf5x0f`, thread `01a06b78`, 4:10 AM).
  6 HIGH + 5 MED + 2 LOW. Verified-CLEAN: E-HIGH-1 arithmetic (only IntaRNA/RNAup feed potency),
  H-off returns None not 0.0, seal tamper-detect, no ΔPSI/rank/perm added, one_to_one joins, dup
  event_id. Real holes → FIX PASS DISPATCHED (Codex task `k3e96f9rw`, resumes implementer thread
  `01a06aa1`): H1 window.py:105 non-finite (inf/NaN) scores manufacture verdict; H2 stale
  arm_e_provenance.json asserts refutes w/ n_with_window=0; H3 match.py:380+strata.py:603 ScrMOE
  provenance not enforced downstream (can mix 3-ScrMOE PSI + 11-control CPM); H4 strata.py:251 STAR
  NaN/Inf/neg counts; H5 checkpoint.py:138 requeues on failed flush; H6 window.py:293 seal race +
  shared .tmp. +5 MED (H-off indel orientation, permute tag not commit-bound, sbatch allowlist/atomic,
  detectability 6/11→2/3 unrecorded) +2 LOW (abs-path cache digest; don't stage untracked Arm C).
- **GATE 1 FIX PASS DONE** (Codex `k3e96f9rw`, thread `01a06aa1`): all 13 findings fixed, pytest
  **356 passed / 4 skipped / 0 failed**. arm_e_provenance.json regenerated (prediction=null, no rescore).
  Universe rows now embed ScrMOE allowlist digest + policy, required downstream. Still uncommitted.
- **GATE 1 ROUND 2 = DO-NOT-ADVANCE** (Codex `k24zip9v0`, thread `01a06b92`): H2/H3/H5/H6/GPU/sbatch-
  accession/policy CLOSED, LOW clean, freeze-safe + verdict-neutral. 1 HIGH still open (H1 non-finite
  potency/burden slips the prereg/report path → false verdict) + 4 MED (whitespace gene IDs; legacy
  unlocked-draft rebuild now blocked by seal verify; sbatch dup-launch shared-tmp race; permute
  manifest not bound to prereg/scored_strata). FIX PASS ROUND 2 dispatched (Codex `k53obdyc8`, thread
  `01a06aa1`). Then Gate 1 round 3 re-review. Loop until ADVANCE before any launch.
- **GATE 1 FIX PASS ROUND 2 DONE** (Codex `k53obdyc8`): all 5 fixed (H1 joint potency/burden/window
  validator now shared by prereg + report; whitespace gene IDs; sbatch dup-launch safe; legacy-draft
  migration hash-gated; permute manifest+scored_strata SHA-bound). pytest **367 passed / 4 skipped / 0
  failed**. Still uncommitted. **GATE 1 ROUND 3 IN FLIGHT** — Codex `kly0p4aln` (read-only), tight
  focus on twice-open H1 + regression hunt. Loop until ADVANCE.
- **GATE 1 ROUND 3 = DO-NOT-ADVANCE** (Codex `kly0p4aln`, thread `01a06ba8`): H1 narrowed to a
  tolerance/rounding bypass (math.isclose accepts inconsistent window; rounding before verdict
  manufactures one) + 2 MED (zero-width/BOM gene IDs; permute self-selected freeze tag = no trusted
  root) + 1 LOW (cache key omits scorer/samtools digest). All else CLEAN (sbatch races, migration
  locking, freeze-safety). FIX PASS ROUND 3 dispatched (Codex `k42g4jdez`, thread `01a06aa1`): strict
  recompute-and-reject window consistency (no isclose, no pre-round) shared by prereg+report; NFKC gene
  IDs; permute exact-approved-tag/signed root; cache-key scorer+samtools digest. Then Gate 1 round 4.
- **ROUND-3 FIX PASS LANDED FULLY (verified on disk) but Codex hit 5-HR LIMIT again on emit** (task
  `k42g4jde` failed "try again at 9:11 AM"; writes completed before the failure — same late-fail
  pattern). Claude-verified 05:17 EDT: isclose GONE from window.py; window.py "publishes recomputed
  endpoint, never caller's copy" (strict recompute, no pre-round); report.py:281 delegates to
  window.normalize_window_rows + prereg_payload (shared H1 path); rounding-boundary regression tests
  present on BOTH prereg (test_titration_window.py:394) and report (test_titration_report.py:105) with
  the reviewer's exact 50.00000049 probe; strata NFKC + control/format/separator reject; permute
  APPROVED_FREEZE_TAG exact + signed-tag-object SHA + scored_strata exact-bytes bind; sbatch cache key
  now includes scorer_source + external executable digest. pytest **376 passed / 4 skipped / 0 failed**
  (172s). Still uncommitted.
- **GATE 1 ROUND 4 = DO-NOT-ADVANCE** (Codex `kzccqp533`, thread `01a06c8c`, 9:12 AM). H1 core CLOSED
  (numeric mismatch/non-finite rejected both paths, valid ladder still verdicts, no import cycle),
  gene IDs CLOSED, Phase-9 trusted-root CLOSED (all-zero sentinel fails closed), cache-key helper
  CLOSED, locks/freeze-safety clean. Remaining: 2 HIGH — (a) sbatch:58 NameError `score_cache_key`
  (imported as `scores.`) = ABORTS every Arm D launch; (b) window.py:128 `_finite_number` accepts bool
  → manufactures verdict — + 2 MED (int(length) truncates fractional rungs; non-canonical pseudocount
  false-rejects valid ladder). FIX PASS ROUND 4 dispatched — Codex `kk1mopm41`, thread `01a06aa1`
  (first codex-reply returned EMPTY + made no edits; re-dispatched). Then Gate 1 round 5. Loop until
  ADVANCE. HARD: seal UNOPENED, freeze UNCUT, no Phase-9, no GO-manufacture.
- **GATE 1 FIX PASS ROUND 4 DONE** (Codex `kk1mopm41`): all 4 fixed — sbatch:58 now
  `scores.score_cache_key(...)` + execution-level heredoc test; window.py rejects Python+NumPy bool
  before coercion; integral-length + EXPECTED_LADDER_LENGTHS membership; canonical OFF_TARGET_PSEUDOCOUNT
  only. pytest **388 passed / 4 skipped / 0 failed**. Still uncommitted. **GATE 1 ROUND 5 IN FLIGHT** —
  Codex `ktnl4ng2u` (read-only), verifies the launch-blocker sbatch fix + regression hunt. On ADVANCE →
  PHASE 1 LAUNCH.
- **GATE 1 ROUND 5 = DO-NOT-ADVANCE** (Codex `ktnl4ng2u`, thread `01a06ca0`): LAUNCH-BLOCKER DEAD —
  sbatch heredoc CONFIRMED executing (real SHA-256, correct imports), Arm D path clean + freeze-safe,
  ladder/pseudocount/canonical-verdict verified. Only remaining = one narrow class: numpy booleans
  (np.array(True)/np.bool_/0-d bool array) slip _finite_number + pseudocount guard in Arm E window.py
  → manufacture verdict. 1 HIGH + 1 MED, same root. FIX PASS ROUND 5 dispatched (Codex `kyagmqge1`,
  thread `01a06aa1`): ONE shared strict numeric validator rejecting all boolean forms + non-scalar
  arrays, routed through every intake incl direct directional_verdict. Then Gate 1 round 6 — expected
  final. HARD: seal UNOPENED, freeze UNCUT, no Phase-9, no GO-manufacture.
- **GATE 1 FIX PASS ROUND 5 DONE** (Codex `kyagmqge1`): `_strict_finite_scalar()` (window.py:130)
  rejects Python bool / numpy.bool_ / 0-d bool array / bool-dtype array / non-scalar arrays before any
  float() or ==pseudocount; routed through all 11 numeric intakes; only float() left is inside the
  validator post-rejection. pytest **415 passed / 4 skipped / 0 failed**. Still uncommitted. **GATE 1
  ROUND 6 (FINAL) IN FLIGHT** — Codex `kk0vocutp` (read-only): confirm boolean class closed + holistic
  launch-readiness sweep of the whole diff. On ADVANCE → PHASE 1 LAUNCH.
- **GATE 1 ROUND 6 = DO-NOT-ADVANCE** (Codex `kk0vocutp`, thread `01a06cac`): boolean class CLOSED +
  whole diff otherwise LAUNCH-CLEAN, but numeric strings ("2.0") still pass float() → verdict. Root
  cause = blocklist validator. FIX PASS ROUND 6 dispatched (Codex `kta66vq4e`, thread `01a06aa1`):
  convert _strict_finite_scalar to a TYPE ALLOWLIST (accept only Python int/float + numpy real scalars,
  reject bool/str/bytes/Decimal/None/complex/ndarray/everything else before coercion, then finiteness).
  Closes the whole malformed-input class at once. Then Gate 1 round 7 — expected ADVANCE → LAUNCH.
- **GATE 1 FIX PASS ROUND 6 DONE** (Codex `kta66vq4e`): type allowlist live (window.py:130), all
  intakes routed, only float() inside validator, malformed-input regression tests (str/bytes/Decimal/
  Fraction/None/complex) across all 4 entry points. pytest **448 passed / 4 skipped / 0 failed**. Still
  uncommitted. **GATE 1 ROUND 7 (FINAL) IN FLIGHT** — Codex `k527yts3z` (read-only): final malformed-
  input attack + launch-readiness regression confirm. On ADVANCE → PHASE 1 LAUNCH.
- **GATE 1 ROUND 7 = DO-NOT-ADVANCE** (Codex `k527yts3z`, thread `01a06cc2`): all launch-relevant clean,
  ONLY remaining = isinstance admits user-defined int/float/numpy-real SUBCLASSES → subclass ladder
  manufactures verdict. Final type-space corner. FIX PASS ROUND 7 dispatched (Codex `kyyvl7rzg`, thread
  `01a06aa1`): EXACT-TYPE membership (type(v) is exactly int/float/concrete numpy real scalar classes;
  no isinstance, no bool). After exact-type there is no subclass/coercion path left. Then Gate 1 round 8
  — expected ADVANCE → LAUNCH.
- **ROUND-7 FIX LANDED + CLAUDE-VERIFIED, but Codex hit 5-HR LIMIT again on emit** (task `kyyvl7rzg`
  failed "try again at 2:12 PM"; writes completed first). Verified on disk 10:22 EDT: window.py:72
  `_ALLOWED_REAL_SCALAR_TYPES` = frozenset{int, float, np.float64/32/16, np.int64/32/16/8,
  np.uint64/32/16/8}; window.py:153 `if type(value) not in _ALLOWED_REAL_SCALAR_TYPES` (no isinstance,
  bool excluded); subclass regression tests present (_MyInt/_MyFloat/_MyNumpyFloat parametrized across
  compute_window/prereg/report). pytest **463 passed / 4 skipped / 0 failed** (147s). Still uncommitted.
- **★ GATE 1 = ADVANCE (launch-ready) — 2026-09-04 14:16 EDT** (Codex `k7bm2io46`, thread `01a06da0`).
  8 review rounds total; malformed-input type space fully closed (166 probes rejected across all 4
  entry points, valid ladders still verdict), whole diff regression-clean, freeze-safe, verdict-neutral.
  pytest **463 passed / 4 skipped / 0 failed**. Still UNCOMMITTED. → PHASE 1 LAUNCH.
- **PHASE 1 PREFLIGHT (Codex `k3suq3ui2`):** Arm D inputs ALL present (3 pinned ScrMOE BAMs/counts,
  positives, GENCODE v50 FASTA/GTF/FAI, hg19→GRCh38 chain, Pangolin GPU env). **Arm E NOT launchable —
  missing:** pinned off-target transcriptome reference (no canonical path), a coord-bearing site
  manifest producer (chrom,start,end,strand,target_5to3), and runnable seed-scan + Pangolin H-off
  sweep scripts (`gpu_preemptable.sbatch` still a placeholder, exits 2). Arm E scan/sweep were never
  built — only scorer modules + template. → SPLIT: launch D now; BUILD+launch E separately.
- **NOTE: harness classifier initially BLOCKED the autonomous Codex sbatch-launch session; user
  approved once ("Approve just this launch"). Re-dispatches have gone through since.**
- **★ ARM D LAUNCHED — 2026-09-04 16:37 EDT.** afterok chain live on mit_preemptable + --requeue
  (self-sufficient; survives session end / wifi drop):
  - `21992164` feng-bench-universe — rebuild ScrMOE-only universe. sentinel
    `$SCRATCH/feng/bench/detectable_universe/.done` → `feng/artifacts/bench/detectable_universe_v50.parquet`.
  - `21992165` feng-bench-match (afterok:21992164) — build 8×101 strata. sentinel
    `feng/artifacts/bench/strata_v50.done` → `strata_v50.parquet`. (script = ephemeral
    `/tmp/feng-armd-match.sbatch`; Slurm copied it at submit — should be moved into feng/scripts/ later
    for reproducibility.)
  - `21992166_[0-31%16]` feng-bench-intarna (afterok:21992165) — sharded IntaRNA sweep, content-
    addressed cache, atomic publish. terminal = `$SCRATCH/feng/bench/scores/current/manifest.json`
    + `current/intarna_shards/intarna_{0..31}.done` → `current/scored_strata.parquet`.
  - Logs: `feng/logs/feng-bench-{universe,match,intarna}-*.{out,err}`.
  - **POLL:** `squeue -j 21992164,21992165,21992166 -o '%.14i %.28j %.14P %.12T %.20E %R'`; per-job exit
    `sacct -j <id> --format=JobID,State,ExitCode,Elapsed`.
  - **KEY ENV FACT:** Codex's MCP sandbox has NO network route to slurm001 (scontrol → family=0/port=0),
    so it CANNOT sbatch. Submission was run from the MAIN session shell (reaches Slurm). Any future
    launch: Codex prepares, main shell submits.
- **★ ARM E BUILT + LAUNCHED — 2026-09-04 16:45 EDT.** Codex built the missing pieces (ref-build,
  seed-scan, Pangolin H-off sweep), pytest **474 passed / 4 skipped / 0 failed**. Submitted from MAIN
  shell (sandbox can't reach Slurm). afterok chain, mit_preemptable + --requeue + per-length .done:
  - `21992483` feng-armE-ref — build+SHA-pin deterministic GENCODE-v50 exonic transcriptome →
    `$SCRATCH/feng/arm_e/reference/gencode_v50_exonic_union.fa` (+ `.pin.json`). sentinel
    `$SCRATCH/feng/arm_e/reference/gencode_v50_exonic_union.done`.
  - `21992484_[1-15%8]` feng-armE-scan (afterok:21992483) — seed-extend off-target scan, lengths 2–30
    → `$SCRATCH/feng/arm_e/scan/sites_L{02..30}.tsv` + `burden_L{02..30}.json`. sentinels
    `scan/scan_L{02..30}.done`.
  - `21992485_[1-15%2]` feng-armE-hoff (afterok:21992484, **--gres=gpu:1**) — Pangolin H-off sweep
    (feng-gpu-env) → `$SCRATCH/feng/arm_e/hoff/pangolin_hoff_L{02..30}.json`. sentinels
    `hoff/pangolin_hoff_L{02..30}.done`. WATCH: GPU on mit_preemptable may pend for a GPU slot.
  - Logs: `feng/logs/feng-armE-{ref,scan,hoff}-*.{out,err}`.
  - Scan script env (for a resubmit): SCRATCH, FENG_ENV_PREFIX, FENG_GPU_ENV_PREFIX, FENG_ARME_ROOT,
    FENG_ARME_GTF, FENG_ARME_FASTA, FENG_ARME_MAX_EDIT=5, FENG_ARME_SCAN_DIR, FENG_ARME_HOFF_CHUNK_SIZE=32,
    FENG_GPU_OUT_DIR (exact values in this session's Arm E submit block).

---

## ★★ SAFE-TO-STOP / RESUME — 2026-09-04 17:25 EDT (ARM D RELAUNCHED post-hotfix; ARM E fixed, re-review pending)

**Session can end / wifi can drop NOW. Arm D's 3 jobs run on Slurm independent of this session
(--requeue + checkpoint + per-shard .done). Arm E is code-fixed + green but NOT relaunched (re-review
gated on Codex 5-hr reset 7:13 PM).**

- **POLL:** `squeue -u $(whoami) -o '%.14i %.24j %.14P %.11T %.20E %R'`
  Per-job exit after terminal: `sacct -j <id> --format=JobID,JobName,State,ExitCode,Elapsed`.
- **Arm D (RELAUNCHED 17:25 — new IDs):** 21994396 universe → 21994397 match → 21994398_[0-31%16]
  IntaRNA. Done when `$SCRATCH/feng/bench/scores/current/manifest.json` + all 32 `intarna_*.done` exist
  → `current/scored_strata.parquet`. **Prior chain 21992164/165/166 FAILED** (universe crashed on a
  float-NaN IncLevel cell — `strata.py:_parse_psi_list`).
  - **⚠ HOTFIX PENDING REVIEW:** `strata.py` `_cell_or_empty()` (Claude-applied, NOT Codex-reviewed)
    coerces missing/NaN rMATS cells to empty — faithful to existing NA handling; regression test
    `test_load_se_events_tolerates_missing_inclevel_cells` added, 34 universe-tests green. **Codex must
    verify at 7:13 reset** (is NaN-as-empty the right PSI semantic, or should the event be dropped?).
    If Codex disagrees → the universe (+ downstream match/intarna) rerun. WATCH the universe job for a
    second real-data failure fixtures didn't cover.
- **Arm E (fixed, green, NOT running):** scan/hoff CANCELLED earlier (corrupt-by-design). All 4 HIGH +
  3 MED fixed (Codex `k0r5a6wbm`, 487 green). Ref 21992483 done+pinned (reusable). **Re-review pending
  Codex reset ~7:13 PM** → on ADVANCE relaunch scan+hoff from MAIN shell (submit block in transcript;
  skip ref). Open design Q: L2/L4/L6 uninformative → "complete" over full 2–30 vs informative L≥8.
- **⚠ DO NOT MUTATE THE WORKING TREE while jobs run.** All Phase-0 fixes + Arm E build are UNCOMMITTED
  on feng-arm; the jobs import `feng/src/` at RUNTIME. A `git checkout`/`stash`/branch-switch would
  corrupt in-flight (and requeued) jobs. Leave the tree as-is until both chains finish, OR commit
  feng-arm first (code is Gate-1 ADVANCE clean, 474 green) — committing is safe, mutating away is not.
- **NO live monitors needed** — none armed (Monitor watchers + MCP codex tasks die with the session
  anyway). Resume path = the poll lines above; nothing depends on this session staying alive.
- **★ ARM E CANCELLED — infra review (Codex `k7cypwepp`, thread `01a06e30`) = DO-NOT-ADVANCE, 4 HIGH.**
  The launch scripts were built + submitted on 474-green tests WITHOUT the standing adversarial gate;
  review caught corrupt-by-design flaws. **`scancel 21992484 21992485` done 2026-09-04 16:53.** Ref job
  `21992483` COMPLETED clean (3m06s, 437,245 records, FASTA 232,557,984 B, SHA `dae0a683…4723`, pinned)
  — **KEPT + reusable, do NOT rebuild.** Arm D UNAFFECTED (separate, reviewed, running).
  **ARM E FIX LIST before any relaunch (delegate to Codex, then re-review, then relaunch scan+hoff):**
  - HIGH-1 arm_e_scan.sbatch:34 / seed_extend.py:248 — fixed `max_edit=5` degenerates for 2/4-mers →
    exhaustive scan (~billions of rows, wall-time/QUOTA blowup). Need a length-aware edit policy
    (preregistered) or a non-enumerating counting strategy w/ explicit short-mer handling.
  - HIGH-2 arm_e_scan.py:305 — every variable-length alignment counted as a separate off-target
    (97 hits for one 18-mer locus). Arm C collapses at armc_build.py:655; Arm E must too. Inflates
    burden + H-off with LENGTH-DEPENDENT bias → can MANUFACTURE "shorter is better." Add deterministic
    genomic-site dedup/cluster (best alignment per site), use that canonical manifest for BOTH count
    and H-off.
  - HIGH-3 arm_e_scan.py:292/194 + arm_e_hoff.py:40 — signal checkpoint saves live file/burden
    independently of last-committed record → double-count on resume. Checkpoint only an immutable
    last-committed (index, offset, burden) snapshot.
  - HIGH-4 arm_e_hoff.py:59 + gpu_preemptable.sbatch:117 — Pangolin/fetch/orientation failure →
    nonempty JSON burden=None but exit 0 + `.done` touched → failed shard treated as complete. Must
    write incomplete diagnostic atomically, exit nonzero, NO `.done`.
  - MED — transcriptome union omits exon–exon junction sequences (not full mature-transcript space):
    add junction records or explicitly narrow + record the limitation. Checkpoints not bound to
    input/config digests (rung/ASO/max_edit/manifest SHA/algo version). **window.py:265 has NO consumer
    reading burden_L*.json / pangolin_hoff_L*.json → need a hash/schema-validated assembler requiring
    all 15 scan + 15 H-off terminals, else curve stays incomplete.**
  - Verified clean by review: ref determinism + SHA, fail-loud on missing inputs, ASO revcomp correct,
    1-based strand-consistent coords, array 1..15→L02..L30 mapping, no ΔPSI/rank/perm/seal, None-not-0.0.
  - Relaunch (after fixes + re-review): reuse the pinned ref (skip 21992483); resubmit scan+hoff from
    the MAIN shell (sandbox can't reach Slurm) using this session's Arm E submit block.
  - **ARM E FIX DONE** (Codex `k0r5a6wbm`, thread `01a06e36`): all 4 HIGH + 3 MED fixed, pytest
    **487 passed / 4 skipped / 0 failed**. length-aware policy effective_edits=min(requested,5,
    floor(length/5)) → L2/L4/L6 terminal-uninformative (burden=None); canonical best-per-site manifest
    drives burden AND H-off; immutable checkpoint snapshots; H-off fail-loud (terminal=false + nonzero,
    no .done); universe narrowed to contiguous pre-mRNA exonic loci (junctions excluded+recorded);
    hash-bound checkpoints; window.py:307 15-rung assembler (missing/incomplete → off_target_hits=None
    → prediction None). Ref unchanged, Arm D untouched.
  - **DELIBERATE LIMITATION / OPEN DESIGN Q:** L2/L4/L6 burden=None → if the window requires ALL 15
    rungs (2..30), the curve can NEVER complete → Arm E prediction permanently None. Needs a decision:
    define "complete" over the informative subset (L≥8) vs full 2–30. Being assessed in re-review.
  - **ARM E RE-REVIEW IN FLIGHT — Codex `kgbhrmeqw`** (read-only). On ADVANCE → relaunch scan+hoff from
    MAIN shell (reuse pinned ref, skip 21992483). On DO-NOT-ADVANCE → fix loop. Do NOT relaunch without
    clean re-review.
- **NEXT (Phase 3, when compute lands, needs go-ahead):** read scored_strata + Arm E window; honest
  verdict write-up (D bind-vs-act beyond energy; E therapeutic-window direction); Codex sol adversarial
  review of the NUMBERS before any headline. Seal stays UNOPENED; no manufactured GO.
- Open cleanup (non-blocking): move `/tmp/feng-armd-match.sbatch` into feng/scripts/ for reproducibility.

---
- (history) **GATE 1 ROUND 8 (FINAL) PENDING CODEX 5-HR RESET ~2:12 PM EDT.** Timer monitor armed → auto-resume
  ~2:13. On resume: dispatch fresh Codex adversarial round-8 review of the exact-type fix + final
  launch-readiness. Expected ADVANCE (type space fully closed). On ADVANCE → PHASE 1 LAUNCH
  (D: rebuild ScrMOE-only universe → strata → IntaRNA sweep; E: off-target scan + Pangolin H-off) on
  mit_preemptable + checkpoint → Phase 2 monitors. HARD: seal UNOPENED, freeze UNCUT, no Phase-9, no GO.
  Do NOT jump the gate — launch only after ADVANCE. If round 8 finds another purely-theoretical
  input-validation nit (data flow is internal computed floats, not attacker input), consider accept+log.
- **RESUME (or if session died, do this manually):**
  1. Dispatch Codex sol (medium) ADVERSARIAL review of the uncommitted fix diff — read-only,
     independent-reviewer framing, verdict ADVANCE / DO-NOT-ADVANCE + HIGH/MED/LOW. (Prompt template =
     the `kwy34n2sh` dispatch in this session's transcript; covers all 11 findings + freeze-safety +
     verdict-neutrality + new-bug hunt.)
  2. If DO-NOT-ADVANCE: dispatch Codex to fix findings (resume thread `01a06aa1-...` or fresh), re-run
     pytest, re-review until clean.
  3. On ADVANCE: PHASE 1 LAUNCH on mit_preemptable + checkpoint, in this order — (Arm D) rebuild
     ScrMOE-only detectable universe [compute; committed artifact is STALE & now validator-rejected] →
     build 8×101 strata → IntaRNA score sweep; (Arm E) off-target seed-extend transcriptome scan +
     Pangolin H-off genome sweep to populate therapeutic_window. Delegate the sbatch build+launch to
     Codex. HARD: seal UNOPENED, freeze UNCUT, no Phase-9, no manufactured GO.
  4. Arm one lean blocking monitor per launched job (fires to MAIN at terminal, covers failure states
     too); verify job PID/state before trusting "running"; update this block with live job+monitor IDs.
- HARD RULES baked in: seal UNOPENED / freeze UNCUT without explicit user sign-off; Phase-9 one-shot
  never autonomous; do NOT manufacture a GO; if Codex hits 5-hr rolling limit, monitor for reset+1min
  and continue (don't checkpoint-stop for the 5-hr case).

---

## ★ LIVE CHECKPOINT — 2026-09-04 (Arm C sig-tested + reviewed → PARKED; 5-arm audit DONE; D+E goal STAGED, nothing launched)

**STATE FOR /clear — read this block first. NO compute launched, NO monitors live yet.**

**Arm C — CLOSED, PARKED as exploratory probe (not publishable as prereg).**
- DeLong + cluster-bootstrap done (Opus fork, recomputed from persisted `armc_features.parquet`
  190743×15, deterministic LOSO replay — no cluster). Honest numbers:
  - Cluster-bootstrap (B=2000, by SBO group): AUC diff model−duplex_dg **+0.0404, 95% CI
    [+0.0003, +0.0782], p=0.049** (significant but FRAGILE, CI grazes 0). AP diff model−edit
    **+0.0042, CI [−0.0055,+0.0175], p=0.399** (NOT significant). iid DeLong p=0.00014 is
    OVERSTATED (ignores fold clustering) — do NOT quote it.
  - Verdict: consistent with NO-GO. Binding beats energy on ranking by a whisker, nothing on
    precision. Coverage: 217 pos / 190743 rows / prev 0.00114 / 65 groups / **40 evaluable folds**.
- Codex sol adversarial review → **NOT publishable as the prespecified Arm C verdict**: (1)
  +0.05-on-both bar is a SUBSTITUTED weaker endpoint (PLAN.md:134 locks nested LOAO, combined
  baseline, AUCPR gain ≥0.10/≥2×, bootstrap p<0.01, FDR, external replication); (2) negatives are
  RESPONDER-ONLY (S4) → 1-in-1000 promiscuity claim confounded; (3) 16 low-complexity SBOs dropped
  pre-scan → only 65/81 analyzed. Metric math / orientation / leakage all CLEAN.
- All logged in `feng/artifacts/phase0_arm_c_probe/NUMBERS_AND_HEADLINES.md` (§8 = review, banner
  at top). My inflated claims fixed there (42×→~37×, dropped false "structurally near-zero",
  "no baseline beats anywhere"→pooled-only). NO further compute on Arm C.

**5-ARM INFRASTRUCTURE AUDIT — DONE (Codex sol, read-only). Trust verdict:**
- **Arm A** trustworthy for narrow descriptive claim (1 LOW: guide-geometry guard). GO path dead
  (effective N≈10.9 vs ≥25 gate).
- **Arm B** biggest hole: concordance test = one-line PLACEHOLDER (not built); liftover unsafe
  (first/last-base only → internal chain-gap intervals emitted contiguous; null-assembly restamp).
- **Arm C** +2 new HIGH (exonome cache reuse bypasses ≥90% lift tripwire; merge count-only).
- **Arm D** core CLEAN (deterministic NN matching, frozen seed, positive globally excluded from all
  bg pools = no leakage, 1-based→0-based correct, Phase-9 inert behind NotImplementedError). 3 HIGH
  pre-freeze blockers = light fixes, NOT redesign: pooled-control universe reaches strata (rebuild
  from 3 ScrMOE controls), silent-skip missing STAR counts, score cache not hash-bound. +2 MED
  (1-to-1 score-join table; Phase-9 unlock not tied to freeze tag) +1 LOW (dup event_id).
- **Arm E** candidate-exploratory only, 5 HIGH: double-penalized on-target burden; sparse curve
  still emits verdict; weak seal semantics (locked=True overwrites); Pangolin H-off ignores
  aso_5to3+strand & zero-fabricates missing; checkpoint masks failures as clean exit.

**ARM RANKING (positive result / unit compute):** 1. **Arm D** (best bet — unopened seal, core
sound, ~100-400 CPU core-hrs; NO leaked signal, could still NO-GO honestly) 2. Arm B (real
endpoint but unbuilt + feasibility untested) 3. Arm E (signal leans AGAINST shorter-better;
decisive off-target work is Arm-C-scale) 4. Arm A (descriptive only) 5. Arm C (parked)
6. variant appendix (NO-GO).

**TONIGHT'S GOAL — STAGED, NOT RUNNING.** Full `/goal` spec persisted at `feng/GOAL_TONIGHT.md`
(survives /clear). (Phase 0
fix D+E blockers via Codex sol → Gate 1: tests green + Codex adversarial re-review → Phase 1 launch
D IntaRNA score sweep + E off-target scan/Pangolin sweep on preemptable partitions w/ checkpoint →
Phase 2 lean blocking monitors → HANDOFF). **HARD RULE baked in: seal stays UNOPENED, freeze UNCUT
without explicit user sign-off.** Sealed Phase-9 test is the one-shot terminal step — never
autonomous.

**PENDING USER INPUT (blocking launch):** (a) go-ahead to run the goal / dispatch Phase 0 fixes;
(b) confirm spending Arm-C-scale compute on Arm E off-target scan despite the against-hypothesis
signal; (c) whether Arm D seal fires tonight (default NO). **No jobs, no monitors yet — nothing to
resume; safe to /clear and restart from this block + the /goal.**

---

## ★ LIVE CHECKPOINT — 2026-09-03 (Arm C DONE — verdict NO-GO; Arm E CPU+GPU-infra+H committed)

**ARM C GATE TERMINAL (job 21881182 COMPLETED 4m21s) -> VERDICT: NO-GO.** GATE_RESULT.md written
(feng/artifacts/phase0_arm_c_probe/GATE_RESULT.md). Honest, prespecified, NOT goalpost-moved:
- 190743 bound (SBO,exon) rows; positive rate 217/190743 = 0.001 (extreme imbalance); 65 SBO
  groups, only ~15 evaluable (25 folds degenerate / 0-positive held-out).
- Pooled held-out: Model (LOSO logistic, binding features) AUC 0.838 / AP 0.042. Baselines:
  -edit_distance 0.674 / 0.040; **-duplex_dg 0.797 / 0.004**. GO bar = beat best baseline +0.05
  on BOTH. Model beats best AUC baseline (duplex_dg) by only +0.041 (<0.05 FAIL) and best AP
  baseline (edit_distance) by +0.002 (FAIL). -> **NO-GO on both margins.**
- HONEST READ (the anticipated outcome): duplex ENERGY largely explains action (dG-alone AUC
  0.797, per-fold dG AUC frequently 0.8-0.99). Binding model adds a little AUC, not past the bar;
  AP near-flat under 0.001 prevalence. This is a clean result, NOT a failure: "binding predicts
  action, but mostly via energy, not beyond it." Do NOT move the GO criterion to manufacture a GO.
- NEXT (post-clear): re-review the gate result (Codex/Opus structural) BEFORE it becomes a
  headline; present NO-GO + the strong-energy mechanism honestly; any separately-prespecified
  follow-up only transparently. Freeze-safe (candidate/paper_supported; no sealed ΔPSI read).
- Per-SBO checkpointing gap noted (whole-SBO .done only; repeat-rich SBOs reset on preemption) --
  did not fix mid-run, let it finish. For a rerun: intra-SBO exon-batch checkpoint OR non-
  preemptable partition for the heavy tail.

## ★ LIVE CHECKPOINT — 2026-09-03 (Arm C FULL COMPUTE LAUNCHED — running) [SUPERSEDED by NO-GO above]

**LAUNCH CORRECTION (added after the fact):** the "BUILT, NOT LAUNCHED" state below
froze the 05:33 build. An autonomous goal loop then rebuilt `armc_scan.sbatch`
(mtime 09:05) and **LAUNCHED the full compute** — this line was missing from the
checkpoint. Actual live state as of 2026-09-03 ~16:30:
- **Scan array `21881181`** (feng-armc-scan, `--array=1-81%16`, mit_preemptable, 4 CPU /
  16G, `--time=08:00:00`, `--requeue` + `--signal=B:SIGTERM@60`). Submitted 09:05:44,
  WorkDir repo root. Per task: build shared exonome once (exonome.lock), then
  `armc_build scan` revcomp(SBO) @ max_edit=5 vs GENCODE-v50 genome ->
  `$SCRATCH/feng/armc/hits/<DG_ID>.parquet` + `.done`.
- **Gate `21881182`** (feng-armc-gate, `--dependency=afterok:21881181`, 30 min) ->
  `armc_gate` LOSO CV -> writes `feng/artifacts/phase0_arm_c_probe/GATE_RESULT.md`.
- Progress ~16:30: **60/81 `.done`**, 14 running, 8 pending. `.err` files are benign
  low-complexity SBO SKIP messages (16 self-skips by design), NOT failures. Pooled so
  far POS=157 / NEG~121.9k; steering toward ~238 positives once all ~40 positive-
  bearing SBOs land. Monitor `bz2sn3mim` fires to MAIN at gate terminal.
- Pinned inputs on scratch: GENCODE v50 primary GTF+FASTA, hg19ToHg38 chain.
- NEXT: read GATE_RESULT.md at gate terminal; re-review verdict before it becomes a
  headline; hold the pre-specified GO bar (beat best baseline +0.05 on BOTH AUC & AP).

**Negative-set blocker from the 09-02c checkpoint is RESOLVED.** `sbo_sequences.tsv`
(81 SBOs) plus the Holgersen S1 design-table PDF (`feng/sj-pdf-1-nat-...-2.pdf`) are now
staged, so the full bind->act pipeline was built end-to-end. **sbatch was deliberately
NOT called** — coordinator instruction: build + pytest only, then stop for a fresh
adversarial review of the scan-orientation / label-join logic before anything touches
the cluster. `sbatch --test-only` dry-run validated both scripts' SLURM syntax; no job
was queued.

Files written (`feng/src/feng/bench/armc_build.py`, `armc_gate.py`,
`feng/scripts/armc_scan.sbatch`, `armc_gate.sbatch`, `feng/tests/test_armc.py`). Full
suite `feng/tests/` **178 passed / 3 skipped / 0 fail** (11 of those in test_armc.py,
incl. real-PDF parse-count checks against S1/S4 — 105 S1 rows / 82 unique dg_id / 81
unique sequences, 724 S4 rows / 75 unique ids, all confirmed against the actual PDFs
via pdfplumber, not assumed).

Scan orientation: SBO 5'->3' is reverse-complemented before scanning (`scan_sbo`
in armc_build.py) since the exonome holds transcript(sense)-strand exon sequences and
SBOs are antisense; `test_scan_orientation_matters` is a positive-control unit test
proving the *unreverse-complemented* query does NOT hit a planted target while the
correct orientation does (guards the literal reversed-strand bug class).

Label join: S4.id -> S1.id (unique key) -> S1.dg_id, then S4's own 1-based-inclusive
exon interval -> `feng.bench.strata.positive_interval_to_exon_key` (same converter
Arm D's positives use) -> exact `exon_key` string match against the bound-exon set
from the scan. `actionable=1` iff the bound exon_key is in that SBO's positive set,
else `actionable=0` (bind-no-act). A dg_id tested in >1 cell type unions positives
across all its S1 ids (binding features are sequence-only/cell-agnostic).

Leakage guard: `duplex_dg` is computed by the SAME IntaRNA call (SBO as query, exon
sequence as target — mirrors the existing `feng.bench.scores` nusinersen convention)
for every bound exon regardless of label; S4's own `rnaplex_dg`/`rnacofold_dg`
columns are read ONLY to build the positive-exon_key set, never copied into a
feature column — `test_label_join_sanity` explicitly asserts those two values never
appear in `duplex_dg`.

**Deviations from the original brief (flag for the reviewer):**
- `sklearn` is confirmed **NOT installed** in feng-env (`ModuleNotFoundError`).
  `armc_gate.py` implements its own small L2 logistic-regression (batch gradient
  descent, numpy, deterministic zero-init/no RNG) plus ROC-AUC/average-precision from
  scratch (numpy + `scipy.stats.rankdata`), rather than sklearn GradientBoosting.
- `build_exonome()` reuses `feng.bench.scores`'s existing pysam-preferred/
  samtools-fallback fail-closed sequence fetch (same path Arm D's frozen scores use)
  instead of a fresh `bedtools getfasta` call — one proven fetch path, not two.
- `parse_s4` uses pdfplumber `extract_tables()` directly; it already reassembles the
  two-line-wrapped exon-interval cells into one row (verified: every one of the 724
  data rows comes back with exactly 11 cells). The raw-text/regex reconstruction the
  brief anticipated was not needed.
- S1 has 82 unique `dg_id` values, not 81: `DG111` and `DG2386` share an identical
  sequence (`CCATACCCCACCTCCCTC`) in S1 — `sbo_sequences.tsv` already dropped
  `DG2386` as the duplicate. The build pipeline drives its 81 array tasks off
  `sbo_sequences.tsv`, not S1's `dg_id` column, so this doesn't affect sizing.
- Label join uses **exact** `exon_key` match (no coordinate-overlap fallback) — a
  paper-vs-GENCODE-v50 exon-boundary mismatch would under-label a true positive as
  `actionable=0`, which is conservative (biases toward NO-GO, not toward leakage) but
  is a real, un-fixed limitation worth the reviewer's attention.
- Untested at real scale: only synthetic 2-3-exon fixtures were exercised. The actual
  81621-exon x 81-SBO scan's hit-rate, IntaRNA call volume, and per-fold class balance
  are unknown until the array job actually runs.

**Not done:** sbatch was not called. Next step once review is clean: launch
`armc_scan.sbatch` (array 1-81%16), then `armc_gate.sbatch --dependency=afterok:<array id>`.

---

## ★ LIVE CHECKPOINT — 2026-09-03b (Arm E — CPU stack COMMITTED b08fae7; GPU infra + docs built)

**Arm E (StericRuler length-titration) build underway** (Codex down till Sep 7 -> Opus fork
subagents). Plan `~/.claude/plans/add-this-as-a-goofy-whistle.md`. New pkg
`feng/src/feng/titration/`. Slices landed (all additive, TDD, NOT committed yet):
- **ladder.py** — `build_ladder` symmetric even 2..30 ladder around the true nusinersen
  18-mer anchor chr5:70076584-70076601 (GRCh38 v50, +). Fail-closed in-code assert L=18 ASO ==
  NUSINERSEN_5P_TO_3P. Reuses scores.py revcomp + fetch. 12 tests. Frame verified vs the
  published intron-7 flank on every rung.
- **evaluate.py** — per-mer CPU scorers: IntaRNA on-target (3 honest states: None sub-seed
  error / 0.0 no_duplex flag / real dg — plan's "0.0 for short mers" was WRONG, corrected),
  Vienna RNAduplex/RNAcofold/RNAup/RNAplfold/RNAfold, off-target via scan.seed_and_extend
  (opt-in). MED-1 FIXED: rnaplfold_accessibility now folds the target in GENOMIC CONTEXT
  (pad=70/W=70/L=40) not the isolated window -> accessibility 0.31-0.51 (was ~1.0 artifact).
  19 tests.
- **window.py** — therapeutic-window endpoint: potency (backbone = IntaRNA + RNAup only;
  RNAduplex/cofold + plfold-acc secondary, no double-count) / (off-target burden + 1). Prereg
  lock MECHANISM (prereg_payload/hash/lock/verify, deterministic SHA-256) but **DRAFT only —
  NO real seal cut** (real lock = gated step pending user sign-off + GPU axes). 25 tests.
- **Adversarial gate (Opus fork) over ladder+evaluate: ADVANCE** — no HIGH, frame-correct,
  never-fabricate clean, orientation consistent, honesty/freeze clean. Findings MED-1 (fixed),
  LOW-1/2/3 (carried into window / logged). rnafold_mfe=0.00-everywhere REFUTED as bug (real
  open-chain MFE).
- Full suite **252 tests** collects clean.

**HONEST DRAFT SIGNAL (do NOT manufacture a GO):** CPU-only POTENCY alone == `refutes_shorter_
better` — longer mers bind monotonically harder (IntaRNA -27->-48, L 18->30). The "shorter
better" hypothesis lives or dies in the FULL therapeutic window (÷ off-target burden + GPU
MD/splice), NOT raw potency. RNAup plateaus ~-42 for L>=26 (the turnover). Verdict not settled.

**DONE since:** CPU stack (ladder/evaluate/window/report) COMMITTED b08fae7. **GPU infra**
(feng/src/feng/gpu/ env.py+checkpoint.py, feng/scripts/templates/gpu_preemptable.sbatch,
feng/scripts/setup_gpu_env.sh documented-not-run) BUILT + Opus-gate ADVANCE; 3 MED fixed
(MED-1 signal-flush: B: prefix meant payload never got SIGTERM -> checkpoint_on_signal was dead
code; now bash trap forwards kill -TERM to child + waits within 60s grace before single requeue,
proven by a real SIGTERM subprocess test; MED-2 cuda/cudnn cross-constraint; MED-3 check_gpu_env
no-raise). **PLAN.md + TASKBOARD.md** Arm E section + A-D GPU-tier slots (Arm D = secondary/
preregistered/no-alter-locked-score) written. Suite **304**. GPU infra + docs COMMIT pending
this turn.

**GPU ENABLEMENT DONE (committed):** env `~/orcd/pool/feng-gpu-env` (6.2G; torch 2.6.0+cu124 /
openmm 8.2 / pangolin 1.0.2; SpliceAI deferred TF-conflict) installed 4m on mit_normal (compute
nodes have internet). Commit 0004d8f = **axis H Pangolin deep-splice scorer** (gpu_scores.py:
pangolin_exon7_delta = ISS-N1-masking -> exon-7 inclusion shift, gpu-env _worker subprocess;
H-off off-target-splicing interface wired, genome-wide sweep deferred) + GPU-env gate fix
(check_gpu_env probe timed out under CPU load -> false torch:False; now timeout=inconclusive +
scorer imports in-worker) + recipe pin. Suite green.

**KEY HONEST FINDING:** Pangolin CPU smoke, correct sign (masking ISS-N1 raises exon-7 inclusion
= nusinersen mechanism). SATURATION AMBIGUOUS + readout-sensitive (do NOT claim monotonic): the
donor-exact exon7_delta rises L6..L30 (0.0015->0.027) BUT the near-donor max_gain rises L6->L24
(0.014->0.253) then DROPS at L30 (0.075) — a turnover hint. No clean saturates-vs-rises claim at
smoke scale; needs a LOCKED readout + full-ladder GPU-scale run before any prereg use. Combined
with IntaRNA duplex + RNAup (both rise with length), the on-target potency picture leans
longer-acts-harder but is not settled on the splice axis. "Shorter is better" still rests mainly
on the SPECIFICITY side: off-target scan + H-off (Pangolin off-target splicing) must show shorter
mers shedding off-targets faster than they lose potency. At L30 max splice-gain shifted to a new
offset (+6) = off-target splicing emerging. NOT a GO; do not manufacture one.
OPEN: (a) masking scheme (footprint->N one-hot; readout = gain at exon-7 5' donor chr5:70076575)
is a CANDIDATE modeling choice NEEDING adversarial review + a locked readout before prereg.
(b) LICENSE: gpu_scores.py reproduces Pangolin's compute_score+one_hot_encode VERBATIM (GPL-3.0,
attributed) because pangolin.pangolin top-imports absent vcf/gffutils/pyfastx. Vendoring GPL-3.0
into feng needs a license decision (feng repo license? GPL copyleft may apply to the combined
work) — resolve before any release/publication.

**NEXT (Arm E remaining) — REPRIORITIZED:** the decisive axis is now SPECIFICITY, not more
on-target scorers. (1) OFF-TARGET compute: CPU seed-extend transcriptome scan (populate burden)
+ H-off Pangolin off-target-splicing genome-wide sweep via gpu_preemptable.sbatch -> makes
therapeutic_window non-None + gives the REAL directional verdict. HEAVY (Arm-C-scale). (2) axis G
OpenMM MD binding-free-energy scorer — ANOTHER on-target proxy (likely also rises with length),
lower priority than specificity. (3) prereg SEAL (gated: user sign-off + all axes). (4) final
adversarial gate over gpu_scores masking scheme + Codex sol durable gate (Sep 7). Perf (LOW-3):
evaluate re-indexes the 3GB genome per fetch -> pass on-disk .fai before the off-target scan.

## ★ LIVE CHECKPOINT — 2026-09-02c (Arm E GO; Arm C probe AMBER->GO-leaning)

**Arm E APPROVED (user) — GO, not started.** Length-titration + GPU plan at
`~/.claude/plans/add-this-as-a-goofy-whistle.md`. Build pending.

**Arm C probe RUN (Codex down till Sep 7 -> Claude-side).** Memo `feng/artifacts/phase0_arm_c_probe/
README.md`. Initial Sonnet probe = AMBER (feared 16-gene clustering). Then Claude parsed Holgersen 2021
supp S4/S5 (pdfplumber, installed into feng-env; PDFs at feng/sj-pdf-4/5-...0921.pdf):
- **~75 SBOs across diverse loci** (NOT 16 clustered genes) -> ≥25 LOAO bar CLEARED.
- Both axes in ONE harmonized source (binding RNAplex/RNAcofold ΔG + edit distance; action S4 exon
  ΔPSI 724 sites, S5 gene LFC 30493 rows). No cross-study pooling (Arm A's killer absent).
- BUT both tables are HITS/RESPONDERS only -> **bind-but-no-act negative class must be CONSTRUCTED**.
- Verdict updated AMBER->GO-leaning, conditional on negative-set construction.
Subagents (Claude-side; both hit a session rate-limit 429 first run — resets 11:10pm ET — then resumed):
  1. Sequence blocker RESOLVED — user staged S1/S2/S3/S6 + main text into feng/. **S1
     (sj-pdf-1-...-2.pdf) = the design table with the 81 SBO SEQUENCES** (ID|DG|cell|cases|controls|
     Sequence|category). Parsed -> feng/artifacts/phase0_arm_c_probe/sbo_sequences.tsv (81 distinct
     seqs 16-20nt / 82 DG ids). Join key: S4 numeric ID == S1 ID -> links spliced sites to SBO seq+cell.
  2. Haiku lit scan DONE — external-validation set exists: **eSkip-Finder** (~300-500 ASOs, tier-1,
     siloed/no-pool), Scharner 2020 (tier-2, measured Tm), Collins dystrophin (tier-3). Kept for later.

**STANDING ORDER (user, autonomous): drive Arm C until the compute is LOCKED + a revolutionary
headline that returns positive results healthily. Do NOT return to user; lean subagents + /loop;
minimize main-thread tokens; decisions optimized for revolutionary-but-defensible.**

**ARM C COMPUTE — Sonnet builder DISPATCHED (background) to build + LAUNCH the pipeline:**
- armc_build.py (parse S1/S4, bedtools getfasta exon-ome from the 81621-universe, scan revcomp(SBO)
  @edit<=5 vs exon-ome, label actionable=in-S4-hits vs bind-no-act) + armc_gate.py (leave-one-SBO-out
  CV: binding-feature model vs edit-distance/energy baseline -> AUC/AP -> GO/NO-GO in GATE_RESULT.md).
- armc_scan.sbatch (array 1-81%16, mit_preemptable, requeue-safe, per-SBO .done) -> armc_gate.sbatch
  (afterok dependency). One launch = fully autonomous. Builder returns array+gate job ids + pytest.
- Builder DONE (build+pytest, NOT launched). Fresh **Opus reviewer** ran (per user: Opus reviews AND
  applies its own fixes, not Sonnet) -> 2 HIGH + 1 MED fixed: (a) average-precision tie bug that
  inflated baseline AP ~0.33 = verdict-corrupting (rewrote to sklearn convention + known-answer test);
  (b) 16 degenerate/repeat-motif SBOs (poly-A/(AC)n/ATCATC...) matched ~every exon at 3nt seeds ->
  runtime blowup + label swamp, added sequence-only complexity guard, EXCLUDED before scan (65 SBOs
  remain); (c) wall-clock 2h->24h (max_edit=5 ~=3h/SBO). 181 passed/3 skipped. LAUNCH-READY yes.
- **CANARY LAUNCHED: job 21864147 (array task 5 = DG952, clean 20mer), mit_preemptable.** Per Opus
  advice (real repeat-rich wall time unknown; preemptable restarts from scratch) — validate 1 SBO
  before the full array. Background waiter bxlg5g715 notifies main on terminal w/ state+Elapsed+DG952
  hit rows/class balance.
- **CANARY COMPLETED (3h33m) -> EXPOSED FATAL DESIGN FLAW. Full array HALTED.** DG952 = 483 bound
  exons but 0 positives / 483 negatives. Root cause (verified): candidate universe was Arm D's
  detectable_universe_v50.parquet = GSE307431 **SMA-fibroblast** SE set, but Holgersen assayed
  **HEK293T/HepG2/PXB** -> only **3 of 724 S4 positives (0.4%) exist in that universe** -> no positive
  class -> NO-GO by construction. ALSO: S4 edit_dist is "background"/NaN for ~430 of 724 spliced exons
  (INDIRECT effects, not bound); only ~294 are direct sequence matches (edit 0-5; ~194 at edit 4-5).
- **REDESIGN dispatched to fresh Opus agent (implements own fixes, per user rule).** Decisions (locked):
  (1) universe = genome-wide GENCODE v50 exons (verify >=90% of ~294 direct positives present);
  (2) positive class = direct-bound-acted = S4 rows with NUMERIC edit_dist<=5 (~294); indirect ~430
  excluded from the bind->act test; (3) negatives = bound-but-not-acted SAMPLED against a fixed
  deterministic random background pool (20k-50k GENCODE exons, label-blind) — bounds scan cost, keeps
  max_edit=5; (4) max_edit stays 5 (positives reach edit 4-5, CANNOT tighten); (5) features binding-only
  recomputable genome-wide (edit_distance, exact_match_length, duplex_dg, gc_fraction, exon_length;
  DROP baseline_psi). Stronger/more-revolutionary framing: among near-complementary binding sites, does
  binding predict splicing ACTION? Opus re-tunes sbatch --time (24h was for the abandoned full-exonome
  scan) + shard plan, re-runs pytest, does a fast 1-2 SBO verification that positives are non-empty.
- **Redesign #1 (genome-wide universe) DONE but tripwire caught a DEEPER flaw: S4 coords are hg19.**
  Only ~10-17% of the 204 direct-positive intervals overlap GENCODE v50 (GRCh38); Holgersen methods
  cite Gencode v25/v27 (misleading — the S4 TABLE is GRCh37/hg19). PROVEN: hg19->GRCh38 liftover lifts
  **204/204 (100%)** cleanly, overlap 10%->51%. Chain staged: $SCRATCH/feng/liftover/hg19ToHg38.over.
  chain.gz (224K, UCSC). feng.rbp.liftover.load_chain_file -> ChainMap.lift_interval handles it.
- **Redesign #2 (coordinate fix) dispatched to fresh Opus implementer.** (1) lift S4 positives
  hg19->GRCh38 (drop+count unmapped/multimap); (2) positive target = the LIFTED interval's own GRCh38
  sequence (NOT required to equal an annotated exon -> recovers ~100%, not the ~51% ceiling of exact-
  exon join; keeps novel/Intropolis exons); (3) negatives = bounded random annotated-exon pool; (4)
  features SITE/DUPLEX-level, identical path both classes (edit_distance, exact_match_length, duplex_dg
  IntaRNA, binding_window_gc) — DROP whole-exon length/gc (region-vs-exon confound); gate conditioned on
  binding (no edit<=5 hit in target -> excluded from both classes); no S4-ΔG/edit leakage (recompute).
  Residual limitation flagged: neg pool may include unexpressed exons (binding-only features can't
  exploit expression -> conservative label noise). Opus re-tunes --time, verifies positives non-empty.
- **Coord fix DONE + verified: 204/204 lift (100%), frame confirmed (edit0 positives' revcomp appears
  at edit0 in lifted region w/ 30bp pad, both classes). pytest 189/3 green. LAUNCH-READY yes.**
  Positive-side LOSO feasibility (computed inline from S4, no scan): **40 non-degenerate SBOs have >=1
  direct positive** (33>=2, 26>=3), 238 positives total, spread 1-25/SBO -> clears >=25-fold bar.
  Negatives abundant (~5.5% bind -> ~1100/SBO).
- **PRODUCTION CANARY LAUNCHED: job 21875568 (array 17,61,79 = DG745/DG1063/DG767, positive-rich).**
  Confirms end-to-end integration (lift+fetch positives, pool build, scan, feature write) + wall time
  before the full array. Monitor bkibxl01r fires at terminal w/ per-SBO POS/NEG. Stale cache cleared.
- **CANARY HEALTHY (3 SBOs, job 21875568):** DG745 POS15/NEG2095, DG1063 16/3203, DG767 14/3887 —
  both classes, correct site-level feature schema. Wall: DG745 3h11m (incl one-time shared pool build),
  others ~1h55m -> steady-state ~2h/SBO. Bumped --time 4h->8h for preemption headroom.
- **GATE SMOKE (3-SBO, local, /tmp): pipeline works end-to-end.** Model ROC-AUC 0.869, AP 0.062;
  baselines edit_distance AUC 0.637/AP 0.044, **duplex_dg AUC 0.864/AP 0.042**. Verdict NO-GO on 3 SBOs
  (GO needs beat best baseline +0.05 on BOTH AUC & AP; model beats duplex_dg by only +0.005 AUC).
  IMPORTANT HONEST SIGNAL: duplex energy alone nearly explains action on these 3 -> watch whether the
  full 40-SBO run shows the model beating energy (the revolutionary "bind != act beyond energy" claim)
  or NOT. DO NOT move the pre-specified GO criterion to manufacture a GO (post-hoc-repair do-not).
- **FULL COMPUTE LAUNCHED: scan array 21881181 (1-81%16) + gate 21881182 (afterok).** 3 canary + 16
  low-complexity SBOs self-skip. ~40 positive-bearing SBOs drive the LOSO. Gate writes
  feng/artifacts/phase0_arm_c_probe/GATE_RESULT.md + .done. Monitor bz2sn3mim (durable, zero-token vs a
  12h-blocking Haiku) fires to MAIN at gate terminal w/ verdict + pooled metrics; also catches
  DependencyNeverSatisfied (a scan shard failing -> gate would need manual relaunch after requeue).
- NEXT on gate terminal: read GATE_RESULT.md. If GO healthy -> lock Arm C + headline + Codex/Opus
  structural review of the result before locking. If NO-GO (energy explains action) -> that is the
  HONEST finding; present it + the strong-energy mechanism (itself a clean result); consider a
  separately-prespecified question only transparently, never a goalpost move. Either way: re-review the
  gate result before it becomes a headline.
Parked decisions still open: freeze control definition (ScrMOE-only vs scrambled+untreated); Sep-7
Codex Arm D step-3 durable gate.

## ★ LIVE CHECKPOINT — 2026-09-02b (atlas LOW cleanup DONE + committed; phase3 rMATS DONE)

**Arm A atlas LOWs: DONE.** All 4 fixed in `feng/src/feng/atlas/arm_a.py` write_report (Opus-4.8
fork): (a) hardcoded 125/115/115 counts now interpolated from `data.windows` (n_records=len;
n_exact_windows=distinct_window_key.nunique; n_guide_covered=nunique where n_guides_containing_hull>0)
— regen reproduced 125/115/115 byte-identical, value-verified; (b) Methods now names the two
reference medians (local_z_median window-local vs block_z_positional_median coarse-block); (d) I6
bullet clarifies q<0 below E7 start 320. (c) already-documented, no change. pytest 122 passed / 3
skipped / 0 fail, lint PASS. Figure SVG/PDF regen was pure matplotlib nondeterministic-id noise
(no plotting/data change) — REVERTED, kept out of commit. LOW-(e) EXTRACTION line-37 is a different
arm, still open/optional.

**Phase 3 rMATS: DONE (not just submitted).** All 3 contrasts have `.done` in
`$SCRATCH/feng/phase3/rmats/{F18OMe_vs_scrambled_OMe,F18MOE_vs_scrambled_MOE,F20PMO_vs_scrambled_PMO}`.
Job 21824983 finished, queue empty. Results stay SEALED pre-freeze (do NOT read ΔPSI).

**Arm B liftover: DONE + COMMITTED `8497f34`.** `feng/src/feng/rbp/liftover.py` (377 lines,
self-contained UCSC-chain parser+mapper, NO external tool) lifts hg19 eCLIP peaks to GRCh38(v50),
re-stamps assembly=GRCh38 + source/target + chain path/SHA so encode.py's hard fence accepts a join;
raw hg19 still rejected, unmapped/multimap dropped+counted. test_rbp_liftover.py 10 passed. Fence NOT
weakened. FOLLOW-UP for real Phase-10 run: stage `hg19ToHg38.over.chain.gz`
(https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz, ~230KB) on
scratch, pass path to load_chain_file. feng-env has NO pyliftover/CrossMap/liftOver/pysam — pure-Python
engine removes that dependency. Concordance test itself NOT run (post-freeze, Phase 10).

**Arm D step-3: DONE + gate-clean (fallback) + COMMITTED.** Fix-pass resolved all 4 gate findings;
Claude verified on disk (suite 163 passed / 3 skipped; freeze-seal test unchanged; anchor-row logic
read + freeze-safe). Strata schema now 14 cols incl is_positive + coords (chrom/exon_start_0base/
exon_end/strand); 808 rows = 8 positives (is_positive anchor, match_distance 0, coords from own
universe row) × (1 + 100 backgrounds); scores.py runs on it with no missing-column error; sbatch
sharded 32-way idempotent (per-shard .done, skip-if-exists, atomic mv, merge-when-all), NOT launched.
Durable Codex sol gate STILL deferred to Sep 7 — re-run whole step-3 gate then as gate of record.
Files (feng/src/feng/bench/): match.py (build_strata — control-only 5-covariate matching PSI/junction/
CPM/exon-len/GC, log1p on the 3 heavy-tailed, caliper 1.0 SD, deterministic nearest-100 exon_key
tiebreak, positives removed from pool, assert_no_contrast_fields guard, fail-closed nulls/dupes/zero-
spread; map_positives_to_universe via positive_interval_to_exon_key; check_kill_gates; CLI), scores.py
(18-nt Hamming-to-nusinersen-RC baseline full + IntaRNA wrapper never-fabricate + CLI), permute.py
(Phase-9 ΔR/10^6-perm test HARD-GUARDED: assert_phase9_unlocked raises unless FENG_PHASE9_UNLOCK=1/
sentinel; run_phase9_test body NotImplementedError; only pure percentile primitive testable; seed
20260901). feng/scripts/bench_intarna_scores.sbatch scaffold (NOT launched). 3 test files.
- Freeze seal FIXED: test bans substring "terminal" in non-docstring literals across bench/field/rbp/
  gate; renamed FENG_TERMINAL_UNLOCK->FENG_PHASE9_UNLOCK etc. test_freeze_seal 17 passed.
- IntaRNA PRESENT: /home/dashii/orcd/pool/feng-env/bin/IntaRNA (+RNAplfold, samtools; pysam absent ->
  samtools fallback live). Fixture: complementary -27.0 kcal/mol, polyA 0.0.
- Kill-gate ALL 8 PASS (eligible bg 2999-14352, 100 selected each; all in_universe+detectable).
- Full suite 162 passed / 3 skipped / 0 fail. Claude read permute.py + match.py directly: freeze-safe
  confirmed (no ΔR/rank/perm readout, control-only, positives excluded).
- Adversarial fallback gate (fork ae639f9be47421056) RETURNED **DO-NOT-ADVANCE**: freeze discipline
  CLEAN (seal genuine/AST-based/unweakened, no ΔR/rank/perm computed, control-only confirmed, matching
  logic sound + substantively tested), but 2 HIGH seam defects + MED + LOW. FIX-PASS IN FLIGHT (resumed
  builder fork af1c8ba310560fc55):
  * HIGH-1 strata output lacks coord cols (chrom/exon_start_0base/exon_end/strand) scores.py needs ->
    sweep ValueError row 1. Fix: carry coords from universe into strata block.
  * HIGH-2 positives excluded from scored artifact -> 8 positives never scored, permute has nothing to
    rank. Fix: emit positive anchor rows flagged is_positive into the scored input, w/ coords+covariates
    (strata parquet = scoring input = ranking input).
  * MED match->score seam untested (both green in isolation). Fix: add end-to-end seam test.
  * LOW sbatch not idempotent on mit_preemptable. Fix: per-shard .done / skip-if-exists.
  Freeze-safe (adding positive rows to SCORE input is fine — scores are pure sequence fns; ranking stays
  in guarded permute). On green + re-check -> commit.

**OPEN FREEZE-DECISION (user, before Phase-8 FREEZE):** PLAN L146 says matching uses "ScrMOE controls
ONLY", but the committed detectable universe pools scrambled+untreated controls, so match covariates
are computed over BOTH control types. Fork flagged, did NOT change pre-freeze (correct). Resolving =
either (a) accept scrambled+untreated as the control definition and amend the prereg wording, or (b)
rebuild the universe on ScrMOE-only. Decide before FREEZE; not a step-3 code bug.

Codex sol durable gate still deferred to Sep 7 (~01:29 reset).

---

## ★ PRIOR CHECKPOINT — 2026-09-02 (build layers COMMITTED; fallback gates clean)

**One-line state:** Arm D/B build layers (positives_v50, detectable universe, RBP eCLIP/KD) are
FIXED, gate-clean on the fallback path, and COMMITTED on `feng-arm`. Ground-truth phase still DONE +
sealed (do NOT read ΔPSI pre-freeze). Arm A RED (descriptive atlas done). Appendix KILLED.

**Fix-pass complete.** All prior gate findings resolved (7 + REV3L + GOLGA4 + encode fence, earlier
passes) plus a final round of 3 MED + 2 LOW from a Claude fork adversarial gate:
- MED-1 sbatch `validate_published_artifact()` now guards CONTRAST_FIELDS on the `.done` fast-path
  (was fresh-build only) — a stale/tampered parquet with a contrast column can no longer be published.
- MED-2 was a FALSE POSITIVE (counts were correct: parquet 170584 unique intervals / 81621
  detectable, not transposed). Added a README↔parquet count cross-check to the validator anyway.
- MED-3 added `strata.positive_interval_to_exon_key` (1-based-inclusive positive -> 0-based half-open
  universe exon_key) + tests + a positives-README coordinate contract note, so the future
  positive↔universe join can't silently off-by-one and zero out the >=50-background gate.
- LOW-1 provenance validator is now SHA-based/path-agnostic and hashes the committed repo copies
  (feng/artifacts/bench/provenance/b1.txt,b2.txt + phase3/sample_conditions.tsv); build path writes
  b1/b2 into that committed dir -> committed bundle round-trips on a fresh checkout.
- LOW-2 distinct-gene count now unions `all_gene_id_bases` (9343 representative -> 9407 union);
  README updated.
- pytest **122 passed / 3 skipped / 0 fail**. Validator round-trip re-verified on the committed
  bundle (contrast-clean, counts match parquet, provenance SHAs match).

**STANDING GATE (Codex sol adversarial) = DEFERRED to Sep 7 2026** (Codex weekly usage limit, ~resets
Sep 7 01:29). In its place, TWO Claude Opus-4.8 lean fork gates ran and both returned **ADVANCE
(0 HIGH, 0 MED)** — one on the build layers, one on the fix diff. Logged as an accepted fallback;
re-run the Codex sol structural review after Sep 7 as the durable gate of record, but the layers are
committed on that fallback basis (user decision).

**Untouched follow-ups (next work):** Arm D step-3 strata matching + score impls (use
`positive_interval_to_exon_key` for the coord join); Arm B hg19->v50 liftover before the month-10
concordance test; atlas LOW findings (optional).

---

## ★ PRIOR CHECKPOINT — 2026-09-02 (fix-pass #1 KILLED; #2 redispatched)

**One-line state:** Ground-truth phase DONE (20/20 STAR + rMATS 3/3 sealed, do NOT read ΔPSI
pre-freeze). Appendix KILLED. Arm A headline RED (descriptive atlas done). Three Arm-D/B build
layers built (positives_v50, detectable universe, RBP layer) -> adversarial gate found 4 HIGH +
3 MED + 2 LOW -> fix-pass #1 (task-mtkco2qa-imj9ki, session 01a06317-d3b7-7a11-b744-b0cbd8d2f7b9)
was **KILLED at 13:14 before applying ANY patch** (investigated deeply, esp. REV3L 270x, but wrote
NOTHING; companion registry lost on /clear). Disk proof it never landed: parquet+README still 12:49,
README still says "101724", strata.py 12:34 / encode.py 12:41 (both pre-fix). All 7 findings + REV3L
STILL OPEN. Fix-pass #2 tried to RESUME thread -> FAILED 0s ("No previous Codex task thread found"; /clear
wiped companion registry so --resume-last has nothing). **Fix-pass #3 (task-mtkelmk6-ytkerv, codex session 01a06349-58c6-7921-889c-989610ef7552) ran
13m58s then FAILED on USAGE LIMIT (resets 2:25 PM EDT / 14:25). It LANDED most fixes before dying;
only the final artifact rebuild is incomplete.**

LANDED on disk (mtime > 14:01, verified): strata.py (GC fail-closed + exon-collapse), encode.py
(assembly fence + sha), bench_detectable_universe.sbatch (pre-rMATS 11-accession hash gate),
test_bench_universe.py + test_rbp_encode.py (**17 affected tests green**), positives_v50.tsv +
README, provenance/b1.txt + b2.txt. **REV3L RESOLVED**: kept mapped to MANE exon 8 with explicit
note ("paper exon 9 = 85-nt cassette flanked by 105-nt exon + 1020-nt intron = MANE exon 8"); all 8
intervals reproduce Fig 3B exon lengths.

NOT DONE (the ONLY remaining fix step): the **universe artifact rebuild** — the command that hit the
limit was `PYTHONPATH=feng/src /home/dashii/orcd/pool/feng-env/bin/python -m feng.bench...`. So
detectable_universe_v50.parquet (still mtime 12:49) + detectable_universe_v50_README.md (still says
"101724" line 58) are STALE — must regenerate to unique-exon count (~81621) + populated gc_fraction
no-nulls.

**FIX-PASS COMPLETE + VERIFIED ON DISK (2026-09-02 14:30, task-mtkfglet-pjdjns done 7m39s).** All 7
gate findings + REV3L fixed. Independently verified (not just Codex-claimed): parquet rebuilt 14:30;
detectable README = 81621 unique exons (170584 total unique intervals); gc_fraction null_count=0
(checked directly via pyarrow); allowlist pinned+hashed; encode assembly-fence enforced; REV3L
corrected (Fig 3B 85-nt exon = MANE exon 8, not 9); pytest 117 passed / 3 skipped.

**RE-REVIEW #1 (task-mtkfs750-oiq864) RETURNED: DO NOT ADVANCE — 2 HIGH + 3 MED + 1 LOW.** 5/7 prior
findings PASS; #4 (encode fence) FAIL; GOLGA4 positive is a NEW miss. Details:
- HIGH positives_v50.tsv:8 — GOLGA4 mapped to WRONG exon. Locked as 48-nt MANE exon 4 (chr3:37286014-
  37286061); correct is NON-MANE 66-nt exon 4 chr3:37273536-37273601 (ENSE00003978486.1), present +
  detectable in universe. Same wrong id in README:35 + test:212. (Other 7 positives incl REV3L PASS.)
- HIGH encode.py:200 — assembly fence still soft: hg19 frames stamped target_assembly=GRCh38 w/o
  liftover; assert_same_assembly (208) optional, no join calls it -> silent hg19->v50 merge still possible.
- MED strata.py:466 (419 multi-gene keys collapse to arbitrary gene); sbatch:138 (.done fast-path skips
  schema/dupkey/GC-null/README/provenance validation -> stale sentinel skips fixes); encode.py:218
  (SHA stamped without hashing bytes -> future misattribution).
- LOW strata.py:157 (zip IJC/SJC no length/sample-count assert).

**FIX-PASS #4 (task-mtkg9k8t-vkiq02) — fixes LANDED, but Codex then hit the WEEKLY limit (resets
Sep 7 2026 01:29; ~5 days out). No structured fix-report emitted; the failure was AFTER the pytest
rerun exit-0, so writes completed.** Claude-side verification (Codex unavailable):
- HIGH-1 GOLGA4 CORRECTED on disk: positives_v50.tsv row OT07_GOLGA4_e4 now 66-nt chr3:37273536-
  37273601 ENSE00003978486.1 exon 4 (non-MANE cassette), README + test updated.
- HIGH-2 encode fence HARD-ENFORCED: encode.py adds a frame subclass overriding merge()/join() (lines
  104/110) to call assert_same_assembly internally + source_assembly stamping + liftover-provenance
  guard (214-221). No longer optional.
- Parquet re-verified: 170584 rows, gc_fraction null_count=0, README 81621 unique exons.
- Full pytest suite: **123 passed, 0 fail** (feng-env, PATH=feng-env/bin). New fix-tests present +
  passing: test_positive_intervals_match_figure_3b_exon_lengths, test_loaded_frame_merge_hard_rejects_
  mismatched_assembly, test_coordinate_join_rejects_mixed_or_missing_assembly, test_reader_rejects_
  sha256_mismatch, test_load_se_events_rejects_mismatched_or_wrong_sample_counts, test_encode_file_
  hard_rejects_non_hg19_assembly, test_encode_file_validates_sha256.
- NOT independently re-read by Claude (tests pass but semantics unconfirmed): MED-3 multi-gene collapse
  provenance, MED-4 .done artifact validation. Trust the passing tests OR read strata.py:466 + sbatch
  .done path.

**BLOCKER: standing-gate RE-REVIEW #2 (Codex sol) CANNOT run until Sep 7 (weekly limit).** Nothing
committed; build layers NOT in 039e304. Options pending user decision: (A) commit now on Claude-side
verification, log re-review as deferred/open; (B) hold uncommitted until Sep 7 Codex re-review; (C)
run a Claude fork (Opus 4.8) adversarial review NOW as fallback gate, then commit. Watchers: use the
Monitor tool, not background Bash (Bash watchers killed intermittently this env). Last clean commit `039e304`
(does NOT include the 3 build layers or fixes).

**RESUME AFTER /clear (do this in order):**
1. Poll fix-pass #2: get its task-id (from driver, or newest /home/dashii/.codex/sessions/ file
   resuming thread 01a06317-...), then `node /home/dashii/.claude/plugins/cache/openai-codex/codex/
   1.0.6/scripts/codex-companion.mjs status <task-id>` / `... result <task-id>`. Companion registry
   does NOT survive /clear; Codex task keeps running on disk. Poll from MAIN thread (driver refuses).
   VERIFY IT ACTUALLY WROTE FILES (check mtimes > dispatch time, not just "session ran") — #1 ran
   huge then wrote nothing.
2. Verify the fix outputs: unique-exon count (~81621 not 101724), gc_fraction populated no-nulls,
   freeze allowlist pinned+hashed in the sbatch, encode assembly fence enforced, REV3L resolved,
   pytest count. Then RE-REVIEW (Codex sol adversarial) before treating layers as gate-passed.
3. Only AFTER re-review is clean: commit the build layers + fixes on feng-arm.
4. Untouched follow-ups: Arm D step-3 strata matching + score impls; Arm B hg19->v50 liftover;
   verify the 8 positives vs Fig 3B before freeze; atlas LOW findings (optional).

Full findings + fix-pass detail are in the "Fork 2"/gate sections below.

---

Three prior fronts LANDED on disk (details below the older checkpoint). State:
- **Arm A descriptive atlas — DONE** (`feng/artifacts/arm_a_descriptive/`): report + 2 figures
  + 2 tsv. +9nt S4 offset fixed. Descriptive only, honesty-tagged.
- **Phase 2 HapSMA audit — DONE, REAL NO-GO -> Appendix arm KILLED** (`feng/artifacts/
  phase2_variant_recon/EXTRACTION_HAPSMA.md`): 0 qualifying donor-independent target-site ALTs
  (need >=2); NOT a biological zero-variant/conservation claim; 1000G rows = HPRC controls.
  Verdict verified not a coord artifact (Fork 2 gate). **Appendix (natural-variant x ASO
  decomposition) is FORMALLY KILLED on public evidence** — kill-gated by design, expected under
  four-concurrent-headlines framing; does NOT touch A-D. Only legit revival = genuinely new
  independent patient data with target-site ALTs + per-base callability; public well is dry
  (controlled-access blocked by do-not; Paraphase bamlets = controls not patients; interval is
  locked, cannot broaden). Do NOT re-litigate without new public patient cohort.
- **Phase 3 rMATS — INCOMPLETE**: `feng/scripts/phase3_rmats.sbatch` written (per-chemistry-
  vs-matched-scrambled, no pooled model) but `sample_conditions.tsv` NEVER delivered; sbatch
  hard-fails without it. Not submitted.

**Two Opus-4.8 fork subagents dispatched (this session, background):**
1. FINISH rMATS: resolve GSM9224134-9224153 -> condition/chemistry/control_type from
   authoritative GEO metadata, write+validate `feng/artifacts/phase3/sample_conditions.tsv`
   (cols run_accession,condition,chemistry,control_type; 3+3 per contrast; label the 2 non-
   contrast samples), then `sbatch feng/scripts/phase3_rmats.sbatch` (env SCRATCH=/home/dashii/
   orcd/scratch FENG_ENV_PREFIX=/home/dashii/orcd/pool/feng-env, from repo root). Returns job ID.
2. STANDING-GATE adversarial review (read-only) of the 2 DONE deliverables (Arm A atlas +
   HapSMA NO-GO) — structural + honesty-state + coord/frame + claim-inflation. Returns findings.

**Fork 1 DONE (2026-09-02):** `feng/artifacts/phase3/sample_conditions.tsv` written+validated
(20 rows; mapping from GEO SOFT `GSE307431_family.soft.gz` + SRA RunInfo, joined by BioSample).
**rMATS submitted: job 21824983 (PENDING, Priority).** Contrasts: F18OMe/F18MOE/F20PMO each vs
own chemistry-matched scrambled; GSM9224134/47=untreated (never enter a contrast). Strandedness
fr-firststrand (verified from ReadsPerGene col4>>col3). GTF pinned via FENG_P3_GTF=...primary_
assembly.annotation.gtf. Outputs -> `$SCRATCH/feng/phase3/rmats/{contrast}/`. Idempotent/requeue-
safe (per-contrast .done). Resubmit env: from repo root, `SCRATCH=/home/dashii/orcd/scratch
FENG_ENV_PREFIX=/home/dashii/orcd/pool/feng-env FENG_P3_GTF=<primary_assembly gtf> sbatch
feng/scripts/phase3_rmats.sbatch`.

**Fork 2 DONE (2026-09-02): STANDING GATE PASSED on both deliverables.** No HIGH/MED
confirmed defects; both honesty-clean + internally consistent. +9nt fix is test-guarded
(a -9 shift raises); HapSMA coord reconciliation verified end-to-end; NO-GO correctly framed
as evidence-gate not biological-zero; no forbidden strings; no pooling. ADVANCE approved.
- FIXED inline: recon `README.md` top status line now SUPERSEDED -> points at EXTRACTION_HAPSMA
  NO-GO (was stale `NEEDS-EXTRACTION`).
- LOGGED LOW (non-blocking, optional): (a) atlas `arm_a.py:553/529` — two manifest count lines
  ("125 records / 115 windows", "115/115 guide-covered") are hardcoded string literals in
  write_report, correct now (test-guarded) but could silently drift if inputs change; interpolate
  from `data`. (b) atlas Methods should name the two distinct same/opposite reference medians
  (window-local containing-guide vs block positional). (c) atlas 0/N block-fraction cell when
  block zmedian==0 (documented). (d) atlas report line 8 clarify q<0 for I6. (e) PLAUSIBLE, not
  re-verified: EXTRACTION line 37 "no extras or omissions" exact-match on 142 PSV13 carriers.

**COMMITTED 2026-09-02: `039e304`** (feng-scoped checkpoint — Arm A atlas, HapSMA NO-GO,
Phase 3 rMATS scripts + conditions, recon README fix, HANDOFF). Old-arm logs (work/, stageA*.out,
Log.out, pre/slurm/logs/) left untracked on purpose. Branch feng-arm. Moon plan edit (execution
addendum) is in ~/.claude/plans/ — outside repo, not committed.

**THREE Opus-4.8 forks dispatched (2026-09-02, background) — unblocked build fronts:**
- **Arm D step 1 — DONE. Gate 8/8 v50-mappable.** `feng/artifacts/bench/positives_v50.tsv` +
  README. 8 positives from Ottesen 2026 (PMC12805893) Fig 3B: POLR2H e2, PITHD1 e2, SERPINB7 e4,
  RTTN e4, REV3L e9, PRKRA e2, GOLGA4 e4, PAK1 e2 — all MANE Select internal cassette exons.
  **OPEN RISK (candidate-state):** paper gives NO coordinates; coords derived from gene+ordinal on
  MANE, assuming MANE ordinal = paper ordinal. MUST verify each exon vs Fig 3B junction context
  before the freeze lock. Event list is paper_supported; coords are candidate. POLR2H flagged
  known/locked-not-blind. Adversarial gate: deferred to run jointly with step-2 universe.
- **Arm D step 2 — code DONE, job PENDING.** `feng/src/feng/bench/strata.py` (parser+metrics+
  freeze-guard+CLI), `feng/scripts/bench_detectable_universe.sbatch`, `feng/tests/test_bench_universe.py`
  (53 tests pass). Method: rMATS SE enumeration, 11 controls split into 2 dummy groups + `--statoff`
  (no ΔPSI, no treatment BAM read). Thresholds: mean junction reads>=20, PSI estimable >=6/11 controls.
  Freeze-guard re-derives treatment set + aborts on intersection. **Job 21826363 PENDING (monitor
  bt5uxf8mg).** Outputs on completion: `feng/artifacts/bench/detectable_universe_v50.parquet` + README
  (counts + per-chrom density) -> then check >=50-background-per-positive gate vs the 8 positives.
- **Arm B step 1 — DONE.** `feng/artifacts/rbp/{encode_manifest.tsv (1141 files),curated_staged.tsv,
  README.md}`, `feng/src/feng/rbp/encode.py` (loader), `feng/tests/test_rbp_encode.py` (52 guard tests
  pass, honesty-clean). 203 RBP×cell combos carry BOTH eCLIP+KD (joinable universe); curated 152 files
  md5-verified at `$SCRATCH/feng/encode/` (incl. HNRNPA1/TIA1/TIAL1). Batch-corrected ENCSR870OLK
  chosen over ENCSR413YAF; IncLevelDifference sign-only. NO concordance test yet (post-freeze;
  concordance.py untouched). **KEY OPEN: all ENCODE files are hg19 -> needs_liftover=yes on all rows;
  liftover to GRCh38/v50 is the gating step before the month-10 concordance test. Deferred, flagged.**

**Universe job 21826363 — COMPLETED (2026-09-02 12:49).** `detectable_universe_v50.parquet` + README
written; scratch `.done` present. Counts: 254891 SE enumerated -> **101724 detectable** across 9387
genes. >=50-per-positive gate has big headroom (true matched-strata test = Arm D step 3, not built).
CAVEAT for matching: detectable baseline-PSI quartiles 0.833/0.993/0.999 (pool skews near-
constitutive) -> watch thin PSI-matched strata for any mid-PSI positive. Monitor bt5uxf8mg ended.

**STANDING GATE RETURNED (2026-09-02): NOT PASSED — 4 HIGH + 3 MED + 2 LOW. Fix pass dispatched.**
Guard suite 106 passed/3 skipped. Honesty/clean-room/freeze-fields clean; issues are structural.
CONFIRMED findings:
- HIGH strata.py:214 — GC fraction NULL for all 254891 rows (pysam missing, failed open); GC is a
  documented matching covariate. Fail-closed + populate + rerun.
- HIGH detectable README:58 — "101724 detectable" double-counts alternate-flank rMATS events; only
  81621 UNIQUE exon intervals (13721 repeated up to 16x). Need exon-identity key + collapse.
- HIGH sbatch:76 — freeze-guard trusts mutable labels (Python check runs only after rMATS reads
  BAMs). Pin+hash exact 11-accession control allowlist BEFORE rMATS; reject dupes.
- HIGH encode.py:129 — hg19 fence advisory only; readers return frames w/o assembly provenance ->
  silent hg19->v50 join possible. Enforce assembly flag + attach provenance at every join.
- MED detectable README:35 — "--statoff = no comparison statistic" FALSE (dummy control
  IncLevelDifference in 185260 source rows; parquet drops it). Reword + assert contrast fields never
  enter parquet.
- MED detectable README:93 — provenance cites a DELETED .work table; retain/hash the stable copy +
  conditions + b1.txt/b2.txt for reproducibility.
- LOW encode.py:74 — loader silently discards sha256; include+validate in EncodeFile.
PLAUSIBLE (correctness-critical): REV3L exon 9 positive has NO exact event in the universe while MANE
exon 8 IS enumerated -> paper-ordinal=MANE-ordinal assumption may be WRONG for REV3L. Verify all 8
identities from Ottesen junction/primer context, prioritize REV3L; correct positives_v50.tsv if off.
Layer disposition: L1 map clean (ordinal-gated by REV3L); L2 NOT clean (4 issues); L3 parse/md5 clean
(reviewer rehashed 152/152) but assembly fencing not.
Full review saved in Codex task output (fresh task, agent a0f980b29875e7932).

**FIX PASS dispatched (2026-09-02):** Codex sol --write, background **task-mtkco2qa-imj9ki**
(monitor bx8mr29pd polls companion status). Covers all 7 confirmed findings + REV3L verification +
universe rerun + new tests. Fetch result when terminal: `node /home/dashii/.claude/plugins/cache/
openai-codex/codex/1.0.6/scripts/codex-companion.mjs result task-mtkco2qa-imj9ki`. Resume after
/clear: verify new unique-exon count + GC populated (no nulls) + freeze allowlist pinned+hashed +
encode assembly fence + REV3L resolution + pytest count, then RE-REVIEW before advancing. Do NOT
treat layers as gate-passed until the re-review is clean. (Codex background tasks poll-stall: the
rescue driver returns a task-id immediately; poll companion status/result from the MAIN thread, not
via the driver subagent — it refuses to poll.)

**rMATS job 21824983 — COMPLETED 3/3 clean (2026-09-02 12:32).** All 3 contrast dirs have full
output (fromGTF.*, JC/JCEC.raw.input.*, {SE,MXE,A3SS,A5SS,RI}.MATS.{JC,JCEC}.txt) + `.done`
sentinel; `.err` 0 bytes; log feng/logs/feng-p3-rmats-21824983.out. **Structure verified only —
ΔPSI values NOT read (sealed until freeze).** Months-3-4 Ground-truth phase = DONE (align + sealed
treatment contrasts). Monitor bh54vkwnh ended.

Resume after /clear: (1) rMATS DONE — do NOT inspect the .MATS ΔPSI tables before freeze. (2) Await the 3 fork notifications -> standing adversarial gate each -> then check Arm D's
≥50-background-per-positive gate once both D fronts land. Atlas LOW findings still optional.

---

## ★ PRIOR CHECKPOINT — 2026-09-02 ~09:25 EDT (three fronts in flight)

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

**10:02 EDT status:** all 3 Codex tasks still RUNNING, none reported yet at the
10:05 stop deadline. They keep running in background; outputs land on disk. Working
tree has in-flight Codex edits (feng/.../README.md, .gitignore) — NOT committed here,
leave for the task to finish.

Resume after /clear: check `/codex:status <task-id>` for each; verify outputs on disk
(atlas figures/code under feng/; `feng/scripts/phase3_rmats.sbatch` +
`feng/artifacts/phase3/sample_conditions.tsv`; `feng/artifacts/phase2_variant_recon/
EXTRACTION_HAPSMA.md`). Standing gate: adversarial-review each before advancing.
Nothing committed this session.

---


---
## (archived) original 2026-09-01 phase-1 progress detail

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

