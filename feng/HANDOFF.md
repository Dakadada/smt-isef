# HANDOFF — feng/ (the new prospective arm)

**Created: 2026-09-01. State: greenfield. Nothing built yet. Plan is frozen.**

This file is the entry point for anyone (Claude or Codex) picking up the `feng/`
arm. The old `src/`/`pre/` program is a separate, committed effort — see
`docs/HANDOFF.md` for it. `feng/` does not touch that code.

---

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
