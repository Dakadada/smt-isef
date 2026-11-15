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
## ⏩ CONFERENCE POSTER — LIVE (2026-09-30)
GOAL (user): get compute out. Keep HANDOFF current so /clear is always safe. Codex MCP still DOWN
(CONNECTION_CLOSED) but the **CLI route WORKS** — see Codex route note below. sol-low rule stands.
All prior monitors/heartbeat STOPPED (user paused, then resumed).
- **★★★ 2026-10-28 — ITEM 2 (PSI efficacy) REANALYSIS: endogenous size series band identity SOLVED by
  blind replicate validation. NOT yet adversarially reviewed (standing gate OPEN).**
  Source data = the two endogenous size-series gels in repo root (untracked): `Singh Lab Quantified
  Data - 7-27 truncated ASOs endogenous (including RT PCR) samples 1-9.csv` (audit "Gel 8") +
  `...7-28 ASO size comparison w samples 1-9 on endogenous gene (+RT PCR).csv` ("Gel 9"). They are a
  REPLICATE PAIR of the same size series → used as the independent (non-circular) scorer.
  METHOD: blind brute force over all 6 band-pair hypotheses (product bands at Rf clusters ~0.40 /
  ~0.48 / ~0.56 / ~0.58; inclusion = lower-Rf = larger amplicon), scored ONLY on cross-gel replicate
  agreement — biology NEVER fed into selection, applied only as an after-the-fact check.
  WINNER: **+E7 (inclusion) = ~0.56 band, ΔE7 (skipping) = ~0.58 band; the big ~0.40 band is a
  non-product (gDNA/heteroduplex), excluded.** Pearson r=0.991 / Spearman ρ=0.967 (NOT ceiling-driven:
  only lanes F18 & 3UP8 peg at 100%, dropping them r=0.995; max lane-to-lane |Δ|=4.0 pp). The lab
  audit's PHOTO band call (+E7=upper0.40 / ΔE7=0.56) ranks DEAD LAST (r=0.030) → this reanalysis
  OVERTURNS the workbook photo assignment in `psi25/SMN2_all_sources_PSI_band_assignment_audit.xlsx`.
  RESULT under the FINAL-PSI loading order (lane2=WT,3=Scr18,4=F18,5=Scr14,6=F14,7=Scr10,8=F10,
  9=Scr8,10=3UP8) — mean inclusion PSI% (G8/G9): F18 100(100/100), 3UP8 98(100/96), F14 95(94/96),
  F10 91(90/93) | Scr18 74, Scr10 71, Scr8 67, WT 59, Scr14 57. All 4 actives (91-100) ≫ all 5
  controls (57-74), clean bimodal gap; **F18>F14>F10 potency order correct** → CORROBORATES item 1 ΔG
  rank (F18 -24.7 > F14 -17.0 > F10 -11.1): binding energy predicts inclusion.
  CAVEATS (do NOT publish as measured): (1) RANK-only, not quantitative — F18/3UP8 peg at 100% ceiling
  (skip band ~absent) + WT=59% baseline is biologically high (untreated SMN2 usually ~10-20%) ⇒
  absolute PSI miscalibrated; safe claim = ORDERING only. (2) Band rule is replicate-VALIDATED but the
  lane→ASO LABELS still assume the FINAL-PSI loading order (the 4-active/5-control split is forced by
  the PSI gap, but which active=F18 vs F14 is not). (3) CONTRADICTS the gel photo → needs ONE
  amplicon-size anchor (ladder MW or junction/seq call) to confirm 0.56=+E7 before publication.
  7-17 POSITIONAL (E3 walk A-F) endogenous pair ATTEMPTED, FAILED validation: trial-1 lanes 5-8 are
  single-band → PSI undefined → replicate anchor collapses (best r=-0.31). Too degenerate to blind-
  validate from this export; NOT presenting a number. Only the size series survives.
  Scripts (session scratchpad, SESSION-SCOPED — rebuildable from this entry if gone):
  `/tmp/claude-264670/-home-dashii-smt-isef/1f0a4778-fd3f-44fb-94ee-1a8a1119f4cd/scratchpad/`
  {brute.py (6-hypothesis blind search), psi_perm2.py (coherence), xlsxdump.py (audit reader)}.
  NEXT: STANDING-GATE Codex sol adversarial review of this reanalysis BEFORE it feeds any artifact;
  lab ask = (a) lane→ASO map + (b) one amplicon-size anchor. No artifact written yet.
- **★★★ 2026-10-05/06 — ITEM 3 JOB 24817368 COMPLETED, B DEFECT CONFIRMED FIXED. REVIEW IN FLIGHT.**
  `sacct`: **COMPLETED exit 0:0**, elapsed **1-11:57:02** (~36h, vs the 4.5h ETA — shards slowed from
  ~1min to ~40min each late in the run), End 2026-10-05T22:07:21, MaxRSS 7.67G of 120G. **241/241**
  accessibility shards `.done`. Tables rewritten 22:07 in `conference/artifacts/item3_spec_{occupancy,
  gated,absolute}.tsv` + `item3_spec_provenance.txt`. NOTE job stdout/stderr logs are **0 bytes** (shard
  loop prints nothing) — shard `.done` mtimes in
  `/orcd/scratch/orcd/011/dashii/conference_smn2_aso/results/item3_spec_v1/accessibility_shards/` are
  the real progress signal, NOT the logs.
  **B IS REAL NOW:** applicable/evaluated/**failed** → F18 10/10/**0**, F14 1751/1751/**0**,
  F10 320202/320202/**0** (was 100% failed). Accessible loci (genes): F18 3 (3), F14 1133 (1114),
  F10 204838 (35268). Gate cuts everyone ~3-4× but **ordering unchanged** — @100nM `sum_theta_off`
  ungated→gated: F18 16→3, F14 4144→1133, F10 577057→148952; burden F10/F18=**49,650.7×**,
  F14/F18=377.7×; specificity ratio F10/F18=**1.995e-5**, F14/F18=2.648e-3; `theta_on_over_sum_theta_off`
  F18 0.3333 / F14 8.826e-4 / F10 6.650e-6. Expression still `unavailable` (expected).
  ⇒ **A, B, C now all agree: length DEGRADES specificity. No regime shows honest parity.** Headline
  stands as "length degrades specificity; occupancy- AND accessibility-weighting confirm it".
  **⚠ STANDING-GATE REVIEW NOT YET RUN — dispatched then CANCELLED by user 2026-10-06 within ~1 min.**
  Task b7gcw4xm1 killed (`TaskStop` confirmed, no codex process left); `item3_review.out` 0 bytes, NO
  verdict produced, sandbox was read-only so repo/artifacts untouched. **The gate is still OPEN — item 3
  must NOT advance to README reframe until it runs.**
  RE-FIRE verbatim (prompt file survives, no rebuild needed):
  `cd /home/dashii/smt-isef && codex exec -m gpt-5.6-sol -c model_reasoning_effort="low" --sandbox
  read-only --skip-git-repo-check "$(cat <scratchpad>/item3_review_prompt.txt)"`
  where `<scratchpad>` = `/tmp/claude-264670/-home-dashii-smt-isef/6ccb40ba-d2c7-4718-8218-723d9b48c985/
  scratchpad` (SESSION-SCOPED — if that dir is gone, the prompt must be rebuilt; its 9 attack points are
  summarized next so it can be).
  Scope was item3_spec.py + test_item3_spec.py + the 3 tables + provenance + sbatch. 9 attack points:
  (1) does the U→T fix merely SILENCE the fail-closed check — could a WRONG locus now pass? (2) F10
  theta_on=0.511@1nM vs ~1.0 for F18/F14 — bug, and is theta_on/sum_theta_off well-posed across ASOs
  with differing theta_on? (3) mass-balance violation — ~149k simultaneously occupied sites from a
  finite pool; does it inflate or deflate the F10-vs-F18 contrast? (4) F18's denominator-of-3
  fragility; (5) is the accessibility gate applied SYMMETRICALLY to on- vs off-target? (6) 3UP8
  saturated-empties not coerced to 0/NaN; (7) shard-merge integrity — header leakage, dupes, does
  anything assert 241? (8) `query_id " [candidate]"` suffix severity + is query_id ever a join key?
  (9) any other "confident number on absent evidence" (the B-was-all-failed class) recurrence.
  **NEXT after review:** fix/log findings → item3 README reframe → then the `query_id " [candidate]"`
  suffix cleanup (cosmetic, severity being re-confirmed by the review).
- **CODEX ROUTE (2026-10-06):** MCP `codex mcp-server` still CONNECTION_CLOSED, but **CLI works** —
  `codex-cli 0.159.2` at `~/.local/bin/codex`, auth OK, smoke-tested `ROUTE_OK`. Dispatch via
  `codex exec -m gpt-5.6-sol -c model_reasoning_effort="low|medium" --sandbox read-only|workspace-write
  --skip-git-repo-check`. ⚠ **`~/.codex/config.toml` default is now `gpt-6.1-sol` effort `low`** — newer
  than the 5.6 named in CLAUDE.md's tiering table, so ALWAYS pass `-m` explicitly or you silently get
  6.1. User DEFERRED deciding whether to move the tiering table to 6.1. Codex sandbox still CANNOT
  reach the Slurm controller → Claude submits all jobs from the main shell.
- **★★ 2026-10-04 — ITEM 3 "B" DEFECT ROOT-CAUSED + FIXED + RERUN LAUNCHED (job 24817368)** (SUPERSEDED
  by the completion entry above; cause/fix detail retained).
  Cause was NOT data: pure ALPHABET mismatch. `feng/src/feng/titration/evaluate.py:313` fail-closed
  compares the rung target against a GENOMIC (DNA) fetch verbatim; conference targets are stored RNA (U).
  genomic `CCAGCATTATGAAAGTGA` vs target `CCAGCAUUAUGAAAGUGA` — identical but for T/U → ViennaError
  "padded window does not contain the rung target" for EVERY U-containing locus. Arm E passes
  DNA-lettered rungs, so Arm E never hit this.
  FIX (conference-local; shared feng/evaluate.py deliberately UNTOUCHED): `item3_spec.py`
  `add_accessibility` now applies `target_5to3.upper().replace("U","T")` before the call, and records
  `failed:<ExcName>` instead of a bare "failed" so the reason survives next time.
  VERIFIED on real loci: 4 available (e.g. accessibility=0.4058) / 36 not_applicable / **0 failed**
  (was 100% failed). `conference/tests/test_item3_spec.py` 6 passed.
  Cleared 241 stale all-failed `*.accessible.tsv`(+`.done`) + the merged accessibility TSV; the 19h
  bowtie/loci stage (`item3_spec_loci.tsv.done`) and `split.done` are REUSED (not recomputed).
  **RESUBMITTED job 24817368** — RUNNING, verified healthy mid-flight: shard 1 gave 288 available /
  9712 not_applicable / **0 failed** (was 100% failed); progress 27/241 shards in ~30min → ETA ~4.5h.
  Loci TSV integrity checked: 0 malformed rows / 2,408,841. Known COSMETIC-only artifact: the
  `query_id` VALUE has " [candidate]" appended inside it (e.g. `EV_F18_a0ae2d60de0b31fa [candidate]`);
  no column shift, dg joins resolved fine — low-severity cleanup, only matters if query_id is used as
  an external join key.
  WATCHER: background task **brafl3zm5** (Bash until-loop, fires ONCE at terminal and prints the B/gated
  accessibility columns). Replaced Monitor bwkpowgyo/byjyq7od4 etc — the Monitor tool caps at 30m here
  and churned re-arms; a backgrounded until-loop is the right tool for a single completion signal.
  ~322k applicable loci
  now actually invoke RNAplfold (previously died instantly on the exception), so expect real runtime;
  2-day limit + per-shard `.done` + requeue make it safe.
  ON COMPLETE: re-read `item3_spec_gated.tsv` — B should finally carry real
  accessible_off_target_locus_count / gated occupancy instead of degenerate 1.0 ratios. Then
  standing-gate review + item3 README reframe.
- **★ 2026-10-04 — ITEM 2: F14/F10 ARE usable, via CENSORED upper bounds.** Their inclusion peak was
  never exported, so no point PSI exists; treat as left-censored: PSI < LOD/(LOD+skip). Gel7 skipping
  volumes F18 239 / F14 217 / F10 159; F18 inclusion 65 → **F18 PSI = 21.4% (detected)**. With LOD=33
  (faintest peak exported anywhere on Gel7): **F14 < 13.2%, F10 < 17.2%**. F18 exceeds BOTH bounds for
  any LOD ≤ 50; at LOD ≥ ~58 the F14 bound crosses F18 and the ordering breaks. So the ordering
  F18 > F14,F10 is real but ASSUMPTION-DEPENDENT — the workbook states no LOD was provided, so this must
  NOT be presented as measured. ROBUST FIX = lab re-quantifies with a manual ROI box at the inclusion
  Rf (~0.55) in the F14/F10 lanes. AWAITING user decision on which item-2 form to publish.
- **ITEM 3 — spec replay job 24031068 COMPLETED** (19h19m, exit 0). Tables in conference/artifacts/
  item3_spec_{occupancy,gated,absolute}.tsv. VERDICT: occupancy-weighting does NOT make shorter
  "about the same" — it makes F10 WORSE. @100nM: off-target occupancy burden F10≈36,066× F18, F14≈259× F18;
  specificity ratio F10≈2.7e-5 of F18. A(occupancy): F18 sum_theta_off=16, F14=4144, F10=577,057 @100nM.
  C(absolute-dG): at loose bars F10≫F18 (−8: F10 49,439 genes vs F18 1,327); at strict −14/−16 F10=0
  (artifact — F10 on-target only −10.3, can't bind anything that tight). No regime shows honest parity.
  ⚠ **B (accessibility gate) DEFECTIVE**: accessibility_evaluated=0 for ALL loci (all "failed") →
  accessible_count=0, gated ratios degenerate=1.0. RNAplfold stage produced nothing usable. NEEDS FIX
  before B is reportable. Expression stayed unavailable (as expected).
  NEXT item3: fix B accessibility stage; standing-gate review; then item3 README reframe (headline =
  "length degrades specificity; occupancy-weighting confirms it" — NOT "about the same").
- **ITEM 1 — REDO with user buffer (50 mM mono + 1.5 mM Mg2+) DONE (code+artifacts+tests+README).**
  thermo.py: new `effective_monovalent_from_mixed_salt` (von Ahsen 2001 [Na+]eq=[mono]+120*sqrt([Mg]mM),
  =0.197 M). item1.py: PRIMARY now the mono+Mg buffer; raw cols renamed nn_tbe_*→nn_buffer_* + added
  monovalent/divalent/effective basis cols; provenance records buffer+equivalence. Artifacts REGENERATED
  (item1_dg_raw/salt_sensitivity/ranks/provenance). NEW NN buffer dG: F18 −24.71, F14 −16.95, F10 −11.11,
  3UP8 −13.09 (raw 1.021M refs unchanged; IntaRNA unchanged −27.0/−18.89/−12.8/−14.4). Ranking unchanged
  F18>F14>3UP8>F10. README updated (method+von Ahsen cite+caveat). test_thermo.py 8 passed (added von Ahsen test).
  CAVEAT: heuristic (DNA-param salt eq + Mg→Na on RNA; blocker #1 compounded). Uncommitted. NO commit (policy).
- **ITEM 2 — band identity RESOLVED by user, but data-limited.** User rule: inclusion=TOP (faint, mostly
  undetected), skipping=BOTTOM consistent line (Gel7 ~Rf0.60), top bands=noise. Inclusion band detected
  ONLY in lanes 3&4 = Scr18, F18. So among F-series ONLY **F18 has quantifiable PSI** (~21% inclusion,
  65/(65+239)); F14(lane6)/F10(lane8) show NO detected inclusion band → PSI below detection, NOT 0 (no LOD).
  3UP8 not in Gel7. ⇒ **ΔG-vs-PSI correlation NOT possible from Gel7 (n=1 F-series point).** Super-minigene
  (7-7/7-9) PSI still blank. NEEDS USER DECISION on how item 2 proceeds. Band source: psi25/
  SMN2_all_sources_PSI_band_assignment_audit.xlsx (Gel 07 sheet).
  Tracked in Task tool: #1 Item1 DONE, #2 Item2 pending(PSI-blocked), #3 Item3 in_progress, #4 Item6 DONE.
- **NEW USER RULE:** dev Codex dispatches = **sol LOW only, never sol high** (memory [[sol-low-only]]).
- **ITEM 3 REFRAME — thermodynamic specificity, tests A+B+C IN FLIGHT.** Raw bowtie -v3 AND energy-filtered
  gene COUNTS both mislead (weight all loci equally; F10=12,532 energy-filtered genes @10x vs F18=2 —
  genomic degeneracy of short seqs, not real bound load). Ordering never flips (F18≪F14≪F10 every
  threshold), so counting can't show "shorter ≈ same". Building predicted-bound-LOAD metrics instead:
  A=occupancy-weighted specificity ratio (theta=[A]/([A]+Kd), Kd=exp(dG/RT), sweep [A] 1nM-1uM);
  B=expression+accessibility gate (EXPRESSION UNAVAILABLE in repo → accessibility-only via RNAplfold);
  C=absolute-dG-threshold recount (sweep -8..-16). Equivalence margin for "F10≈F18" DEFERRED (compute
  raw ratios first, decide later — NO pass/fail baked in). 3UP8 stays saturated.
  **Codex task kjqzfaxmj DONE** (thread 01a0df31, gpt-5.6-sol LOW): built conference/aso_pipeline/item3_spec.py
  + tests/test_item3_spec.py (6 passed) + conference/scripts/item3_spec_replay.sbatch + README/provenance +
  3 tables (item3_spec_{occupancy,gated,absolute}.tsv) marked `pending_slurm` until the job runs. 3UP8=saturated.
  **SUBMITTED replay: job 24031068** (mit_preemptable 32c/120G/2d, requeue+SIGTERM@180, preflights inputs
  exit-2). MONITOR byjyq7od4 (30m — RE-ARM on expiry; job may queue/run long). Job atomically replaces the
  pending tables with real F18/F14/F10 occupancy(A)+gated(B,accessibility-only)+absolute-dG(C) numbers.
  NEXT: on COMPLETE → read A/B/C + F10/F18 ratios → user sets "about the same" margin → standing-gate
  adversarial review (sol LOW, fresh thread) of item3_spec code+results → then item3 README reframe + blockers.
- Energy-filter artifacts already landed (prior session): conference/artifacts/item3_energy_*.
- (below: prior 03:16 checkpoint, superseded where it overlaps)

## ⏩ CONFERENCE POSTER ITEMS 1/3/6 — earlier checkpoint (2026-09-26 03:16)
GOAL: 4-ASO F-series poster deliverables (conference/PLAN.md scope). Items 1,3,6 dispatched; 2 blocked on PSI.
- **ASO set** (5'->3', RNA:RNA @37°C assumed): F18=TCACTTTCATAATGCTGG, F14=TTTCATAATGCTGG,
  F10=ATAATGCTGG, 3UP8=GCTGGCAG (3UP8 target site UNCERTAIN — flagged for confirm).
- **★ ITEMS 1 + 6 DONE (k2n0udpyv, 22/0/0 tests).** Code in conference/aso_pipeline/, artifacts in
  conference/artifacts/ (item1_dg_raw.tsv, item1_salt_sensitivity.tsv, item6_thermo_null_summary.tsv +
  READMEs). ΔG both models rank F18>F14>3UP8>F10; NN@0.089M-TBE: F18 -23.17, F14 -15.86, 3UP8 -12.45,
  F10 -10.29 kcal/mol. Item6 null: F18/F14 real enrichment (p=0.0099); F10/3UP8 nulls degenerate (only
  8/6 unique dinuc shuffles) → not interpretable. RNG seed 20260926.
- **★ ITEM 3 — job 23993733 OOM-KILLED at 40min** (bowtie -a -v3 -p16 on F10's 21.6M hits blew 64G).
  FIX IN FLIGHT: **Codex task kithorcst** (thread 01a0dc93…) — bounded-RAM streaming gene-id dedup
  (disk sort -u, not in-RAM set), drop -p to 4-8, per-ASO steps so F10 gets full headroom, --mem to
  180-250G safety, keep requeue+ckpt. FIX DONE (kithorcst, 27/0/0, bounded-RAM proven: 100k dup hits <5MB).
  **★ RESUBMITTED — job 23996546** (mit_preemptable, 8cpu/220G, disk sort -u dedup, per-ASO steps, F10 split
  into 8 singleton passes -p1). **MONITOR bfwna3jmi** (re-arm on expiry). NOTE: Codex sandbox CANNOT reach
  Slurm controller — Claude submits from main shell every time. Old job 23993733 (OOM) + monitors dead.
  Preflight: F18 (2109), F10 HEADLINE (21.6M), F14 (246k), 3UP8 SATURATED (168M, skip).
- **★ job 23996546 COMPLETED but OUTPUT DEFECTIVE (2 bugs).** (1) bowtie -v3 counts ALL ZERO despite
  bowtie logging 21.8M/23.8M/24.9M alignments → counter drops all bowtie hits (blast path works, bug is
  bowtie-specific: --suppress column layout / stdout parse mismatch). bowtie -v3 is THE headline metric.
  (2) All rows keyed by opaque Q_<hash>, no map to (aso_id, real-vs-scramble) → summary unreadable. Also
  blast E=1000 too permissive (thousands of genes, not ~445 ballpark) — flagged for user, not retuned.
  FIX IN FLIGHT: **Codex task kvf7mmpd2** (thread 01a0dc93…) — fix bowtie→gene mapper w/ known-answer
  regression test, persist Q_hash→ASO map, emit clean labeled summary (real F18/F14/F10 both aligners +
  scramble-null dist separately), flag blast E-value. Will resubmit (Claude, main shell). Monitor bfwna3jmi ended.
- **★ BOWTIE BUG ROOT CAUSE + FIX (kvf7mmpd2, 29/0/0):** trailing `-` in bowtie cmd = literal FILENAME
  (dumped 846MiB to ./- , cleaned) not stdout → counter got nothing. Fixed bowtie_stream.sh; bowtie stage
  versioned v2 so bad zero .done ignored; persists Q_hash→ASO map; emits item3_headline_f18_vs_f10.tsv +
  item3_real_aso_gene_counts.tsv + item6_alignment_scramble_summary.tsv + item3_query_identity_map.tsv.
  blast E=1000 flagged maximally-permissive in README; Codex rec = bowtie -v3 as headline, tighten blast
  only w/ user sign-off. **★ RESUBMITTED job 23997871, MONITOR bbkmysq5v** (re-arm on expiry).
- **★★ ITEM 3 COMPLETE (job 23997871 COMPLETED, bowtie counts now real).** Outputs in
  conference/artifacts/: item3_headline_f18_vs_f10.tsv, item3_real_aso_gene_counts.tsv,
  item3_query_identity_map.tsv, item6_alignment_scramble_summary.tsv. HEADLINE (bowtie -v3, distinct
  transcript-compat genes): **F18=1518, F14=37713, F10=71148** → SHORTER mer hits MORE genes (F10 ~71k
  ≈ whole genome). blastn-short E1000: F18=1588, F14=684, F10=0 (10-mer below blast-short threshold).
  ⚠ **RESULT CONTRADICTS the poster's naive item-3 framing** ("18-mer hits ~445, 10-mer fewer"). At FIXED
  ≤3 mismatches a 10-mer (7/10 fixed) matches ~everywhere; 18-mer (15/18) is specific. The 10-mer
  specificity story is THERMODYNAMIC (item1: F10 ΔG -10 vs F18 -23 → only near-perfect sites bind stably
  @37), NOT fixed-mismatch alignment count. Likely need an ENERGY-FILTERED off-target count to make the
  point honestly. DECISION for user: reframe / add energy-filtered count / proceed to standing-gate review.
- **STANDING-GATE ADVERSARIAL REVIEW DONE (sol high, task khrhb5be7) → DO-NOT-ADVANCE, 6 blockers**
  (mostly honesty/framing, NOT broken compute). Tests 21/0/0 read-only (8 need writable tmp; live=29/0/0):
  1. HIGH thermo.py:75 — SantaLucia DNA salt eq on Turner RNA energies + Tris-as-0.089M-monovalent = not
     validated for RNA/TBE → F18 -23.17/F10 -10.29 are heuristic, not defensible quant TBE ΔG. RAW
     IntaRNA/RNAduplex values still usable if separated from the salt interpretation. (Energy-filter task
     uses IntaRNA not NN-salt → sidesteps.)
  2. HIGH item3_README:44 — calls bowtie -v3 "defensible/directly interpretable" but 30% leniency (F10)
     vs 16.7% (F18) means it CANNOT support specificity; must reframe as negative/saturation result +
     point to energy-filtered analysis. (← kklvndnxf provides exactly this.)
  3. HIGH offtarget.py:195 — F10 BLAST "0" is short-query non-detection, not zero off-targets; needs
     below-threshold/unsupported state, not numeric 0 (false-negative fabrication risk).
  4. MED item6 F10/3UP8 nulls degenerate (8/6 unique seqs, real appears 12/13x) but presented w/ precise
     p-values as if valid → label non-interpretable. F18/F14 fine.
  5. MED "off-target genes" = gene-span intersection incl whole introns, ignores expr/accessibility →
     must be NAMED as annotated gene-span intersections.
  6. MED sbatch reuses indexes by filename+.done, no FASTA/GTF/index HASHES in provenance → stale-index risk.
  CLEAN: no_duplex flagged, raw/derived separated, no post-hoc tuning, seed/versions/tags recorded, 3UP8
  saturation→missing not fabricated, psi_stub refuses, bowtie regression fixed w/ known-answer test.
- **ENERGY-FILTERED off-target (item 3-thermo) IN FLIGHT — Codex task kklvndnxf** (thread 01a0de28, sol
  high). Threshold = <=10x weaker equil binding vs on-target (ΔG<=ΔG_on+RT·ln10≈+1.42), sweep 1x-100x,
  IntaRNA scorer, enumerate <=3-mm variants→score→viable-only exact-align→count genes. Resolves blockers
  #1/#2 scientific core. NOTE: old Codex threads die on session restart — start FRESH each session.
- **FIX PLAN for the 6 blockers:** ONE Codex fix task AFTER energy-filter lands (so item3 README reframe
  can cite it). Mostly README honesty + code: false-neg-zero state (#3), degenerate-null labels (#4),
  metric rename (#5), reference-hash provenance (#6), salt-model caveat (#1). Then terra re-review.
- **⚠ TWO USER CONFIRMS NEEDED:** (1) 3UP8 target = CUGCCAGC (ctx pos 7-14, 3nt upstream of F18) — derived,
  UNCONFIRMED. (2) chem model = RNA:RNA @37 assumed; DNA:RNA rerun needs a validated hybrid param file.
- **CODEX thread 01a0dc93… — design APPROVED, IMPLEMENTING as task k2n0udpyv** (was k0c0duufr for the
  design pass). Guardrail added: 8-mer 3UP8 at -v3 is combinatorially saturated → mark uninformative, don't
  grind; HEADLINE = F18 vs F10 exhaustive. 3UP8 target derived = context pos ~7-14 CUGCCAGC (FLAG confirm).
  NN model = RNAduplex Turner-2004 @37 + entropy salt correction, 1x TBE ≈ 0.089M monovalent + sensitivity grid.
  No commits (repo policy). Direct MCP notifies on completion.
- **(design pass) CODEX TASK k0c0duufr** (gpt-5.6-sol high, workspace-write, cwd repo root) — item 1 (ΔG two models:
  IntaRNA-default + NN@37°C w/ 1x-TBE monovalent correction, side by side), item 3 (bowtie -v3 + BLASTn-short
  gene-count 18 vs 10, sbatch on GRCh38 v50 FASTA @ ~/orcd/scratch/feng/gencode/v50/), item 6 (dinuc-preserving
  scramble ×100 through BOTH pipelines). Direct MCP → notifies on completion. Item-3 alignment is a sub-sbatch.
- **ITEM 2 (ΔG-vs-PSI correlation+regression): BLOCKED on user PSI table** (endogenous 7-27 truncated series,
  replicate-level, 4 ASOs + no-ASO control). Codex leaving a plug-in stub. Fires when PSI lands.
- **NEXT AFTER k0c0duufr:** standing-gate Codex sol adversarial review of the 1/3/6 deliverables before advancing.
- User setting alarm was one-shot cron 0786df47 (fired). Guard rails checked: quota healthy, motd clean.

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
- **#3 MINOR FIX DONE (thread 01a09070…):** M1 strip 16-byte .pyc header (benign recompile inert,
  body-tamper caught, truncated fails closed); M2 committed feng/src/feng/pangolin_package_seed_
  provenance.json. NEW pin PINNED_PANGOLIN_PACKAGE_SHA256 =
  0af9857403226c5e72719e42d52fa5babfcda8fde9957ce7824df277e0051234. Tests per-file GREEN
  (test_source_auth 12/0/0, test_arm_e_launch 61/0/0, test_titration_gpu_scores 23/1/0).
- **★ commit-1 LANDED = `83415d94a829b979c2ccbd4282a5d1569788bf1d`** (feng-arm, 2026-09-11). CODE
  trust root: source_auth.py + pangolin_package_seed_provenance.json + titration/{arm_e_hoff,
  arm_e_scan,gpu_scores} + bench/{scores,strata,match,permute} + gpu/checkpoint + scripts + tests +
  PREREG (DRAFT). report/window/armc_* + artifacts deliberately EXCLUDED (not in H-off closure).
- **FINAL RE-VERIFY DONE (task on thread 01a09154…, terra med) → LAUNCH-READY: YES.** B1+B2 both
  CLOSED as committed; pin+manifest+count6+version1.0.2 match; launcher SHA-bound; digest binds
  tree+commit+Pangolin identity; closure clean (no report/window dep); prereg DRAFT + TOCTOU retained.
- **⏳ LAUNCH PENDING USER:** harness classifier blocked the sbatch submit from Claude. User fires:
  `SCRATCH=/orcd/scratch/orcd/011/dashii RUN_ID=$(date -u +%Y%m%dT%H%M%SZ); bash feng/scripts/
  submit_arm_e_hoff.sh 83415d94a829b979c2ccbd4282a5d1569788bf1d "$RUN_ID"`. MOTD clean, SCRATCH 81%
  used (watch). On submit: capture RUN_ID + JOBID → re-arm state monitor (L10 ~20h long pole) →
  verify old L02-L10 payloads superseded into run dir → on 15/15 VALIDATE (no burden read) →
  **commit-2** flip prereg DRAFT→SEALED stamping commit-1 SHA + recompute job id + manifest sha.
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
