J5-FIX2 implemented; production DO-NOT-ADVANCE pending orchestrator lookup refresh, review, and intended-config Slurm profile.
conference/model_v2/jobs/j5_encode_lookup.py:32 — enforce collection paths, canonical @id/@type, file→experiment dataset binding, fetched/embedded label/term agreement, and independent biosample-term consistency; rereview repro now lookup_failed.
conference/model_v2/j5_rbp_context.py:183 — refuse explicitly unverified identity_status during ingestion; existing pinned legacy fixtures remain supported.
conference/model_v2/jobs/j5_benchmark.py:117 — deterministic top-12 density + 4 median windows from manifest-named J6 .done row counts only; no directory scans.
conference/model_v2/jobs/j5_benchmark.py:28 — reference-bound timing checkpoints; per-window feature wall, fresh-worker peak RSS, final logical/allocated output disk and index bytes; worker startup/index wall recorded separately.
conference/model_v2/jobs/j5_benchmark.sbatch:17 — forward production J5_FASTA/J5_GENE_BOUNDS, pinned J6 manifest, and tunable J5_PROFILE_HIGH/J5_PROFILE_MEDIAN; same production FASTA reader/GC provenance and gene checks.
conference/model_v2/jobs/j5_full.sbatch:11 — tune SBATCH_ARRAY_INX/SBATCH_TIMELIMIT/SBATCH_MEM_PER_NODE; shard count follows array size or J5_SHARD_COUNT, with contiguous-array validation.
conference/model_v2/tests/test_j5_rbp_context.py:705 — 32 added cases cover malformed/missing/conflicting identity, embedded objects, density selection, GC output/reference reuse, FASTA reader, and ingestion refusal; all 58 original tests retained.
Requested root pytest: 90 passed, 0 skipped, 0 failed (5.68 s); bash -n and scoped git diff --check passed. Full feng suite not run: task-scoped J5 verification.
Remaining evidence/limits: no network lookup, production benchmark, submissions, commits, wet-lab reads or /orcd traversal; fresh workers do not force cold OS/filesystem caches; disk is final artifact size, not temporary peak usage.
