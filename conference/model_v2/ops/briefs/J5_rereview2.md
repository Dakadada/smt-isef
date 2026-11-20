# TASK J5-REREVIEW2 (Codex gpt-6.1-terra, medium, targeted, READ-ONLY except one report).
Verify the findings in conference/model_v2/ops/briefs/J5_REREVIEW_RESULT.md (identity bypass in jobs/j5_encode_lookup.py; benchmark config) are fixed by the J5 fix2 (J5_FIX2_RESULT.md; do not trust it): rerun the identity-spoof repro (fake /biosample-types/ target path, conflicting ontology term) against current code; check path/@id/@type/dataset binding for bypasses; confirm j5_rbp_context.py refuses unverified identity_status.

Then judge the REAL benchmark (Slurm 29280509; read only ops/logs/mv2-j5-profile-29280509.out and small json in ~/orcd/scratch/feng/model_v2/j5_profile_v3): per-window wall 11-14 s even for 217-site windows (fixed per-window cost, fresh worker each), ~70-88 sites/s at ~1000 sites/window, peak RSS 2.1 GB, output ~7 KB/site (~77 GB for 11.0M sites; scratch has ~250 GB free). Is jobs/j5_full.sbatch (64 shards, contiguous windows, 2h, 16 GiB, one persistent worker per shard) sized safely, and is output/requeue/gate handling sound?

Run from repo root: PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 /home/dashii/orcd/pool/smt-isef-env/bin/python -m pytest conference/model_v2/tests/test_j5_rbp_context.py -q

No source edits, submissions, wet-lab files, or /orcd traversal.

Write <=15 lines to conference/model_v2/ops/briefs/J5_REREVIEW2_RESULT.md: verdict ADVANCE/DO-NOT-ADVANCE, and if ADVANCE the exact recommended sbatch override env (SBATCH_ARRAY_INX / SBATCH_TIMELIMIT / SBATCH_MEM_PER_NODE) and the exact gate check to run afterward.
