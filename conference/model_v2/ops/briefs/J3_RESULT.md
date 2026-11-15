# J3_RESULT — 2026-11-15
J3 implemented and tested; repository gate remains DO-NOT-ADVANCE for the existing feng failure below.
- j3_scoring.py:76 — antiparallel WC/GU/mismatch descriptors, ASO-oriented 1-based positions.
- j3_scoring.py:226 — independent status/reason/evidence channels; failures stay explicit nulls.
- jobs/j3_vienna.py:9 — ViennaRNA 2.7.2 duplexfold and exact whole-footprint RNAplfold probability.
- j3_scoring.py:69 — normalized (ASO,target,left/right flank) SHA256; config/tool namespaced SQLite cache.
- j3_scoring.py:156 — transactional per-channel checkpoints; completed channels survive interruption.
- j3_scoring.py:379 — deterministic site-ID sharding, input SHA validation, atomic parquet + .done SHA256.
- jobs/j3_sentinel.sbatch:1 — mit_quicktest, 00:10:00, one CPU, --requeue and --requeue-safe; bash -n passed.
- tests/test_j3_scoring.py:60 — duplex oracle, failure/null preservation, independence, cache, shard union, recovery, CLI.
J3 pytest: 33 passed, 0 skipped, 0 failed (OMP_NUM_THREADS=1; smt-isef-env python).
model_v2 pytest: 151 passed, 1 skipped, 0 failed (OMP_NUM_THREADS=1 PYTHONPATH=feng/src; smt-isef-env python).
feng gate: 581 passed, 5 skipped, 1 failed (required command, feng-env python; run once).
Unfixed/outside ownership: test_cleanroom.py:23 rejects existing 'src/' strings in unchanged source_auth.py:31,314,400.
One unrestricted model_v2 rerun aborted in PyArrow (exit 134); the one-thread rerun above passed.
Exact J3 test command: OMP_NUM_THREADS=1 /home/dashii/orcd/pool/smt-isef-env/bin/python -m pytest conference/model_v2/tests/test_j3_scoring.py -q
Exact submit command (NOT executed): sbatch /home/dashii/smt-isef/conference/model_v2/jobs/j3_sentinel.sbatch
Input: ~/orcd/scratch/feng/model_v2/j6_sentinel (22 committed shards / 6,145 rows at metadata inspection).
Output: ~/orcd/scratch/feng/model_v2/j3_sentinel; cache j3_cache.sqlite; parquet commit marker has channel counts.
Assumption: "unmodified RNA; chemistry NOT modelled"; energy is not potency or occupancy.
Limitation: duplex minimum can lie in the flank; no footprint anchor. Accessibility is a local-sequence proxy.
RNAplfold W=120/L=60 clipped to context; 80 nt flanks clipped at contig boundaries; no J2 grid approximation/join.
Transient failures are uncached; remove the affected J3 .done marker to rebuild a committed failed shard and retry.
Slurm/shared-filesystem recovery and full sentinel runtime remain untested because submission was prohibited.
Scope: task conference paths/interpreter override repo feng defaults; no feng changes, commit, or job submission.
