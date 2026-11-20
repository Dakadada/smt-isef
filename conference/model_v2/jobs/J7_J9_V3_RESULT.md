# J7/J9 v3 — UNREVIEWED; production rerun BLOCKED

## Changes
- `conference/model_v2/j7_context_features.py:140`: canonical exon-map provenance includes the same whole-implementation SHA256 used by shard identity.
- `conference/model_v2/j7_context_features.py:159,345`: shard reads and merge independently reject missing/stale builder SHA; existing map hashes remain bound into shard and merge commitments.
- `conference/model_v2/j9_smn2_prediction.py:1335,1370`: stage a complete generation, verify parquet hashes/bytes/rows, write its manifest last, fsync files, then atomically replace `.current`; serialize writers with flock. Prior generations are retained. Existing flat datasets require a fresh root.
- `conference/model_v2/j9_smn2_prediction.py:1349`: `read_dataset` pins one generation and verifies manifest hashes before returning tables. README documents the reader contract. Compatibility aliases remain, but separate opens of them across a publication are NOT a snapshot.
- J9 allowlist adds exactly `j9_smn2_v3`; no broadened source roots. No sbatch scripts changed.
- `tests/test_j7_context_features.py:102` and `tests/test_j7_j9_j10_regressions.py:32`: synthetic map fixtures now authenticate the builder; the implementation-change test rebuilds its synthetic map for each simulated implementation.
- `tests/test_j7_j9_v3_regressions.py:19`: builder provenance, independent shard/merge rejection, interrupted staging/pointer replacement, concurrent publication, snapshot retention, integrity rejection, and v3 allowlist coverage.

## Tests
Before implementation: **8 failed, 0 passed, 0 skipped**, as expected (missing builder SHA/checks, interrupted publication, missing snapshot reader). The v3 allowlist regression separately failed before its fix (**1 failed, 8 passed, 0 skipped**). Two additional concurrency/fault-injection checks were added after implementation.
Final focused per-file results: **141 passed, 0 skipped, 0 failed**:
- `test_j7_context_features.py`: 29 passed.
- `test_j9_smn2_prediction.py`: 58 passed.
- `test_j7_j9_j10_regressions.py`: 43 passed.
- `test_j7_j9_v3_regressions.py`: 11 passed.
Command for each file: `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 /home/dashii/orcd/pool/smt-isef-env/bin/python -m pytest conference/model_v2/tests/<file> -q -p no:cacheprovider`.
No full suite run, per explicit OOM constraint. Shell syntax and scoped `git diff --check` passed. Initial edit command used unavailable `python`; it made no edits and was rerun with `python3`.

## Implementation identities
- J7 SHA256: `a419389b637251f7ef34a1ff28b24f3e241f84e8735427c5563568125f3d2c79`.
- J9 SHA256: `c4ae2025fcbafefd0fabb58745804c2cd4c0573ff4ad3d9b7749a94bf163ecaf`.
These are source hashes, **not production artifact certification**.

## Submission attempts and blocker
Both roots were absent at preflight:
- `/home/dashii/orcd/scratch/feng/model_v2/j7_context_v3`
- `/home/dashii/orcd/scratch/feng/model_v2/j9_smn2_v3`
Attempted J7 builder once with `J7_ROOT` above, `J3_ROOT=/home/dashii/orcd/scratch/feng/model_v2/j3_genome_v1`, `J7_EXPECTED_ROWS=11002701`, `sbatch --parsable --partition=mit_preemptable --export=ALL jobs/j7_build_exons.sbatch` (full repository-relative script path used).
Attempted J9 once with `J9_OUT_DIR` above, `J9_ATLAS_DIR=/home/dashii/orcd/scratch/feng/model_v2/j4_atlas_v1/SMN2_exon7_issn1_fine`, `sbatch --parsable --partition=mit_preemptable --export=ALL jobs/j9_smn2.sbatch`.
Both bounded submission commands timed out (exit **124**) without a job ID. **Acceptance is unconfirmed; check accounting before any retry.** No retries, array submission, or merge submission were performed.
- `squeue`: `Unable to contact slurm controller (connect failure)`; another bounded query reported `Invalid user: dashii`.
- `sacct`: address family 0 unsupported; persistent connection to `slurm003:7819` failed; `Problem talking to the database: Resource temporarily unavailable`.
- At 2026-11-16T22:19:05Z, exact-path stat checks found no v3 canonical exon map, J7 merge marker, or J9 manifest.
- No Slurm IDs or completion states can be reported. No production logs can be attributed without IDs. Completion polling and output-SHA verification are blocked by unavailable scheduler access.

## Resume / limitations
First reconcile both timed-out requests in accounting to avoid duplicate jobs. If neither was accepted, submit the builder and J9 into the named fresh roots; submit J7 array with `--array=0-255%3 --dependency=afterok:<builder>` and merge with `--dependency=afterok:<array>`, all on `mit_preemptable` with `--export=ALL`. The array cap allows at most 6G alongside J9's 4G; account for other active jobs before launching. Existing scripts already use `--requeue`, and J7 shards are restartable.
Poll bounded squeue/sacct every 3–5 minutes; inspect failed log tails and retry at most once. After completion, authenticate canonical-map builder SHA, all merge commitments and 11,002,701 rows; verify J9 through `read_dataset` and compare its recorded source SHA against the current implementation. No v2 root was overwritten.
J7 v3 and J9 v3 remain **UNREVIEWED**. Production outputs and their SHAs are **unverified**. Separate adversarial review remains required. Generation retention consumes disk over repeated J9 publications; no automatic cleanup was introduced. File fsync is used, but power-loss directory durability is not claimed. No staging, commit, push, or unrelated edits performed.
