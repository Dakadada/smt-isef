J4-FIX — all three requested code fixes implemented; real GPU rerun remains unverified.
Paths below are relative to conference/model_v2/.
gpu_perturbation_atlas.py:232 — resolve source paths and repository root consistently before relative-path provenance hashing; file hashing/model/artifact paths canonicalized too.
gpu_perturbation_atlas.py:65 — require all 25 GRCh38 primary chromosome lengths; reject wrong/missing/duplicate/malformed index entries before constructing Pangolin, for default/override FASTA and smoke/config modes.
gpu_perturbation_atlas.py:560 — bind canonical FASTA/index paths, index SHA-256, and primary contig signature into scorer provenance, cache identity, manifests, and Parquet metadata.
jobs/j4_sentinel.sbatch:31 — fresh readiness marker per invocation; check after child exit, allow pre-handler termination to requeue, refuse post-readiness flush failures, defer signals across spawn/PID assignment.
tests/test_gpu_perturbation_atlas.py:251 — real delayed-numpy-import termination, successful flush, failed flush, and readiness/TERM race exercise the actual shell with mocked scontrol.
tests/test_gpu_perturbation_atlas.py:435 — symlinked-root constructor regression; :483 wrong-reference matrix; :510 provenance binding, alias stability, and changed-index cache refusal.
TDD: initial targeted run reproduced 11 failures / 2 passes; final requested gate: 47 passed, 1 skipped, 0 failed (1.60s). bash -n passed.
Command: PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 /home/dashii/orcd/pool/smt-isef-env/bin/python -m pytest conference/model_v2/tests/test_gpu_perturbation_atlas.py -q
Live default .fai validated read-only; SHA-256 a2c323ea4cff34d7123ace4578f7e122b2d2f5a22f40dc23eb8b97d17723d169.
Length source: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
Unfixed requested defects: none. No commit, submission, real model inference, or unrelated edits; conference scope/interpreter/test command override AGENTS.md feng defaults.
Limitations: contig-length authentication is not a whole-genome sequence checksum; a forged/stale matching index is not excluded. Existing sequence-slice hashes remain in run identity.
Skipped: real CPU/GPU tolerance requires a Slurm GPU allocation. Sentinel 29165682 was not rerun; real weight loading and GPU behavior remain unverified.
