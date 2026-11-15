# J4/F2 result — implementation ready; broader gates blocked
- `gpu_perturbation_atlas.py:142`: batched deterministic Pangolin adapter; authenticated package, model/weight/code hashes; existing channels 1/4/7/10 (splice-site probabilities).
- `gpu_perturbation_atlas.py:227`: bounded batches, recursive OOM backoff, explicit null failure rows; no partial ensemble averages.
- `gpu_perturbation_atlas.py:336`: content-addressed reference cache; strand-correct N-ablation grid; raw_model_delta, never ΔPSI/occupancy.
- `gpu_perturbation_atlas.py:369`: deterministic sharding, immutable reference snapshots on resume, atomic parquet, checksummed manifests/checkpoints, SIGTERM flush.
- `gpu_perturbation_atlas.py:470`: optional real CPU/GPU tolerance comparison; no real inference outside Slurm.
- `tests/test_gpu_perturbation_atlas.py:1`: stub-only CPU tests; adapter channels, coordinates, shards, corruption, atomic failure, OOM, SIGTERM and requeue ordering.
- `jobs/j4_sentinel.sbatch:1`: 15-minute mit_preemptable GPU sentinel; shell alone owns requeue; Python flushes with requeue=False.
- Scope exception: this explicitly requested result brief is additional to the owner file list. No other source files edited; no commit or job submitted.

## Exact final pytest results
- J4 (`smt-isef-env`): **26 passed, 1 skipped, 0 failed** (1.40s).
- Skip: real CPU/GPU tolerance requires a Slurm GPU allocation; no login-node inference.
- Full `conference/model_v2/tests/` (`smt-isef-env`): **129 passed, 1 skipped, 1 failed** (2.96s).
- Unfixed: `tests/test_j3_scoring.py:64::test_antiparallel_complementarity_and_wobble`; expected one mismatch, implementation returns two. Concurrent J3 work, outside scope.
- Full `feng/tests/` (standing gate's `feng-env`, its binaries on PATH): **584 passed, 2 skipped, 1 failed** (222.94s).
- Unfixed: `feng/tests/test_cleanroom.py:23::test_implementation_has_no_legacy_reference[source_file42]`; existing `feng/src/feng/source_auth.py` contains forbidden `src/`. feng edits prohibited.
- Earlier combined run in requested smt-isef-env: 667 passed, 6 skipped, 15 failed; missing feng dependencies/binaries. Authoritative feng-env rerun above supersedes that environment check.
- Earlier conference runs hit OpenBLAS thread exhaustion and temporarily missing concurrent J3 module; child test now caps threads, final rerun above supersedes them.
- Test commands: `PYTHONPATH=feng/src <ENV>/bin/python -m pytest <TEST_PATH> -q -rfs`; set OPENBLAS_NUM_THREADS/OMP_NUM_THREADS/MKL_NUM_THREADS/NUMEXPR_NUM_THREADS=1.

## Sentinel (human submits, from repository root)
```bash
sbatch conference/model_v2/jobs/j4_sentinel.sbatch
# Optional real CPU/GPU comparison, in the same GPU allocation:
sbatch --export=ALL,J4_CHECK_TOLERANCE=1,J4_OUT_DIR=$HOME/orcd/scratch/feng/model_v2/j4_tolerance conference/model_v2/jobs/j4_sentinel.sbatch
```
- Default output: `~/orcd/scratch/feng/model_v2/j4_sentinel/{reference_scores,ablation_maps}/`; `shard_0.done` is a checksummed JSON manifest, not a touch marker.
- Expected smoke: 4 reference rows + 20 ablation rows (5 footprints × 4 tissues); F18 70076584–70076601 included; donor readout 70076574.
- Verified on disk: FASTA smoke context = 10,101 bases; F18 = 18 bases; GPU python resolves to `/orcd/pool/007/dashii/feng-gpu-env/bin/python3.11`; package hash matches committed seed.
- `bash -n` passed; no real model scores produced. GPU env has torch/Pangolin/pyarrow but no pytest; sentinel invokes the module directly.

## Open questions / limitations
- GPU smoke, numerical tolerance and accelerator OOM behavior remain unverified on hardware; CPU stub results are not model evidence.
- Config supplies prespecified splice sites and one absolute tiling interval; automatic GTF site discovery and genome-wide scheduling are outside this engine.
- Startup dependency/weight failures abort without .done; per-record inference failures are null rows and yield complete_with_failures plus nonzero CLI exit.
- Resume requires unchanged input/model/source/batch plan; use a new output directory for changed configurations or explicit retries of terminal failed rows.
- Materialized task plans and one checkpoint part per batch need sizing before a genome-wide run. No sealed data, outcome-based sampling, or scientific gate changes.
