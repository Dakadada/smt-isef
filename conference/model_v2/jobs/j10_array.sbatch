#!/bin/bash
#SBATCH --job-name=mv2-j10-offtarget
#SBATCH --partition=mit_normal
#SBATCH --array=0-255%64
#SBATCH --time=00:20:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%A_%a.out
# Recorded J3 size test: 10,837 rows, 642.72s, 1,231,396 KiB RSS,
# 9,814,128 filesystem-input blocks (includes scanning J6 + Vienna computation).
# J10 reads one ~42,979-row scored shard twice for SHA plus one streaming pass.
# Estimate 1-5 min / <1 GiB per task; request 20 min / 2 GiB conservatively.
# No isolated J3 parquet read timing available: estimate, not a new measurement.
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
exec /usr/bin/time -v /home/dashii/orcd/pool/smt-isef-env/bin/python -m conference.model_v2.j10_offtarget_profile partial \
  --shard-index "${SLURM_ARRAY_TASK_ID}" --shard-count 256 \
  --output-dir ~/orcd/scratch/feng/model_v2/j10_offtarget_v2
