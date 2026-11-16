#!/bin/bash
#SBATCH --job-name=mv2-j8-bounded
#SBATCH --partition=mit_preemptable
#SBATCH --time=00:20:00
#SBATCH --mem=6G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# No inference. Each verified source is a resumable SQLite transaction;
# each output shard is committed atomically with a SHA-verified .done marker.
# Historical pre-fix bounded footprint (login-node probe, 2026-11-15): ~55 s wall,
# peak RSS 1.6 GB (the J3 parquet row group dominates), ~0.18 GB index.
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
exec /home/dashii/orcd/pool/smt-isef-env/bin/python -m conference.model_v2.j8_feature_join \
  --root "$HOME/orcd/scratch/feng/model_v2" \
  --out "$HOME/orcd/scratch/feng/model_v2/j8_bounded_v2" --bounded --max-index-gib 4
