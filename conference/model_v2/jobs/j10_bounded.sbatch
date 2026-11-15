#!/bin/bash
#SBATCH --job-name=mv2-j10-bounded
#SBATCH --partition=mit_preemptable
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# 300 already-scored J6/J3 sites plus J9 F18/F14/F10; 3UP8 excluded.
# Final atomic JSON is an idempotent checkpoint. An interrupted bounded pass reruns.
# No folding backend here: only the pure-sequence complementarity channel is
# rescored for decoys; duplex/accessibility nulls are reported blocked.
# Measured (login-node probe, 2026-11-15): 2.9 s wall, peak RSS 1.6 GB, 152 KB JSON.
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
exec /home/dashii/orcd/pool/smt-isef-env/bin/python -m conference.model_v2.j10_robustness \
  --j9 "$HOME/orcd/scratch/feng/model_v2/j9_smn2_v1/features.parquet" \
  --j3 "$HOME/orcd/scratch/feng/model_v2/j3_genome_v1/j3_shard0000_of0256.parquet" \
  --out "$HOME/orcd/scratch/feng/model_v2/j10_bounded_v1/report.json" --limit 300 --seed 1729
