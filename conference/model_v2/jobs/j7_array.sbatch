#!/bin/bash
#SBATCH --job-name=mv2-j7-array
#SBATCH --partition=mit_normal
#SBATCH --array=0-255%16
#SBATCH --time=02:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%A_%a.out
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
J7_ROOT=${J7_ROOT:-$HOME/orcd/scratch/feng/model_v2/j7_context_v2}
J3_ROOT=${J3_ROOT:-$HOME/orcd/scratch/feng/model_v2/j3_genome_v1}
exec /home/dashii/orcd/pool/smt-isef-env/bin/python -m conference.model_v2.j7_context_features shard \
  --input-dir "$J3_ROOT" --exons "$J7_ROOT/canonical_exons.json" \
  --output-dir "$J7_ROOT/shards" --index "${SLURM_ARRAY_TASK_ID}" --count 256
