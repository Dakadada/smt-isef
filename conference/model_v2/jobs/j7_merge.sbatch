#!/bin/bash
#SBATCH --job-name=mv2-j7-merge
#SBATCH --partition=mit_normal
#SBATCH --time=02:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
J7_ROOT=${J7_ROOT:-$HOME/orcd/scratch/feng/model_v2/j7_context_v2}
J3_ROOT=${J3_ROOT:-$HOME/orcd/scratch/feng/model_v2/j3_genome_v1}
# Explicit expected count is required; obtain it from the accepted J3 gate.
: "${J7_EXPECTED_ROWS:?set J7_EXPECTED_ROWS to accepted J3 row count}"
exec /home/dashii/orcd/pool/smt-isef-env/bin/python -m conference.model_v2.j7_context_features merge \
  --input-dir "$J3_ROOT" --exons "$J7_ROOT/canonical_exons.json" \
  --output-dir "$J7_ROOT/shards" --count 256 --expected-rows "$J7_EXPECTED_ROWS"
