#!/bin/bash
#SBATCH --job-name=mv2-j10-offtarget-merge
#SBATCH --partition=mit_normal
#SBATCH --time=00:20:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# Submit with --dependency=afterok:ARRAY_JOB_ID (exact commands in j10_RESULT.md).
# All 256 partials, both checksums, common J3 config and 11,002,701 raw rows required.
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
exec /usr/bin/time -v /home/dashii/orcd/pool/smt-isef-env/bin/python -m conference.model_v2.j10_offtarget_profile merge \
  --shard-count 256 --expected-rows 11002701
