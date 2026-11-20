#!/bin/bash
#SBATCH --job-name=mv2-j5-sentinel
#SBATCH --partition=mit_preemptable
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# Builds a regional index in this job; no full index dependency.
# Four queries: F18/F14/F10 at chr5:70076584 and the +/-1 Mb exon-7 anchor window.
# This is a regulatory-context smoke test; absent evidence is not a biological NO-GO.
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
: "${J5_LOOKUP:?Set refreshed experiment-verified lookup TSV}"
: "${J5_LOOKUP_SHA256:?Pin refreshed lookup SHA256}"
exec /usr/bin/time -v /home/dashii/orcd/pool/smt-isef-env/bin/python \
  -m conference.model_v2.j5_rbp_context \
  --encode-root "$HOME/orcd/scratch/feng/encode" \
  --assembly-lookup "$J5_LOOKUP" --lookup-sha256 "$J5_LOOKUP_SHA256" \
  --out "${J5_SENTINEL_OUT:-$HOME/orcd/scratch/feng/model_v2/j5_sentinel_fix_v2}" \
  --sentinel --require-query-done
