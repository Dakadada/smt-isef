#!/bin/bash
#SBATCH --job-name=mv2-j5-gate
#SBATCH --partition=mit_normal
#SBATCH --time=12:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# Submit with --dependency=afterok:<FULL_ARRAY_JOB_ID>. Streams all output rows.
# Small CPU/RAM request; hashing ~77GB and counting JSON statuses can take hours.
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
: "${J5_LOOKUP:?Set refreshed lookup TSV}"
: "${J5_LOOKUP_SHA256:?Pin refreshed TSV SHA256}"
args=(--index-dir "${J5_INDEX:-$HOME/orcd/scratch/feng/model_v2/j5_index_v4}"
  --assembly-lookup "$J5_LOOKUP" --lookup-sha256 "$J5_LOOKUP_SHA256"
  --j6-genome-dir "${J5_J6_DIR:-$HOME/orcd/scratch/feng/model_v2/j6_genome_v1}"
  --j6-manifest-sha256 "${J5_J6_MANIFEST_SHA256:-b7b64df4609781af7dbb7aebe9f32091f7e22fca5480bd2371506585a246cdd5}"
  --full-out "${J5_FULL_OUT:-$HOME/orcd/scratch/feng/model_v2/j5_j6_v4}"
  --shard-count "${J5_SHARD_COUNT:-64}"
  --report "${J5_GATE_REPORT:-$HOME/orcd/scratch/feng/model_v2/j5_gate_v4.json}")
if [[ -n ${J5_FASTA:-} ]]; then args+=(--fasta "$J5_FASTA"); fi
if [[ -n ${J5_GENE_BOUNDS:-} ]]; then args+=(--genes "$J5_GENE_BOUNDS"); fi
exec /home/dashii/orcd/pool/smt-isef-env/bin/python -m conference.model_v2.jobs.j5_gate "${args[@]}"
