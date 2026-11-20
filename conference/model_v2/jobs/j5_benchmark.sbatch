#!/bin/bash
#SBATCH --job-name=mv2-j5-profile
#SBATCH --partition=mit_preemptable
#SBATCH --time=01:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# One persistent worker, 40 consecutive windows of one production shard.
# Export the same J5_FASTA/J5_GENE_BOUNDS as production. No production outputs.
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
: "${J5_LOOKUP:?Set refreshed ENCODE lookup TSV}"
: "${J5_LOOKUP_SHA256:?Pin refreshed table SHA256}"
args=(--index-dir "${J5_INDEX:-$HOME/orcd/scratch/feng/model_v2/j5_index_v4}"
  --out "${J5_PROFILE_OUT:-$HOME/orcd/scratch/feng/model_v2/j5_profile_v4}"
  --assembly-lookup "$J5_LOOKUP" --lookup-sha256 "$J5_LOOKUP_SHA256"
  --j6-genome-dir "${J5_J6_DIR:-$HOME/orcd/scratch/feng/model_v2/j6_genome_v1}"
  --j6-manifest-sha256 "${J5_J6_MANIFEST_SHA256:-b7b64df4609781af7dbb7aebe9f32091f7e22fca5480bd2371506585a246cdd5}"
  --persistent --windows "${J5_PROFILE_WINDOWS:-40}"
  --shard-id "${J5_PROFILE_SHARD:-0}" --shard-count "${J5_SHARD_COUNT:-64}"
  --offset "${J5_PROFILE_OFFSET:-0}")
if [[ -n ${J5_FASTA:-} ]]; then args+=(--fasta "$J5_FASTA"); fi
if [[ -n ${J5_GENE_BOUNDS:-} ]]; then args+=(--genes "$J5_GENE_BOUNDS"); fi
exec /usr/bin/time -v /home/dashii/orcd/pool/smt-isef-env/bin/python \
  -m conference.model_v2.jobs.j5_benchmark "${args[@]}"
