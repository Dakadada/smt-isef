#!/bin/bash
#SBATCH --job-name=mv2-j5-full
#SBATCH --partition=mit_preemptable
#SBATCH --array=0-63%4
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%A_%a.out
# Conservative defaults exceed 1.5x prior 3.2h / 2059MiB measurements.
# Review persistent_profile.json before submission; request at least its
# recommended_seconds/recommended_mem_mib (already include 1.5x margin).
# The job refuses an undersized or mismatched profile before processing windows.
# Unmeasured chromosome peaks/contention remain a limit of sampled projections.
# Tune submission via SBATCH_ARRAY_INX (e.g. 0-127%4), SBATCH_TIMELIMIT and
# SBATCH_MEM_PER_NODE; Slurm's SBATCH_* environment overrides these directives.
# J5_SHARD_COUNT defaults to SLURM_ARRAY_TASK_COUNT; use a complete 0-based array.
# Run j5_benchmark.sbatch with the SAME J5_FASTA/J5_GENE_BOUNDS/index first.
# Requeue verifies existing .done hashes and resumes at completed J6 windows.
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
: "${J5_LOOKUP:?Set the refreshed experiment-verified lookup TSV}"
: "${J5_LOOKUP_SHA256:?Pin the refreshed table SHA256}"
shard_count=${J5_SHARD_COUNT:-${SLURM_ARRAY_TASK_COUNT:-64}}
if [[ ! $shard_count =~ ^[1-9][0-9]*$ ]] ||
   (( ${SLURM_ARRAY_TASK_MIN:-0} != 0 || ${SLURM_ARRAY_TASK_STEP:-1} != 1 ||
      ${SLURM_ARRAY_TASK_COUNT:-$shard_count} != shard_count ||
      ${SLURM_ARRAY_TASK_MAX:-$((shard_count-1))} != shard_count-1 )); then
  echo 'J5 requires a complete contiguous 0-based array matching J5_SHARD_COUNT' >&2; exit 2
fi
args=(--index-dir "${J5_INDEX:-$HOME/orcd/scratch/feng/model_v2/j5_index_v4}"
      --out "${J5_FULL_OUT:-$HOME/orcd/scratch/feng/model_v2/j5_j6_v4}/shard-${SLURM_ARRAY_TASK_ID}"
      --assembly-lookup "$J5_LOOKUP" --lookup-sha256 "$J5_LOOKUP_SHA256"
      --j6-genome-dir "${J5_J6_DIR:-$HOME/orcd/scratch/feng/model_v2/j6_genome_v1}"
      --j6-manifest-sha256 "${J5_J6_MANIFEST_SHA256:-b7b64df4609781af7dbb7aebe9f32091f7e22fca5480bd2371506585a246cdd5}"
      --shard-id "$SLURM_ARRAY_TASK_ID" --shard-count "$shard_count" --require-query-done)
if [[ -n ${J5_FASTA:-} ]]; then args+=(--fasta "$J5_FASTA"); fi
if [[ -n ${J5_GENE_BOUNDS:-} ]]; then args+=(--genes "$J5_GENE_BOUNDS"); fi
# Read the actual time limit, including sbatch command-line/environment overrides.
job_details=$(scontrol show job "${SLURM_JOB_ID:?Slurm job ID required}" -o)
if [[ $job_details =~ (^|[[:space:]])TimeLimit=([^[:space:]]+) ]]; then
  allocation_time=${BASH_REMATCH[2]}
else
  echo 'Cannot determine allocated Slurm TimeLimit' >&2; exit 2
fi
profile_args=(--index-dir "${J5_INDEX:-$HOME/orcd/scratch/feng/model_v2/j5_index_v4}"
  --assembly-lookup "$J5_LOOKUP" --lookup-sha256 "$J5_LOOKUP_SHA256"
  --j6-manifest-sha256 "${J5_J6_MANIFEST_SHA256:-b7b64df4609781af7dbb7aebe9f32091f7e22fca5480bd2371506585a246cdd5}"
  --shard-count "$shard_count"
  --check-allocation "${J5_PROFILE_OUT:-$HOME/orcd/scratch/feng/model_v2/j5_profile_v4}/persistent_profile.json"
  --allocation-time "$allocation_time"
  --allocation-mib "${SLURM_MEM_PER_NODE:?Slurm memory allocation required}")
if [[ -n ${J5_FASTA:-} ]]; then profile_args+=(--fasta "$J5_FASTA"); fi
if [[ -n ${J5_GENE_BOUNDS:-} ]]; then profile_args+=(--genes "$J5_GENE_BOUNDS"); fi
/home/dashii/orcd/pool/smt-isef-env/bin/python -m conference.model_v2.jobs.j5_benchmark "${profile_args[@]}"
exec /usr/bin/time -v /home/dashii/orcd/pool/smt-isef-env/bin/python \
  -m conference.model_v2.j5_rbp_context "${args[@]}"
