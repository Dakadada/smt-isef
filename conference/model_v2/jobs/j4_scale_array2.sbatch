#!/usr/bin/env bash
# J4 scale-out atlas array runner, v2: ONE engine process per array task, which
# loads Pangolin ONCE and drives every config in the shard in-process
# (jobs/j4_scale_driver.py). v1 (jobs/j4_scale_array.sbatch) spawned one process
# per config and paid the measured 23.1 s model load 20,395 times (~130 GPU-h of
# startup); this runner pays it 319 times (~2 GPU-h).
#
# NO --array directive on purpose: the array size is set at submit time from the
# shard ledger, e.g.
#   NSHARDS=$(grep -vc '^#' ~/orcd/scratch/feng/model_v2/j4_configs/j4_config_shards.tsv)
#   sbatch --array=0-$((NSHARDS-1))%8 conference/model_v2/jobs/j4_scale_array2.sbatch
# --time IS fixed here: 64 configs x 3.16 s/config (batch 32 bench, Slurm
# 29169335) + 23 s load = ~4 min of work, so 00:45:00 is ~10x headroom for a
# slow preemptable node and FS contention.
# Variables: J4_SHARD_DIR, J4_ATLAS_DIR, J4_BATCH_SIZE (default 32), J4_GPU_ENV,
#            J4_FASTA.
# Masked rows keep evidence_class sequence_ablation_proxy (set by the engine).
#SBATCH --job-name=mv2-j4-scale-array2
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:45:00
#SBATCH --requeue
#SBATCH --signal=B:SIGTERM@60
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%A_%a.out
#SBATCH --error=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%A_%a.err
set -euo pipefail
cd /home/dashii/smt-isef
: "${SLURM_ARRAY_TASK_ID:?submit with --array=0-<NSHARDS-1>; see header}"
module load cuda/12.9.1
module load cudnn/9.8.0.87-cuda12
module load miniforge/25.11.0-0
GPU_ENV="${J4_GPU_ENV:-/home/dashii/orcd/pool/feng-gpu-env}"
PYTHON_BIN="${GPU_ENV}/bin/python"
test -x "${PYTHON_BIN}"
export PATH="${GPU_ENV}/bin:${PATH}"
export PYTHONPATH="${PWD}/feng/src:${PWD}"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4
export CUBLAS_WORKSPACE_CONFIG=:4096:8
SHARD_DIR="${J4_SHARD_DIR:-${HOME}/orcd/scratch/feng/model_v2/j4_configs/shards}"
ATLAS_DIR="${J4_ATLAS_DIR:-${HOME}/orcd/scratch/feng/model_v2/j4_atlas_v1}"
BATCH_SIZE="${J4_BATCH_SIZE:-32}"
SHARD_FILE="${SHARD_DIR}/$(printf 'shard_%04d.jsonl' "${SLURM_ARRAY_TASK_ID}")"
test -s "${SHARD_FILE}"
mkdir -p "${ATLAS_DIR}"
FASTA_ARGS=()
if [[ -n "${J4_FASTA:-}" ]]; then FASTA_ARGS=(--fasta "${J4_FASTA}"); fi
# A fresh, private marker dir per invocation prevents stale readiness on requeue.
READY_DIR=$(mktemp -d "${ATLAS_DIR}/.j4-ready.XXXXXX")
export J4_READY_FILE="${READY_DIR}/ready"
trap 'rm -f -- "${J4_READY_FILE}"; rmdir -- "${READY_DIR}" 2>/dev/null || true' EXIT

# Sole requeue owner, following feng/scripts/templates/gpu_preemptable.sbatch
# and jobs/j4_sentinel.sbatch. The driver installs
# checkpoint_on_signal(requeue=False) (and run_atlas nests its own per-config
# flush), writes its final checkpoint, and exits 0.
PAYLOAD_PID=""
requeue_handler() {
    trap '' SIGTERM SIGUSR1
    if [[ -n "${PAYLOAD_PID}" ]]; then
        kill -TERM "${PAYLOAD_PID}" 2>/dev/null || true
        if ! wait "${PAYLOAD_PID}"; then
            # Read after the child exits: it may have installed its handler
            # concurrently with TERM. Before readiness it cannot have written
            # atlas work, so a default TERM exit is safe to requeue.
            if [[ -e "${J4_READY_FILE}" ]]; then
                echo 'ERROR: payload checkpoint flush failed; refusing requeue' >&2
                exit 1
            fi
        fi
    fi
    if [[ -z "${SLURM_JOB_ID:-}" ]] || ! scontrol requeue "${SLURM_JOB_ID}"; then
        echo 'ERROR: checkpoint saved but requeue failed' >&2
        exit 1
    fi
    exit 0
}
# Defer handling across spawn/$! assignment, so even an immediate TERM has
# a known child to stop and reap before requeue. Completed configs are already
# durable (.done + .j4_config_sha256), so the requeued task skips them.
PREEMPT_REQUESTED=0
trap 'PREEMPT_REQUESTED=1' SIGTERM SIGUSR1
"${PYTHON_BIN}" conference/model_v2/jobs/j4_scale_driver.py \
    --shard-dir "${SHARD_DIR}" --atlas-dir "${ATLAS_DIR}" \
    --task-id "${SLURM_ARRAY_TASK_ID}" --batch-size "${BATCH_SIZE}" \
    "${FASTA_ARGS[@]+"${FASTA_ARGS[@]}"}" &
PAYLOAD_PID=$!
trap requeue_handler SIGTERM SIGUSR1
if (( PREEMPT_REQUESTED )); then requeue_handler; fi
wait "${PAYLOAD_PID}"
trap - SIGTERM SIGUSR1
