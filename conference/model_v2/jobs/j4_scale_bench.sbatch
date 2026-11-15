#!/usr/bin/env bash
# Submit from repo root: sbatch conference/model_v2/jobs/j4_scale_bench.sbatch
# Short real-Pangolin throughput bench over ~50 real J4 configs, --batch-size {4,16,32}.
# Prints inferences/s and peak GPU memory per sweep; feeds array shard size and --time.
#SBATCH --job-name=mv2-j4-scale-bench
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:40:00
#SBATCH --requeue
#SBATCH --signal=B:SIGTERM@60
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
#SBATCH --error=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.err
set -euo pipefail
cd /home/dashii/smt-isef
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
OUT="${J4_BENCH_DIR:-${HOME}/orcd/scratch/feng/model_v2/j4_scale_bench}"
BENCH_CONFIGS="${J4_BENCH_CONFIGS:-50}"
BENCH_BATCH_SIZES="${J4_BENCH_BATCH_SIZES:-4,16,32}"
BENCH_BUDGET_S="${J4_BENCH_BUDGET_S:-300}"
test -d "${SHARD_DIR}"
mkdir -p "${OUT}"
# A fresh, private marker per invocation prevents stale readiness on requeue.
READY_DIR=$(mktemp -d "${OUT}/.j4-ready.XXXXXX")
export J4_READY_FILE="${READY_DIR}/ready"
trap 'rm -f -- "${J4_READY_FILE}"; rmdir -- "${READY_DIR}"' EXIT

# Sole requeue owner, following feng/scripts/templates/gpu_preemptable.sbatch.
# Python installs checkpoint_on_signal(requeue=False), flushes, and exits 0.
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
# a known child to stop and reap before requeue.
PREEMPT_REQUESTED=0
trap 'PREEMPT_REQUESTED=1' SIGTERM SIGUSR1
"${PYTHON_BIN}" conference/model_v2/jobs/j4_scale_bench.py \
    --shard-dir "${SHARD_DIR}" --out-dir "${OUT}" \
    --configs "${BENCH_CONFIGS}" --batch-sizes "${BENCH_BATCH_SIZES}" \
    --budget-seconds "${BENCH_BUDGET_S}" &
PAYLOAD_PID=$!
trap requeue_handler SIGTERM SIGUSR1
if (( PREEMPT_REQUESTED )); then requeue_handler; fi
wait "${PAYLOAD_PID}"
trap - SIGTERM SIGUSR1
