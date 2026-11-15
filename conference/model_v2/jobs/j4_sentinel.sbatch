#!/usr/bin/env bash
# Submit from repo root: sbatch conference/model_v2/jobs/j4_sentinel.sbatch
# Real Pangolin only; five 18-base footprints, including F18 70076584..70076601.
#SBATCH --job-name=mv2-j4-sentinel
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:15:00
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
OUT="${J4_OUT_DIR:-${HOME}/orcd/scratch/feng/model_v2/j4_sentinel}"
mkdir -p "${OUT}"

# Sole requeue owner, following feng/scripts/templates/gpu_preemptable.sbatch.
# Python installs checkpoint_on_signal(requeue=False), flushes, and exits 0.
PAYLOAD_PID=""
requeue_handler() {
    trap '' SIGTERM SIGUSR1
    if [[ -n "${PAYLOAD_PID}" ]]; then
        kill -TERM "${PAYLOAD_PID}" 2>/dev/null || true
        if ! wait "${PAYLOAD_PID}"; then
            echo 'ERROR: payload checkpoint flush failed; refusing requeue' >&2
            exit 1
        fi
    fi
    if [[ -z "${SLURM_JOB_ID:-}" ]] || ! scontrol requeue "${SLURM_JOB_ID}"; then
        echo 'ERROR: checkpoint saved but requeue failed' >&2
        exit 1
    fi
    exit 0
}
trap requeue_handler SIGTERM SIGUSR1
CHECK_ARGS=()
if [[ "${J4_CHECK_TOLERANCE:-0}" == 1 ]]; then
    CHECK_ARGS+=(--check-tolerance)
fi
"${PYTHON_BIN}" -m conference.model_v2.gpu_perturbation_atlas \
    --smoke --output-dir "${OUT}" --batch-size 2 "${CHECK_ARGS[@]}" &
PAYLOAD_PID=$!
wait "${PAYLOAD_PID}"
trap - SIGTERM SIGUSR1
