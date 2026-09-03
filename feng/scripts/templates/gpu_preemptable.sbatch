#!/usr/bin/env bash
# =============================================================================
# TEMPLATE — shared feng GPU preemptable sbatch. NOT a runnable job as-is.
#
# Copy this file, then fill every __PLACEHOLDER__ for your arm's payload:
#   __JOB_NAME__     e.g. feng-armE-md
#   __ARRAY_RANGE__  e.g. 1-15%2   (one shard per mer; %2 stays under any cap)
#   __WALLTIME__     e.g. 12:00:00 (<= 2-00:00:00 on mit_preemptable)
#   __PAYLOAD_CMD__  the per-shard GPU command (see the payload block below)
#
# Design (mirrors feng/scripts/armc_scan.sbatch + bench_intarna_scores.sbatch):
#   * mit_preemptable + --requeue + --signal=B:SIGTERM@60 : evictable partition
#     with real GPU headroom; SLURM warns ~60s before eviction.
#   * per-shard .done + skip-if-exists + atomic mv : idempotent across requeues.
#   * checkpoint/requeue trap wired to feng.gpu.checkpoint : a preempted shard
#     flushes its checkpoint and requeues, resuming instead of restarting.
#   * GPU conda env at ~/orcd/pool/feng-gpu-env (created by setup_gpu_env.sh).
#
# GPU partition facts (2026-09-03): mit_normal_gpu = 6h, QOS cap 2 GPU/user;
# mit_preemptable = 2-day, effectively unlimited GPU but REQUIRES this
# checkpoint+requeue discipline. Heavy sweeps -> this template (preemptable).
# =============================================================================
#SBATCH --job-name=__JOB_NAME__
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=__WALLTIME__
#SBATCH --requeue
#SBATCH --signal=B:SIGTERM@60
#SBATCH --array=__ARRAY_RANGE__
#SBATCH --output=feng/logs/%x-%A_%a.out
#SBATCH --error=feng/logs/%x-%A_%a.err

set -euo pipefail

# --- resolved module combo (keep in sync with feng.gpu.env.resolve_modules) ---
module load cuda/12.9.1
module load cudnn/9.8.0.87-cuda12
module load gromacs/2024.1          # drop for deep-splice-only (no MD) jobs
module load miniforge/25.11.0-0

# --- activate the GPU compute env (NOT the CPU feng-env) ----------------------
GPU_ENV_PREFIX="${FENG_GPU_ENV_PREFIX:-${HOME}/orcd/pool/feng-gpu-env}"
if [[ ! -d "${GPU_ENV_PREFIX}" ]]; then
    echo "ERROR: GPU env missing: ${GPU_ENV_PREFIX}" >&2
    echo "       create it once via feng/scripts/setup_gpu_env.sh (review first)" >&2
    exit 1
fi
CONDA_BASE="$(conda info --base)"
# shellcheck disable=SC1091
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${GPU_ENV_PREFIX}"

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
PYTHON_BIN="${GPU_ENV_PREFIX}/bin/python"
export PYTHONPATH="${SUBMIT_DIR}/feng/src:${PYTHONPATH:-}"

# --- output layout: per-shard result + sentinel + checkpoint ------------------
TASK="${SLURM_ARRAY_TASK_ID:?ERROR: run me as a SLURM array job}"
OUT_DIR="${FENG_GPU_OUT_DIR:-${SCRATCH:-}/feng/__JOB_NAME__}"
mkdir -p "${OUT_DIR}" feng/logs
OUT="${OUT_DIR}/shard_${TASK}.parquet"
DONE="${OUT_DIR}/shard_${TASK}.done"
CKPT="${OUT_DIR}/shard_${TASK}.ckpt.json"

# --- idempotent skip: a completed shard is never recomputed (requeue-safe) -----
if [[ -s "${OUT}" && -f "${DONE}" ]]; then
    echo "shard ${TASK} already done; skipping"
    exit 0
fi

# --- preemption trap: SLURM's --signal=B:SIGTERM@60 delivers to THIS batch shell
# ONLY (the `B:` prefix), never to the python payload. So the trap must FORWARD
# the signal to the payload child, WAIT for it to flush its final checkpoint
# within the ~60s grace, and THEN requeue exactly once. The bash trap is the
# SOLE requeue owner: the payload's checkpoint_on_signal MUST use requeue=False
# (flush-only) so requeue does not double-fire.
#
# The payload is responsible for: (1) reading ${CKPT} on start and resuming/
# skipping finished work; (2) writing ${CKPT} on the feng.gpu.checkpoint
# should_checkpoint cadence; (3) installing feng.gpu.checkpoint.checkpoint_on_
# signal(..., requeue=False) so a forwarded SIGTERM flushes one final checkpoint.
PAYLOAD_PID=""
requeue_handler() {
    echo "caught eviction signal; forwarding to payload ${PAYLOAD_PID:-none} to flush" >&2
    if [[ -n "${PAYLOAD_PID}" ]]; then
        kill -TERM "${PAYLOAD_PID}" 2>/dev/null || true
        # bounded by SLURM's grace window; wait for the child's final flush+exit
        wait "${PAYLOAD_PID}" 2>/dev/null || true
    fi
    echo "requeuing job ${SLURM_JOB_ID:-unknown} (single requeue owner: sbatch trap)" >&2
    [[ -n "${SLURM_JOB_ID:-}" ]] && scontrol requeue "${SLURM_JOB_ID}" || true
    exit 0
}
trap requeue_handler SIGTERM SIGUSR1

# --- payload -----------------------------------------------------------------
# Write results to a temp file and atomically mv into place so a preemption
# mid-write never leaves a half-written ${OUT} that the skip check would trust.
# The partial is keyed on the STABLE array task id (NOT ${SLURM_JOB_ID}, which
# changes every requeue and would leak orphaned .partial.<oldid> files); clear
# any stale partial from a prior preempted attempt before writing.
TMP_OUT="${OUT}.partial.${TASK}"
rm -f "${TMP_OUT}"
echo "shard ${TASK}: starting payload -> ${OUT} (ckpt ${CKPT})"

# __PAYLOAD_CMD__ — replace with the arm's GPU scorer, e.g.:
#   "${PYTHON_BIN}" -m feng.titration.gpu_md \
#       --shard "${TASK}" --checkpoint "${CKPT}" --out "${TMP_OUT}"
# The module MUST: (1) load ${CKPT} if present and resume; (2) checkpoint on the
# feng.gpu.checkpoint cadence; (3) install checkpoint_on_signal(requeue=False)
# so a forwarded SIGTERM flushes a final checkpoint; (4) write ${TMP_OUT}.
# RUN IT IN THE BACKGROUND and `wait`, so the trap above can forward the signal
# to it — a foreground payload would never see the forwarded SIGTERM:
#
#   "${PYTHON_BIN}" -m __PAYLOAD_CMD__ \
#       --shard "${TASK}" --checkpoint "${CKPT}" --out "${TMP_OUT}" &
#   PAYLOAD_PID=$!
#   wait "${PAYLOAD_PID}"
echo "ERROR: __PAYLOAD_CMD__ not filled in — this is a template, not a live job" >&2
exit 2

# --- commit: atomic publish + sentinel (unreachable until payload filled) -----
mv -f "${TMP_OUT}" "${OUT}"
touch "${DONE}"
trap - SIGTERM SIGUSR1
echo "shard ${TASK}: done -> ${OUT}"
