#!/usr/bin/env bash
# =============================================================================
# Arm E Pangolin H-off sweep — concrete preemptable implementation.
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
#SBATCH --job-name=feng-armE-hoff
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --requeue
#SBATCH --signal=B:SIGTERM@60
#SBATCH --array=1-15%2
#SBATCH --output=feng/logs/%x-%A_%a.out
#SBATCH --error=feng/logs/%x-%A_%a.err

set -euo pipefail

# --- authenticate the committed launcher and select the sealed source tree ----
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
: "${SCRATCH:=$(readlink -f "${HOME}/orcd/scratch" 2>/dev/null || true)}"
export SCRATCH
SEALED_COMMIT="${FENG_ARME_SEALED_COMMIT:?full committed source SHA is required}"
EXPECTED_LAUNCHER_SHA="${FENG_ARME_LAUNCHER_BLOB_SHA256:?launcher blob SHA-256 is required}"
RUN_ID="${FENG_ARME_RUN_ID:?UTC run id like 20260910T120000Z is required}"
[[ "${SEALED_COMMIT}" =~ ^([0-9a-f]{40}|[0-9a-f]{64})$ ]] || {
    echo "ERROR: FENG_ARME_SEALED_COMMIT must be a full lowercase commit SHA" >&2
    exit 1
}
[[ "${EXPECTED_LAUNCHER_SHA}" =~ ^[0-9a-f]{64}$ ]] || {
    echo "ERROR: expected launcher SHA-256 must be 64 lowercase hex characters" >&2
    exit 1
}
[[ "${RUN_ID}" =~ ^[0-9]{8}T[0-9]{6}Z$ ]] || {
    echo "ERROR: FENG_ARME_RUN_ID must be a UTC basic timestamp" >&2
    exit 1
}
ACTUAL_LAUNCHER_SHA="$(sha256sum "$0" | awk '{print $1}')"
[[ "${ACTUAL_LAUNCHER_SHA}" == "${EXPECTED_LAUNCHER_SHA}" ]] || {
    echo "ERROR: launcher SHA-256 mismatch; refusing unauthenticated launcher" >&2
    exit 1
}
RESOLVED_COMMIT="$(git -C "${SUBMIT_DIR}" rev-parse --verify "${SEALED_COMMIT}^{commit}")"
[[ "${RESOLVED_COMMIT}" == "${SEALED_COMMIT}" ]] || {
    echo "ERROR: sealed commit does not resolve exactly" >&2
    exit 1
}
COMMITTED_LAUNCHER_SHA="$(
    git -C "${SUBMIT_DIR}" cat-file blob \
        "${SEALED_COMMIT}:feng/scripts/arm_e_hoff.sbatch" |
        sha256sum | awk '{print $1}'
)"
[[ "${EXPECTED_LAUNCHER_SHA}" == "${COMMITTED_LAUNCHER_SHA}" ]] || {
    echo "ERROR: launcher hash is not the sealed commit blob" >&2
    exit 1
}

TASK="${SLURM_ARRAY_TASK_ID:?ERROR: run me as a SLURM array job}"
ARRAY_JOB_ID="${SLURM_ARRAY_JOB_ID:?SLURM_ARRAY_JOB_ID is required}"
LENGTH="$(printf '%02d' "$((2*TASK))")"
OUT_DIR="${FENG_GPU_OUT_DIR:-${SCRATCH:?SCRATCH is required}/feng/arm_e/hoff}"
mkdir -p "${OUT_DIR}" "${OUT_DIR}/.locks" "${SUBMIT_DIR}/feng/logs"

# Hold this lock until task exit, not merely through preparation.  This prevents
# duplicate array submissions for one length from superseding or writing while
# an authenticated task is scoring.
LOCK_PATH="${OUT_DIR}/.locks/pangolin_hoff_L${LENGTH}.lock"
exec {LENGTH_LOCK_FD}>"${LOCK_PATH}"
if ! flock -n "${LENGTH_LOCK_FD}"; then
    echo "ERROR: another H-off task owns length L${LENGTH}" >&2
    exit 1
fi

SOURCE_PARENT="${SCRATCH:?SCRATCH is required}/feng/arm_e/source/${SEALED_COMMIT}/${ARRAY_JOB_ID}/${TASK}"
mkdir -p "${SOURCE_PARENT}"
SOURCE_ROOT="$(mktemp -d "${SOURCE_PARENT}/snapshot.XXXXXX")"
git -C "${SUBMIT_DIR}" archive --format=tar "${SEALED_COMMIT}" feng/src \
    | tar -xf - -C "${SOURCE_ROOT}"
AUTH_BOOTSTRAP="${SOURCE_ROOT}/feng/src/feng/source_auth.py"
[[ -f "${AUTH_BOOTSTRAP}" ]] || {
    echo "ERROR: committed source-auth bootstrap missing from source snapshot" >&2
    exit 1
}
EXPECTED_BOOTSTRAP_SHA="$(
    git -C "${SUBMIT_DIR}" cat-file blob "${SEALED_COMMIT}:feng/src/feng/source_auth.py" |
        sha256sum | awk '{print $1}'
)"
ACTUAL_BOOTSTRAP_SHA="$(sha256sum "${AUTH_BOOTSTRAP}" | awk '{print $1}')"
[[ "${ACTUAL_BOOTSTRAP_SHA}" == "${EXPECTED_BOOTSTRAP_SHA}" ]] || {
    echo "ERROR: source-auth bootstrap SHA-256 mismatch" >&2
    exit 1
}

# --- module environment; no feng import is allowed before bootstrap verify ----
module load cuda/12.9.1
module load cudnn/9.8.0.87-cuda12
module load miniforge/25.11.0-0
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
PYTHON_BIN="${GPU_ENV_PREFIX}/bin/python"

# Build the manifest from Git's sealed object store, then re-hash every staged
# file before any feng import.  All later parent/worker processes repeat this
# same full-tree check in-process immediately before runpy imports feng.
SOURCE_MANIFEST="${SOURCE_ROOT}/source-manifest.json"
SOURCE_TREE_SHA="$(
    "${PYTHON_BIN}" "${AUTH_BOOTSTRAP}" create-manifest \
        --git-repo "${SUBMIT_DIR}" --commit "${SEALED_COMMIT}" \
        --source-root "${SOURCE_ROOT}" --manifest "${SOURCE_MANIFEST}"
)"
[[ "${SOURCE_TREE_SHA}" =~ ^[0-9a-f]{64}$ ]] || {
    echo "ERROR: authenticated source-tree identity is malformed" >&2
    exit 1
}
export FENG_ARME_SOURCE_ROOT="${SOURCE_ROOT}"
export FENG_ARME_SOURCE_MANIFEST="${SOURCE_MANIFEST}"
export FENG_ARME_SOURCE_TREE_SHA256="${SOURCE_TREE_SHA}"
export FENG_ARME_SEALED_COMMIT="${SEALED_COMMIT}"
export FENG_ARME_LAUNCHER_BLOB_SHA256="${ACTUAL_LAUNCHER_SHA}"
export FENG_REQUIRE_SOURCE_AUTH=1
export PYTHONPATH="${SOURCE_ROOT}/feng/src"
export PYTHONSAFEPATH=1
export PYTHONNOUSERSITE=1
export PYTHONDONTWRITEBYTECODE=1

# --- output layout: per-shard result + sentinel + checkpoint ------------------
SCAN_DIR="${FENG_ARME_SCAN_DIR:-${SCRATCH}/feng/arm_e/scan}"
MANIFEST="${SCAN_DIR}/sites_L${LENGTH}.tsv"
SCAN_BURDEN="${SCAN_DIR}/burden_L${LENGTH}.json"
GENOME_FASTA="${FENG_ARME_FASTA:-${SCRATCH}/feng/gencode/v50/GRCh38.primary_assembly.genome.fa}"
GENOME_PIN="${FENG_ARME_PIN:-${SCRATCH}/feng/arm_e/reference/gencode_v50_exonic_union.pin.json}"
OUT="${OUT_DIR}/pangolin_hoff_L${LENGTH}.json"
DONE="${OUT_DIR}/pangolin_hoff_L${LENGTH}.done"
CKPT="${OUT_DIR}/pangolin_hoff_L${LENGTH}.ckpt.json"

# --- crash-recoverable Path-A supersession transaction ------------------------
PREPARED="${OUT_DIR}/.prepared-${RUN_ID}-L${LENGTH}"
SUPERSEDED_DIR="${OUT_DIR}/superseded/${RUN_ID}/L${LENGTH}"
if [[ ! -f "${PREPARED}" ]]; then
    mkdir -p "${SUPERSEDED_DIR}"
    for prior in "${DONE}" "${OUT}" "${CKPT}"; do
        destination="${SUPERSEDED_DIR}/$(basename "${prior}")"
        if [[ -e "${prior}" && -e "${destination}" ]]; then
            echo "ERROR: both live and superseded copies exist: ${prior}" >&2
            exit 1
        elif [[ -e "${prior}" ]]; then
            mv -- "${prior}" "${destination}"
        elif [[ -e "${destination}" ]]; then
            # Recovery after preemption between an earlier move and marker publish.
            :
        fi
    done
    PREPARED_TMP="$(mktemp "${OUT_DIR}/.${RUN_ID}-L${LENGTH}.prepared.XXXXXX")"
    printf '%s\n' "${SEALED_COMMIT}" >"${PREPARED_TMP}"
    mv -- "${PREPARED_TMP}" "${PREPARED}"
fi

# The marker must predate all current-run output.  Because the length lock is
# held for this process's lifetime, a duplicate launch cannot cross this point.
[[ -f "${PREPARED}" ]] || {
    echo "ERROR: supersession transaction did not publish ${PREPARED}" >&2
    exit 1
}
cd "${SOURCE_ROOT}"
HOFF_COMMAND=(
    "${PYTHON_BIN}" "${AUTH_BOOTSTRAP}" run-module
    feng.titration.arm_e_hoff --
)

# --- idempotent skip after this authenticated run has been prepared -----------
if [[ -s "${OUT}" && -f "${DONE}" ]] &&
   "${HOFF_COMMAND[@]}" \
        --manifest "${MANIFEST}" --checkpoint "${CKPT}" --out "${OUT}" \
        --genome-fasta "${GENOME_FASTA}" --gpu-python "${PYTHON_BIN}" \
        --genome-pin "${GENOME_PIN}" \
        --scan-burden "${SCAN_BURDEN}" \
        --chunk-size "${FENG_ARME_HOFF_CHUNK_SIZE:-32}" --validate-output; then
    echo "shard ${TASK} already done and authenticated; skipping"
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
        if ! wait "${PAYLOAD_PID}"; then
            echo "ERROR: payload failed while handling checkpoint signal" >&2
            exit 1
        fi
    fi
    echo "requeuing job ${SLURM_JOB_ID:-unknown} (single requeue owner: sbatch trap)" >&2
    if [[ -z "${SLURM_JOB_ID:-}" ]] || ! scontrol requeue "${SLURM_JOB_ID}"; then
        echo "ERROR: checkpoint saved but requeue failed" >&2
        exit 1
    fi
    exit 0
}
trap requeue_handler SIGTERM SIGUSR1

# --- payload -----------------------------------------------------------------
# Write results to a temp file and atomically mv into place so a preemption
# mid-write never leaves a half-written ${OUT} that the skip check would trust.
# The partial is keyed on the STABLE array task id (NOT ${SLURM_JOB_ID}, which
# changes every requeue and would leak orphaned .partial.<oldid> files); clear
# any stale partial from a prior preempted attempt before writing.
echo "shard ${TASK}: starting payload -> ${OUT} (ckpt ${CKPT})"

[[ -s "${MANIFEST}" ]] || { echo "ERROR: site manifest missing: ${MANIFEST}" >&2; exit 1; }
[[ -s "${SCAN_BURDEN}" ]] || { echo "ERROR: scan burden missing: ${SCAN_BURDEN}" >&2; exit 1; }
[[ -s "${GENOME_FASTA}" && -s "${GENOME_FASTA}.fai" ]] || { echo "ERROR: genome FASTA/FAI missing" >&2; exit 1; }
[[ -s "${GENOME_PIN}" ]] || { echo "ERROR: genome source pin missing" >&2; exit 1; }
"${HOFF_COMMAND[@]}" \
    --manifest "${MANIFEST}" --checkpoint "${CKPT}" --out "${OUT}" \
    --genome-fasta "${GENOME_FASTA}" --gpu-python "${PYTHON_BIN}" \
    --genome-pin "${GENOME_PIN}" \
    --scan-burden "${SCAN_BURDEN}" \
    --chunk-size "${FENG_ARME_HOFF_CHUNK_SIZE:-32}" &
PAYLOAD_PID=$!
wait "${PAYLOAD_PID}"
[[ -s "${OUT}" ]] || { echo "ERROR: H-off output missing: ${OUT}" >&2; exit 1; }
"${HOFF_COMMAND[@]}" \
    --manifest "${MANIFEST}" --checkpoint "${CKPT}" --out "${OUT}" \
    --genome-fasta "${GENOME_FASTA}" --gpu-python "${PYTHON_BIN}" \
    --genome-pin "${GENOME_PIN}" \
    --scan-burden "${SCAN_BURDEN}" \
    --chunk-size "${FENG_ARME_HOFF_CHUNK_SIZE:-32}" --validate-output

# --- commit: atomic publish + sentinel (unreachable until payload filled) -----
touch "${DONE}"
trap - SIGTERM SIGUSR1
echo "shard ${TASK}: done -> ${OUT}"
