#!/usr/bin/env bash
#SBATCH --job-name=feng-armE-ref
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --requeue
#SBATCH --signal=B:SIGTERM@60
#SBATCH --output=feng/logs/%x-%j.out
#SBATCH --error=feng/logs/%x-%j.err
set -euo pipefail
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
PYTHON_BIN="${FENG_ENV_PREFIX:-/home/dashii/orcd/pool/feng-env}/bin/python"
ROOT="${FENG_ARME_ROOT:-${SCRATCH:?SCRATCH is required}/feng/arm_e}"
GTF="${FENG_ARME_GTF:-${SCRATCH}/feng/gencode/v50/gencode.v50.primary_assembly.annotation.gtf}"
FASTA="${FENG_ARME_FASTA:-${SCRATCH}/feng/gencode/v50/GRCh38.primary_assembly.genome.fa}"
REF="${ROOT}/reference/gencode_v50_exonic_union.fa"
PIN="${ROOT}/reference/gencode_v50_exonic_union.pin.json"
CKPT="${ROOT}/reference/gencode_v50_exonic_union.ckpt.json"
DONE="${ROOT}/reference/gencode_v50_exonic_union.done"
mkdir -p "${ROOT}/reference" "${SUBMIT_DIR}/feng/logs"
[[ -s "${GTF}" ]] || { echo "ERROR: GTF missing: ${GTF}" >&2; exit 1; }
[[ -s "${FASTA}" && -s "${FASTA}.fai" ]] || { echo "ERROR: FASTA/FAI missing: ${FASTA}" >&2; exit 1; }
export PYTHONPATH="${SUBMIT_DIR}/feng/src:${PYTHONPATH:-}"
if [[ -s "${REF}" && -s "${PIN}" ]]; then
  "${PYTHON_BIN}" -m feng.titration.arm_e_scan validate-reference --reference "${REF}" --pin "${PIN}"
  touch "${DONE}"
  exit 0
fi
PAYLOAD_PID=""
requeue_handler() { [[ -n "${PAYLOAD_PID}" ]] && kill -TERM "${PAYLOAD_PID}" 2>/dev/null || true; [[ -z "${PAYLOAD_PID}" ]] || wait "${PAYLOAD_PID}" || exit 1; scontrol requeue "${SLURM_JOB_ID:?}" || exit 1; exit 0; }
trap requeue_handler SIGTERM
"${PYTHON_BIN}" -m feng.titration.arm_e_scan build-reference --gtf "${GTF}" --genome-fasta "${FASTA}" --out "${REF}" --pin "${PIN}" --checkpoint "${CKPT}" &
PAYLOAD_PID=$!
wait "${PAYLOAD_PID}"
[[ -s "${REF}" && -s "${PIN}" ]] || { echo "ERROR: reference outputs missing" >&2; exit 1; }
touch "${DONE}"
trap - SIGTERM
