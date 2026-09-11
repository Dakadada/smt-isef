#!/usr/bin/env bash
#SBATCH --job-name=feng-armE-scan
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --requeue
#SBATCH --signal=B:SIGTERM@60
#SBATCH --array=1-15%8
#SBATCH --output=feng/logs/%x-%A_%a.out
#SBATCH --error=feng/logs/%x-%A_%a.err
set -euo pipefail
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
PYTHON_BIN="${FENG_ENV_PREFIX:-/home/dashii/orcd/pool/feng-env}/bin/python"
ROOT="${FENG_ARME_ROOT:-${SCRATCH:?SCRATCH is required}/feng/arm_e}"
FASTA="${FENG_ARME_FASTA:-${SCRATCH}/feng/gencode/v50/GRCh38.primary_assembly.genome.fa}"
REF="${ROOT}/reference/gencode_v50_exonic_union.fa"
PIN="${ROOT}/reference/gencode_v50_exonic_union.pin.json"
TASK="${SLURM_ARRAY_TASK_ID:?array task required}"
LENGTH="$(printf '%02d' "$((2*TASK))")"
MANIFEST="${ROOT}/scan/sites_L${LENGTH}.tsv"
BURDEN="${ROOT}/scan/burden_L${LENGTH}.json"
CKPT="${ROOT}/scan/scan_L${LENGTH}.ckpt.json"
DONE="${ROOT}/scan/scan_L${LENGTH}.done"
mkdir -p "${ROOT}/scan" "${SUBMIT_DIR}/feng/logs"
[[ -s "${REF}" && -s "${PIN}" ]] || { echo "ERROR: pinned transcriptome reference missing" >&2; exit 1; }
[[ -s "${FASTA}" && -s "${FASTA}.fai" ]] || { echo "ERROR: genome FASTA/FAI missing" >&2; exit 1; }
export PYTHONPATH="${SUBMIT_DIR}/feng/src:${PYTHONPATH:-}"
# Candidate length-aware policy is enforced by arm_e_scan: the requested ceiling
# is capped at 20% of rung length (and five), while L2/L4/L6 are terminal-
# uninformative with burden=None. This avoids exhaustive short-mer enumeration.
if [[ -s "${MANIFEST}" && -s "${BURDEN}" && -f "${DONE}" ]]; then
  "${PYTHON_BIN}" -m feng.titration.arm_e_scan validate-burden --manifest "${MANIFEST}" --burden-out "${BURDEN}"
  exit 0
fi
PAYLOAD_PID=""
requeue_handler() { [[ -n "${PAYLOAD_PID}" ]] && kill -TERM "${PAYLOAD_PID}" 2>/dev/null || true; [[ -z "${PAYLOAD_PID}" ]] || wait "${PAYLOAD_PID}" || exit 1; scontrol requeue "${SLURM_JOB_ID:?}" || exit 1; exit 0; }
trap requeue_handler SIGTERM
"${PYTHON_BIN}" -m feng.titration.arm_e_scan scan --task "${TASK}" --reference "${REF}" --pin "${PIN}" --genome-fasta "${FASTA}" --manifest "${MANIFEST}" --burden-out "${BURDEN}" --checkpoint "${CKPT}" --max-edit "${FENG_ARME_MAX_EDIT:-5}" &
PAYLOAD_PID=$!
wait "${PAYLOAD_PID}"
[[ -s "${MANIFEST}" && -s "${BURDEN}" ]] || { echo "ERROR: scan outputs missing" >&2; exit 1; }
"${PYTHON_BIN}" -m feng.titration.arm_e_scan validate-burden --manifest "${MANIFEST}" --burden-out "${BURDEN}"
touch "${DONE}"
trap - SIGTERM
