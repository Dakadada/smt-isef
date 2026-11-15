#!/usr/bin/env bash
#SBATCH --job-name=feng-armc-gate
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --requeue
#SBATCH --output=feng/logs/%x-%j.out
#SBATCH --error=feng/logs/%x-%j.err

# Arm C bind->act gate, stage 2: merge the 81 per-SBO hits parquets written by
# feng/scripts/armc_scan.sbatch into one armc_features.parquet, then run the
# leave-one-SBO-out gate and write GATE_RESULT.md. Submit with
# --dependency=afterok:<armc_scan array job id> so this only runs once every
# scan task has succeeded. Idempotent (skip if .done already present).

set -euo pipefail

module load miniforge/25.11.0-0

if [[ -n "${FENG_ENV_PREFIX:-}" ]]; then
    [[ -d "${FENG_ENV_PREFIX}" ]] || { echo "ERROR: FENG_ENV_PREFIX missing: ${FENG_ENV_PREFIX}" >&2; exit 1; }
    CONDA_BASE="$(conda info --base)"
    # shellcheck disable=SC1091
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "${FENG_ENV_PREFIX}"
fi

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
ENV_PREFIX="${FENG_ENV_PREFIX:-/home/dashii/orcd/pool/feng-env}"
PYTHON_BIN="${FENG_BENCH_PYTHON:-${ENV_PREFIX}/bin/python}"

ARMC_DIR="${FENG_ARMC_DIR:-${SCRATCH:-}/feng/armc}"
HITS_DIR="${ARMC_DIR}/hits"
FEATURES_OUT="${FENG_ARMC_FEATURES:-${SUBMIT_DIR}/feng/artifacts/phase0_arm_c_probe/armc_features.parquet}"
GATE_OUT="${FENG_ARMC_GATE_RESULT:-${SUBMIT_DIR}/feng/artifacts/phase0_arm_c_probe/GATE_RESULT.md}"
DONE="${GATE_OUT}.done"

SBO_TSV="${FENG_ARMC_SBO_TSV:-${SUBMIT_DIR}/feng/artifacts/phase0_arm_c_probe/sbo_sequences.tsv}"

if [[ -s "${GATE_OUT}" && -f "${DONE}" ]]; then
    echo "gate result already done; skipping ($GATE_OUT)"
    exit 0
fi

[[ -d "${HITS_DIR}" ]] || { echo "ERROR: hits dir missing: ${HITS_DIR}" >&2; exit 1; }

n_hits="$(find "${HITS_DIR}" -maxdepth 1 -name '*.parquet' | wc -l)"
n_expected="$(($(tail -n +2 "${SBO_TSV}" | wc -l)))"
if [[ "${n_hits}" -lt "${n_expected}" ]]; then
    echo "ERROR: only ${n_hits}/${n_expected} per-SBO hits parquets present in ${HITS_DIR}" >&2
    echo "       (did the armc_scan array job finish? check its .done sentinels)" >&2
    exit 1
fi

export PYTHONPATH="${SUBMIT_DIR}/feng/src:${PYTHONPATH:-}"
mkdir -p feng/logs "$(dirname "${FEATURES_OUT}")"

echo "merging ${n_hits} per-SBO hits parquets -> ${FEATURES_OUT}"
"${PYTHON_BIN}" -m feng.bench.armc_build merge \
    --hits-dir "${HITS_DIR}" --out "${FEATURES_OUT}"

echo "running LOSO gate -> ${GATE_OUT}"
"${PYTHON_BIN}" -m feng.bench.armc_gate \
    --features "${FEATURES_OUT}" --out "${GATE_OUT}"

touch "${DONE}"
echo "done: ${GATE_OUT}"
