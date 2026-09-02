#!/usr/bin/env bash
#SBATCH --job-name=feng-p3-star
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --requeue
#SBATCH --array=1-20%10
#SBATCH --output=feng/logs/%x-%A_%a.out
#SBATCH --error=feng/logs/%x-%A_%a.err

# Phase 3 alignment: one STAR job per GSE307431 sample (array task = manifest
# data row). Paired-end, gzipped input read via zcat; sorted BAM + GeneCounts.
# Requeue-safe: a per-sample .done sentinel is written only after an atomic
# publish of the output dir, so preemption/requeue skips finished samples and
# never leaves a half-written result in place. %10 caps concurrency to be a
# considerate guest on the shared preemptable pool.

set -euo pipefail

module load miniforge/25.11.0-0

if [[ -n "${FENG_ENV_PREFIX:-}" ]]; then
    [[ -d "${FENG_ENV_PREFIX}" ]] || { echo "ERROR: FENG_ENV_PREFIX missing: ${FENG_ENV_PREFIX}" >&2; exit 1; }
    CONDA_BASE="$(conda info --base)"
    # shellcheck disable=SC1091
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "${FENG_ENV_PREFIX}"
fi

# --- inputs / paths -------------------------------------------------------
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
MANIFEST="${FENG_P3_MANIFEST:-${SUBMIT_DIR}/feng/artifacts/phase3/ena_manifest.tsv}"

if [[ -z "${SCRATCH:-}" && ( -z "${FENG_P3_FASTQ_DIR:-}" || -z "${FENG_P3_STAR_DIR:-}" ) ]]; then
    echo "ERROR: set SCRATCH, or set both FENG_P3_FASTQ_DIR and FENG_P3_STAR_DIR" >&2
    exit 1
fi
FASTQ_DIR="${FENG_P3_FASTQ_DIR:-${SCRATCH}/feng/phase3/fastq}"
STAR_OUT="${FENG_P3_STAR_DIR:-${SCRATCH}/feng/phase3/star}"
STAR_INDEX="${FENG_STAR_INDEX_DIR:-${SCRATCH}/feng/star/gencode-v50-sjdb149}"
THREADS="${SLURM_CPUS_PER_TASK:-16}"
STAR_BIN="${FENG_STAR_BIN:-$(command -v STAR || true)}"

[[ -s "${MANIFEST}" ]] || { echo "ERROR: manifest missing: ${MANIFEST}" >&2; exit 1; }
[[ -n "${STAR_BIN}" && -x "${STAR_BIN}" ]] || { echo "ERROR: STAR not on PATH; activate feng env or set FENG_STAR_BIN" >&2; exit 1; }
[[ -s "${STAR_INDEX}/SAindex" ]] || { echo "ERROR: STAR index incomplete at ${STAR_INDEX}" >&2; exit 1; }
mkdir -p "${STAR_OUT}"

TASK="${SLURM_ARRAY_TASK_ID:?ERROR: this script must run as a SLURM array job}"

# --- resolve this task's run + paired FASTQs from the manifest ------------
# Data row N (1-based) = manifest line N+1 (header skipped).
row="$(tail -n +2 "${MANIFEST}" | sed -n "${TASK}p")"
[[ -n "${row}" ]] || { echo "ERROR: no manifest data row for array task ${TASK}" >&2; exit 1; }

IFS=$'\t' read -r run sample layout rc bc bytes md5s ftps <<< "${row}"
[[ -n "${run}" ]] || { echo "ERROR: empty run accession at row ${TASK}" >&2; exit 1; }
if [[ "${layout}" != "PAIRED" ]]; then
    echo "ERROR: ${run} layout is '${layout}', expected PAIRED" >&2
    exit 1
fi

R1="${FASTQ_DIR}/${run}_1.fastq.gz"
R2="${FASTQ_DIR}/${run}_2.fastq.gz"
[[ -s "${R1}" && -s "${R2}" ]] || { echo "ERROR: staged FASTQs missing for ${run} (run phase3_stage first)" >&2; exit 1; }

FINAL_DIR="${STAR_OUT}/${run}"
DONE="${FINAL_DIR}/.done"
if [[ -f "${DONE}" ]]; then
    echo "SKIP ${run}: already complete (${DONE})"
    exit 0
fi

# --- align into a task-scoped work dir, then publish atomically ----------
WORK="${STAR_OUT}/.work-${run}-${SLURM_JOB_ID:-manual}"
rm -rf "${WORK}"
mkdir -p "${WORK}"
trap 'rm -rf "${WORK}"' EXIT

echo "STAR ${run}: threads=${THREADS} index=${STAR_INDEX}"
"${STAR_BIN}" \
    --runThreadN "${THREADS}" \
    --genomeDir "${STAR_INDEX}" \
    --readFilesIn "${R1}" "${R2}" \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix "${WORK}/${run}." \
    --outTmpDir "${WORK}/_tmp"

# Require the two outputs the downstream ΔPSI / count cross-check depends on.
[[ -s "${WORK}/${run}.Aligned.sortedByCoord.out.bam" ]] || { echo "ERROR: STAR produced no sorted BAM for ${run}" >&2; exit 1; }
[[ -s "${WORK}/${run}.ReadsPerGene.out.tab" ]] || { echo "ERROR: STAR produced no gene counts for ${run}" >&2; exit 1; }

# Atomic publish: swap the finished work dir into place, then mark done.
rm -rf "${FINAL_DIR}.replacing"
mv "${WORK}" "${FINAL_DIR}.replacing"
trap - EXIT
rm -rf "${FINAL_DIR}"
mv "${FINAL_DIR}.replacing" "${FINAL_DIR}"
touch "${DONE}"
echo "DONE ${run}: ${FINAL_DIR}"
