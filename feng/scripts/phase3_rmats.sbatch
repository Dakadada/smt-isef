#!/usr/bin/env bash
#SBATCH --job-name=feng-p3-rmats
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --requeue
#SBATCH --output=feng/logs/%x-%j.out
#SBATCH --error=feng/logs/%x-%j.err

# Phase 3 rMATS 4.3.0 contrasts for GSE307431. Each targeted ASO is compared
# only with its own chemistry-matched scrambled controls. Untreated samples are
# retained in the metadata but never added to a chemistry contrast.
#
# Requeue safety: each contrast runs in a job-scoped work directory and is
# atomically published with its .done sentinel only after required outputs pass
# validation. A requeue skips already published contrasts.

set -euo pipefail

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
CONDITIONS="${FENG_P3_CONDITIONS:-${SUBMIT_DIR}/feng/artifacts/phase3/sample_conditions.tsv}"
ENV_PREFIX="${FENG_ENV_PREFIX:-/home/dashii/orcd/pool/feng-env}"
PYTHON_BIN="${FENG_P3_PYTHON:-${ENV_PREFIX}/bin/python}"
RMATS="${FENG_P3_RMATS:-${ENV_PREFIX}/bin/rmats.py}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
READ_LENGTH="${FENG_P3_READ_LENGTH:-100}"
LIB_TYPE="${FENG_P3_LIBTYPE:-fr-firststrand}"

if [[ -z "${SCRATCH:-}" && ( -z "${FENG_P3_STAR_DIR:-}" || -z "${FENG_P3_RMATS_DIR:-}" || -z "${FENG_P3_GTF:-}" ) ]]; then
    echo "ERROR: set SCRATCH, or explicitly set FENG_P3_STAR_DIR, FENG_P3_RMATS_DIR, and FENG_P3_GTF" >&2
    exit 1
fi
STAR_DIR="${FENG_P3_STAR_DIR:-${SCRATCH}/feng/phase3/star}"
OUT_ROOT="${FENG_P3_RMATS_DIR:-${SCRATCH}/feng/phase3/rmats}"

# The requested canonical path is checked first. The existing Phase 3 GENCODE
# v50 staging uses the upstream primary-assembly filename, accepted as a
# documented fallback rather than silently referencing a missing file.
GTF="${FENG_P3_GTF:-${SCRATCH}/feng/gencode/v50/annotation.gtf}"
if [[ ! -s "${GTF}" && -z "${FENG_P3_GTF:-}" ]]; then
    GTF_FALLBACK="${SCRATCH}/feng/gencode/v50/gencode.v50.primary_assembly.annotation.gtf"
    if [[ -s "${GTF_FALLBACK}" ]]; then
        echo "WARNING: ${GTF} is absent; using ${GTF_FALLBACK}" >&2
        GTF="${GTF_FALLBACK}"
    fi
fi

[[ -s "${CONDITIONS}" ]] || {
    echo "ERROR: validated conditions table missing: ${CONDITIONS}" >&2
    echo "Run: bash feng/scripts/fetch_geo_metadata.sh" >&2
    exit 1
}
[[ -s "${GTF}" ]] || { echo "ERROR: GENCODE v50 GTF missing: ${GTF}" >&2; exit 1; }
[[ -x "${PYTHON_BIN}" ]] || { echo "ERROR: Python missing: ${PYTHON_BIN}" >&2; exit 1; }
[[ -s "${RMATS}" ]] || { echo "ERROR: rMATS missing: ${RMATS}" >&2; exit 1; }
[[ "${READ_LENGTH}" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid read length: ${READ_LENGTH}" >&2; exit 1; }

for required in run_accession condition chemistry control_type; do
    if ! awk -F '\t' -v wanted="${required}" 'NR==1 {for (i=1;i<=NF;i++) if ($i==wanted) found=1} END {exit !found}' "${CONDITIONS}"; then
        echo "ERROR: conditions table lacks column: ${required}" >&2
        exit 1
    fi
done

mkdir -p "${OUT_ROOT}"

column_number() {
    awk -F '\t' -v wanted="$1" 'NR==1 {for (i=1;i<=NF;i++) if ($i==wanted) {print i; exit}}' "${CONDITIONS}"
}
RUN_COL="$(column_number run_accession)"
CONDITION_COL="$(column_number condition)"
CHEMISTRY_COL="$(column_number chemistry)"
CONTROL_COL="$(column_number control_type)"

collect_runs() {
    local wanted_condition="$1"
    local wanted_chemistry="$2"
    local wanted_control="$3"
    awk -F '\t' \
        -v rc="${RUN_COL}" -v cc="${CONDITION_COL}" -v hc="${CHEMISTRY_COL}" -v tc="${CONTROL_COL}" \
        -v condition="${wanted_condition}" -v chemistry="${wanted_chemistry}" -v control="${wanted_control}" \
        'NR > 1 && $cc == condition && $hc == chemistry && $tc == control {print $rc}' \
        "${CONDITIONS}"
}

cleanup_work=""
cleanup() {
    if [[ -n "${cleanup_work}" && "${cleanup_work}" == "${OUT_ROOT}"/.work-* ]]; then
        rm -rf -- "${cleanup_work}"
    fi
}
trap cleanup EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

for specification in "F18OMe:OMe" "F18MOE:MOE" "F20PMO:PMO"; do
    condition="${specification%%:*}"
    chemistry="${specification##*:}"
    contrast="${condition}_vs_scrambled_${chemistry}"
    final_dir="${OUT_ROOT}/${contrast}"
    done_file="${final_dir}/.done"

    if [[ -f "${done_file}" ]]; then
        echo "SKIP ${contrast}: already complete (${done_file})"
        continue
    fi
    if [[ -e "${final_dir}" ]]; then
        echo "ERROR: incomplete final directory exists without sentinel: ${final_dir}" >&2
        exit 1
    fi

    mapfile -t target_runs < <(collect_runs "${condition}" "${chemistry}" "targeted")
    mapfile -t control_runs < <(collect_runs "control" "${chemistry}" "scrambled")
    if [[ "${#target_runs[@]}" -ne 3 || "${#control_runs[@]}" -ne 3 ]]; then
        echo "ERROR: ${contrast} requires exactly 3 targeted and 3 matched scrambled samples; got ${#target_runs[@]} and ${#control_runs[@]}" >&2
        exit 1
    fi

    target_bams=()
    control_bams=()
    for run in "${target_runs[@]}"; do
        [[ -f "${STAR_DIR}/${run}/.done" ]] || { echo "ERROR: alignment sentinel missing for ${run}" >&2; exit 1; }
        bam="${STAR_DIR}/${run}/${run}.Aligned.sortedByCoord.out.bam"
        [[ -s "${bam}" ]] || { echo "ERROR: BAM missing for ${run}: ${bam}" >&2; exit 1; }
        target_bams+=("${bam}")
    done
    for run in "${control_runs[@]}"; do
        [[ -f "${STAR_DIR}/${run}/.done" ]] || { echo "ERROR: alignment sentinel missing for ${run}" >&2; exit 1; }
        bam="${STAR_DIR}/${run}/${run}.Aligned.sortedByCoord.out.bam"
        [[ -s "${bam}" ]] || { echo "ERROR: BAM missing for ${run}: ${bam}" >&2; exit 1; }
        control_bams+=("${bam}")
    done

    cleanup_work="${OUT_ROOT}/.work-${contrast}-${SLURM_JOB_ID:-manual}"
    cleanup
    mkdir -p "${cleanup_work}/output" "${cleanup_work}/tmp"
    (IFS=,; printf '%s\n' "${target_bams[*]}") > "${cleanup_work}/target_bams.txt"
    (IFS=,; printf '%s\n' "${control_bams[*]}") > "${cleanup_work}/control_bams.txt"

    echo "RUN ${contrast}: 3 targeted vs 3 chemistry-matched scrambled; readLength=${READ_LENGTH}; libType=${LIB_TYPE}"
    "${PYTHON_BIN}" "${RMATS}" \
        --b1 "${cleanup_work}/target_bams.txt" \
        --b2 "${cleanup_work}/control_bams.txt" \
        --gtf "${GTF}" \
        --od "${cleanup_work}/output" \
        --tmp "${cleanup_work}/tmp" \
        --task both \
        -t paired \
        --readLength "${READ_LENGTH}" \
        --libType "${LIB_TYPE}" \
        --nthread "${THREADS}"

    for event in SE MXE A5SS A3SS RI; do
        [[ -s "${cleanup_work}/output/${event}.MATS.JC.txt" ]] || {
            echo "ERROR: rMATS output missing for ${contrast}: ${event}.MATS.JC.txt" >&2
            exit 1
        }
    done
    cp "${cleanup_work}/target_bams.txt" "${cleanup_work}/output/target_bams.txt"
    cp "${cleanup_work}/control_bams.txt" "${cleanup_work}/output/control_bams.txt"
    printf 'condition\t%s\nchemistry\t%s\ncontrol_type\tscrambled\nread_length\t%s\nlib_type\t%s\n' \
        "${condition}" "${chemistry}" "${READ_LENGTH}" "${LIB_TYPE}" \
        > "${cleanup_work}/output/contrast_metadata.tsv"
    touch "${cleanup_work}/output/.done"
    mv "${cleanup_work}/output" "${final_dir}"
    cleanup
    cleanup_work=""
    echo "DONE ${contrast}: ${final_dir}"
done

echo "All three chemistry-matched rMATS contrasts are complete."
