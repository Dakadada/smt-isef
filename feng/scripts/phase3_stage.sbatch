#!/usr/bin/env bash
#SBATCH --job-name=feng-p3-stage
#SBATCH --account=mit_general
#SBATCH --partition=mit_data_transfer
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --requeue
#SBATCH --output=feng/logs/%x-%j.out
#SBATCH --error=feng/logs/%x-%j.err

# Phase 3 staging: download GSE307431 / PRJNA1321660 paired FASTQs from ENA and
# MD5-verify every file. Runs on mit_data_transfer (has outbound internet; keeps
# the load off compute nodes). Downloads run FENG_P3_MAXPAR-at-a-time because a
# single ENA stream is ~1 MB/s and 60 GB would not fit the 12 h wall otherwise;
# 4 concurrent streams stay polite to ENA and the shared node.
#
# Requeue/resume-safe: partial files resume with curl -C -, verified files carry
# a .ok sentinel and are skipped, so preemption or a re-submit never re-downloads
# or corrupts a completed file. The job's pass/fail is decided by a final
# whole-manifest verification pass, not by the background workers' own status.

set -euo pipefail

# --- inputs / paths -------------------------------------------------------
SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
MANIFEST="${FENG_P3_MANIFEST:-${SUBMIT_DIR}/feng/artifacts/phase3/ena_manifest.tsv}"
MAXPAR="${FENG_P3_MAXPAR:-4}"

if [[ -z "${FENG_P3_FASTQ_DIR:-}" ]]; then
    if [[ -z "${SCRATCH:-}" ]]; then
        echo "ERROR: set SCRATCH, or set FENG_P3_FASTQ_DIR" >&2
        exit 1
    fi
    FENG_P3_FASTQ_DIR="${SCRATCH}/feng/phase3/fastq"
fi
OUT="${FENG_P3_FASTQ_DIR}"

[[ -s "${MANIFEST}" ]] || { echo "ERROR: manifest not found or empty: ${MANIFEST}" >&2; exit 1; }
[[ "${MAXPAR}" =~ ^[0-9]+$ && "${MAXPAR}" -ge 1 ]] || { echo "ERROR: FENG_P3_MAXPAR must be a positive integer" >&2; exit 1; }
mkdir -p "${OUT}"
[[ -w "${OUT}" ]] || { echo "ERROR: FASTQ output dir not writable: ${OUT}" >&2; exit 1; }
command -v curl >/dev/null || { echo "ERROR: curl not found" >&2; exit 1; }
command -v md5sum >/dev/null || { echo "ERROR: md5sum not found" >&2; exit 1; }

# --- verify one file against an expected md5 (0 = match) ------------------
verify_md5() {
    local path="$1" want="$2" got
    [[ -s "${path}" ]] || return 1
    got="$(md5sum "${path}" | awk '{print $1}')"
    [[ "${got}" == "${want}" ]]
}

# --- download one file with retries, then md5-verify ----------------------
# Args: url md5. Writes ${OUT}/<basename>.ok on success; returns nonzero on
# permanent failure. Safe to run many at once (distinct dest files).
download_one() {
    local url="$1" want="$2" dest fname attempt
    case "${url}" in
        http://*|https://*) : ;;
        *) url="https://${url}" ;;
    esac
    fname="$(basename "${url}")"
    dest="${OUT}/${fname}"

    if [[ -f "${dest}.ok" ]] && verify_md5 "${dest}" "${want}"; then
        echo "SKIP verified ${fname}"
        return 0
    fi

    local max_attempts="${FENG_P3_ATTEMPTS:-12}"
    for ((attempt=1; attempt<=max_attempts; attempt++)); do
        echo "GET ${fname} (attempt ${attempt}/${max_attempts})"
        if curl -L --fail -sS --connect-timeout 60 \
                --retry 4 --retry-delay 10 --retry-connrefused \
                -C - -o "${dest}" "${url}"; then
            if verify_md5 "${dest}" "${want}"; then
                touch "${dest}.ok"
                echo "OK md5 ${fname}"
                return 0
            fi
            echo "WARN ${fname}: md5 mismatch after attempt ${attempt}; refetching clean" >&2
            rm -f "${dest}"
        else
            # rc=18 (partial) is the common ENA drop: keep the partial so the
            # next attempt resumes it with -C - instead of restarting.
            echo "WARN ${fname}: curl attempt ${attempt} failed (rc=$?); resuming" >&2
        fi
        sleep $(( (attempt < 10 ? attempt : 10) * 15 ))
    done
    echo "ERROR: ${fname} failed md5 after all retries (want ${want}); removing" >&2
    rm -f "${dest}" "${dest}.ok"
    return 1
}

# --- flatten manifest to a url<TAB>md5 worklist ---------------------------
WORKLIST="$(mktemp)"
trap 'rm -f "${WORKLIST}"' EXIT
tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r run sample layout rc bc bytes md5s ftps; do
    [[ -n "${run}" ]] || continue
    IFS=';' read -r -a URL_ARR <<< "${ftps}"
    IFS=';' read -r -a MD5_ARR <<< "${md5s}"
    if [[ "${#URL_ARR[@]}" -ne "${#MD5_ARR[@]}" || "${#URL_ARR[@]}" -eq 0 ]]; then
        echo "ERROR: url/md5 count mismatch for ${run}" >&2
        exit 1
    fi
    for i in "${!URL_ARR[@]}"; do
        printf '%s\t%s\n' "${URL_ARR[$i]}" "${MD5_ARR[$i]}"
    done
done > "${WORKLIST}"

[[ -s "${WORKLIST}" ]] || { echo "ERROR: no files parsed from manifest" >&2; exit 1; }
echo "STAGING: $(wc -l < "${WORKLIST}") files, up to ${MAXPAR} concurrent, into ${OUT}"

# --- parallel download with a concurrency gate ----------------------------
while IFS=$'\t' read -r url want; do
    [[ -n "${url}" ]] || continue
    download_one "${url}" "${want}" &
    # Throttle: block while MAXPAR workers are already running.
    while [[ "$(jobs -rp | wc -l)" -ge "${MAXPAR}" ]]; do
        wait -n 2>/dev/null || true
    done
done < "${WORKLIST}"
wait   # let the last batch finish

# --- final gate: every manifest file must be present and md5-verified -----
tail -n +2 "${MANIFEST}" | while IFS=$'\t' read -r run sample layout rc bc bytes md5s ftps; do
    [[ -n "${run}" ]] || continue
    IFS=';' read -r -a URL_ARR <<< "${ftps}"
    IFS=';' read -r -a MD5_ARR <<< "${md5s}"
    for i in "${!URL_ARR[@]}"; do
        url="${URL_ARR[$i]}"; want="${MD5_ARR[$i]}"
        fname="$(basename "${url}")"; dest="${OUT}/${fname}"
        if ! { [[ -f "${dest}.ok" ]] && verify_md5 "${dest}" "${want}"; }; then
            echo "INCOMPLETE ${fname}" >&2
            exit 3
        fi
    done
done

echo "STAGING COMPLETE: all manifest files present and md5-verified in ${OUT}"
