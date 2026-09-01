#!/usr/bin/env bash
set -euo pipefail

FASTA_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/GRCh38.primary_assembly.genome.fa.gz"
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.primary_assembly.annotation.gtf.gz"
FASTA_MD5="da1a11258be075cfa7af718162c894e7"
GTF_MD5="289b91e5e95e8b0450d223246f10a12e"

download_resume() {
    local url="$1"
    local partial="$2"

    if command -v curl >/dev/null 2>&1; then
        curl --fail --location --retry 3 --retry-delay 5 \
            --continue-at - --output "${partial}" "${url}"
    elif command -v wget >/dev/null 2>&1; then
        wget --continue --output-document="${partial}" "${url}"
    else
        echo "ERROR: curl or wget is required" >&2
        return 1
    fi
}

verify_pinned_md5() {
    local artifact="$1"
    local expected_md5="$2"
    local actual_md5

    if [[ ! "${expected_md5}" =~ ^[0-9a-fA-F]{32}$ || "${expected_md5}" == "00000000000000000000000000000000" ]]; then
        echo "ERROR: authoritative MD5 is empty, invalid, or all-zero for ${artifact}" >&2
        return 1
    fi
    actual_md5="$(md5sum "${artifact}" | awk '{print $1}')"
    if [[ "${actual_md5}" != "${expected_md5,,}" ]]; then
        echo "ERROR: authoritative MD5 mismatch for ${artifact}" >&2
        echo "expected ${expected_md5,,}; observed ${actual_md5}" >&2
        return 1
    fi
    echo "verified authoritative MD5  ${actual_md5}  ${artifact}"
}

record_secondary_hashes() {
    local archive="$1"
    local authoritative_md5="$2"
    local sha_sidecar="${archive}.sha256"
    local md5_sidecar="${archive}.md5"
    local actual_sha256
    local recorded_sha256=""
    local temporary

    actual_sha256="$(sha256sum "${archive}" | awk '{print $1}')"
    if [[ -s "${sha_sidecar}" ]]; then
        recorded_sha256="$(awk 'NR == 1 {print $1}' "${sha_sidecar}")"
        if [[ "${actual_sha256}" != "${recorded_sha256}" ]]; then
            echo "ERROR: recorded secondary SHA-256 mismatch for ${archive}" >&2
            echo "expected ${recorded_sha256}; observed ${actual_sha256}" >&2
            return 1
        fi
    else
        temporary="${sha_sidecar}.tmp.$$"
        printf '%s  %s\n' "${actual_sha256}" "$(basename "${archive}")" > "${temporary}"
        mv -- "${temporary}" "${sha_sidecar}"
    fi

    temporary="${md5_sidecar}.tmp.$$"
    printf '%s  %s\n' "${authoritative_md5,,}" "$(basename "${archive}")" > "${temporary}"
    mv -- "${temporary}" "${md5_sidecar}"
    echo "recorded secondary SHA-256  ${actual_sha256}  ${archive}"
}

fetch_archive() {
    local url="$1"
    local archive="$2"
    local expected_md5="$3"
    local partial="${archive}.part"

    if [[ ! -f "${archive}" ]]; then
        echo "fetching ${url}"
        download_resume "${url}" "${partial}"
        verify_pinned_md5 "${partial}" "${expected_md5}"
        mv -- "${partial}" "${archive}"
    else
        echo "found ${archive}; verifying before reuse"
        verify_pinned_md5 "${archive}" "${expected_md5}"
    fi
    record_secondary_hashes "${archive}" "${expected_md5}"
}

validate_fasta() {
    local path="$1"
    local first_line=""
    IFS= read -r first_line < "${path}" || true
    [[ "${first_line}" == '>'* ]]
}

validate_gtf() {
    local path="$1"
    awk -F '\t' '
        BEGIN { found = 0 }
        /^[[:space:]]*#/ || /^[[:space:]]*$/ { next }
        { found = 1; if (NF != 9) exit 1; exit 0 }
        END { if (!found) exit 1 }
    ' "${path}"
}

validate_format() {
    local output="$1"
    local format="$2"
    if [[ "${format}" == "fasta" ]]; then
        validate_fasta "${output}"
    elif [[ "${format}" == "gtf" ]]; then
        validate_gtf "${output}"
    else
        echo "ERROR: unknown decompressed format: ${format}" >&2
        return 1
    fi
}

decompress_atomic() {
    local archive="$1"
    local output="$2"
    local archive_md5="$3"
    local format="$4"
    local provenance="${output}.provenance"
    local temporary="${output}.tmp.$$"
    local provenance_temporary="${provenance}.tmp.$$"
    local recorded_archive_md5=""
    local recorded_output_md5=""
    local actual_output_md5=""

    if [[ -s "${output}" && -s "${provenance}" ]]; then
        recorded_archive_md5="$(awk -F= '$1 == "archive_md5" {print $2}' "${provenance}")"
        recorded_output_md5="$(awk -F= '$1 == "output_md5" {print $2}' "${provenance}")"
        actual_output_md5="$(md5sum "${output}" | awk '{print $1}')"
        if [[ "${recorded_archive_md5}" == "${archive_md5,,}" \
            && -n "${recorded_output_md5}" \
            && "${actual_output_md5}" == "${recorded_output_md5}" ]] \
            && validate_format "${output}" "${format}"; then
            echo "verified decompressed artifact and provenance: ${output}"
            return 0
        fi
        echo "decompressed provenance or content mismatch; regenerating ${output}"
    fi

    if ! gzip -dc "${archive}" > "${temporary}"; then
        rm -f -- "${temporary}"
        return 1
    fi
    if ! validate_format "${temporary}" "${format}"; then
        echo "ERROR: decompressed ${format} validation failed for ${archive}" >&2
        rm -f -- "${temporary}"
        return 1
    fi
    actual_output_md5="$(md5sum "${temporary}" | awk '{print $1}')"
    printf 'archive_md5=%s\noutput_md5=%s\n' \
        "${archive_md5,,}" "${actual_output_md5}" > "${provenance_temporary}"
    mv -- "${temporary}" "${output}"
    mv -- "${provenance_temporary}" "${provenance}"
}

main() {
    if [[ -z "${FENG_DATA_ROOT:-}" ]]; then
        if [[ -z "${SCRATCH:-}" ]]; then
            echo "ERROR: set FENG_DATA_ROOT or SCRATCH to writable cluster storage" >&2
            return 1
        fi
        FENG_DATA_ROOT="${SCRATCH}/feng/gencode/v50"
    fi

    if ! mkdir -p "${FENG_DATA_ROOT}" \
        || [[ ! -d "${FENG_DATA_ROOT}" || ! -w "${FENG_DATA_ROOT}" ]]; then
        echo "ERROR: GENCODE data root is not writable: ${FENG_DATA_ROOT}" >&2
        return 1
    fi

    local fasta_archive="${FENG_DATA_ROOT}/GRCh38.primary_assembly.genome.fa.gz"
    local gtf_archive="${FENG_DATA_ROOT}/gencode.v50.primary_assembly.annotation.gtf.gz"
    local fasta_path="${fasta_archive%.gz}"
    local gtf_path="${gtf_archive%.gz}"

    fetch_archive "${FASTA_URL}" "${fasta_archive}" "${FASTA_MD5}"
    fetch_archive "${GTF_URL}" "${gtf_archive}" "${GTF_MD5}"
    decompress_atomic "${fasta_archive}" "${fasta_path}" "${FASTA_MD5}" "fasta"
    decompress_atomic "${gtf_archive}" "${gtf_path}" "${GTF_MD5}" "gtf"

    echo "GENCODE v50 primary-assembly inputs ready under ${FENG_DATA_ROOT}"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi
