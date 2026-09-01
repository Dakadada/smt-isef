#!/usr/bin/env bash
#SBATCH --job-name=feng-star-v50
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=02:00:00
#SBATCH --requeue
#SBATCH --output=feng/logs/%x-%j.out
#SBATCH --error=feng/logs/%x-%j.err

# account=mit_general, partition=mit_preemptable filled from the cluster's SLURM
# convention. This script is requeue-safe (idempotent complete-index skip + atomic
# publish), so preemption on mit_preemptable resumes cleanly on requeue.

set -euo pipefail

module load miniforge/25.11.0-0

if [[ -n "${FENG_ENV_PREFIX:-}" ]]; then
    if [[ ! -d "${FENG_ENV_PREFIX}" ]]; then
        echo "ERROR: FENG_ENV_PREFIX does not exist: ${FENG_ENV_PREFIX}" >&2
        exit 1
    fi
    CONDA_BASE="$(conda info --base)"
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "${FENG_ENV_PREFIX}"
fi

if [[ -z "${FENG_DATA_ROOT:-}" || -z "${FENG_STAR_INDEX_DIR:-}" ]]; then
    if [[ -z "${SCRATCH:-}" ]]; then
        echo "ERROR: set SCRATCH, or set both FENG_DATA_ROOT and FENG_STAR_INDEX_DIR" >&2
        exit 1
    fi
fi

FENG_DATA_ROOT="${FENG_DATA_ROOT:-${SCRATCH}/feng/gencode/v50}"
FENG_SJDB_OVERHANG="${FENG_SJDB_OVERHANG:-149}"
FENG_STAR_INDEX_DIR="${FENG_STAR_INDEX_DIR:-${SCRATCH}/feng/star/gencode-v50-sjdb${FENG_SJDB_OVERHANG}}"
FENG_STAR_RAM_BYTES="${FENG_STAR_RAM_BYTES:-55000000000}"
THREADS="${SLURM_CPUS_ON_NODE:-16}"
STAR_BIN="${FENG_STAR_BIN:-$(command -v STAR || true)}"
PYTHON_BIN="${FENG_PYTHON_BIN:-$(command -v python || true)}"

FASTA_PATH="${FENG_DATA_ROOT}/GRCh38.primary_assembly.genome.fa"
GTF_PATH="${FENG_DATA_ROOT}/gencode.v50.primary_assembly.annotation.gtf"
INDEX_PARENT="$(dirname "${FENG_STAR_INDEX_DIR}")"
BUILD_DIR="${FENG_STAR_INDEX_DIR}.building-${SLURM_JOB_ID:-manual}"

if ! [[ "${FENG_SJDB_OVERHANG}" =~ ^[0-9]+$ \
    && "${FENG_STAR_RAM_BYTES}" =~ ^[0-9]+$ \
    && "${THREADS}" =~ ^[0-9]+$ ]]; then
    echo "ERROR: sjdbOverhang, RAM bytes, and thread count must be nonnegative integers" >&2
    exit 1
fi
if [[ -z "${STAR_BIN}" || ! -x "${STAR_BIN}" ]]; then
    echo "ERROR: STAR is not on PATH; activate the feng env or set FENG_STAR_BIN" >&2
    exit 1
fi
if [[ -z "${PYTHON_BIN}" || ! -x "${PYTHON_BIN}" ]]; then
    echo "ERROR: python is not on PATH; activate the feng env or set FENG_PYTHON_BIN" >&2
    exit 1
fi
if [[ ! -s "${FASTA_PATH}" || ! -s "${GTF_PATH}" ]]; then
    echo "ERROR: fetched FASTA/GTF inputs are missing under ${FENG_DATA_ROOT}" >&2
    exit 1
fi
if ! mkdir -p "${INDEX_PARENT}" || [[ ! -d "${INDEX_PARENT}" || ! -w "${INDEX_PARENT}" ]]; then
    echo "ERROR: STAR index parent is not writable: ${INDEX_PARENT}" >&2
    exit 1
fi

case "${BUILD_DIR}" in
    ""|"/")
        echo "ERROR: unsafe STAR build directory: ${BUILD_DIR}" >&2
        exit 1
        ;;
esac

STAR_VERSION="$("${STAR_BIN}" --version 2>&1 | head -n 1)"
FASTA_MD5="$(md5sum "${FASTA_PATH}" | awk '{print $1}')"
GTF_MD5="$(md5sum "${GTF_PATH}" | awk '{print $1}')"

REQUIRED_INDEX_FILES=(
    Genome
    SA
    SAindex
    chrLength.txt
    chrName.txt
    chrNameLength.txt
    chrStart.txt
    exonGeTrInfo.tab
    exonInfo.tab
    geneInfo.tab
    genomeParameters.txt
    sjdbInfo.txt
    sjdbList.fromGTF.out.tab
    sjdbList.out.tab
    transcriptInfo.tab
)

index_files_complete() {
    local index_dir="$1"
    local required_file
    for required_file in "${REQUIRED_INDEX_FILES[@]}"; do
        if [[ ! -s "${index_dir}/${required_file}" ]]; then
            return 1
        fi
    done
}

marker_matches() {
    local index_dir="$1"
    "${PYTHON_BIN}" - \
        "${index_dir}/feng_index_manifest.json" \
        "${FASTA_MD5}" "${GTF_MD5}" "${STAR_VERSION}" \
        "${FENG_SJDB_OVERHANG}" "${FENG_STAR_RAM_BYTES}" "${THREADS}" <<'PY'
import json
from pathlib import Path
import sys

marker_path = Path(sys.argv[1])
expected = {
    "schema_version": 1,
    "inputs": {
        "fasta_md5": sys.argv[2],
        "gtf_md5": sys.argv[3],
    },
    "star_version": sys.argv[4],
    "generation_params": {
        "run_mode": "genomeGenerate",
        "run_threads": int(sys.argv[7]),
        "sjdb_gtf": True,
        "sjdb_overhang": int(sys.argv[5]),
        "limit_genome_generate_ram": int(sys.argv[6]),
    },
}
try:
    observed = json.loads(marker_path.read_text(encoding="utf-8"))
except (OSError, json.JSONDecodeError):
    raise SystemExit(1)
raise SystemExit(0 if observed == expected else 1)
PY
}

write_marker() {
    local index_dir="$1"
    "${PYTHON_BIN}" - \
        "${index_dir}/feng_index_manifest.json" \
        "${FASTA_MD5}" "${GTF_MD5}" "${STAR_VERSION}" \
        "${FENG_SJDB_OVERHANG}" "${FENG_STAR_RAM_BYTES}" "${THREADS}" <<'PY'
import json
from pathlib import Path
import sys

marker = {
    "schema_version": 1,
    "inputs": {
        "fasta_md5": sys.argv[2],
        "gtf_md5": sys.argv[3],
    },
    "star_version": sys.argv[4],
    "generation_params": {
        "run_mode": "genomeGenerate",
        "run_threads": int(sys.argv[7]),
        "sjdb_gtf": True,
        "sjdb_overhang": int(sys.argv[5]),
        "limit_genome_generate_ram": int(sys.argv[6]),
    },
}
path = Path(sys.argv[1])
temporary = path.with_name(path.name + ".tmp")
temporary.write_text(json.dumps(marker, indent=2, sort_keys=True) + "\n", encoding="utf-8")
temporary.replace(path)
PY
}

if index_files_complete "${FENG_STAR_INDEX_DIR}" && marker_matches "${FENG_STAR_INDEX_DIR}"; then
    echo "verified matching STAR index at ${FENG_STAR_INDEX_DIR}; nothing to do"
    exit 0
fi

if [[ -e "${FENG_STAR_INDEX_DIR}" ]]; then
    echo "existing STAR index is incomplete or does not match requested inputs/parameters; rebuilding"
    rm -rf -- "${FENG_STAR_INDEX_DIR}"
fi

cleanup() {
    if [[ -n "${BUILD_DIR:-}" && -d "${BUILD_DIR}" ]]; then
        rm -rf -- "${BUILD_DIR}"
    fi
}
trap cleanup EXIT INT TERM

# A requeued job gets the same job id. Remove only its incomplete staging directory,
# then restart STAR; the completed index is published atomically below.
if [[ -d "${BUILD_DIR}" ]]; then
    rm -rf -- "${BUILD_DIR}"
fi
mkdir -p "${BUILD_DIR}"

"${STAR_BIN}" \
    --runMode genomeGenerate \
    --runThreadN "${THREADS}" \
    --genomeDir "${BUILD_DIR}" \
    --genomeFastaFiles "${FASTA_PATH}" \
    --sjdbGTFfile "${GTF_PATH}" \
    --sjdbOverhang "${FENG_SJDB_OVERHANG}" \
    --limitGenomeGenerateRAM "${FENG_STAR_RAM_BYTES}"

if ! index_files_complete "${BUILD_DIR}"; then
    echo "ERROR: STAR exited without a complete index" >&2
    exit 1
fi
write_marker "${BUILD_DIR}"
if ! marker_matches "${BUILD_DIR}"; then
    echo "ERROR: STAR completion marker does not match requested build" >&2
    exit 1
fi
if [[ -e "${FENG_STAR_INDEX_DIR}" ]]; then
    echo "ERROR: final index path appeared during build: ${FENG_STAR_INDEX_DIR}" >&2
    exit 1
fi

mv -- "${BUILD_DIR}" "${FENG_STAR_INDEX_DIR}"
BUILD_DIR=""
echo "STAR index complete: ${FENG_STAR_INDEX_DIR}"
