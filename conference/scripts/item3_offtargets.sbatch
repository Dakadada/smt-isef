#!/usr/bin/env bash
#SBATCH --job-name=conf-aso-ot
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G
#SBATCH --time=1-00:00:00
#SBATCH --requeue
#SBATCH --signal=B:SIGTERM@120
#SBATCH --output=conference/logs/%x-%j.out
#SBATCH --error=conference/logs/%x-%j.err

set -euo pipefail

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
PROJECT_PYTHON="/home/dashii/orcd/pool/feng-env/bin/python"
REFERENCE_DIR="/home/dashii/orcd/scratch/feng/gencode/v50"
GENOME_FASTA="${REFERENCE_DIR}/GRCh38.primary_assembly.genome.fa"
GENCODE_GTF="${REFERENCE_DIR}/gencode.v50.primary_assembly.annotation.gtf.gz"
SCRATCH_ROOT="/orcd/scratch/orcd/011/dashii"
WORK_DIR="${SCRATCH_ROOT}/conference_smn2_aso"
TOOL_ENV="${WORK_DIR}/tools"
INDEX_DIR="${WORK_DIR}/indexes"
RESULT_DIR="${WORK_DIR}/results"
GENE_TABLE="${WORK_DIR}/gencode_v50_gene_spans.tsv"
QUERY_FASTA="${SUBMIT_DIR}/conference/artifacts/item3_unique_queries.fa"
QUERY_MANIFEST="${SUBMIT_DIR}/conference/artifacts/item3_query_manifest.tsv"
ARTIFACT_DIR="${SUBMIT_DIR}/conference/artifacts"
BOWTIE_THREADS=4
BLAST_THREADS=4
SORT_MEMORY=2G
BOWTIE_COUNT_VERSION=v2

export PYTHONPATH="${SUBMIT_DIR}/feng/src:${SUBMIT_DIR}:${PYTHONPATH:-}"
export TMPDIR="/dev/shm"
export SCRATCH="${SCRATCH_ROOT}"

mkdir -p "${WORK_DIR}" "${INDEX_DIR}" "${RESULT_DIR}"

stage_tools() {
    if [[ -x "${TOOL_ENV}/bin/bowtie" && -x "${TOOL_ENV}/bin/blastn" ]]; then
        return
    fi
    module load miniforge/25.11.0-0
    conda create --yes --prefix "${TOOL_ENV}" --channel conda-forge --channel bioconda \
        bowtie=1.3.1 blast=2.16.0
    [[ -x "${TOOL_ENV}/bin/bowtie" && -x "${TOOL_ENV}/bin/blastn" ]]
}

stage_genes() {
    if [[ -s "${GENE_TABLE}" && -f "${GENE_TABLE}.done" ]]; then
        return
    fi
    "${PROJECT_PYTHON}" -m conference.aso_pipeline.item3 genes \
        --gtf "${GENCODE_GTF}" --output "${GENE_TABLE}"
    touch "${GENE_TABLE}.done"
}

stage_bowtie_index() {
    local prefix="${INDEX_DIR}/grch38_v50"
    if [[ -f "${prefix}.1.ebwt" && -f "${prefix}.rev.2.ebwt" && -f "${prefix}.done" ]]; then
        return
    fi
    "${TOOL_ENV}/bin/bowtie-build" --threads "${SLURM_CPUS_PER_TASK}" \
        "${GENOME_FASTA}" "${prefix}"
    touch "${prefix}.done"
}

stage_blast_db() {
    local prefix="${INDEX_DIR}/grch38_v50_blast"
    if [[ -f "${prefix}.ndb" && -f "${prefix}.done" ]]; then
        return
    fi
    "${TOOL_ENV}/bin/makeblastdb" -in "${GENOME_FASTA}" -dbtype nucl \
        -parse_seqids -out "${prefix}"
    touch "${prefix}.done"
}

run_bowtie_partition() {
    local label="$1"
    local query_fasta="$2"
    local threads="$3"
    local output="${RESULT_DIR}/bowtie1_v3_${BOWTIE_COUNT_VERSION}_${label}.tsv"
    local raw_counts="${RESULT_DIR}/bowtie1_v3_${BOWTIE_COUNT_VERSION}_${label}_raw.tsv"
    local sort_dir="${WORK_DIR}/sort/bowtie1_v3_${BOWTIE_COUNT_VERSION}_${label}"
    if [[ -s "${output}" && -f "${output}.done" ]]; then
        return
    fi
    case "${sort_dir}" in
        "${WORK_DIR}"/sort/bowtie1_v3_*) ;;
        *) echo "ERROR: unsafe sort temp path: ${sort_dir}" >&2; exit 2 ;;
    esac
    rm -rf -- "${sort_dir}"
    mkdir -p "${sort_dir}"
    bash "${SUBMIT_DIR}/conference/scripts/bowtie_stream.sh" \
        "${TOOL_ENV}/bin/bowtie" "${INDEX_DIR}/grch38_v50" \
        "${query_fasta}" "${threads}" 3 \
        | "${PROJECT_PYTHON}" -m conference.aso_pipeline.item3 map \
            --aligner bowtie1-v3 --gene-table "${GENE_TABLE}" \
            --query-fasta "${query_fasta}" --raw-counts "${raw_counts}" \
        | LC_ALL=C sort -u -S "${SORT_MEMORY}" --parallel=2 -T "${sort_dir}" \
        | "${PROJECT_PYTHON}" -m conference.aso_pipeline.item3 reduce \
            --aligner bowtie1-v3 --query-fasta "${query_fasta}" \
            --raw-counts "${raw_counts}" --output "${output}"
    touch "${output}.done"
    rm -rf -- "${sort_dir}"
}

stage_bowtie_counts() {
    local parts=(
        "${RESULT_DIR}/bowtie1_v3_${BOWTIE_COUNT_VERSION}_F18.tsv"
        "${RESULT_DIR}/bowtie1_v3_${BOWTIE_COUNT_VERSION}_F14.tsv"
    )
    run_bowtie_partition F18 "${ARTIFACT_DIR}/item3_unique_queries_F18.fa" \
        "${BOWTIE_THREADS}"
    run_bowtie_partition F14 "${ARTIFACT_DIR}/item3_unique_queries_F14.fa" \
        "${BOWTIE_THREADS}"

    local f10_fastas=("${ARTIFACT_DIR}/item3_unique_queries_F10"/*.fa)
    if [[ ! -e "${f10_fastas[0]}" ]]; then
        echo "ERROR: no per-query F10 FASTAs found" >&2
        exit 2
    fi
    local query_fasta label
    for query_fasta in "${f10_fastas[@]}"; do
        label="F10_$(basename "${query_fasta}" .fa)"
        run_bowtie_partition "${label}" "${query_fasta}" 1
        parts+=("${RESULT_DIR}/bowtie1_v3_${BOWTIE_COUNT_VERSION}_${label}.tsv")
    done

    local output="${RESULT_DIR}/bowtie1_v3_counts.tsv"
    "${PROJECT_PYTHON}" -m conference.aso_pipeline.item3 merge \
        --inputs "${parts[@]}" --output "${output}"
    touch "${output}.done"
}

run_blast_partition() {
    local family="$1"
    local query_fasta="${ARTIFACT_DIR}/item3_unique_queries_${family}.fa"
    local output="${RESULT_DIR}/blastn_short_${family}.tsv"
    local raw_counts="${RESULT_DIR}/blastn_short_${family}_raw.tsv"
    local sort_dir="${WORK_DIR}/sort/blastn_short_${family}"
    if [[ -s "${output}" && -f "${output}.done" ]]; then
        return
    fi
    case "${sort_dir}" in
        "${WORK_DIR}"/sort/blastn_short_*) ;;
        *) echo "ERROR: unsafe sort temp path: ${sort_dir}" >&2; exit 2 ;;
    esac
    rm -rf -- "${sort_dir}"
    mkdir -p "${sort_dir}"
    "${TOOL_ENV}/bin/blastn" -task blastn-short -query "${query_fasta}" \
        -db "${INDEX_DIR}/grch38_v50_blast" -evalue 1000 -dust no -strand both \
        -num_threads "${BLAST_THREADS}" \
        -outfmt '6 qseqid sseqid sstart send' \
        | "${PROJECT_PYTHON}" -m conference.aso_pipeline.item3 map \
            --aligner blastn-short --gene-table "${GENE_TABLE}" \
            --query-fasta "${query_fasta}" --raw-counts "${raw_counts}" \
        | LC_ALL=C sort -u -S "${SORT_MEMORY}" --parallel=2 -T "${sort_dir}" \
        | "${PROJECT_PYTHON}" -m conference.aso_pipeline.item3 reduce \
            --aligner blastn-short --query-fasta "${query_fasta}" \
            --raw-counts "${raw_counts}" --output "${output}"
    touch "${output}.done"
    rm -rf -- "${sort_dir}"
}

stage_blast_counts() {
    local parts=()
    local family
    for family in F18 F14 F10; do
        run_blast_partition "${family}"
        parts+=("${RESULT_DIR}/blastn_short_${family}.tsv")
    done
    local output="${RESULT_DIR}/blastn_short_counts.tsv"
    "${PROJECT_PYTHON}" -m conference.aso_pipeline.item3 merge \
        --inputs "${parts[@]}" --output "${output}"
    touch "${output}.done"
}

stage_finalize() {
    "${PROJECT_PYTHON}" -m conference.aso_pipeline.item3 finalize \
        --manifest "${QUERY_MANIFEST}" \
        --bowtie-counts "${RESULT_DIR}/bowtie1_v3_counts.tsv" \
        --blast-counts "${RESULT_DIR}/blastn_short_counts.tsv" \
        --output-dir "${ARTIFACT_DIR}"
    {
        printf '[candidate]\n'
        "${TOOL_ENV}/bin/bowtie" --version | sed -n '1,2p'
        "${TOOL_ENV}/bin/blastn" -version | sed -n '1,2p'
        printf 'GENCODE v50 GRCh38 primary assembly\n'
        printf 'Bowtie: -f -a -v 3; F18/F14 -p %s; each F10 unique query -p 1\n' "${BOWTIE_THREADS}"
        printf 'Bowtie count pipeline version: %s (stdout boundary regression)\n' "${BOWTIE_COUNT_VERSION}"
        printf 'BLAST: -task blastn-short -evalue 1000 -dust no -strand both\n'
        printf 'Gene dedup: disk-backed GNU sort -u, %s RAM ceiling per partition\n' "${SORT_MEMORY}"
    } > "${ARTIFACT_DIR}/item3_software_versions.txt.tmp"
    mv "${ARTIFACT_DIR}/item3_software_versions.txt.tmp" \
        "${ARTIFACT_DIR}/item3_software_versions.txt"
}

run_payload() {
    stage_tools
    stage_genes
    stage_bowtie_index
    stage_blast_db
    stage_bowtie_counts
    stage_blast_counts
    stage_finalize
}

if [[ "${1:-}" == "--payload" ]]; then
    run_payload
    exit 0
fi

PAYLOAD_PID=""
requeue_handler() {
    trap - SIGTERM SIGUSR1
    echo "checkpoint signal received; stopping payload process group" >&2
    if [[ -n "${PAYLOAD_PID}" ]]; then
        kill -TERM -- "-${PAYLOAD_PID}" 2>/dev/null || true
        wait "${PAYLOAD_PID}" 2>/dev/null || true
    fi
    if [[ -z "${SLURM_JOB_ID:-}" ]] || ! scontrol requeue "${SLURM_JOB_ID}"; then
        echo "ERROR: completed stage markers preserved but requeue failed" >&2
        exit 1
    fi
    exit 0
}
trap requeue_handler SIGTERM SIGUSR1

setsid bash "$0" --payload &
PAYLOAD_PID=$!
wait "${PAYLOAD_PID}"
trap - SIGTERM SIGUSR1
