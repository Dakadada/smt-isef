#!/usr/bin/env bash
# Candidate-only ENERGY-FILTERED extension of conference item 3.
set -euo pipefail

SUBMIT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PROJECT_PYTHON="/home/dashii/orcd/pool/feng-env/bin/python"
INTARNA="/home/dashii/orcd/pool/feng-env/bin/IntaRNA"
SCRATCH_ROOT="/orcd/scratch/orcd/011/dashii"
WORK_DIR="${SCRATCH_ROOT}/conference_smn2_aso"
BOWTIE="${WORK_DIR}/tools/bin/bowtie"
BOWTIE_INDEX="${WORK_DIR}/indexes/grch38_v50"
GENE_TABLE="${WORK_DIR}/gencode_v50_gene_spans.tsv"
RUN_DIR="${ITEM3_ENERGY_RUN_DIR:-${WORK_DIR}/results/item3_energy_v1}"
RESULT_DIR="${RUN_DIR}/results"
SORT_DIR="${RUN_DIR}/sort"
ARTIFACT_DIR="${SUBMIT_DIR}/conference/artifacts"
SCORE_TABLE="${ARTIFACT_DIR}/item3_energy_variant_scores.tsv"
QUERY_FASTA="${ARTIFACT_DIR}/item3_energy_viable_100x.fa"
RAW_COUNTS="${RESULT_DIR}/exact_v0_raw_counts.tsv"
EXACT_COUNTS="${ARTIFACT_DIR}/item3_energy_exact_v0_counts.tsv"
RAW_GENE_COUNTS="${ARTIFACT_DIR}/item3_real_aso_gene_counts.tsv"
BOWTIE_THREADS="${BOWTIE_THREADS:-4}"
INTARNA_THREADS="${INTARNA_THREADS:-4}"
SORT_MEMORY="${SORT_MEMORY:-1G}"

export PYTHONPATH="${SUBMIT_DIR}/feng/src:${SUBMIT_DIR}:${PYTHONPATH:-}"
export SCRATCH="${SCRATCH_ROOT}"

for required in \
    "${PROJECT_PYTHON}" "${INTARNA}" "${BOWTIE}" \
    "${BOWTIE_INDEX}.1.ebwt" "${BOWTIE_INDEX}.rev.2.ebwt" \
    "${GENE_TABLE}" "${RAW_GENE_COUNTS}"; do
    if [[ ! -f "${required}" ]]; then
        echo "ERROR: required input missing: ${required}" >&2
        exit 2
    fi
done

mkdir -p "${RESULT_DIR}" "${ARTIFACT_DIR}"
case "${SORT_DIR}" in
    "${WORK_DIR}"/results/item3_energy_v1/sort) ;;
    /tmp/conference_smn2_aso_item3_energy_v1/sort) ;;
    *) echo "ERROR: unsafe sort path: ${SORT_DIR}" >&2; exit 2 ;;
esac
rm -rf -- "${SORT_DIR}"
mkdir -p "${SORT_DIR}"
trap 'rm -rf -- "${SORT_DIR}"' EXIT

"${PROJECT_PYTHON}" -m conference.aso_pipeline.item3_energy score \
    --output-dir "${ARTIFACT_DIR}" --intarna "${INTARNA}" \
    --threads "${INTARNA_THREADS}"

bash "${SUBMIT_DIR}/conference/scripts/bowtie_stream.sh" \
    "${BOWTIE}" "${BOWTIE_INDEX}" "${QUERY_FASTA}" \
    "${BOWTIE_THREADS}" 0 \
    | "${PROJECT_PYTHON}" -m conference.aso_pipeline.item3_energy map \
        --gene-table "${GENE_TABLE}" --score-table "${SCORE_TABLE}" \
        --raw-counts "${RAW_COUNTS}" \
    | LC_ALL=C sort -u -S "${SORT_MEMORY}" --parallel=2 -T "${SORT_DIR}" \
    | "${PROJECT_PYTHON}" -m conference.aso_pipeline.item3_energy reduce \
        --raw-counts "${RAW_COUNTS}" --output "${EXACT_COUNTS}"

"${PROJECT_PYTHON}" -m conference.aso_pipeline.item3_energy finalize \
    --score-table "${SCORE_TABLE}" --exact-counts "${EXACT_COUNTS}" \
    --raw-gene-counts "${RAW_GENE_COUNTS}" --output-dir "${ARTIFACT_DIR}"

{
    printf '[candidate]\n'
    "${INTARNA}" --version
    "${BOWTIE}" --version | sed -n '1,2p'
    printf 'Temperature: 37 C (310.15 K); RNA:RNA; IntaRNA default model\n'
    printf 'IntaRNA columns/settings: E, one best interaction; ASO=query, isolated target variant=target\n'
    printf 'Enumeration: complete unique Hamming ball <=3 substitutions; no indels/bulges\n'
    printf 'Threshold: variant_dG <= on_target_dG + RT*ln(fold), R=0.00198720425864083 kcal mol^-1 K^-1\n'
    printf 'Sweep folds: 1,3,10,30,100; primary=10; margins computed at 310.15 K\n'
    printf 'No-duplex convention: IntaRNA 0.0 plus explicit flag; always non-viable\n'
    printf 'Alignment: Bowtie 1 -f -a -v 0 -p %s; viable-100x union only\n' "${BOWTIE_THREADS}"
    printf 'Query orientation: reverse complement of genomic target variant\n'
    printf 'Gene count: union across viable variants; GENCODE v50 gene spans; transcript-compatible strand\n'
    printf 'Genome: /home/dashii/orcd/scratch/feng/gencode/v50/GRCh38.primary_assembly.genome.fa\n'
    printf 'GTF: /home/dashii/orcd/scratch/feng/gencode/v50/gencode.v50.primary_assembly.annotation.gtf.gz\n'
    printf 'Bowtie index: %s\n' "${BOWTIE_INDEX}"
    printf 'Raw comparator: %s (Bowtie -a -v 3, previously completed)\n' "${RAW_GENE_COUNTS}"
    printf '3UP8: excluded because the 8-mer <=3-mismatch search is saturated\n'
    printf 'Prediction predeclared before scoring: at 10x, F10 has fewer viable variants and fewer filtered genes than F18\n'
} > "${ARTIFACT_DIR}/item3_energy_provenance.txt.tmp"
mv "${ARTIFACT_DIR}/item3_energy_provenance.txt.tmp" \
    "${ARTIFACT_DIR}/item3_energy_provenance.txt"
