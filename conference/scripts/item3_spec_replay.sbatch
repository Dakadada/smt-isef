#!/usr/bin/env bash
#SBATCH --job-name=conf-item3-spec
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=2-00:00:00
#SBATCH --requeue
#SBATCH --signal=B:SIGTERM@180
#SBATCH --output=conference/logs/%x-%j.out
#SBATCH --error=conference/logs/%x-%j.err

set -euo pipefail

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
PYTHON="/home/dashii/orcd/pool/feng-env/bin/python"
RNAPLFOLD="/home/dashii/orcd/pool/feng-env/bin/RNAplfold"
SAMTOOLS="/home/dashii/orcd/pool/feng-env/bin/samtools"
GENOME="/home/dashii/orcd/scratch/feng/gencode/v50/GRCh38.primary_assembly.genome.fa"
WORK_DIR="/orcd/scratch/orcd/011/dashii/conference_smn2_aso"
BOWTIE="${WORK_DIR}/tools/bin/bowtie"
INDEX="${WORK_DIR}/indexes/grch38_v50"
GENES="${WORK_DIR}/gencode_v50_gene_spans.tsv"
ARTIFACTS="${SUBMIT_DIR}/conference/artifacts"
RUN_DIR="${WORK_DIR}/results/item3_spec_v1"
QUERY_FASTA="${RUN_DIR}/item3_spec_queries.fa"
LOCI="${RUN_DIR}/item3_spec_loci.tsv"
SHARDS="${RUN_DIR}/accessibility_shards"
ACCESSIBLE_LOCI="${RUN_DIR}/item3_spec_loci_accessibility.tsv"
SCORES="${ARTIFACTS}/item3_energy_variant_scores.tsv"

export PYTHONPATH="${SUBMIT_DIR}/feng/src:${SUBMIT_DIR}:${PYTHONPATH:-}"
export SCRATCH="/orcd/scratch/orcd/011/dashii"
mkdir -p "${RUN_DIR}" "${SHARDS}" "${ARTIFACTS}"

for required in "${PYTHON}" "${RNAPLFOLD}" "${SAMTOOLS}" "${GENOME}" \
    "${BOWTIE}" "${INDEX}.1.ebwt" "${INDEX}.rev.2.ebwt" "${GENES}" \
    "${SCORES}" "${ARTIFACTS}/item1_dg_raw.tsv"; do
    [[ -f "${required}" ]] || { echo "missing input: ${required}" >&2; exit 2; }
done

if [[ ! -f "${LOCI}.done" ]]; then
    "${PYTHON}" -m conference.aso_pipeline.item3_spec prepare \
        --scores "${SCORES}" --output "${QUERY_FASTA}"
    bash "${SUBMIT_DIR}/conference/scripts/bowtie_stream.sh" \
        "${BOWTIE}" "${INDEX}" "${QUERY_FASTA}" "${SLURM_CPUS_PER_TASK}" 0 \
        | "${PYTHON}" -m conference.aso_pipeline.item3_spec map \
            --scores "${SCORES}" --gene-table "${GENES}" > "${LOCI}.tmp"
    mv "${LOCI}.tmp" "${LOCI}"
    touch "${LOCI}.done"
fi

# Checkpoint accessibility in bounded shards. Header is copied into every shard.
if [[ ! -f "${SHARDS}/split.done" ]]; then
    rm -f "${SHARDS}"/part-*
    tail -n +2 "${LOCI}" | split -l 10000 - "${SHARDS}/part-"
    header="$(head -n 1 "${LOCI}")"
    for part in "${SHARDS}"/part-*; do
        { printf '%s\n' "${header}"; sed -n '1,$p' "${part}"; } > "${part}.tsv"
        rm -f "${part}"
    done
    touch "${SHARDS}/split.done"
fi

for shard in "${SHARDS}"/part-*.tsv; do
    output="${shard%.tsv}.accessible.tsv"
    [[ -f "${output}.done" ]] && continue
    "${PYTHON}" -m conference.aso_pipeline.item3_spec accessibility \
        --loci "${shard}" --output "${output}.tmp" \
        --rnaplfold "${RNAPLFOLD}" --genome-fasta "${GENOME}" \
        --samtools "${SAMTOOLS}" --workers "${SLURM_CPUS_PER_TASK}"
    mv "${output}.tmp" "${output}"
    touch "${output}.done"
done

accessible_shards=("${SHARDS}"/part-*.accessible.tsv)
[[ -e "${accessible_shards[0]}" ]] || { echo "no accessibility shards" >&2; exit 2; }
head -n 1 "${accessible_shards[0]}" > "${ACCESSIBLE_LOCI}.tmp"
for shard in "${accessible_shards[@]}"; do
    tail -n +2 "${shard}" >> "${ACCESSIBLE_LOCI}.tmp"
done
mv "${ACCESSIBLE_LOCI}.tmp" "${ACCESSIBLE_LOCI}"

"${PYTHON}" -m conference.aso_pipeline.item3_spec finalize \
    --loci "${ACCESSIBLE_LOCI}" --item1 "${ARTIFACTS}/item1_dg_raw.tsv" \
    --output-dir "${ARTIFACTS}"

{
    printf '[candidate]\n'
    printf 'command: sbatch conference/scripts/item3_spec_replay.sbatch\n'
    printf 'R=1.987204e-3 kcal mol^-1 K^-1; T=310.15 K\n'
    printf 'concentrations_M=1e-9,1e-8,1e-7,1e-6\n'
    printf 'absolute_dG_cutoffs_kcal_mol=-8,-10,-12,-14,-16\n'
    printf 'accessibility_cutoff=0.5 mean per-base unpaired probability\n'
    printf 'RNAplfold defaults reused from feng.titration.evaluate: pad=70,W=80,L=40\n'
    printf 'expression_state=unavailable; cited evidence.json contains event PSI, not transcriptome-wide expression\n'
    printf '3UP8=saturated; no numeric metric fabricated\n'
    printf 'random_seed=none; deterministic enumeration, exact alignment, and formulas\n'
    "${BOWTIE}" --version | sed -n '1,2p'
    "${RNAPLFOLD}" --version | sed -n '1p'
} > "${ARTIFACTS}/item3_spec_provenance.txt.tmp"
mv "${ARTIFACTS}/item3_spec_provenance.txt.tmp" "${ARTIFACTS}/item3_spec_provenance.txt"
