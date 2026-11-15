#!/usr/bin/env bash
#SBATCH --job-name=feng-armc-scan
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
# Wall clock (REDESIGN, bounded pool): each SBO is scanned against a fixed
# ~20k-exon background pool (+ direct positives), NOT the whole ~1.1M-exon
# GENCODE annotation. Measured on random genome-wide GENCODE v50 exons at
# max_edit=5: ~255 ms/exon (the balanced pigeonhole seeds are only 3 nt at
# budget 5 and prune little; genome-wide exons run ~2x the canary's
# detectable-subset rate because they include long, repeat-rich exons). So
# pool=20000 projects to ~85 min/SBO on average; 4 h leaves ~3x headroom for
# repeat-rich outliers and lets a preemption-requeued task (the whole SBO scan
# restarts from scratch -- no intra-scan checkpoint) finish one uninterrupted
# pass. The site-level redesign does NOT raise this: the dominant cost is the
# seed_and_extend scan over the ~20k-exon pool (unchanged), and IntaRNA now runs
# on the short best-binding WINDOW (~18-25 nt target) instead of the whole exon,
# so it is if anything cheaper; --time stays 04:00:00. RUN A 1-SBO CANARY FIRST.
# Retune FENG_ARMC_POOL_SIZE and this --time together if the pool size changes
# (time scales linearly in pool size).
#SBATCH --time=08:00:00
#SBATCH --requeue
#SBATCH --signal=B:SIGTERM@60
#SBATCH --array=1-81%16
#SBATCH --output=feng/logs/%x-%A_%a.out
#SBATCH --error=feng/logs/%x-%A_%a.err

# Arm C bind->act gate, stage 1: one SLURM array task per SBO (81 total).
# Each task: (a) build the shared TARGET exonome FASTA cache exactly once, behind
# a flock so concurrent tasks don't race the first build -- the target is a
# fixed, deterministic random pool of GENCODE v50 exons UNION the LIFTED
# (hg19->GRCh38) direct-bound S4 positive loci (see feng.bench.armc_build); the
# build ALSO runs the >=90% positive-LIFT tripwire and hard-fails (nonzero) if
# the positives did not lift (wrong/absent chain), aborting the pipeline before
# any scan is spent; (b) scan its own SBO (revcomp'd) against the cached exonome
# at max-edit 5; (c) label bound targets against Holgersen S1/S4 (positive =
# direct-bound-and-acted, numeric edit<=5) with SITE-LEVEL features and write ONE
# per-SBO features parquet + .done sentinel. Idempotent: a re-run
# (incl. after a preemption requeue) skips any dg_id whose output + sentinel
# already exist. See feng/scripts/armc_gate.sbatch for the merge + LOSO gate step
# (dependency=afterok on this array).

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
SAMTOOLS="${FENG_BENCH_SAMTOOLS:-${ENV_PREFIX}/bin/samtools}"
INTARNA="${FENG_BENCH_INTARNA:-${ENV_PREFIX}/bin/IntaRNA}"
MAX_EDIT="${FENG_ARMC_MAX_EDIT:-5}"

SBO_TSV="${FENG_ARMC_SBO_TSV:-${SUBMIT_DIR}/feng/artifacts/phase0_arm_c_probe/sbo_sequences.tsv}"
S1_PDF="${FENG_ARMC_S1_PDF:-${SUBMIT_DIR}/feng/sj-pdf-1-nat-10.1089_nat.2020.0921-2.pdf}"
S4_PDF="${FENG_ARMC_S4_PDF:-${SUBMIT_DIR}/feng/sj-pdf-4-nat-10.1089_nat.2020.0921.pdf}"
GTF="${FENG_ARMC_GTF:-${SCRATCH:-}/feng/gencode/v50/gencode.v50.primary_assembly.annotation.gtf}"
GENOME_FASTA="${FENG_ARMC_FASTA:-${SCRATCH:-}/feng/gencode/v50/GRCh38.primary_assembly.genome.fa}"
CHAIN="${FENG_ARMC_CHAIN:-${SCRATCH:-}/feng/liftover/hg19ToHg38.over.chain.gz}"
POOL_SIZE="${FENG_ARMC_POOL_SIZE:-20000}"
POOL_SEED="${FENG_ARMC_POOL_SEED:-20260903}"
TARGET_PAD="${FENG_ARMC_PAD:-30}"
ARMC_DIR="${FENG_ARMC_DIR:-${SCRATCH:-}/feng/armc}"
EXONOME_FASTA="${ARMC_DIR}/exonome.fa"
HITS_DIR="${ARMC_DIR}/hits"

[[ -s "${SBO_TSV}" ]] || { echo "ERROR: sbo_sequences.tsv missing: ${SBO_TSV}" >&2; exit 1; }
[[ -s "${S1_PDF}" ]] || { echo "ERROR: S1 PDF missing: ${S1_PDF}" >&2; exit 1; }
[[ -s "${S4_PDF}" ]] || { echo "ERROR: S4 PDF missing: ${S4_PDF}" >&2; exit 1; }
[[ -s "${GTF}" ]] || { echo "ERROR: GENCODE GTF missing: ${GTF}" >&2; exit 1; }
[[ -s "${GENOME_FASTA}" ]] || { echo "ERROR: reference FASTA missing: ${GENOME_FASTA}" >&2; exit 1; }
[[ -s "${CHAIN}" ]] || { echo "ERROR: hg19->GRCh38 chain missing: ${CHAIN}" >&2; exit 1; }
export FENG_ARMC_CHAIN="${CHAIN}"
mkdir -p "${ARMC_DIR}" "${HITS_DIR}" feng/logs

export PYTHONPATH="${SUBMIT_DIR}/feng/src:${PYTHONPATH:-}"

TASK="${SLURM_ARRAY_TASK_ID:?ERROR: this script must run as a SLURM array job}"

# --- requeue-safe: on SIGTERM (preemption warning) or SIGUSR1, ask SLURM to
# requeue THIS array task rather than letting it die uncleanly. Safe to call
# repeatedly / on an already-finished task: everything below is idempotent.
requeue_handler() {
    echo "caught termination signal; requesting requeue of job ${SLURM_JOB_ID:-unknown}" >&2
    if [[ -n "${SLURM_JOB_ID:-}" ]]; then
        scontrol requeue "${SLURM_JOB_ID}" || true
    fi
    exit 0
}
trap requeue_handler SIGTERM SIGUSR1

# --- resolve this task's dg_id (row N of sbo_sequences.tsv, header skipped) --
DG_ID="$(tail -n +2 "${SBO_TSV}" | sed -n "${TASK}p" | cut -f1)"
[[ -n "${DG_ID}" ]] || { echo "ERROR: no dg_id for array task ${TASK}" >&2; exit 1; }

OUT="${HITS_DIR}/${DG_ID}.parquet"
DONE="${HITS_DIR}/${DG_ID}.done"
if [[ -s "${OUT}" && -f "${DONE}" ]]; then
    echo "dg_id ${DG_ID} (task ${TASK}) already done; skipping"
    exit 0
fi

# --- build (or verify) the shared target exonome cache, serialized via flock.
# The build runs the >=90% positive-coverage tripwire and exits nonzero if the
# positive class is empty; `set -e` then aborts this task (and, since every task
# hits the same tripwire, the whole array) before any scan compute is spent.
LOCK="${ARMC_DIR}/exonome.lock"
(
    flock -x 200
    if [[ ! -s "${EXONOME_FASTA}" ]]; then
        echo "task ${TASK} (dg_id ${DG_ID}) building target exonome (pool=${POOL_SIZE}) -> ${EXONOME_FASTA}"
        "${PYTHON_BIN}" -m feng.bench.armc_build build-exonome \
            --gtf "${GTF}" --genome-fasta "${GENOME_FASTA}" \
            --cache-path "${EXONOME_FASTA}" \
            --s1-pdf "${S1_PDF}" --s4-pdf "${S4_PDF}" --chain "${CHAIN}" \
            --pool-size "${POOL_SIZE}" --pool-seed "${POOL_SEED}" \
            --max-edit "${MAX_EDIT}" --pad "${TARGET_PAD}" --samtools "${SAMTOOLS}"
    fi
) 200>"${LOCK}"

# --- scan this SBO against the cached exonome, write its labeled feature rows -
echo "task ${TASK}: scanning dg_id ${DG_ID}"
"${PYTHON_BIN}" -m feng.bench.armc_build scan \
    --dg-id "${DG_ID}" \
    --sbo-tsv "${SBO_TSV}" --s1-pdf "${S1_PDF}" --s4-pdf "${S4_PDF}" \
    --exonome-cache "${EXONOME_FASTA}" --chain "${CHAIN}" \
    --intarna "${INTARNA}" --max-edit "${MAX_EDIT}" \
    --out "${OUT}"
touch "${DONE}"
echo "task ${TASK}: dg_id ${DG_ID} done -> ${OUT}"

trap - SIGTERM SIGUSR1
