#!/usr/bin/env bash
#SBATCH --job-name=feng-bench-intarna
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --requeue
#SBATCH --output=feng/logs/%x-%j.out
#SBATCH --error=feng/logs/%x-%j.err

# Arm D frozen sequence scores (SCAFFOLD -- do not launch without sign-off).
#
# Computes the two preregistered ranking scores for every exon in the matched
# strata (positives + backgrounds):
#   * baseline: min 18-nt Hamming distance to nusinersen's reverse complement
#   * test:     best accessibility-aware IntaRNA interaction energy of
#               nusinersen against each exon's +/-100-nt transcript-strand window
#
# The IntaRNA sweep is a 100-400 CPU-hour job over the full universe/strata. It
# is a deterministic sequence function of the genome and the fixed nusinersen
# oligo -- it does NOT read any treatment alignment, contrast table, or sealed
# ΔPSI, and it does NOT compute the benchmark outcome (delta_r / permutation
# p-value). Those live in feng.bench.permute, sealed until the Phase 9 unlock.
#
# Freeze note: run this only to materialise the frozen score COLUMNS that the
# preregistration will hash. Do not run permute.run_phase9_test here.

set -euo pipefail

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
ENV_PREFIX="${FENG_ENV_PREFIX:-/home/dashii/orcd/pool/feng-env}"
PYTHON_BIN="${FENG_BENCH_PYTHON:-${ENV_PREFIX}/bin/python}"
SAMTOOLS="${FENG_BENCH_SAMTOOLS:-${ENV_PREFIX}/bin/samtools}"
INTARNA="${FENG_BENCH_INTARNA:-${ENV_PREFIX}/bin/IntaRNA}"
WINDOW="${FENG_BENCH_WINDOW:-100}"
N_SHARDS="${FENG_BENCH_N_SHARDS:-32}"

# Inputs: the matched strata (positives + backgrounds) and the reference FASTA.
STRATA="${FENG_BENCH_STRATA:-${SUBMIT_DIR}/feng/artifacts/bench/strata_v50.parquet}"
GENOME_FASTA="${FENG_BENCH_FASTA:-${SCRATCH:-}/feng/gencode/v50/GRCh38.primary_assembly.genome.fa}"
OUT_DIR="${FENG_BENCH_OUT_DIR:-${SCRATCH:-}/feng/bench/scores}"
SHARD_DIR="${OUT_DIR}/intarna_shards"

[[ -s "${STRATA}" ]] || { echo "ERROR: strata parquet not found: ${STRATA}" >&2; exit 1; }
[[ -s "${GENOME_FASTA}" ]] || { echo "ERROR: reference FASTA missing: ${GENOME_FASTA}" >&2; exit 1; }
[[ -x "${INTARNA}" || -f "${INTARNA}" ]] || { echo "ERROR: IntaRNA not found: ${INTARNA}" >&2; exit 1; }
mkdir -p "${OUT_DIR}" "${SHARD_DIR}" feng/logs

export PYTHONPATH="${SUBMIT_DIR}/feng/src:${PYTHONPATH:-}"

# Baseline (cheap) -- skip if already published.
if [[ -s "${OUT_DIR}/strata_baseline.parquet" && -f "${OUT_DIR}/strata_baseline.done" ]]; then
    echo "baseline already done; skipping"
else
    "${PYTHON_BIN}" -m feng.bench.scores \
        --events "${STRATA}" --genome-fasta "${GENOME_FASTA}" --samtools "${SAMTOOLS}" \
        --window "${WINDOW}" --score baseline \
        --out "${OUT_DIR}/strata_baseline.parquet"
    touch "${OUT_DIR}/strata_baseline.done"
fi

# IntaRNA (expensive, 100-400 CPU-h). Sharded with per-shard .done sentinels so
# a --requeue after preemption resumes at the shard boundary instead of
# restarting the whole sweep. Each shard writes atomically then drops a sentinel.
for (( shard=0; shard<N_SHARDS; shard++ )); do
    shard_out="${SHARD_DIR}/intarna_${shard}.parquet"
    shard_done="${SHARD_DIR}/intarna_${shard}.done"
    if [[ -s "${shard_out}" && -f "${shard_done}" ]]; then
        echo "intarna shard ${shard}/${N_SHARDS} already done; skipping"
        continue
    fi
    tmp_out="${shard_out}.tmp"
    "${PYTHON_BIN}" -m feng.bench.scores \
        --events "${STRATA}" --genome-fasta "${GENOME_FASTA}" --samtools "${SAMTOOLS}" \
        --window "${WINDOW}" --score intarna --intarna "${INTARNA}" \
        --n-shards "${N_SHARDS}" --shard-index "${shard}" \
        --out "${tmp_out}"
    mv -f "${tmp_out}" "${shard_out}"
    touch "${shard_done}"
    echo "intarna shard ${shard}/${N_SHARDS} done"
done

# Merge shards only when all are present.
"${PYTHON_BIN}" - "${SHARD_DIR}" "${N_SHARDS}" "${OUT_DIR}/strata_intarna.parquet" <<'PYEOF'
import sys, pandas as pd
from pathlib import Path
shard_dir, n_shards, out = Path(sys.argv[1]), int(sys.argv[2]), Path(sys.argv[3])
parts = []
for i in range(n_shards):
    p = shard_dir / f"intarna_{i}.parquet"
    if not p.is_file():
        sys.exit(f"ERROR: shard {i} missing at {p}; not merging")
    parts.append(pd.read_parquet(p))
pd.concat(parts, ignore_index=True).to_parquet(out, index=False)
print(f"merged {n_shards} shards -> {out}")
PYEOF

echo "Arm D score columns written to ${OUT_DIR}"
