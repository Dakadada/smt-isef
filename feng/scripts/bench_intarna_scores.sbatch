#!/usr/bin/env bash
#SBATCH --job-name=feng-bench-intarna
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --array=0-31%16
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

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    [[ "${SLURM_ARRAY_TASK_COUNT}" -eq "${N_SHARDS}" ]] || {
        echo "ERROR: array task count ${SLURM_ARRAY_TASK_COUNT} must equal N_SHARDS ${N_SHARDS}" >&2
        exit 1
    }
    SHARD_FIRST="${SLURM_ARRAY_TASK_ID}"
    SHARD_STEP="${SLURM_ARRAY_TASK_COUNT}"
else
    # Manual/non-array execution remains a supported recovery path.
    SHARD_FIRST=0
    SHARD_STEP=1
fi

# Inputs: the matched strata (positives + backgrounds) and the reference FASTA.
STRATA="${FENG_BENCH_STRATA:-${SUBMIT_DIR}/feng/artifacts/bench/strata_v50.parquet}"
GENOME_FASTA="${FENG_BENCH_FASTA:-${SCRATCH:-}/feng/gencode/v50/GRCh38.primary_assembly.genome.fa}"
OUT_DIR="${FENG_BENCH_OUT_DIR:-${SCRATCH:-}/feng/bench/scores}"

[[ -s "${STRATA}" ]] || { echo "ERROR: strata parquet not found: ${STRATA}" >&2; exit 1; }
[[ -s "${GENOME_FASTA}" ]] || { echo "ERROR: reference FASTA missing: ${GENOME_FASTA}" >&2; exit 1; }
[[ -x "${INTARNA}" || -f "${INTARNA}" ]] || { echo "ERROR: IntaRNA not found: ${INTARNA}" >&2; exit 1; }
mkdir -p "${OUT_DIR}" feng/logs

export PYTHONPATH="${SUBMIT_DIR}/feng/src:${PYTHONPATH:-}"

# Content-address every cache entry to all scientific inputs and score settings.
# A changed strata table, reference, external executable, Python scorer source,
# setting, or oligo therefore selects a fresh namespace instead of stale .done files.
CACHE_KEY="$("${PYTHON_BIN}" - "${STRATA}" "${GENOME_FASTA}" "${INTARNA}" "${SAMTOOLS}" "${WINDOW}" "${N_SHARDS}" <<'PYEOF'
from pathlib import Path
import sys
import feng.bench.scores as scores
print(scores.score_cache_key(
    {"strata": Path(sys.argv[1]), "reference_fasta": Path(sys.argv[2]),
     "intarna_executable": Path(sys.argv[3]),
     "samtools_executable": Path(sys.argv[4]),
     "scorer_source": Path(scores.__file__)},
    {"window": int(sys.argv[5]), "n_shards": int(sys.argv[6]),
     "query_5to3": scores.NUSINERSEN_5P_TO_3P, "cache_schema": 2},
))
PYEOF
)"
CACHE_DIR="${OUT_DIR}/cache/${CACHE_KEY}"
SHARD_DIR="${CACHE_DIR}/intarna_shards"
mkdir -p "${CACHE_DIR}" "${SHARD_DIR}"
echo "score cache key: ${CACHE_KEY}"

# Baseline (cheap) -- skip if already published.
if [[ -s "${CACHE_DIR}/strata_baseline.parquet" && -f "${CACHE_DIR}/strata_baseline.done" ]]; then
    echo "baseline already done; skipping"
else
    baseline_tmp="$(mktemp "${CACHE_DIR}/.strata_baseline.tmp.XXXXXXXX.parquet")"
    "${PYTHON_BIN}" -m feng.bench.scores \
        --events "${STRATA}" --genome-fasta "${GENOME_FASTA}" --samtools "${SAMTOOLS}" \
        --window "${WINDOW}" --score baseline \
        --out "${baseline_tmp}"
    (
        flock -x 9
        if [[ -s "${CACHE_DIR}/strata_baseline.parquet" && -f "${CACHE_DIR}/strata_baseline.done" ]]; then
            rm -f "${baseline_tmp}"
        else
            mv -f "${baseline_tmp}" "${CACHE_DIR}/strata_baseline.parquet"
            baseline_done_tmp="$(mktemp "${CACHE_DIR}/.strata_baseline.done.XXXXXXXX")"
            mv -f "${baseline_done_tmp}" "${CACHE_DIR}/strata_baseline.done"
        fi
    ) 9>"${CACHE_DIR}/strata_baseline.lock"
fi

# IntaRNA (expensive, 100-400 CPU-h). A completed shard is accepted only after
# its full source slice and finite score column validate. Unique candidate names
# plus a per-shard lock make duplicate launches safe.
validate_shard() {
    "${PYTHON_BIN}" - "$1" "${STRATA}" "${N_SHARDS}" "$2" <<'PYEOF'
import math, sys
from pathlib import Path
import pandas as pd

candidate, source = Path(sys.argv[1]), Path(sys.argv[2])
n_shards, shard = int(sys.argv[3]), int(sys.argv[4])
expected = pd.read_parquet(source).iloc[shard::n_shards].reset_index(drop=True)
actual = pd.read_parquet(candidate).reset_index(drop=True)
missing = set(expected.columns) - set(actual.columns)
if missing or "intarna_min_energy" not in actual:
    raise SystemExit(f"invalid shard schema; missing={sorted(missing)}")
pd.testing.assert_frame_equal(
    actual[list(expected.columns)], expected, check_dtype=False, check_like=False
)
scores = pd.to_numeric(actual["intarna_min_energy"], errors="coerce")
if scores.isna().any() or not scores.map(math.isfinite).all():
    raise SystemExit("invalid shard: IntaRNA scores must all be finite")
PYEOF
}

for (( shard=SHARD_FIRST; shard<N_SHARDS; shard+=SHARD_STEP )); do
    shard_out="${SHARD_DIR}/intarna_${shard}.parquet"
    shard_done="${SHARD_DIR}/intarna_${shard}.done"
    if [[ -s "${shard_out}" && -f "${shard_done}" ]] && validate_shard "${shard_out}" "${shard}"; then
        echo "intarna shard ${shard}/${N_SHARDS} already done; skipping"
        continue
    fi
    tmp_out="$(mktemp "${SHARD_DIR}/.intarna_${shard}.tmp.XXXXXXXX.parquet")"
    if ! "${PYTHON_BIN}" -m feng.bench.scores \
        --events "${STRATA}" --genome-fasta "${GENOME_FASTA}" --samtools "${SAMTOOLS}" \
        --window "${WINDOW}" --score intarna --intarna "${INTARNA}" \
        --n-shards "${N_SHARDS}" --shard-index "${shard}" \
        --out "${tmp_out}"; then
        rm -f "${tmp_out}"
        exit 1
    fi
    validate_shard "${tmp_out}" "${shard}"
    (
        flock -x 9
        if [[ -s "${shard_out}" && -f "${shard_done}" ]] && validate_shard "${shard_out}" "${shard}"; then
            rm -f "${tmp_out}"
        else
            mv -f "${tmp_out}" "${shard_out}"
            done_tmp="$(mktemp "${SHARD_DIR}/.intarna_${shard}.done.XXXXXXXX")"
            mv -f "${done_tmp}" "${shard_done}"
        fi
    ) 9>"${SHARD_DIR}/intarna_${shard}.lock"
    echo "intarna shard ${shard}/${N_SHARDS} done"
done

# Array tasks that finish before the final shard publish checkpoint and exit.
# The last arriving task observes every validated shard sentinel and performs
# the merge; duplicate finalizers are harmless because publication is locked.
for (( shard=0; shard<N_SHARDS; shard++ )); do
    shard_out="${SHARD_DIR}/intarna_${shard}.parquet"
    shard_done="${SHARD_DIR}/intarna_${shard}.done"
    if [[ ! -s "${shard_out}" || ! -f "${shard_done}" ]] || ! validate_shard "${shard_out}" "${shard}"; then
        echo "checkpoint complete; waiting for remaining shards"
        exit 0
    fi
done

# Merge shards to a per-writer candidate only when all are present and valid.
scored_tmp="$(mktemp "${CACHE_DIR}/.scored_strata.tmp.XXXXXXXX.parquet")"
"${PYTHON_BIN}" - "${SHARD_DIR}" "${N_SHARDS}" \
    "${CACHE_DIR}/strata_baseline.parquet" "${scored_tmp}" <<'PYEOF'
import sys, pandas as pd
from pathlib import Path
from feng.bench.scores import join_score_tables
shard_dir, n_shards = Path(sys.argv[1]), int(sys.argv[2])
baseline_path, out = Path(sys.argv[3]), Path(sys.argv[4])
parts = []
for i in range(n_shards):
    p = shard_dir / f"intarna_{i}.parquet"
    if not p.is_file():
        sys.exit(f"ERROR: shard {i} missing at {p}; not merging")
    parts.append(pd.read_parquet(p))
test = pd.concat(parts, ignore_index=True)
joined = join_score_tables(pd.read_parquet(baseline_path), test)
joined.to_parquet(out, index=False)
print(f"merged {n_shards} shards and joined {len(joined)} rows 1-to-1 -> {out}")
PYEOF

# Serialize creation of the immutable parquet+manifest cache package. A second
# launcher validates and reuses the completed package instead of overwriting it.
(
flock -x 9
"${PYTHON_BIN}" - "${CACHE_DIR}" "${CACHE_KEY}" "${scored_tmp}" <<'PYEOF'
from pathlib import Path
import hashlib, json, os, sys, tempfile
cache_dir, cache_key, candidate = Path(sys.argv[1]), sys.argv[2], Path(sys.argv[3])
parquet = cache_dir / "scored_strata.parquet"
manifest_path = cache_dir / "manifest.json"
if parquet.is_file() and manifest_path.is_file():
    existing = json.loads(manifest_path.read_text(encoding="utf-8"))
    digest = hashlib.sha256(parquet.read_bytes()).hexdigest()
    if existing == {"cache_key": cache_key, "parquet": parquet.name,
                    "parquet_sha256": digest}:
        candidate.unlink()
        raise SystemExit(0)
os.replace(candidate, parquet)
manifest = {
    "cache_key": cache_key,
    "parquet": parquet.name,
    "parquet_sha256": hashlib.sha256(parquet.read_bytes()).hexdigest(),
}
fd, tmp_name = tempfile.mkstemp(prefix="manifest.", suffix=".tmp", dir=cache_dir)
try:
    with os.fdopen(fd, "w", encoding="utf-8") as handle:
        json.dump(manifest, handle, sort_keys=True)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp_name, manifest_path)
finally:
    if os.path.exists(tmp_name):
        os.unlink(tmp_name)
PYEOF
) 9>"${CACHE_DIR}/publish.lock"

# Publish by one atomic symlink switch. The immutable cache directory contains
# both parquet and manifest, so readers can never pair new scores with an old key.
"${PYTHON_BIN}" - "${OUT_DIR}" "${CACHE_KEY}" <<'PYEOF'
from pathlib import Path
import os, secrets, sys

out_dir, cache_key = Path(sys.argv[1]), sys.argv[2]
current = out_dir / "current"
tmp = out_dir / f".current.tmp.{os.getpid()}.{secrets.token_hex(8)}"
try:
    os.symlink(f"cache/{cache_key}", tmp)
    os.replace(tmp, current)
finally:
    if tmp.is_symlink():
        tmp.unlink()
PYEOF

echo "Arm D score columns written to ${OUT_DIR}"
