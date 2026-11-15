#!/usr/bin/env bash
#SBATCH --job-name=feng-bench-universe
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --requeue
#SBATCH --output=feng/logs/%x-%j.out
#SBATCH --error=feng/logs/%x-%j.err

# Arm D control-only detectable-exon universe.
#
# Enumerates all transcript-internal exons from GENCODE v50 and scores each for
# baseline junction depth and inclusion level from STAR junction tables;
# rMATS-turbo SE (skipped-exon) events remain source-event provenance but do not
# gate candidate identity. Counts use ONLY the three scrambled-MOE (ScrMOE)
# alignments from GSE307431. Both
# rMATS sample groups are drawn from the control pool and the run uses
# --statoff, which disables statistical testing. rMATS may still emit a
# dummy-control IncLevelDifference; the builder refuses all contrast columns.
# The two groups exist only to satisfy the rMATS two-argument interface.
#
# Freeze discipline: no F18OMe / F18MOE / F20PMO alignment is read, and no
# treatment effect table is consumed. The Python builder re-checks the control
# selection and aborts on any treatment leakage.
#
# Requeue safety: work happens in a job-scoped directory and is atomically
# published with a .done sentinel only after outputs validate; a requeue that
# finds the sentinel exits immediately.

set -euo pipefail

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
CONDITIONS="${FENG_BENCH_CONDITIONS:-${SUBMIT_DIR}/feng/artifacts/phase3/sample_conditions.tsv}"
ENV_PREFIX="${FENG_ENV_PREFIX:-/home/dashii/orcd/pool/feng-env}"
PYTHON_BIN="${FENG_BENCH_PYTHON:-${ENV_PREFIX}/bin/python}"
RMATS="${FENG_BENCH_RMATS:-${ENV_PREFIX}/bin/rmats.py}"
SAMTOOLS="${FENG_BENCH_SAMTOOLS:-${ENV_PREFIX}/bin/samtools}"
THREADS="${SLURM_CPUS_PER_TASK:-8}"
READ_LENGTH="${FENG_BENCH_READ_LENGTH:-100}"
LIB_TYPE="${FENG_BENCH_LIBTYPE:-fr-firststrand}"
PKG_PATH="${SUBMIT_DIR}/feng/src"

if [[ -z "${SCRATCH:-}" && ( -z "${FENG_BENCH_STAR_DIR:-}" || -z "${FENG_BENCH_OUT_DIR:-}" || -z "${FENG_BENCH_GTF:-}" ) ]]; then
    echo "ERROR: set SCRATCH, or explicitly set FENG_BENCH_STAR_DIR, FENG_BENCH_OUT_DIR, and FENG_BENCH_GTF" >&2
    exit 1
fi
STAR_DIR="${FENG_BENCH_STAR_DIR:-${SCRATCH}/feng/phase3/star}"
OUT_DIR="${FENG_BENCH_OUT_DIR:-${SCRATCH}/feng/bench/detectable_universe}"
GENOME_FASTA="${FENG_BENCH_FASTA:-${SCRATCH}/feng/gencode/v50/GRCh38.primary_assembly.genome.fa}"

GTF="${FENG_BENCH_GTF:-${SCRATCH}/feng/gencode/v50/annotation.gtf}"
if [[ ! -s "${GTF}" && -z "${FENG_BENCH_GTF:-}" ]]; then
    GTF_FALLBACK="${SCRATCH}/feng/gencode/v50/gencode.v50.primary_assembly.annotation.gtf"
    if [[ -s "${GTF_FALLBACK}" ]]; then
        echo "WARNING: ${GTF} is absent; using ${GTF_FALLBACK}" >&2
        GTF="${GTF_FALLBACK}"
    fi
fi
[[ -s "${GTF}" ]] || { echo "ERROR: GTF not found: ${GTF}" >&2; exit 1; }
[[ -s "${CONDITIONS}" ]] || { echo "ERROR: conditions table missing: ${CONDITIONS}" >&2; exit 1; }
[[ -s "${GENOME_FASTA}" ]] || { echo "ERROR: reference FASTA missing: ${GENOME_FASTA}" >&2; exit 1; }
[[ -x "${SAMTOOLS}" ]] || { echo "ERROR: samtools missing: ${SAMTOOLS}" >&2; exit 1; }
[[ "${READ_LENGTH}" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: invalid read length: ${READ_LENGTH}" >&2; exit 1; }

# Repo-relative artifact outputs (committed).
ARTIFACT_DIR="${SUBMIT_DIR}/feng/artifacts/bench"
BUNDLE_ROOT="${ARTIFACT_DIR}/detectable_universe_v50.bundles"
CURRENT_POINTER="${ARTIFACT_DIR}/detectable_universe_v50.current"
OUT_PARQUET="${CURRENT_POINTER}/detectable_universe_v50.parquet"
OUT_README="${CURRENT_POINTER}/detectable_universe_v50_README.md"

# Small provenance inputs are committed into the repo (the multi-gigabyte SE
# table stays on scratch). The published README and the .done validator both
# hash these committed copies so a fresh checkout round-trips without scratch.
PROV_DIR="${ARTIFACT_DIR}/provenance"
COMMIT_CONDITIONS="${SUBMIT_DIR}/feng/artifacts/phase3/sample_conditions.tsv"
COMMIT_B1="${PROV_DIR}/b1.txt"
COMMIT_B2="${PROV_DIR}/b2.txt"

mkdir -p "${OUT_DIR}" "${ARTIFACT_DIR}" "${BUNDLE_ROOT}" "${PROV_DIR}" "${SUBMIT_DIR}/feng/logs"

# --- verify the immutable control allowlist before any BAM is opened ---------
PINNED_CONTROL_SHA256="295328ba55a16d2d009c0ff69657a2d4ab8547cbe254917bcc50b4344b4a412e"
CONTROL_RUNS=(
    SRR35298704
    SRR35298699
    SRR35298715
)
actual_allowlist_hash="$(printf '%s\n' "${CONTROL_RUNS[@]}" | sha256sum | awk '{print $1}')"
[[ "${actual_allowlist_hash}" == "${PINNED_CONTROL_SHA256}" ]] || {
    echo "ERROR: pinned control allowlist hash mismatch: ${actual_allowlist_hash}" >&2
    exit 1
}
[[ "${#CONTROL_RUNS[@]}" -eq 3 ]] || {
    echo "ERROR: pinned control allowlist must contain exactly 3 ScrMOE runs" >&2
    exit 1
}
[[ "$(printf '%s\n' "${CONTROL_RUNS[@]}" | sort | uniq -d | wc -l)" -eq 0 ]] || {
    echo "ERROR: duplicate accession in pinned control allowlist" >&2
    exit 1
}

# Derive the metadata selection only as a cross-check against the pinned IDs.
readarray -t METADATA_CONTROLS < <(
    awk -F '\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) { col[$i] = i }
            next
        }
        $col["condition"] == "control" &&
        $col["control_type"] == "scrambled" && $col["chemistry"] == "MOE" {
            print $col["run_accession"]
        }
    ' "${CONDITIONS}"
)
if [[ "$(printf '%s\n' "${METADATA_CONTROLS[@]}" | sort | uniq -d | wc -l)" -ne 0 ]]; then
    echo "ERROR: duplicate control accession in conditions table" >&2
    exit 1
fi
if ! diff -u \
    <(printf '%s\n' "${CONTROL_RUNS[@]}" | sort) \
    <(printf '%s\n' "${METADATA_CONTROLS[@]}" | sort); then
    echo "ERROR: conditions ScrMOE controls differ from pinned 3-accession allowlist" >&2
    exit 1
fi
echo "Pinned control samples (${#CONTROL_RUNS[@]}, sha256=${PINNED_CONTROL_SHA256}): ${CONTROL_RUNS[*]}"

# Guard: refuse any run that is not classed control in the table.
for run in "${CONTROL_RUNS[@]}"; do
    klass="$(awk -F '\t' -v r="${run}" '
        NR == 1 { for (i=1;i<=NF;i++) col[$i]=i; next }
        $col["run_accession"] == r { print $col["condition"] "/" $col["control_type"] }
    ' "${CONDITIONS}")"
    case "${klass}" in
        control/scrambled) : ;;
        *) echo "ERROR: refusing non-control run ${run} (${klass})" >&2; exit 1 ;;
    esac
done

# Stable SE table (large; stays on scratch) named in the README and hashed by
# the .done guard. The small conditions/b1/b2 inputs are hashed from their
# committed repo copies (COMMIT_*) instead.
STABLE_SE_TABLE="${OUT_DIR}/SE.MATS.JC.control-only.txt"

validate_published_artifact() {
    local parquet_path="${1:-${OUT_PARQUET}}"
    local readme_path="${2:-${OUT_README}}"
    PYTHONPATH="${PKG_PATH}${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON_BIN}" - \
        "${parquet_path}" "${readme_path}" "${STABLE_SE_TABLE}" \
        "${COMMIT_CONDITIONS}" "${COMMIT_B1}" "${COMMIT_B2}" \
        "${GTF}" "${GENOME_FASTA}" "${STAR_DIR}" "${READ_LENGTH}" \
        "${CONTROL_RUNS[@]}" <<'PY'
from pathlib import Path
import hashlib
import json
import re
import sys

import pandas as pd
from feng.bench.strata import (
    assert_scrmoe_universe,
    validate_bundle_manifest,
    validate_input_provenance,
)

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()

# Contrast/statistic fields must never reach the published universe (freeze
# discipline): kept in sync with feng.bench.strata.CONTRAST_FIELDS.
CONTRAST_FIELDS = {"PValue", "FDR", "IncLevelDifference"}

parquet, readme, se_table, conditions, b1, b2, gtf, fasta, star_dir = map(Path, sys.argv[1:10])
read_length = int(sys.argv[10])
expected_runs = set(sys.argv[11:])
sj_paths = {
    run: star_dir / run / f"{run}.SJ.out.tab" for run in expected_runs
}
count_paths = {
    run: star_dir / run / f"{run}.ReadsPerGene.out.tab" for run in expected_runs
}
required_files = (
    parquet, readme, se_table, conditions, b1, b2, gtf, fasta,
    *sj_paths.values(), *count_paths.values(),
)
missing = [str(path) for path in required_files if not path.is_file() or path.stat().st_size == 0]
if missing:
    raise SystemExit("missing/empty published artifact input: " + ", ".join(missing))

frame = pd.read_parquet(parquet)
validate_bundle_manifest(parquet.parent)
assert_scrmoe_universe(frame)
required_columns = {
    "exon_key", "chrom", "exon_start_0base", "exon_end", "strand",
    "detectable", "gc_fraction", "rmats_event_ids", "n_rmats_events",
    "all_gene_ids", "all_gene_id_bases", "all_gene_symbols", "n_gene_ids",
    "multi_gene", "control_mean_cpm_by_gene",
    "annotation_event_ids", "n_annotation_events",
    "gencode_gtf_sha256", "scrmoe_sj_sha256", "junction_read_length",
    "junction_read_length_sha256", "genome_fasta_sha256",
    "scrmoe_reads_per_gene_sha256",
}
missing_columns = sorted(required_columns.difference(frame.columns))
if missing_columns:
    raise SystemExit(f"published parquet missing schema columns: {missing_columns}")
# Freeze guard: the .done fast-path must reject a stale/tampered parquet that
# carries any treatment contrast column, exactly as the fresh build path does.
leaked = sorted(CONTRAST_FIELDS.intersection(frame.columns))
if leaked:
    raise SystemExit(f"published parquet contains forbidden contrast fields: {leaked}")
if frame.empty:
    raise SystemExit("published parquet is empty")
if frame["exon_key"].duplicated().any():
    raise SystemExit("published parquet contains duplicate exon_key values")
if frame["gc_fraction"].isna().any():
    raise SystemExit("published parquet contains null gc_fraction values")
try:
    input_provenance = validate_input_provenance(
        frame,
        gtf,
        list(sj_paths.values()),
        fasta,
        list(count_paths.values()),
        read_length,
    )
except Exception as exc:
    raise SystemExit(str(exc)) from exc
required_annotation_only = {
    "chr18:63793160-63793277:+",  # SERPINB7_e4
    "chr18:70201893-70201983:-",  # RTTN_e4
}
detectable_keys = set(frame.loc[frame["detectable"], "exon_key"])
missing_annotation_only = sorted(required_annotation_only.difference(detectable_keys))
if missing_annotation_only:
    raise SystemExit(
        "published parquet is still gated by rMATS SE enumeration; missing detectable "
        f"annotated exons: {missing_annotation_only}"
    )

text = readme.read_text(encoding="utf-8")
if "101724" in text:
    raise SystemExit("published README contains stale pre-collapse count 101724")
if "three scrambled MOE (ScrMOE) control alignments only" not in text:
    raise SystemExit("published README is not bound to the ScrMOE-only universe")

observed_runs = []
for group_file in (b1, b2):
    bam_paths = [item for item in group_file.read_text().strip().split(",") if item]
    if not bam_paths:
        raise SystemExit(f"empty rMATS control group: {group_file}")
    for bam_path in bam_paths:
        match = re.fullmatch(r"(SRR[0-9]+)\.Aligned\.sortedByCoord\.out\.bam", Path(bam_path).name)
        if match is None:
            raise SystemExit(f"unrecognized BAM path in {group_file}: {bam_path}")
        observed_runs.append(match.group(1))
if len(observed_runs) != len(set(observed_runs)):
    raise SystemExit("duplicate BAM accession across rMATS b1/b2")
if set(observed_runs) != expected_runs:
    raise SystemExit(
        f"rMATS b1/b2 accessions differ from pinned ScrMOE allowlist: "
        f"observed={sorted(observed_runs)} expected={sorted(expected_runs)}"
    )
# Tie the README headline counts to the parquet reality so a mislabelled or
# stale README (e.g. unique/detectable transposed) cannot be published.
n_rows = len(frame)
n_detectable = int(frame["detectable"].sum())
if f"after alternate-flank collapse: {n_rows}" not in text:
    raise SystemExit(f"README unique-interval count does not match parquet ({n_rows})")
if f"Detectable background exons: {n_detectable}" not in text:
    raise SystemExit(f"README detectable count does not match parquet ({n_detectable})")
# Provenance is verified by content hash, not path string, so a committed
# bundle round-trips regardless of absolute vs repo-relative path rendering.
sj_digests = json.loads(str(input_provenance["scrmoe_sj_sha256"]))
count_digests = json.loads(str(input_provenance["scrmoe_reads_per_gene_sha256"]))
authoritative_digests = {
    gtf: str(input_provenance["gencode_gtf_sha256"]),
    fasta: str(input_provenance["genome_fasta_sha256"]),
    **{path: sj_digests[run] for run, path in sj_paths.items()},
    **{path: count_digests[run] for run, path in count_paths.items()},
}
for label, path in (
    ("SE source", se_table),
    ("GENCODE GTF", gtf),
    ("Reference FASTA", fasta),
    *[(f"ScrMOE STAR junction {run}", path) for run, path in sj_paths.items()],
    *[(f"ScrMOE STAR gene count {run}", path) for run, path in count_paths.items()],
    ("Conditions table", conditions),
    ("rMATS b1 control list", b1),
    ("rMATS b2 control list", b2),
):
    digest = (
        authoritative_digests[path]
        if path in authoritative_digests
        else sha256(path)
    )
    if f"{label} SHA-256: `{digest}`" not in text:
        raise SystemExit(f"{label} provenance SHA-256 mismatch in published README")
if f"Junction read length: `{read_length}`" not in text:
    raise SystemExit("junction read length provenance mismatch in published README")
PY
}

DONE_FILE="${CURRENT_POINTER}/.done"
LEGACY_DONE_FILE="${OUT_DIR}/.done"
if [[ -f "${DONE_FILE}" ]]; then
    if validate_published_artifact; then
        echo "SKIP: validated detectable universe already complete (${DONE_FILE})"
        exit 0
    fi
    echo "STALE: ${DONE_FILE} exists but the published bundle failed validation; rebuilding" >&2
    rm -f "${DONE_FILE}"
fi

# --- build the two control-only rMATS group files (arbitrary baseline split) ---
WORK="${OUT_DIR}/.work-${SLURM_JOB_ID:-manual}"
rm -rf "${WORK}"
mkdir -p "${WORK}/output" "${WORK}/tmp"
BUILD_ID="${SLURM_JOB_ID:-manual}-$$"
STAGED_BUNDLE="${BUNDLE_ROOT}/.${BUILD_ID}.tmp"
VERSIONED_BUNDLE="${BUNDLE_ROOT}/${BUILD_ID}"
STAGED_PARQUET="${STAGED_BUNDLE}/detectable_universe_v50.parquet"
STAGED_README="${STAGED_BUNDLE}/detectable_universe_v50_README.md"
cleanup_staged() {
    rm -rf -- "${STAGED_BUNDLE}"
}
trap cleanup_staged EXIT
mkdir -p "${STAGED_BUNDLE}"

group1=()
group2=()
idx=0
for run in "${CONTROL_RUNS[@]}"; do
    sentinel="${STAR_DIR}/${run}/.done"
    bam="${STAR_DIR}/${run}/${run}.Aligned.sortedByCoord.out.bam"
    [[ -f "${sentinel}" ]] || { echo "ERROR: alignment sentinel missing for ${run}" >&2; exit 1; }
    [[ -s "${bam}" ]] || { echo "ERROR: BAM missing for ${run}: ${bam}" >&2; exit 1; }
    if (( idx % 2 == 0 )); then group1+=("${bam}"); else group2+=("${bam}"); fi
    idx=$((idx + 1))
done
# rMATS requires both groups non-empty; a single control lands in group2.
if [[ "${#group2[@]}" -eq 0 ]]; then
    group2+=("${group1[-1]}")
    unset 'group1[${#group1[@]}-1]'
fi
(IFS=,; printf '%s\n' "${group1[*]}") > "${WORK}/b1.txt"
(IFS=,; printf '%s\n' "${group2[*]}") > "${WORK}/b2.txt"

echo "rMATS control-only enumeration: b1=${#group1[@]} b2=${#group2[@]} (both baseline); --statoff"
"${PYTHON_BIN}" "${RMATS}" \
    --b1 "${WORK}/b1.txt" \
    --b2 "${WORK}/b2.txt" \
    --gtf "${GTF}" \
    --od "${WORK}/output" \
    --tmp "${WORK}/tmp" \
    --task both \
    -t paired \
    --readLength "${READ_LENGTH}" \
    --libType "${LIB_TYPE}" \
    --statoff \
    --nthread "${THREADS}"

SE_TABLE="${WORK}/output/SE.MATS.JC.txt"
[[ -s "${SE_TABLE}" ]] || { echo "ERROR: rMATS SE output missing: ${SE_TABLE}" >&2; exit 1; }

# Preserve build inputs: the SE table on scratch, the small control lists as
# committed repo provenance (conditions is already committed at COMMIT_CONDITIONS).
cp "${SE_TABLE}" "${STABLE_SE_TABLE}"
cp "${WORK}/b1.txt" "${COMMIT_B1}"
cp "${WORK}/b2.txt" "${COMMIT_B2}"

# --- derive the detectable universe (parse only; no treatment data touched) ---
PYTHONPATH="${PKG_PATH}${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON_BIN}" -m feng.bench.strata \
    --conditions "${COMMIT_CONDITIONS}" \
    --se-table "${STABLE_SE_TABLE}" \
    --gtf "${GTF}" \
    --star-dir "${STAR_DIR}" \
    --genome-fasta "${GENOME_FASTA}" \
    --samtools "${SAMTOOLS}" \
    --read-length "${READ_LENGTH}" \
    --b1 "${COMMIT_B1}" \
    --b2 "${COMMIT_B2}" \
    --out-parquet "${STAGED_PARQUET}" \
    --out-readme "${STAGED_README}"

[[ -s "${STAGED_PARQUET}" ]] || { echo "ERROR: parquet not staged: ${STAGED_PARQUET}" >&2; exit 1; }
PYTHONPATH="${PKG_PATH}${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON_BIN}" - \
    "${STAGED_BUNDLE}" <<'PY'
from pathlib import Path
import sys

from feng.bench.strata import write_bundle_manifest

write_bundle_manifest(Path(sys.argv[1]))
PY
validate_published_artifact "${STAGED_PARQUET}" "${STAGED_README}" || {
    echo "ERROR: staged universe failed validation" >&2
    exit 1
}

PYTHONPATH="${PKG_PATH}${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON_BIN}" - \
    "${STAGED_BUNDLE}" "${VERSIONED_BUNDLE}" "${CURRENT_POINTER}" \
    "${LEGACY_DONE_FILE}" <<'PY'
from pathlib import Path
import sys

from feng.bench.strata import publish_staged_bundle

publish_staged_bundle(*(Path(value) for value in sys.argv[1:]))
PY
rm -rf "${WORK}"
echo "DONE: detectable universe published -> ${OUT_PARQUET}"
