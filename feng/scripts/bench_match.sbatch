#!/usr/bin/env bash
#SBATCH --job-name=feng-bench-match
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --requeue
#SBATCH --output=feng/logs/%x-%j.out
#SBATCH --error=feng/logs/%x-%j.err

# Arm D control-only matched-strata stage. The matching CLI writes a staged
# pair beside the final artifacts. An independent validator checks the full
# frozen 8 x (1 positive + 100 backgrounds) contract before the existing
# atomic bundle publisher replaces either final file or creates .done.

set -euo pipefail

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
ENV_PREFIX="${FENG_ENV_PREFIX:-/home/dashii/orcd/pool/feng-env}"
PYTHON_BIN="${FENG_BENCH_PYTHON:-${ENV_PREFIX}/bin/python}"
PKG_PATH="${SUBMIT_DIR}/feng/src"

if [[ -z "${SCRATCH:-}" && ( -z "${FENG_BENCH_UNIVERSE_DIR:-}" || -z "${FENG_BENCH_MATCH_DIR:-}" ) ]]; then
    echo "ERROR: set SCRATCH, or explicitly set FENG_BENCH_UNIVERSE_DIR and FENG_BENCH_MATCH_DIR" >&2
    exit 1
fi

ARTIFACT_DIR="${FENG_BENCH_ARTIFACT_DIR:-${SUBMIT_DIR}/feng/artifacts/bench}"
POSITIVES="${FENG_BENCH_POSITIVES:-${ARTIFACT_DIR}/positives_v50.tsv}"
UNIVERSE_DIR="${FENG_BENCH_UNIVERSE_DIR:-${SCRATCH}/feng/bench/detectable_universe}"
UNIVERSE_POINTER="${FENG_BENCH_UNIVERSE_POINTER:-${ARTIFACT_DIR}/detectable_universe_v50.current}"
UNIVERSE_BUNDLES="${ARTIFACT_DIR}/detectable_universe_v50.bundles"
[[ -L "${UNIVERSE_POINTER}" ]] || {
    echo "ERROR: detectable universe current pointer missing: ${UNIVERSE_POINTER}" >&2
    exit 1
}
PINNED_UNIVERSE_BUNDLE="$(readlink -f -- "${UNIVERSE_POINTER}")" || {
    echo "ERROR: detectable universe current pointer is dangling: ${UNIVERSE_POINTER}" >&2
    exit 1
}
PINNED_BUNDLE_ROOT="$(readlink -f -- "${UNIVERSE_BUNDLES}")" || {
    echo "ERROR: detectable universe bundle root missing: ${UNIVERSE_BUNDLES}" >&2
    exit 1
}
case "${PINNED_UNIVERSE_BUNDLE}" in
    "${PINNED_BUNDLE_ROOT}"/*) ;;
    *)
        echo "ERROR: detectable universe current pointer escapes bundle root" >&2
        exit 1
        ;;
esac
UNIVERSE="${PINNED_UNIVERSE_BUNDLE}/detectable_universe_v50.parquet"
UNIVERSE_DONE="${PINNED_UNIVERSE_BUNDLE}/.done"
UNIVERSE_MANIFEST="${PINNED_UNIVERSE_BUNDLE}/manifest.json"
if [[ -n "${FENG_BENCH_UNIVERSE:-}" ]]; then
    OVERRIDE_UNIVERSE="$(readlink -f -- "${FENG_BENCH_UNIVERSE}")" || true
    if [[ "${OVERRIDE_UNIVERSE}" != "${UNIVERSE}" ]]; then
        echo "ERROR: FENG_BENCH_UNIVERSE may only name the universe pinned by current" >&2
        exit 1
    fi
fi
if [[ -n "${FENG_BENCH_UNIVERSE_DONE:-}" ]]; then
    OVERRIDE_UNIVERSE_DONE="$(readlink -f -- "${FENG_BENCH_UNIVERSE_DONE}")" || true
    if [[ "${OVERRIDE_UNIVERSE_DONE}" != "${UNIVERSE_DONE}" ]]; then
        echo "ERROR: FENG_BENCH_UNIVERSE_DONE may only name the sentinel in the pinned bundle" >&2
        exit 1
    fi
fi
MATCH_DIR="${FENG_BENCH_MATCH_DIR:-${SCRATCH}/feng/bench/matched_strata}"
OUT_STRATA="${FENG_BENCH_STRATA:-${ARTIFACT_DIR}/strata_v50.parquet}"
OUT_SUMMARY="${FENG_BENCH_SUMMARY:-${ARTIFACT_DIR}/strata_summary_v50.tsv}"
DONE_FILE="${MATCH_DIR}/.done"
MATCH_SOURCE="${FENG_BENCH_MATCH_SOURCE:-${PKG_PATH}/feng/bench/match.py}"

# Frozen matching settings. Publication under alternate thresholds is not
# exposed through environment overrides.
N_BACKGROUND=100
CALIPER_SD=1.0
MIN_BACKGROUND=50
SEED=20260901

[[ -x "${PYTHON_BIN}" ]] || { echo "ERROR: Python interpreter missing: ${PYTHON_BIN}" >&2; exit 1; }
[[ -s "${UNIVERSE}" ]] || { echo "ERROR: detectable universe missing: ${UNIVERSE}" >&2; exit 1; }
[[ -s "${POSITIVES}" ]] || { echo "ERROR: positives table missing: ${POSITIVES}" >&2; exit 1; }
[[ -s "${MATCH_SOURCE}" ]] || { echo "ERROR: matcher source missing: ${MATCH_SOURCE}" >&2; exit 1; }
[[ -f "${UNIVERSE_DONE}" ]] || {
    echo "ERROR: detectable universe completion sentinel missing: ${UNIVERSE_DONE}" >&2
    exit 1
}
[[ -s "${UNIVERSE_MANIFEST}" ]] || {
    echo "ERROR: detectable universe bundle manifest missing: ${UNIVERSE_MANIFEST}" >&2
    exit 1
}

PYTHONPATH="${PKG_PATH}${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON_BIN}" - \
    "${PINNED_UNIVERSE_BUNDLE}" "${UNIVERSE_DONE}" <<'PY'
from pathlib import Path
import hashlib
import sys

from feng.bench.strata import validate_bundle_manifest

bundle, done = map(Path, sys.argv[1:3])
try:
    validate_bundle_manifest(bundle)
    expected = hashlib.sha256((bundle / "manifest.json").read_bytes()).hexdigest()
    if done.read_text(encoding="ascii").strip() != expected:
        raise ValueError("completion sentinel does not bind the bundle manifest")
except (OSError, ValueError) as exc:
    raise SystemExit(f"invalid pinned detectable-universe bundle: {exc}") from exc
PY

mkdir -p "${MATCH_DIR}" "$(dirname "${OUT_STRATA}")" "$(dirname "${OUT_SUMMARY}")" "${SUBMIT_DIR}/feng/logs"

validate_matched_outputs() {
    local strata_path="${1:-${OUT_STRATA}}"
    local summary_path="${2:-${OUT_SUMMARY}}"
    PYTHONPATH="${PKG_PATH}${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON_BIN}" - \
        "${strata_path}" "${summary_path}" "${POSITIVES}" "${UNIVERSE}" \
        "${N_BACKGROUND}" "${MIN_BACKGROUND}" <<'PY'
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from pandas.api.types import is_bool_dtype

from feng.bench.match import COORD_COLUMNS, COVARIATES
from feng.bench.strata import (
    INPUT_PROVENANCE_COLUMNS,
    UNIVERSE_PROVENANCE,
    UNIVERSE_PROVENANCE_COLUMNS,
)

strata_path, summary_path, positives_path, universe_path = map(Path, sys.argv[1:5])
n_background = int(sys.argv[5])
min_background = int(sys.argv[6])

missing_files = [
    str(path) for path in (strata_path, summary_path, positives_path, universe_path)
    if not path.is_file() or path.stat().st_size == 0
]
if missing_files:
    raise SystemExit("missing/empty matched-strata input: " + ", ".join(missing_files))

try:
    strata = pd.read_parquet(strata_path)
    summary = pd.read_csv(summary_path, sep="\t")
    positives = pd.read_csv(positives_path, sep="\t")
    universe = pd.read_parquet(universe_path)
except Exception as exc:
    raise SystemExit(f"could not read matched-strata bundle: {exc}") from exc

canonical_positive_ids = {
    "OT01_POLR2H_e2",
    "OT02_PITHD1_e2",
    "OT03_SERPINB7_e4",
    "OT04_RTTN_e4",
    "OT05_REV3L_e9",
    "OT06_PRKRA_e2",
    "OT07_GOLGA4_e4",
    "OT08_PAK1_e2",
}
if "event_id" not in positives.columns:
    raise SystemExit("positives table lacks event_id")
positive_ids = positives["event_id"].astype(str)
if len(positive_ids) != 8 or positive_ids.duplicated().any() or set(positive_ids) != canonical_positive_ids:
    raise SystemExit("positives table is not the canonical eight-event Arm D panel")

required_summary = {
    "event_id", "positive_exon_key", "in_universe", "detectable",
    "eligible_background", "n_selected", "kill_gate_pass",
}
missing_summary = sorted(required_summary.difference(summary.columns))
if missing_summary:
    raise SystemExit(f"strata summary missing columns: {missing_summary}")
summary_ids = summary["event_id"].astype(str)
if len(summary) != 8 or summary_ids.duplicated().any() or set(summary_ids) != canonical_positive_ids:
    raise SystemExit("each canonical positive must occur exactly once in the strata summary")
for column in ("in_universe", "detectable", "kill_gate_pass"):
    if not summary[column].eq(True).all():
        raise SystemExit(f"summary gate is not uniformly true: {column}")
eligible = pd.to_numeric(summary["eligible_background"], errors="coerce")
selected = pd.to_numeric(summary["n_selected"], errors="coerce")
if (
    not np.isfinite(eligible).all()
    or (eligible < 0).any()
    or not np.equal(eligible, np.floor(eligible)).all()
):
    raise SystemExit("eligible background counts must be nonnegative integers")
if (
    not np.isfinite(selected).all()
    or (selected < 0).any()
    or not np.equal(selected, np.floor(selected)).all()
):
    raise SystemExit("selected background counts must be nonnegative integers")
if not (eligible >= min_background).all() or not (eligible >= selected).all():
    raise SystemExit(f"every positive must have at least {min_background} eligible backgrounds")
if not np.isfinite(selected).all() or not selected.eq(n_background).all():
    raise SystemExit(f"every positive summary must select exactly {n_background} backgrounds")
if summary["positive_exon_key"].isna().any() or summary["positive_exon_key"].duplicated().any():
    raise SystemExit("summary positive exon keys must be present and unique")

required_strata = {
    "positive_event_id", "positive_exon_key", "is_positive", "match_distance",
    "exon_key", *COORD_COLUMNS, *COVARIATES, *UNIVERSE_PROVENANCE_COLUMNS,
    *INPUT_PROVENANCE_COLUMNS,
}
missing_strata = sorted(required_strata.difference(strata.columns))
if missing_strata:
    raise SystemExit(f"strata parquet missing columns: {missing_strata}")
if (
    strata.empty
    or not is_bool_dtype(strata["is_positive"])
    or not strata["is_positive"].notna().all()
):
    raise SystemExit("strata is_positive must be a non-empty boolean column")
if len(strata) != len(canonical_positive_ids) * (n_background + 1):
    raise SystemExit("strata row count does not equal the frozen 8 x (1 + N) contract")
if strata["positive_event_id"].isna().any() or set(strata["positive_event_id"].astype(str)) != canonical_positive_ids:
    raise SystemExit("strata does not contain exactly the canonical positive IDs")

numeric_columns = ["match_distance", "exon_start_0base", "exon_end", *COVARIATES]
for column in numeric_columns:
    values = pd.to_numeric(strata[column], errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise SystemExit(f"strata column contains non-finite values: {column}")
if (pd.to_numeric(strata["match_distance"]) < 0).any():
    raise SystemExit("match distances must be non-negative")

summary_keys = summary.set_index(summary_ids)["positive_exon_key"].astype(str)
positive_exon_keys = set(summary_keys)
for event_id, group in strata.groupby(strata["positive_event_id"].astype(str), sort=False):
    if len(group) != n_background + 1:
        raise SystemExit(f"{event_id} row count does not equal one positive plus {n_background} backgrounds")
    anchors = group.loc[group["is_positive"]]
    backgrounds = group.loc[~group["is_positive"]]
    if len(anchors) != 1 or len(backgrounds) != n_background:
        raise SystemExit(
            f"{event_id} must contain one positive and {n_background} backgrounds"
        )
    expected_key = summary_keys.loc[event_id]
    if group["positive_exon_key"].astype(str).nunique() != 1 or str(group["positive_exon_key"].iloc[0]) != expected_key:
        raise SystemExit(f"{event_id} carries an inconsistent positive_exon_key")
    if str(anchors["exon_key"].iloc[0]) != expected_key or float(anchors["match_distance"].iloc[0]) != 0.0:
        raise SystemExit(f"{event_id} positive anchor identity or distance is invalid")
    if backgrounds["exon_key"].astype(str).duplicated().any():
        raise SystemExit(f"{event_id} contains duplicate background exons")

background_keys = set(strata.loc[~strata["is_positive"], "exon_key"].astype(str))
if positive_exon_keys.intersection(background_keys):
    raise SystemExit("a positive exon leaked into the matched backgrounds")

for column in UNIVERSE_PROVENANCE_COLUMNS:
    values = strata[column]
    if values.isna().any() or values.nunique(dropna=False) != 1:
        raise SystemExit(f"universe provenance is missing or internally inconsistent: {column}")
    if values.iloc[0] != UNIVERSE_PROVENANCE[column]:
        raise SystemExit(f"universe provenance does not match the frozen policy: {column}")

required_universe = {"exon_key", *COORD_COLUMNS, *COVARIATES, *INPUT_PROVENANCE_COLUMNS}
missing_universe = sorted(required_universe.difference(universe.columns))
if missing_universe:
    raise SystemExit(f"current universe lacks provenance/source columns: {missing_universe}")
if universe["exon_key"].isna().any() or universe["exon_key"].duplicated().any():
    raise SystemExit("current universe has null or duplicate exon identities")
for column in INPUT_PROVENANCE_COLUMNS:
    source_values = universe[column]
    output_values = strata[column]
    if (
        source_values.isna().any()
        or source_values.nunique(dropna=False) != 1
        or output_values.isna().any()
        or output_values.nunique(dropna=False) != 1
        or output_values.iloc[0] != source_values.iloc[0]
    ):
        raise SystemExit(f"strata input provenance does not match current universe: {column}")

source_fields = [*COORD_COLUMNS, *COVARIATES]
source = universe.set_index(universe["exon_key"].astype(str))[source_fields]
observed = source.reindex(strata["exon_key"].astype(str))
if observed.isna().any().any():
    raise SystemExit("strata contains an exon absent from the current universe")
for column in source_fields:
    left = strata[column]
    right = observed[column]
    if column in COVARIATES or column in {"exon_start_0base", "exon_end"}:
        if not np.array_equal(
            pd.to_numeric(left, errors="coerce").to_numpy(dtype=float),
            pd.to_numeric(right, errors="coerce").to_numpy(dtype=float),
        ):
            raise SystemExit(f"strata {column} does not match the current universe")
    elif not np.array_equal(
        left.astype(str).to_numpy(), right.astype(str).to_numpy()
    ):
        raise SystemExit(f"strata {column} does not match the current universe")
PY
}

manifest_for() {
    local strata_path="$1"
    local summary_path="$2"
    PYTHONPATH="${PKG_PATH}${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON_BIN}" - \
        "${UNIVERSE}" "${UNIVERSE_MANIFEST}" "${POSITIVES}" "${MATCH_SOURCE}" "${strata_path}" "${summary_path}" \
        "${N_BACKGROUND}" "${CALIPER_SD}" "${MIN_BACKGROUND}" "${SEED}" <<'PY'
import hashlib
import json
from pathlib import Path
import sys

universe, universe_manifest, positives, matcher, strata, summary = map(Path, sys.argv[1:7])
n_background, caliper_sd, min_background, seed = sys.argv[7:11]

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()

print(json.dumps({
    "schema": "feng-bench-match-v1",
    "inputs": {
        "universe_sha256": sha256(universe),
        "universe_bundle_manifest_sha256": sha256(universe_manifest),
        "positives_sha256": sha256(positives),
        "matcher_source_sha256": sha256(matcher),
    },
    "settings": {
        "n_background": int(n_background),
        "caliper_sd": float(caliper_sd),
        "min_background": int(min_background),
        "seed": int(seed),
    },
    "outputs": {
        "strata_sha256": sha256(strata),
        "summary_sha256": sha256(summary),
    },
}, sort_keys=True, separators=(",", ":")))
PY
}

done_matches_current_bundle() {
    local expected
    expected="$(manifest_for "${OUT_STRATA}" "${OUT_SUMMARY}")"
    [[ "$(<"${DONE_FILE}")" == "${expected}" ]]
}

if [[ -f "${DONE_FILE}" ]]; then
    if validate_matched_outputs && done_matches_current_bundle; then
        echo "SKIP: validated matched strata already complete (${DONE_FILE})"
        exit 0
    fi
    echo "STALE: ${DONE_FILE} exists but the published pair failed validation; rebuilding" >&2
    rm -f "${DONE_FILE}"
fi

BUILD_ID="${SLURM_JOB_ID:-manual}-$$"
STRATA_PARENT="$(dirname "${OUT_STRATA}")"
SUMMARY_PARENT="$(dirname "${OUT_SUMMARY}")"
STAGED_STRATA="${STRATA_PARENT}/.$(basename "${OUT_STRATA}").${BUILD_ID}.tmp"
STAGED_SUMMARY="${SUMMARY_PARENT}/.$(basename "${OUT_SUMMARY}").${BUILD_ID}.tmp"
cleanup_staged() {
    rm -f "${STAGED_STRATA}" "${STAGED_SUMMARY}"
}
trap cleanup_staged EXIT

PYTHONPATH="${PKG_PATH}${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON_BIN}" -m feng.bench.match \
    --universe "${UNIVERSE}" \
    --positives "${POSITIVES}" \
    --out-strata "${STAGED_STRATA}" \
    --out-summary "${STAGED_SUMMARY}" \
    --n-background "${N_BACKGROUND}" \
    --caliper-sd "${CALIPER_SD}" \
    --min-background "${MIN_BACKGROUND}" \
    --seed "${SEED}"

if ! validate_matched_outputs "${STAGED_STRATA}" "${STAGED_SUMMARY}"; then
    echo "ERROR: staged matched strata failed validation" >&2
    exit 1
fi

DONE_CONTENTS="$(manifest_for "${STAGED_STRATA}" "${STAGED_SUMMARY}")"

PYTHONPATH="${PKG_PATH}${PYTHONPATH:+:${PYTHONPATH}}" "${PYTHON_BIN}" - \
    "${STAGED_STRATA}" "${STAGED_SUMMARY}" "${OUT_STRATA}" \
    "${OUT_SUMMARY}" "${DONE_FILE}" "${DONE_CONTENTS}" <<'PY'
from pathlib import Path
import sys

from feng.bench.strata import publish_staged_pair

publish_staged_pair(
    *(Path(value) for value in sys.argv[1:6]), done_contents=sys.argv[6]
)
PY

echo "DONE: matched strata published -> ${OUT_STRATA}, ${OUT_SUMMARY}"
