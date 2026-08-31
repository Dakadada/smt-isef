#!/usr/bin/env bash
#SBATCH --job-name=smt-stageA-rescue
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=36:00:00
#SBATCH --partition=mit_data_transfer
#SBATCH --requeue
#SBATCH --signal=USR1@300
#SBATCH --output=stageA-rescue-%j.out
set -euo pipefail

REPO=${REPO:-/home/dashii/smt-isef}
RUN_DIR=${RUN_DIR:?set RUN_DIR}
RESCUE_MANIFEST=${RESCUE_MANIFEST:?set RESCUE_MANIFEST to the scoped frozen-manifest subset}
CONFIG="$RUN_DIR/stageA.targeted.json"
FULL_MANIFEST="$RUN_DIR/all_runs.frozen.tsv"
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

requeue_self() { touch "$RUN_DIR/RESCUE_REQUEUE_REQUESTED.${SLURM_JOB_ID}"; scontrol requeue "${SLURM_JOB_ID}"; exit 0; }
trap requeue_self USR1

/usr/bin/python3.12 -m rawread_reanalysis process --config "$CONFIG" --manifest "$RESCUE_MANIFEST"
/usr/bin/python3.12 - "$CONFIG" "$FULL_MANIFEST" <<'PY'
import sys
from rawread_reanalysis.common import validate_config
from rawread_reanalysis.quantify import rebuild_combined
cfg = validate_config(sys.argv[1])
rebuild_combined(cfg, sys.argv[2])
PY
/usr/bin/python3.12 -m rawread_reanalysis compute --config "$CONFIG"
/usr/bin/python3.12 -m rawread_reanalysis harmonize --config "$CONFIG"
/usr/bin/python3.12 -m rawread_reanalysis mark-completion --config "$CONFIG" --run-dir "$RUN_DIR"
