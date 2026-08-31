#!/usr/bin/env bash
#SBATCH --job-name=smt-stageA-real
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=36:00:00
#SBATCH --partition=mit_data_transfer
#SBATCH --requeue
#SBATCH --signal=USR1@300
#SBATCH --output=stageA-%j.out
set -euo pipefail

REPO=${REPO:-/home/dashii/smt-isef}
RUN_DIR=${RUN_DIR:?set RUN_DIR to the directory made by build_manifest.py}
CONFIG="$RUN_DIR/stageA.targeted.json"
MANIFEST="$RUN_DIR/all_runs.frozen.tsv"
EVENTS="$RUN_DIR/event_junctions.GRCh38.run.tsv"
export PYTHONPATH="$REPO/src${PYTHONPATH:+:$PYTHONPATH}"

requeue_self() { touch "$RUN_DIR/REQUEUE_REQUESTED.${SLURM_JOB_ID}"; scontrol requeue "${SLURM_JOB_ID}"; exit 0; }
trap requeue_self USR1
(while true; do date -u +%FT%TZ > "$RUN_DIR/heartbeat.${SLURM_JOB_ID}"; sleep 600; done) &
HEARTBEAT_PID=$!
trap 'kill "$HEARTBEAT_PID" 2>/dev/null || true' EXIT

(cd "$RUN_DIR" && sha256sum -c "$(basename "$CONFIG.sha256")" && sha256sum -c "$(basename "$MANIFEST.sha256")" && sha256sum -c "$(basename "$EVENTS.sha256")")
/usr/bin/python3.12 -c 'import sys; from rawread_reanalysis.common import validate_config; c=validate_config(sys.argv[1]); assert c["alignment"]["tier"]=="targeted"' "$CONFIG"

while IFS= read -r family; do
  (cd "$RUN_DIR" && sha256sum -c "${family}.frozen.tsv.sha256")
  /usr/bin/python3.12 -m rawread_reanalysis process --config "$CONFIG" --manifest "$RUN_DIR/${family}.frozen.tsv"
done < <(/usr/bin/python3.12 -c 'import json,sys; print(*json.load(open(sys.argv[1]))["studies"], sep="\n")' "$CONFIG")
/usr/bin/python3.12 - "$CONFIG" "$MANIFEST" <<'PY'
import sys
from rawread_reanalysis.common import validate_config
from rawread_reanalysis.quantify import rebuild_combined
cfg=validate_config(sys.argv[1]); rebuild_combined(cfg,sys.argv[2])
PY
/usr/bin/python3.12 -m rawread_reanalysis compute --config "$CONFIG"
/usr/bin/python3.12 -m rawread_reanalysis harmonize --config "$CONFIG"
/usr/bin/python3.12 -m rawread_reanalysis mark-completion --config "$CONFIG" --run-dir "$RUN_DIR"
