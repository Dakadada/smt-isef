#!/usr/bin/env bash
#SBATCH --job-name=rawread-preempt
#SBATCH --partition=PREEMPTABLE
#SBATCH --cpus-per-task=24
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:0
#SBATCH --signal=B:USR1@120
#SBATCH --requeue
#SBATCH --output=slurm/logs/%x-%j.out
set -euo pipefail
: "${CONFIG:=config/pipeline.yaml}"
: "${STAGE:=all}"
mkdir -p slurm/logs data/reanalysis/heartbeats
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
heartbeat="data/reanalysis/heartbeats/${SLURM_JOB_ID}.json"
write_heartbeat() {
  tmp="${heartbeat}.tmp.$$"
  printf '{"job_id":"%s","time":%s,"stage":"%s"}\n' "$SLURM_JOB_ID" "$(date +%s)" "$STAGE" > "$tmp"
  mv "$tmp" "$heartbeat"
}
on_preempt() { write_heartbeat; exit 99; }
trap on_preempt USR1 TERM
(while true; do write_heartbeat; sleep 600; done) &
heartbeat_pid=$!
trap 'kill "$heartbeat_pid" 2>/dev/null || true' EXIT
python3 -m rawread_reanalysis "$STAGE" --config "$CONFIG"
