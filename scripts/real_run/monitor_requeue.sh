#!/usr/bin/env bash
set -euo pipefail
job_id=${1:?usage: monitor_requeue.sh JOB_ID [RUN_DIR]}
run_dir=${2:-.}
while squeue -h -j "$job_id" -o '%T' | grep -q .; do
  squeue -h -j "$job_id" -o '%i %T %M %R'
  test -f "$run_dir/STAGE_A_COMPLETE" && exit 0
  sleep 300
done
test -f "$run_dir/STAGE_A_COMPLETE" || { echo "job left queue without completion marker; inspect sacct -j $job_id and re-submit the same command" >&2; exit 1; }
