#!/bin/bash
#SBATCH --job-name=cryptic-sweep
#SBATCH --account=mit_general
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:45:00
#SBATCH --array=0-799
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/cryptic-sweep-%A_%a.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/cryptic-sweep-%A_%a.err

# A7: tile effector footprints around the STMN2 cryptic 3' splice site and
# the UNC13A cryptic donor + both alternative acceptors, score each with the
# generalized signed OCCLUDE readout (smt.thermo.occlusion) over the
# max_bp_span ladder, exactly mirroring slurm/sweep_array.sbatch's SMN2 /
# FREE-regime sweep (smt.validate.length_sweep) in structure.
#
# 17,820 footprints (4 targets x 4,455) x 5-span ladder = 89,100 scored rows.
# 800 shards gives 22-23 footprints/shard.
#
# TIMING, re-measured 2026-08-23 before submission. The original note here
# said ~35 s/footprint -> ~13 min/shard under a 20 min limit. Both halves of
# that were a problem:
#   * A smoke run of shard 0 did NOT finish in 10 min. Measured cost was
#     35-41 s/footprint, so ~14 min/shard against a 20 min wall -- only ~40%
#     headroom on a shared, variably-loaded node.
#   * --requeue does NOT save a shard that hits its wall clock. It covers
#     preemption; a TIMEOUT is a plain failure, leaves no output file, and the
#     skip-if-done check then has nothing to skip. That is the A3 540-failed-
#     task pattern (see docs/HANDOFF.md), which is why this is not left alone.
# Fixed the cost rather than only the limit: the unoccluded baseline is a
# function of (construct, span) alone, so it is now computed 20 times per
# shard instead of once per footprint (smt.validate.cryptic_sweep.
# build_baselines). Measured 40.7 s -> 20.4 s per footprint, a verified 1.99x,
# with rows bit-identical to the uncached path (tests/test_cryptic.py).
# Now ~8 min/shard (23 x 20.4 s + ~70 s of baselines) under a 45 min wall.
#
# --requeue + --open-mode=append: a preempted/killed shard is requeued and
# its log is appended to, not clobbered. Each shard skips if its output file
# already exists and writes via a .tmp -> rename, so a shard killed mid-write
# never leaves a half-file the skip-check would trust (see
# smt.validate.cryptic_sweep.run_shard).
#
# NOT SUBMITTED by the agent that wrote this file, per instructions -- left
# for the main session to review and submit.

set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
conda activate /home/dashii/orcd/pool/smt-isef-env

cd /home/dashii/smt-isef
export PYTHONPATH=/home/dashii/smt-isef/src

python -m smt.validate.cryptic_sweep \
    --shard "${SLURM_ARRAY_TASK_ID}" \
    --nshards 800 \
    --out /home/dashii/orcd/pool/smt-isef-data/cryptic_sweep
