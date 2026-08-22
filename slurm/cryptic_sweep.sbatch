#!/bin/bash
#SBATCH --job-name=cryptic-sweep
#SBATCH --account=mit_general
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
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
# 17,820 footprints x 4 targets-worth already folded into that count x 5-span
# ladder = 89,100 scored rows. Sized from a direct timing measurement: ~35 s
# per footprint (5 spans x baseline+blocked constrained partition functions),
# so 800 shards gives ~22 footprints/shard, ~13 min/shard -- comfortably
# under the preemption-horizon guidance in docs/HANDOFF.md (<15 min/unit).
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
