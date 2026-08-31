#!/bin/bash
#SBATCH --job-name=minigene
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --array=0-2
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --signal=B:USR1@120
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/minigene-%A_%a.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/minigene-%A_%a.err

# A6: run the minigene compression algorithm (smt.minigene.report) for one
# gene per array task. NOT SUBMITTED as part of this build -- the direct
# interactive run already completes for all three genes within this
# session's budget (see docs/HANDOFF.md for the actual numbers). Kept here,
# per project convention (see slurm/sweep_array.sbatch, slurm/cryptic_sweep.
# sbatch), for reproducibility and for a future re-run at a wider threshold
# sweep or additional genes where the interactive path may no longer be
# tractable.
#
# Idempotent: smt.minigene.report.run_one_gene_idempotent skips a gene whose
# <out>/<GENE>.done marker already exists, and writes that marker via
# .done.tmp + atomic rename, so --requeue on preemption costs at most the
# one in-flight gene, never a partial/corrupt marker.

set -euo pipefail
trap 'echo "preemption signal received, exiting for clean requeue"; exit 0' USR1

source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
conda activate /home/dashii/orcd/pool/smt-isef-env

cd /home/dashii/smt-isef
export PYTHONPATH=/home/dashii/smt-isef/src

GENES=(SMN2 STMN2 UNC13A)
GENE="${GENES[$SLURM_ARRAY_TASK_ID]}"

python -m smt.minigene.report \
    --gene "${GENE}" \
    --out /home/dashii/orcd/pool/smt-isef-data/minigene
