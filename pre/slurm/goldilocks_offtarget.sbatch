#!/bin/bash
#SBATCH --job-name=g1-offtarget
#SBATCH --account=mit_general
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --array=0-10
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/g1-%A_%a.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/g1-%A_%a.err
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
conda activate /home/dashii/orcd/pool/smt-isef-env
cd /home/dashii/smt-isef
export PYTHONPATH=/home/dashii/smt-isef/src
# One length per task. Each task skips if its output file already exists and
# publishes via .tmp -> os.replace, so a requeue costs at most one length.
python -m smt.validate.goldilocks_offtarget \
    --task "$SLURM_ARRAY_TASK_ID" \
    --out /home/dashii/orcd/pool/smt-isef-data/goldilocks_offtarget
