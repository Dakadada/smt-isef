#!/bin/bash
#SBATCH --job-name=alu-identity2
#SBATCH --account=mit_general
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/alu2-%j.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/alu2-%j.err
# P, re-run orientation-aware. The flank scan is ~990 RNAduplex calls (~3 min)
# and the null is now 2x the work, since every draw is scored in both
# orientations to match how the observation is scored.
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
conda activate /home/dashii/orcd/pool/smt-isef-env
cd /home/dashii/smt-isef
export PYTHONPATH=/home/dashii/smt-isef/src
python -m smt.kinetics.alu_identity --top 3
