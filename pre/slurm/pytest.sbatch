#!/bin/bash
#SBATCH --job-name=pytest
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/pytest-%j.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/pytest-%j.err
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
conda activate /home/dashii/orcd/pool/smt-isef-env
cd /home/dashii/smt-isef
export PYTHONPATH=/home/dashii/smt-isef/src
python -m pytest -q --durations=10
