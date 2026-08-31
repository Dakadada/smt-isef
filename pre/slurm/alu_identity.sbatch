#!/bin/bash
#SBATCH --job-name=alu-identity
#SBATCH --account=mit_general
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/alu-identity-%j.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/alu-identity-%j.err

# P: are the C2B34 flanking inverted repeats actually Alu elements?
# The flank scan is ~3 min (990 RNAduplex calls) plus null trials, which is
# over the interactive limit and not a considerate login-node job anyway.
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
conda activate /home/dashii/orcd/pool/smt-isef-env
cd /home/dashii/smt-isef
export PYTHONPATH=/home/dashii/smt-isef/src
python -m smt.kinetics.alu_identity --top 3 > docs/p-alu-identity.txt 2>&1
cat docs/p-alu-identity.txt
