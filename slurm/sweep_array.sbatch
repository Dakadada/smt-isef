#!/bin/bash
#SBATCH --job-name=aso-sweep
#SBATCH --account=mit_general
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --array=0-399
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/sweep-%A_%a.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/sweep-%A_%a.err

# A3 falsification test: tile ASOs of length 6-20 across and around ISS-N1,
# score each by delta A_5ss, and check the model independently recovers the
# published 8-mer GC-rich-core result without being told the answer.

set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
conda activate /home/dashii/orcd/pool/smt-isef-env

cd /home/dashii/smt-isef
export PYTHONPATH=/home/dashii/smt-isef/src

python -m smt.validate.length_sweep \
    --shard "${SLURM_ARRAY_TASK_ID}" \
    --nshards 400 \
    --out /home/dashii/orcd/pool/smt-isef-data/sweep
