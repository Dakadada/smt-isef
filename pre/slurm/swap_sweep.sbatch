#!/bin/bash
#SBATCH --job-name=swap-sweep
#SBATCH --account=mit_general
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --array=0-399%60
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/swap-%A_%a.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/swap-%A_%a.err

# N4: does "occlusion drives the 5'ss partner downstream" predict the sign of
# delta_A? Pre-registered in docs/preregistration.md; the 90% threshold comes
# from N3's observed separation and is NOT re-fit on this sweep.
#
# 1065 footprints x 5 spans = 5325 partition functions, ~3 PF/footprint-span at
# 798 nt. Smoke shard: 3 footprints in well under a minute.
# Skip-if-done + atomic .tmp rename, so requeue costs at most one shard.
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
conda activate /home/dashii/orcd/pool/smt-isef-env
cd /home/dashii/smt-isef
export PYTHONPATH=/home/dashii/smt-isef/src
python -m smt.validate.swap_sweep --shard "$SLURM_ARRAY_TASK_ID" --nshards 400
