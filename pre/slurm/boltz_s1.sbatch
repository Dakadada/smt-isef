#!/bin/bash
#SBATCH --job-name=boltz-s1
#SBATCH --account=mit_general
#SBATCH --partition=mit_normal_gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --array=0-10
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/boltz-s1-%A_%a.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/boltz-s1-%A_%a.err

# S1: 44 inputs (4 unbound controls + 40 prereg oligos) x 5 diffusion samples.
#
# --diffusion_samples 5 rather than 5 separate --seed runs: featurization is
# shared, so five structures cost far less than five predictions. Seed spread
# is what makes the within-pair comparison statistical instead of anecdotal.
#
# Sizing measured, not guessed: one 356 nt input took 5m45s (job 21071541,
# COMPLETED, no OOM), which is what closed the open window-size risk. Four
# inputs x 5 samples per task against a 3 h wall.
#
# Boltz consumes a DIRECTORY, so inputs are pre-split into part00..part10.
# Skip-if-done is on the part's output dir, so a requeued task costs at most
# that part.
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
module load cuda/13.0.1
conda activate /home/dashii/orcd/pool/boltz-env

export BOLTZ_CACHE=/home/dashii/orcd/pool/boltz-cache
export HF_HOME=/home/dashii/orcd/pool/hf-cache
S=/home/dashii/orcd/pool/smt-isef-data/boltz_s
PART=$(printf "part%02d" "${SLURM_ARRAY_TASK_ID}")
OUT="$S/s1_out/$PART"

if [ -f "$OUT/.complete" ]; then echo "$PART already done"; exit 0; fi
mkdir -p "$OUT"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

boltz predict "$S/s1_split/$PART" \
    --out_dir "$OUT" \
    --cache "$BOLTZ_CACHE" \
    --accelerator gpu --devices 1 \
    --diffusion_samples 5 \
    --output_format mmcif \
    --use_msa_server --override
touch "$OUT/.complete"
echo "S1 $PART complete"
