#!/bin/bash
#SBATCH --job-name=boltz-bp
#SBATCH --account=mit_general
#SBATCH --partition=mit_normal_gpu
#SBATCH --gres=gpu:l40s:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=05:45:00
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/boltz-bp-%j.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/boltz-bp-%j.err

# A4 validation gate: Boltz-2 must reproduce known Cas13 ternary complexes
# before any novel dCas13x-SMN2 prediction is trusted. Inputs are generated
# from the deposited structures; output is RMSD against the deposited coords.

set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
module load cuda/13.0.1
conda activate /home/dashii/orcd/pool/boltz-env

export BOLTZ_CACHE=/home/dashii/orcd/pool/boltz-cache
export HF_HOME=/home/dashii/orcd/pool/hf-cache
mkdir -p "$BOLTZ_CACHE" "$HF_HOME"

IN_DIR=${1:?usage: sbatch boltz_backpredict.sbatch <input_dir> <out_dir>}
OUT_DIR=${2:?usage: sbatch boltz_backpredict.sbatch <input_dir> <out_dir>}
mkdir -p "$OUT_DIR"

nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

boltz predict "$IN_DIR" \
    --out_dir "$OUT_DIR" \
    --cache "$BOLTZ_CACHE" \
    --accelerator gpu \
    --devices 1 \
    --output_format mmcif \
    --override

echo "boltz back-prediction complete: $OUT_DIR"
