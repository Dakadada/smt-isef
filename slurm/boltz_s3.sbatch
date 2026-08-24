#!/bin/bash
#SBATCH --job-name=boltz-s3
#SBATCH --account=mit_general
#SBATCH --partition=mit_normal_gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --array=0-1
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/boltz-s3-%A_%a.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/boltz-s3-%A_%a.err

# S3: hnRNP A1 steric parity. Task 0 = the GATE (deposited hnRNP A1 : RNA
# complexes), task 1 = the four SMN2 systems.
#
# The gate is graded BEFORE any S3 number is interpreted -- running both here
# saves a queue round-trip, it does not license reading S3 first. If the gate
# fails, S3 is reported as ungated and no structural claim is made.
#
# Largest input is 117 nt RNA + 182 aa protein; a 356 nt RNA took 5m45s on this
# partition (job 21071541), so a 3 h wall is generous.
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
module load cuda/13.0.1
conda activate /home/dashii/orcd/pool/boltz-env

export BOLTZ_CACHE=/home/dashii/orcd/pool/boltz-cache
export HF_HOME=/home/dashii/orcd/pool/hf-cache
S=/home/dashii/orcd/pool/smt-isef-data/boltz_s

case "$SLURM_ARRAY_TASK_ID" in
  0) IN="$S/s3gate_in"; OUT="$S/s3gate_out" ;;
  1) IN="$S/s3_in";     OUT="$S/s3_out" ;;
  *) echo "bad task"; exit 1 ;;
esac

if [ -f "$OUT/.complete" ]; then echo "$OUT already done"; exit 0; fi
mkdir -p "$OUT"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader

boltz predict "$IN" \
    --out_dir "$OUT" \
    --cache "$BOLTZ_CACHE" \
    --accelerator gpu --devices 1 \
    --diffusion_samples 5 \
    --output_format mmcif \
    --use_msa_server --override
touch "$OUT/.complete"
echo "S3 task $SLURM_ARRAY_TASK_ID complete"
