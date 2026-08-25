#!/bin/bash
#SBATCH --job-name=boltz-s5-smn2
#SBATCH --account=mit_general
#SBATCH --partition=mit_normal_gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/boltz-s5-smn2-%j.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/boltz-s5-smn2-%j.err

# S5 HELD-OUT GATE ONLY. Four deposited protein-free RNA : small-molecule
# complexes (8QMH, 7QSH, 1F27, 6XRQ), all found by querying RCSB and screened
# to exclude every risdiplam / SMN-C5 / branaplam-class entry -- gating on the
# Campagne deposition and then predicting it would be a memorisation test.
#
# THE SMN2 PREDICTION IS DELIBERATELY NOT IN THIS JOB. The gate is graded, and
# read, before any target prediction exists. If it fails, S5 is reported
# ungated and makes no claim -- exactly as S4 was declined and S3 was reported
# UNRESOLVED.
#
# Sizes are tiny by this project's standards (24-30 nt RNA + one ligand)
# against the 356 nt RNA that took 5m45s on this partition, so 2 h is ample.
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
module load cuda/13.0.1
conda activate /home/dashii/orcd/pool/boltz-env

export BOLTZ_CACHE=/home/dashii/orcd/pool/boltz-cache
export HF_HOME=/home/dashii/orcd/pool/hf-cache
S=/home/dashii/orcd/pool/smt-isef-data/boltz_s5
IN="$S/smn2_in"
OUT="$S/smn2_out"

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
echo "S5 SMN2 UNGATED complete"
