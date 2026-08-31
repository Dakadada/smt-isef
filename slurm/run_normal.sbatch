#!/usr/bin/env bash
#SBATCH --job-name=rawread-normal
#SBATCH --cpus-per-task=24
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:0
#SBATCH --output=slurm/logs/%x-%j.out
set -euo pipefail
: "${CONFIG:=config/pipeline.yaml}"
: "${STAGE:=all}"
mkdir -p slurm/logs
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
python3 -m rawread_reanalysis "$STAGE" --config "$CONFIG"
