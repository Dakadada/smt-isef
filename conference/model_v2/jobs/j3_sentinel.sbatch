#!/bin/bash
#SBATCH --job-name=mv2-j3-sentinel
#SBATCH --partition=mit_quicktest
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src${PYTHONPATH:+:$PYTHONPATH}
export OMP_NUM_THREADS=1
exec /home/dashii/orcd/pool/smt-isef-env/bin/python \
  -m conference.model_v2.j3_scoring \
  --input-dir "$HOME/orcd/scratch/feng/model_v2/j6_sentinel" \
  --output-dir "$HOME/orcd/scratch/feng/model_v2/j3_sentinel" \
  --fasta "$HOME/orcd/scratch/feng/gencode/v50/GRCh38.primary_assembly.genome.fa" \
  --vienna-python /home/dashii/orcd/pool/feng-env/bin/python \
  --flank 80 --window 120 --span 60 --requeue-safe
