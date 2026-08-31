#!/bin/bash
#SBATCH --job-name=rmats
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/rmats-%j.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/rmats-%j.err
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
conda activate /home/dashii/orcd/pool/rna-env

SC=/home/dashii/orcd/scratch/smt-rna
T=/home/dashii/orcd/pool/smt-isef-data/transcriptome
OUT=$SC/rmats/out

# Idempotent: marker written only after rmats returns 0.
if [ -f "$OUT/.complete" ]; then echo "rmats already done"; exit 0; fi
rm -rf "$OUT" "$SC/rmats/tmp"; mkdir -p "$OUT" "$SC/rmats/tmp"

# Read length 100 (2x100 paired, avgLength 200 per SRA runinfo).
# --variable-read-length tolerates trimmed reads without discarding them.
rmats.py --b1 "$SC/rmats/b1.txt" --b2 "$SC/rmats/b2.txt" \
         --gtf "$T/gencode.v44.basic.annotation.gtf" \
         -t paired --readLength 100 --variable-read-length \
         --nthread 8 --od "$OUT" --tmp "$SC/rmats/tmp"
touch "$OUT/.complete"
echo "rmats complete"
