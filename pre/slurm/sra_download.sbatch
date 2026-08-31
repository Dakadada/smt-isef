#!/bin/bash
#SBATCH --job-name=sra-dl
#SBATCH --account=mit_general
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --array=0-19
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/sra-%A_%a.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/sra-%A_%a.err
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
conda activate /home/dashii/orcd/pool/rna-env

SC=/home/dashii/orcd/scratch/smt-rna
RUNS=($(cat /home/dashii/smt-isef/slurm/sra_runs.txt))
ACC=${RUNS[$SLURM_ARRAY_TASK_ID]}

# Idempotent: a .done marker is written only after both mates land intact.
if [ -f "$SC/fastq/$ACC.done" ]; then echo "$ACC already downloaded"; exit 0; fi

cd "$SC/fastq"
rm -f "${ACC}_1.fastq.gz.tmp" "${ACC}_2.fastq.gz.tmp"
fasterq-dump "$ACC" --split-files --threads 6 --temp "$SC/fastq/_tmp_$ACC" -O "$SC/fastq"
pigz -p 6 -f "${ACC}_1.fastq" "${ACC}_2.fastq"
rm -rf "$SC/fastq/_tmp_$ACC"
touch "$SC/fastq/$ACC.done"
echo "$ACC downloaded"
