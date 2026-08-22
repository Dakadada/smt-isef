#!/bin/bash
#SBATCH --job-name=star-aln
#SBATCH --account=mit_general
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --array=0-19
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/aln-%A_%a.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/aln-%A_%a.err
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
conda activate /home/dashii/orcd/pool/rna-env

SC=/home/dashii/orcd/scratch/smt-rna
RUNS=($(cat /home/dashii/smt-isef/slurm/sra_runs.txt))
ACC=${RUNS[$SLURM_ARRAY_TASK_ID]}
OUT=$SC/bam/$ACC

if [ -f "$OUT.done" ]; then echo "$ACC already aligned"; exit 0; fi
[ -f "$SC/fastq/$ACC.done" ] || { echo "fastq for $ACC not ready"; exit 1; }

STAR --genomeDir "$SC/star_index" \
     --readFilesIn "$SC/fastq/${ACC}_1.fastq.gz" "$SC/fastq/${ACC}_2.fastq.gz" \
     --readFilesCommand zcat \
     --runThreadN 8 \
     --outFileNamePrefix "${OUT}." \
     --outSAMtype BAM SortedByCoordinate \
     --outTmpDir "$SC/_alntmp_${ACC}_$SLURM_JOB_ID"
samtools index "${OUT}.Aligned.sortedByCoord.out.bam"
touch "$OUT.done"
echo "$ACC aligned"
