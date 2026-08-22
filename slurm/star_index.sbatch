#!/bin/bash
#SBATCH --job-name=star-idx
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/star-idx-%j.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/star-idx-%j.err
set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
conda activate /home/dashii/orcd/pool/rna-env

# The toolchain may still be installing when this job starts. Wait, bounded,
# rather than failing and burning a queue slot.
for i in $(seq 1 60); do
  [ -x /home/dashii/orcd/pool/rna-env/bin/STAR ] && break
  echo "waiting for STAR ($i/60)"; sleep 60
done
command -v STAR >/dev/null || { echo "STAR never appeared"; exit 1; }

T=/home/dashii/orcd/pool/smt-isef-data/transcriptome
SC=/home/dashii/orcd/scratch/smt-rna
IDX=$SC/star_index

# Idempotent: STAR writes SA last, so its presence means the index is complete.
# Idempotency marker written only after STAR returns 0. Checking SA alone is
# unsafe under preemption: a killed build leaves a partial SA that looks done.
if [ -f "$IDX/.complete" ]; then echo "index already built, skipping"; exit 0; fi
rm -rf "$IDX"; mkdir -p "$IDX"

[ -s "$T/gencode.v44.basic.annotation.gtf" ] || \
  pigz -dc "$T/gencode.v44.basic.annotation.gtf.gz" > "$T/gencode.v44.basic.annotation.gtf"

STAR --runMode genomeGenerate \
     --runThreadN 16 \
     --genomeDir "$IDX" \
     --genomeFastaFiles "$T/GRCh38.primary_assembly.genome.fa" \
     --sjdbGTFfile "$T/gencode.v44.basic.annotation.gtf" \
     --sjdbOverhang 99 \
     --outTmpDir "$SC/_startmp_$SLURM_JOB_ID"
touch "$IDX/.complete"
echo "STAR index complete"
