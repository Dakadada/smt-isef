#!/bin/bash
#SBATCH --job-name=mv2-j6-array
#SBATCH --partition=mit_normal
#SBATCH --time=04:00:00
#SBATCH --mem=4G
#SBATCH --array=0-127
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%A_%a.out
# J6 genome-wide pre-mRNA (gene-body, exons+introns) site universe: F18/F14/F10, <=3 mismatches. 3UP8 excluded (placement unconfirmed).
set -euo pipefail
cd /home/dashii/smt-isef; export PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
S=$HOME/orcd/scratch/feng
exec /home/dashii/orcd/pool/smt-isef-env/bin/python -m conference.model_v2.j6_site_universe \
  --fasta $S/gencode/v50/GRCh38.primary_assembly.genome.fa --gtf $S/model_v2/gtf_genes_v50.gtf \
  --output-dir $S/model_v2/j6_genome_v1 --asos F18 F14 F10 --max-mismatch 3 \
  --shard-index ${SLURM_ARRAY_TASK_ID} --shard-count 128
