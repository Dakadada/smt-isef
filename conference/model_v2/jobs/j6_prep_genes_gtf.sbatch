#!/bin/bash
#SBATCH --job-name=mv2-j6-prep
#SBATCH --partition=mit_normal
#SBATCH --time=00:30:00
#SBATCH --mem=1G
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# one-time: header + gene rows only (avoids every array task rescanning 4.7GB GTF)
set -e
S=$HOME/orcd/scratch/feng; G=$S/gencode/v50/gencode.v50.primary_assembly.annotation.gtf; O=$S/model_v2/gtf_genes_v50.gtf
mkdir -p $S/model_v2
awk -F'\t' '/^#/ || $3=="gene"' $G > $O.tmp && mv $O.tmp $O
wc -l $O; sha256sum $O; head -3 $O | cut -c1-120
