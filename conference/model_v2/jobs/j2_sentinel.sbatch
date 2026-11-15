#!/bin/bash
#SBATCH --job-name=mv2-j2-sentinel
#SBATCH --partition=mit_quicktest
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
D=/home/dashii/smt-isef/conference/model_v2
FA=$HOME/orcd/scratch/feng/gencode/v50/GRCh38.primary_assembly.genome.fa
OUT=$HOME/orcd/scratch/feng/model_v2/j2_sentinel
$HOME/orcd/pool/feng-env/bin/python $D/jobs/j2_fold.py $D/manifests/j2_sentinel_loci.tsv $FA $OUT 0 1
