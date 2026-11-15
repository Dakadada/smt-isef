#!/bin/bash
#SBATCH --job-name=mv2-j2-array-p2
#SBATCH --partition=mit_normal
#SBATCH --time=03:00:00
#SBATCH --mem=3G
#SBATCH --array=0-127
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%A_%a.out
# idempotent: j2_fold.py skips shards with .done; sha-verified deterministic (see idem test 29162809)
D=/home/dashii/smt-isef/conference/model_v2
FA=$HOME/orcd/scratch/feng/gencode/v50/GRCh38.primary_assembly.genome.fa
exec $HOME/orcd/pool/feng-env/bin/python $D/jobs/j2_fold.py $D/manifests/j2_loci_phase2.tsv $FA $HOME/orcd/scratch/feng/model_v2/j2_phase2 ${SLURM_ARRAY_TASK_ID} 128
