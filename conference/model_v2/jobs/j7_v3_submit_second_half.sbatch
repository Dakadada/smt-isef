#!/bin/bash
#SBATCH --job-name=mv2-j7v3-submit
#SBATCH --partition=mit_normal
#SBATCH --time=00:05:00
#SBATCH --mem=256M
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# per-user submit cap blocks 256-task array: tasks 0-127 submitted first; this submits 128-255 + merge after 0-127 succeeds.
export J7_ROOT="$HOME/orcd/scratch/feng/model_v2/j7_context_v3" J3_ROOT="$HOME/orcd/scratch/feng/model_v2/j3_genome_v1" J7_EXPECTED_ROWS=11002701
cd /home/dashii/smt-isef && J=conference/model_v2/jobs && A=$(sbatch --parsable --partition=mit_preemptable --export=ALL --array=128-255%3 $J/j7_array.sbatch) && sbatch --partition=mit_preemptable --export=ALL --dependency=afterok:${A%%;*} $J/j7_merge.sbatch
