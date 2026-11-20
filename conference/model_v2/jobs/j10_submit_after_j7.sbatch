#!/bin/bash
#SBATCH --job-name=mv2-j10-submit
#SBATCH --partition=mit_normal
#SBATCH --time=00:05:00
#SBATCH --mem=256M
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# association submit cap 500: J10 array (256) submitted after J7 (258) finishes.
cd /home/dashii/smt-isef && A=$(sbatch --parsable conference/model_v2/jobs/j10_array.sbatch) && sbatch --dependency=afterok:${A%%;*} conference/model_v2/jobs/j10_merge.sbatch
