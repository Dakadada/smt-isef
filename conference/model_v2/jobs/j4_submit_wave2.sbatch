#!/bin/bash
#SBATCH --job-name=mv2-j4-wave2-submit
#SBATCH --partition=mit_normal
#SBATCH --time=00:05:00
#SBATCH --mem=256M
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# QOS submit cap (mit_preemptable 256 jobs/user): wave 2 (tasks 200-318) submitted after wave 1 ends (afterany; gate reports gaps).
cd /home/dashii/smt-isef && sbatch --array=200-318%16 conference/model_v2/jobs/j4_scale_array2.sbatch
