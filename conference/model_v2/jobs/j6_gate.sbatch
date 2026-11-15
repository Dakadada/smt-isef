#!/bin/bash
#SBATCH --job-name=mv2-j6-gate
#SBATCH --partition=mit_normal
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
/home/dashii/orcd/pool/smt-isef-env/bin/python /home/dashii/smt-isef/conference/model_v2/jobs/j6_gate.py $HOME/orcd/scratch/feng/model_v2/j6_genome_v1
