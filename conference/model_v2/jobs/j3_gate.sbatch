#!/bin/bash
#SBATCH --job-name=mv2-j3-gate
#SBATCH --partition=mit_normal
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
/home/dashii/orcd/pool/smt-isef-env/bin/python /home/dashii/smt-isef/conference/model_v2/jobs/j3_gate.py $HOME/orcd/scratch/feng/model_v2/j3_genome_v1 256 11002701
