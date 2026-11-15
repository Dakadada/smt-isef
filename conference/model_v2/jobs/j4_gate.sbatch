#!/bin/bash
#SBATCH --job-name=mv2-j4-gate
#SBATCH --partition=mit_normal
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
cd /home/dashii/smt-isef; export PYTHONPATH=feng/src:.
/home/dashii/orcd/pool/smt-isef-env/bin/python conference/model_v2/jobs/j4_scale_gate.py --deep --report conference/model_v2/ops/briefs/j4_gate.json
