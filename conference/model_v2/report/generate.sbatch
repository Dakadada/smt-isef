#!/bin/bash
#SBATCH --job-name=mv2-j11-report
#SBATCH --partition=mit_preemptable
#SBATCH --requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/report/j11-%j.out
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
/home/dashii/orcd/pool/smt-isef-env/bin/python conference/model_v2/j11_report.py generate
/home/dashii/orcd/pool/smt-isef-env/bin/python conference/model_v2/j11_report.py verify
