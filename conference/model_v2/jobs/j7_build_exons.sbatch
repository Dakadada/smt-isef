#!/bin/bash
#SBATCH --job-name=mv2-j7-exons
#SBATCH --partition=mit_normal
#SBATCH --time=00:30:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
D=conference/model_v2
J7_ROOT=${J7_ROOT:-$HOME/orcd/scratch/feng/model_v2/j7_context_v2}
mkdir -p "$J7_ROOT"
# Exactly one full GTF scan; never run this on the login node.
exec /home/dashii/orcd/pool/smt-isef-env/bin/python -m conference.model_v2.j7_context_features build \
  --gtf "$HOME/orcd/scratch/feng/gencode/v50/gencode.v50.primary_assembly.annotation.gtf" \
  --manifests "$D/manifests/j2_loci_phase1.tsv" "$D/manifests/j2_loci_phase2.tsv" \
  --output "$J7_ROOT/canonical_exons.json"
