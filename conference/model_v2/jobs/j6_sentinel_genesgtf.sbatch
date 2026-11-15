#!/bin/bash
#SBATCH --job-name=mv2-j6-sentinel-genesgtf
#SBATCH --partition=mit_quicktest
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src${PYTHONPATH:+:$PYTHONPATH}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
# The previous sentinel directory predates immutable run manifests; preserve it.
exec /usr/bin/time -v /home/dashii/orcd/pool/smt-isef-env/bin/python \
  -m conference.model_v2.j6_site_universe \
  --fasta "$HOME/orcd/scratch/feng/gencode/v50/GRCh38.primary_assembly.genome.fa" \
  --gtf "$HOME/orcd/scratch/feng/model_v2/gtf_genes_v50.gtf" \
  --output-dir "$HOME/orcd/scratch/feng/model_v2/j6_sentinel_genesgtf" \
  --region chr5:70000000-71100000 --asos F18 F14 F10 --max-mismatch 3 \
  --check-smn-anchors --profile-memory
