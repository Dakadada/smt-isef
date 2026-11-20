#!/bin/bash
#SBATCH --job-name=mv2-j9-smn2
#SBATCH --partition=mit_quicktest
#SBATCH --time=00:10:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# J9: build the SMN2 ASO feature/ordering dataset WITH real ViennaRNA channels.
# CPU only, 3 ASOs x 3 prespecified fold alternatives = 9 duplexfold + 9 plfold calls.
# Never run the ViennaRNA worker on the login node; this script is the only route.
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src${PYTHONPATH:+:$PYTHONPATH}
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export ARROW_NUM_THREADS=1

OUT_DIR="${J9_OUT_DIR:-$HOME/orcd/scratch/feng/model_v2/j9_smn2_v2}"
ATLAS_DIR="${J9_ATLAS_DIR:-$HOME/orcd/scratch/feng/model_v2/j4_sentinel_v2}"

mkdir -p "$OUT_DIR"
exec /home/dashii/orcd/pool/smt-isef-env/bin/python \
  -m conference.model_v2.j9_smn2_prediction build \
  --out-dir "$OUT_DIR" \
  --atlas-dir "$ATLAS_DIR" \
  --fasta "$HOME/orcd/scratch/feng/gencode/v50/GRCh38.primary_assembly.genome.fa" \
  --vienna-python /home/dashii/orcd/pool/feng-env/bin/python \
  --cache "$OUT_DIR"
