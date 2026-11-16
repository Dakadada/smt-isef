#!/bin/bash
#SBATCH --job-name=mv2-j5-index
#SBATCH --partition=mit_preemptable
#SBATCH --time=04:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# Index build only. Checkpoints each source file; index.done appears last.
# No implicit assembly override/liftover. Without accession-bound ENCODE metadata
# or an explicit file header, files are recorded as refused, not as negatives.
# J5_METADATA: JSON accession -> original ENCODE File object (including md5sum).
# J5_GENE_BOUNDS: optional GRCh38/GENCODE v50 TSV described in the module docstring.
# Resource requests are provisional; production size/RSS has not been measured.
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
args=(--encode-root "$HOME/orcd/scratch/feng/encode"
      --out "${J5_INDEX:-$HOME/orcd/scratch/feng/model_v2/j5_index_v1}")
if [[ -n ${J5_METADATA:-} ]]; then args+=(--metadata "$J5_METADATA"); fi
if [[ -n ${J5_GENE_BOUNDS:-} ]]; then args+=(--genes "$J5_GENE_BOUNDS"); fi
exec /usr/bin/time -v /home/dashii/orcd/pool/smt-isef-env/bin/python \
  -m conference.model_v2.j5_rbp_context "${args[@]}"
