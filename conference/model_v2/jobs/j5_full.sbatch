#!/bin/bash
#SBATCH --job-name=mv2-j5-full
#SBATCH --partition=mit_preemptable
#SBATCH --time=24:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# SKELETON: no submission made. Full J6 shard list is supplied explicitly.
# J5_QUERY_LIST is an ordered newline list of committed J6 *.parquet paths,
# normally from ~/orcd/scratch/feng/model_v2/j6_genome_v1. Each .done hash/count
# is verified before reading. No recursive discovery of scratch.
# Outputs checkpoint by query shard. Requeue skips matching verified .done files.
# Optional GC controls: set J5_FASTA to a pinned GRCh38 primary FASTA + .fai;
# build index with J5_GENE_BOUNDS first. Query-supplied gene bounds also work.
# Without reference/bounds, controls remain explicitly not_available.
# Runtime, full index disk size, and production RSS remain unmeasured.
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
: "${J5_QUERY_LIST:?Set an ordered newline list of committed J6 Parquet paths}"
args=(--index-dir "${J5_INDEX:-$HOME/orcd/scratch/feng/model_v2/j5_index_v1}"
      --out "${J5_FULL_OUT:-$HOME/orcd/scratch/feng/model_v2/j5_j6_v1}"
      --query-list "$J5_QUERY_LIST" --require-query-done)
if [[ -n ${J5_FASTA:-} ]]; then args+=(--fasta "$J5_FASTA"); fi
if [[ -n ${J5_GENE_BOUNDS:-} ]]; then args+=(--genes "$J5_GENE_BOUNDS"); fi
exec /usr/bin/time -v /home/dashii/orcd/pool/smt-isef-env/bin/python \
  -m conference.model_v2.j5_rbp_context "${args[@]}"
