#!/bin/bash
#SBATCH --job-name=mv2-j8-full
#SBATCH --partition=mit_preemptable
#SBATCH --time=08:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# Full-genome J8 join. No inference, no GPU, single core, single writer.
# Resume: every ingested input and every output shard is an atomic committed
# checkpoint, so a preemption + requeue re-enters and skips completed work.
#
# Budget extrapolated from measured single-shard costs (login-node probe, 2026-11-15):
#   J3   256 shards x 14.1 s /  92.7 MB index  = ~60 min / 23.7 GB
#   J2p1  32 shards x  5.7 s /  14.4 MB index  = ~ 3 min /  0.5 GB
#   J2p2 128 shards x 25.5 s /  64.9 MB index  = ~54 min /  8.3 GB
#   J4  20395 configs x 216 rows (~0.2 MB each) = ~10 min /  4.1 GB
#   join 11,002,701 rows at 1542 rows/s         = ~119 min / 1.4 GB parquet
# Totals: ~4 h wall, ~37 GB SQLite index, ~1.4 GB output over 52,424 shards.
# --max-index-gib must therefore exceed the 20 GiB default; 48 leaves headroom.
# Check free SCRATCH in ~/orcd/.quota before submitting: this needs ~40 GB.
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
exec /home/dashii/orcd/pool/smt-isef-env/bin/python -m conference.model_v2.j8_feature_join \
  --root "$HOME/orcd/scratch/feng/model_v2" \
  --out "$HOME/orcd/scratch/feng/model_v2/j8_genome_v1" --max-index-gib 48
