#!/bin/bash
#SBATCH --job-name=mv2-j8-full
#SBATCH --partition=mit_preemptable
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# Full-genome J8 join. No inference, no GPU, single core, single writer.
# Resume: every ingested input and every output shard is an atomic committed
# checkpoint, so a preemption + requeue re-enters and skips completed work.
#
# Revised J8J10-FIX sizing (synthetic density bench; not a production measurement):
#   1,000 configs x 216 cells = 216,000 indexed rows, all 18nt, clustered per locus.
#   Miss work: 58 -> 41 VM steps as the prefix grows 100x; 216-hit lookup: 1,556.
#   Annotation: 329..340 dense rows/s, 9,534..10,161 miss rows/s; synthetic peak RSS ~168 MiB.
#   Dense slice index: 67.1 MB; actual atlas payloads differ, so retain prior disk
#   projection below. The added width index is small and does not enter RAM in full.
# Prior per-source measurements (2026-11-15), retained as extrapolation inputs:
#   J3   256 x 14.1 s / 92.7 MB = 60 min / 23.7 GB
#   J2   32 x 5.7 s + 128 x 25.5 s = 57 min / 8.8 GB
#   J4   20,395 x 216 cells = 4,405,320 rows; ~10 min / 4.1 GB
#   Prior other join/IO allowance: 119 min for 11,002,701 rows / 52,424 shards.
# Add 19..557 min for new miss..dense annotation => ~4.4..13.4 h total;
# this conservatively counts some join overhead twice. Request 24 h for shared IO.
# This is a density scenario envelope, NOT a bound for arbitrarily overlapping loci.
# RAM estimate ~2 GiB from prior 1.6 GB J3 row-group peak + 32 MiB SQLite cache
# + bounded Arrow/JSON batches; retain 8G allocation. Full production RSS unmeasured.
# Disk estimate ~37 GB index (+ ~0.05 GB width index); 48 GiB cap retained.
# Output estimate ~1.4 GB is prior sparse-overlap extrapolation, NOT a dense bound.
# Recheck cached SCRATCH quota before submission; dense overlaps can enlarge output.
# Fresh v2 directory required: old catalogs/contracts lack authenticated marker paths.
set -euo pipefail
cd /home/dashii/smt-isef
export PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
exec /home/dashii/orcd/pool/smt-isef-env/bin/python -m conference.model_v2.j8_feature_join \
  --root "$HOME/orcd/scratch/feng/model_v2" \
  --out "$HOME/orcd/scratch/feng/model_v2/j8_genome_v2" --max-index-gib 48
