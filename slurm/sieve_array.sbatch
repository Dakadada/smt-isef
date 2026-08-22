#!/bin/bash
#SBATCH --job-name=aso-sieve
#SBATCH --account=mit_general
#SBATCH --partition=mit_normal
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --array=0-199
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=/home/dashii/smt-isef/slurm/logs/sieve-%A_%a.out
#SBATCH --error=/home/dashii/smt-isef/slurm/logs/sieve-%A_%a.err

# A5: transcriptome-wide off-target hybridization site scan for one oligo
# (default nusinersen) against the pinned GENCODE v44 canonical-transcript
# index (smt.sieve.transcriptome). Stage 1 (seed+extend) is cheap; stage 2
# (IntaRNA dG_hybrid + RNAplfold dG_unfold) runs only on survivors, so shard
# runtime is dominated by however many stage-1 candidates happen to fall in
# that shard's slice of transcripts -- most shards will be seconds, a few
# with a real hit will take longer. 200 shards over ~20,000 transcripts is
# ~100 transcripts/shard, matching slurm/sweep_array.sbatch's sizing.
#
# Idempotent + shardable exactly like smt.validate.length_sweep: each shard
# writes its own scan_XXXX_of_XXXX.csv, skips if that file already exists,
# and writes atomically via a .tmp -> os.replace rename, so --requeue on
# preemption/failure costs at most one shard's work.

set -euo pipefail
source /usr/share/lmod/lmod/init/bash
module load miniforge/25.11.0-0
conda activate /home/dashii/orcd/pool/smt-isef-env

cd /home/dashii/smt-isef
export PYTHONPATH=/home/dashii/smt-isef/src

python -m smt.sieve.scan \
    --shard "${SLURM_ARRAY_TASK_ID}" \
    --nshards 200 \
    --out /home/dashii/orcd/pool/smt-isef-data/sieve/scan

# Once ALL 200 shards have completed (check with `ls .../sieve/scan | wc -l`),
# merge + classify with:
#   python -m smt.sieve.polarity \
#       --scan-dir /home/dashii/orcd/pool/smt-isef-data/sieve/scan \
#       --out /home/dashii/orcd/pool/smt-isef-data/sieve/annotated.csv
