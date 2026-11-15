#!/bin/bash
#SBATCH --job-name=mv2-j3-array
#SBATCH --partition=mit_normal
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%A_%a.out
# J3 channels over full J6 genome universe; 256 shards (size test 1/1024 shard = 10m42s, 1.2GB RSS => ~43 min each)
set -euo pipefail
cd /home/dashii/smt-isef; export PYTHONPATH=feng/src OMP_NUM_THREADS=1
S=$HOME/orcd/scratch/feng; O=$S/model_v2/j3_genome_v1
exec /home/dashii/orcd/pool/smt-isef-env/bin/python -m conference.model_v2.j3_scoring \
  --input-dir $S/model_v2/j6_genome_v1 --output-dir $O --cache $O/cache_${SLURM_ARRAY_TASK_ID}.sqlite \
  --fasta $S/gencode/v50/GRCh38.primary_assembly.genome.fa --vienna-python /home/dashii/orcd/pool/feng-env/bin/python \
  --flank 80 --window 120 --span 60 --shard-index ${SLURM_ARRAY_TASK_ID} --shard-count 256 --requeue-safe
