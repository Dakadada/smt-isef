#!/bin/bash
#SBATCH --job-name=mv2-j3-size
#SBATCH --partition=mit_normal
#SBATCH --time=01:00:00
#SBATCH --mem=6G
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# sizing: 1 of 1024 shards of the full J6 genome universe; reports wall + RSS
set -euo pipefail
cd /home/dashii/smt-isef; export PYTHONPATH=feng/src OMP_NUM_THREADS=1
S=$HOME/orcd/scratch/feng
exec /usr/bin/time -v /home/dashii/orcd/pool/smt-isef-env/bin/python -m conference.model_v2.j3_scoring \
  --input-dir $S/model_v2/j6_genome_v1 --output-dir $S/model_v2/j3_size_test --cache $S/model_v2/j3_size_test/cache.sqlite \
  --fasta $S/gencode/v50/GRCh38.primary_assembly.genome.fa --vienna-python /home/dashii/orcd/pool/feng-env/bin/python \
  --flank 80 --window 120 --span 60 --shard-index 0 --shard-count 1024 --requeue-safe
