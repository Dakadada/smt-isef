#!/bin/bash
#SBATCH --job-name=mv2-j4-equiv
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --requeue
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%j.out
# Real-model driver-vs-CLI equivalence on 3 configs. exit 1 unless run_id and ablation sha match for each.
set -euo pipefail
cd /home/dashii/smt-isef
module load cuda/12.9.1 cudnn/9.8.0.87-cuda12 miniforge/25.11.0-0
PY=/home/dashii/orcd/pool/feng-gpu-env/bin/python
export PYTHONPATH=$PWD/feng/src:$PWD PYTHONNOUSERSITE=1 CUBLAS_WORKSPACE_CONFIG=:4096:8
P=$HOME/orcd/scratch/feng/model_v2/j4_equiv_probe; rm -rf $P; mkdir -p $P/shards; S=$P/shards/shard_9999.jsonl
head -3 $HOME/orcd/scratch/feng/model_v2/j4_configs/shards/shard_0000.jsonl > $S
$PY conference/model_v2/jobs/j4_scale_driver.py --shard-file $S --atlas-dir $P/driver --task-id 9999
fail=0
while read -r L; do
  I=$(jq -r .config_id <<<"$L"); printf '%s' "$L" > $P/$I.json
  $PY -m conference.model_v2.gpu_perturbation_atlas --config $P/$I.json --batch-size 32 --output-dir $P/cli/$I
  a=$(jq -r '.run_id+" "+.ablation.sha256' $P/driver/$I/shard_0.done); b=$(jq -r '.run_id+" "+.ablation.sha256' $P/cli/$I/shard_0.done)
  echo "$I driver=[$a] cli=[$b]"; [ "$a" = "$b" ] || fail=1
done < $S
[ $fail = 0 ] && echo EQUIVALENCE_OK || { echo EQUIVALENCE_MISMATCH; exit 1; }
