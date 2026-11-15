# J4-DRIVER result (no submit, no GPU run, no commit)
Tests (brief's exact command): `39 passed in 2.01s`, was `4 failed, 35 passed`. 1 code bug, 3 bad tests:
- CODE `drive()`: skip branch never flushed -> an all-skip resume left the zeroed initial checkpoint (`skipped=0` while stdout said 3). Flushes per skip now.
- TEST flush-hook: captured the LAST `checkpoint_on_signal` hook = `run_atlas`'s nested per-config flush (shared `ck` module), not the task-level one. Takes `hooks[1]`, asserts `2+len(TOY)` hooks, all `requeue=False`.
- TEST sigterm-mid-config: killed at `len(calls)==3`, which lands in config A, so A never completed. Now fires once A is `.done` and B has 1 part -> real partial resume (`next==len(parts)==1`).
- TEST array2-sbatch: counted `j4_scale_driver.py` incl. the header comment; counts non-comment lines only now.
Spec verified: one scorer load/task; driver == per-config CLI on stub (same `run_id`, `ablation.sha256`, parquet rows); idempotent skip (done+parquet+sha); SIGTERM flush with the shell as sole requeue owner; `j4_scale_gate.py`; `--time=00:45:00`.
Submit, NOT executed (319 shards + 20,395 emitted configs confirmed; batch 32 = bench 29169335 winner):
    cd /home/dashii/smt-isef && sbatch --array=0-318%16 conference/model_v2/jobs/j4_scale_array2.sbatch
    python conference/model_v2/jobs/j4_scale_gate.py --deep --report conference/model_v2/ops/briefs/j4_gate.json
Walltime: 64 x 3.16 s + 23 s load = ~3.8 min/task (45 min cap = 12x headroom); 20 waves at %16 = ~1.3 h wall clean, ~20 GPU-h total vs the ~130 GPU-h startup waste of v1.
Post-submit real-model equivalence, 3 configs (GPU node, real Pangolin; PY=/home/dashii/orcd/pool/feng-gpu-env/bin/python):
    P=~/orcd/scratch/feng/model_v2/j4_equiv_probe; mkdir -p $P/shards; S=$P/shards/shard_9999.jsonl
    head -3 ~/orcd/scratch/feng/model_v2/j4_configs/shards/shard_0000.jsonl > $S
    srun -p mit_preemptable --gres=gpu:1 --mem=32G -c 4 -t 30 bash -c "export PYTHONPATH=\$PWD/feng/src:\$PWD; \
      \$PY conference/model_v2/jobs/j4_scale_driver.py --shard-file $S --atlas-dir $P/driver --task-id 9999; \
      while read -r L; do I=\$(jq -r .config_id <<<\"\$L\"); printf '%s' \"\$L\" > $P/\$I.json; \
      \$PY -m conference.model_v2.gpu_perturbation_atlas --config $P/\$I.json --batch-size 32 --output-dir $P/cli/\$I; \
      jq -r '.run_id,.ablation.sha256' $P/driver/\$I/shard_0.done $P/cli/\$I/shard_0.done; done < $S"
    # PASS iff each config prints its run_id and ablation sha256 twice (driver line == CLI line).
