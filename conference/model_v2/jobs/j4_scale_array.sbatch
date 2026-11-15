#!/usr/bin/env bash
# J4 scale-out atlas array runner: one task per config shard, one engine
# invocation per config, idempotent and requeue-safe.
#
# NO --array / --time directive on purpose: both are set at submit time from the
# bench numbers (jobs/j4_scale_bench.sbatch), e.g.
#   NSHARDS=$(grep -vc '^#' ~/orcd/scratch/feng/model_v2/j4_configs/j4_config_shards.tsv)
#   sbatch --array=0-$((NSHARDS-1))%8 --time=06:00:00 \
#          conference/model_v2/jobs/j4_scale_array.sbatch
# Variables: J4_SHARD_DIR, J4_ATLAS_DIR, J4_BATCH_SIZE (winning bench sweep),
#            J4_GPU_ENV, J4_FASTA.
# Masked rows keep evidence_class sequence_ablation_proxy (set by the engine).
#SBATCH --job-name=mv2-j4-scale-array
#SBATCH --account=mit_general
#SBATCH --partition=mit_preemptable
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --requeue
#SBATCH --signal=B:SIGTERM@60
#SBATCH --output=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%A_%a.out
#SBATCH --error=/home/dashii/smt-isef/conference/model_v2/ops/logs/%x-%A_%a.err
set -euo pipefail
cd /home/dashii/smt-isef
: "${SLURM_ARRAY_TASK_ID:?submit with --array=0-<NSHARDS-1>; see header}"
module load cuda/12.9.1
module load cudnn/9.8.0.87-cuda12
module load miniforge/25.11.0-0
GPU_ENV="${J4_GPU_ENV:-/home/dashii/orcd/pool/feng-gpu-env}"
PYTHON_BIN="${GPU_ENV}/bin/python"
test -x "${PYTHON_BIN}"
export PATH="${GPU_ENV}/bin:${PATH}"
export PYTHONPATH="${PWD}/feng/src:${PWD}"
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 NUMEXPR_NUM_THREADS=4
export CUBLAS_WORKSPACE_CONFIG=:4096:8
SHARD_DIR="${J4_SHARD_DIR:-${HOME}/orcd/scratch/feng/model_v2/j4_configs/shards}"
ATLAS_DIR="${J4_ATLAS_DIR:-${HOME}/orcd/scratch/feng/model_v2/j4_atlas_v1}"
BATCH_SIZE="${J4_BATCH_SIZE:-16}"
SHARD_FILE="${SHARD_DIR}/$(printf 'shard_%04d.jsonl' "${SLURM_ARRAY_TASK_ID}")"
test -s "${SHARD_FILE}"
mkdir -p "${ATLAS_DIR}"
FASTA_ARGS=()
if [[ -n "${J4_FASTA:-}" ]]; then FASTA_ARGS=(--fasta "${J4_FASTA}"); fi
# A fresh, private marker dir per invocation prevents stale readiness on requeue.
READY_DIR=$(mktemp -d "${ATLAS_DIR}/.j4-ready.XXXXXX")
export J4_READY_FILE="${READY_DIR}/ready"
trap 'rm -f -- "${J4_READY_FILE}" "${READY_DIR}/plan.tsv"; rmdir -- "${READY_DIR}" 2>/dev/null || true' EXIT

# Sole requeue owner, following feng/scripts/templates/gpu_preemptable.sbatch.
# The engine installs checkpoint_on_signal(requeue=False), flushes, and exits 0.
PAYLOAD_PID=""
requeue_handler() {
    trap '' SIGTERM SIGUSR1
    if [[ -n "${PAYLOAD_PID}" ]]; then
        kill -TERM "${PAYLOAD_PID}" 2>/dev/null || true
        if ! wait "${PAYLOAD_PID}"; then
            # Read after the child exits: it may have installed its handler
            # concurrently with TERM. Before readiness it cannot have written
            # atlas work, so a default TERM exit is safe to requeue.
            if [[ -e "${J4_READY_FILE}" ]]; then
                echo 'ERROR: payload checkpoint flush failed; refusing requeue' >&2
                exit 1
            fi
        fi
    fi
    if [[ -z "${SLURM_JOB_ID:-}" ]] || ! scontrol requeue "${SLURM_JOB_ID}"; then
        echo 'ERROR: checkpoint saved but requeue failed' >&2
        exit 1
    fi
    exit 0
}
# Installed before any work and kept installed BETWEEN configs, so preemption in
# a gap requeues cleanly instead of killing the task.
trap requeue_handler SIGTERM SIGUSR1

# Materialize one single-config JSON per engine invocation (the engine takes one
# absolute genomic tiling interval per config) and emit a plan: id, path, sha256.
PLAN="${READY_DIR}/plan.tsv"
"${PYTHON_BIN}" - "${SHARD_FILE}" "${ATLAS_DIR}" "${PLAN}" <<'PY'
import hashlib, json, os, sys
from pathlib import Path

shard, atlas_dir, plan_path = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
plan = []
for line in shard.read_text().splitlines():
    if not line.strip():
        continue
    config = json.loads(line)
    config_id = config["config_id"]
    if "/" in config_id or config_id.startswith("."):
        raise ValueError(f"unsafe config_id {config_id!r}")
    directory = atlas_dir / config_id
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / "config.json"
    temp = target.with_name(f"config.json.tmp.{os.getpid()}")
    temp.write_text(line + "\n")
    os.replace(temp, target)
    # sha256 of the canonical shard line == config_sha256 in j4_configs_manifest.tsv
    plan.append(f"{config_id}\t{target}\t{hashlib.sha256(line.encode()).hexdigest()}")
plan_path.write_text("\n".join(plan) + "\n")
print(f"PLAN shard={shard.name} configs={len(plan)}")
PY

FAILED=0 DONE=0 SKIPPED=0
while IFS=$'\t' read -r CONFIG_ID CONFIG_PATH CONFIG_SHA; do
    [[ -n "${CONFIG_ID}" ]] || continue
    CFGOUT="${ATLAS_DIR}/${CONFIG_ID}"
    # Idempotent skip: completed marker + published parquet + matching config hash.
    if [[ -s "${CFGOUT}/shard_0.done" && -s "${CFGOUT}/ablation_maps/shard_0.parquet" \
          && -f "${CFGOUT}/.j4_config_sha256" \
          && "$(cat "${CFGOUT}/.j4_config_sha256")" == "${CONFIG_SHA}" ]]; then
        SKIPPED=$((SKIPPED + 1))
        continue
    fi
    rm -f -- "${J4_READY_FILE}"
    # Defer handling across spawn/$! assignment, so even an immediate TERM has
    # a known child to stop and reap before requeue.
    PREEMPT_REQUESTED=0
    trap 'PREEMPT_REQUESTED=1' SIGTERM SIGUSR1
    "${PYTHON_BIN}" -m conference.model_v2.gpu_perturbation_atlas \
        --config "${CONFIG_PATH}" --output-dir "${CFGOUT}" \
        --batch-size "${BATCH_SIZE}" "${FASTA_ARGS[@]+"${FASTA_ARGS[@]}"}" &
    PAYLOAD_PID=$!
    trap requeue_handler SIGTERM SIGUSR1
    if (( PREEMPT_REQUESTED )); then requeue_handler; fi
    set +e
    wait "${PAYLOAD_PID}"
    STATUS=$?
    set -e
    PAYLOAD_PID=""
    if (( STATUS == 0 )); then
        printf '%s' "${CONFIG_SHA}" > "${CFGOUT}/.j4_config_sha256.tmp.$$"
        mv -f -- "${CFGOUT}/.j4_config_sha256.tmp.$$" "${CFGOUT}/.j4_config_sha256"
        DONE=$((DONE + 1))
    else
        # Keep going: one bad locus must not strand the rest of the shard.
        echo "ERROR: config ${CONFIG_ID} exited ${STATUS}" >&2
        printf '%s\t%s\n' "${CONFIG_ID}" "${STATUS}" \
            >> "${ATLAS_DIR}/failures_task_${SLURM_ARRAY_TASK_ID}.tsv"
        FAILED=$((FAILED + 1))
    fi
done < "${PLAN}"

echo "SHARD ${SLURM_ARRAY_TASK_ID} done=${DONE} skipped=${SKIPPED} failed=${FAILED}"
exit $(( FAILED > 0 ))
