#!/usr/bin/env python
"""J4 scale-out in-process atlas driver: ONE model load per array task.

The per-config CLI path (jobs/j4_scale_array.sbatch) spawns one
``python -m conference.model_v2.gpu_perturbation_atlas`` per config, so each of
the 20,395 configs pays the measured 23.1 s Pangolin load => ~130 GPU-h of pure
startup. This driver reuses the in-process mechanism already proven by
jobs/j4_scale_bench.py: authenticate the reference once, construct
``PangolinScorer`` once, then call the engine's public ``run_atlas`` for every
config in the shard.

Semantics are deliberately identical to the per-config CLI path, config by
config:
  * output tree ``<atlas-dir>/<config_id>/`` with the engine's own
    ``shard_0.done`` / ``shard_0.checkpoint.json`` / ``parts/`` /
    ``ablation_maps/shard_0.parquet`` / ``reference_scores/`` layout (the
    reference cache stays per config, exactly as the CLI default does, so
    manifest-cited artifact paths match),
  * ``config.json`` written atomically into the config dir (the canonical shard
    line verbatim, as the shell plan step wrote it),
  * sha-verified config: the sha256 of the canonical shard line is checked
    against a re-canonicalization of the parsed JSON, and published to
    ``.j4_config_sha256`` only AFTER a clean run,
  * idempotent skip when ``shard_0.done`` + published parquet + matching
    ``.j4_config_sha256`` are all present,
  * a config whose run raises, or completes with failed rows, is an explicit
    failure row (``failures_task_<id>.tsv``, same two columns the shell runner
    wrote, plus ``failures_task_<id>.jsonl`` with the reason) and never blocks
    the rest of the shard,
  * SIGTERM/SIGUSR1 flushes a checkpoint and returns; the batch shell
    (jobs/j4_scale_array2.sbatch, following jobs/j4_sentinel.sbatch) stays the
    SOLE requeue owner. ``run_atlas`` installs its own nested flush, so a
    preemption mid-config resumes that config from its engine checkpoint.

Masked rows keep evidence_class sequence_ablation_proxy (set by the engine).
No stub scorer path: this driver only ever runs real Pangolin.

Usage (single invocation per array task):
  python conference/model_v2/jobs/j4_scale_driver.py \
      --shard-dir ~/orcd/scratch/feng/model_v2/j4_configs/shards \
      --atlas-dir ~/orcd/scratch/feng/model_v2/j4_atlas_v1 \
      --task-id ${SLURM_ARRAY_TASK_ID} --batch-size 32
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

from conference.model_v2 import gpu_perturbation_atlas as atlas
from conference.model_v2 import trans_intervention as ti
from feng.gpu import checkpoint as ck

DEFAULT_SHARD_DIR = Path("~/orcd/scratch/feng/model_v2/j4_configs/shards")
DEFAULT_ATLAS_DIR = Path("~/orcd/scratch/feng/model_v2/j4_atlas_v1")
DEFAULT_BATCH_SIZE = 32          # winning bench sweep: 17.4 inf/s, 0.28 GiB peak
SHA_MARKER = ".j4_config_sha256"
DONE_NAME = "shard_0.done"
PARQUET_NAME = "ablation_maps/shard_0.parquet"


def canonical(value) -> str:
    """Byte-for-byte the builder's canonical form (j4_build_configs.canonical)."""
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def shard_file_for(shard_dir, task_id: int) -> Path:
    if type(task_id) is not int or task_id < 0:
        raise ValueError("task id must be a non-negative integer")
    return Path(shard_dir).expanduser() / f"shard_{task_id:04d}.jsonl"


def safe_config_id(config_id) -> str:
    if (not isinstance(config_id, str) or not config_id or "/" in config_id
            or config_id.startswith(".") or config_id in (".", "..")):
        raise ValueError(f"unsafe config_id {config_id!r}")
    return config_id


def read_plan(shard_file) -> list[dict]:
    """Parse a shard jsonl into per-config plan records, fail-closed on duplicates.

    ``sha256`` is the sha256 of the canonical shard line (no newline), i.e. the
    ``config_sha256`` column of j4_configs_manifest.tsv.
    """
    path = Path(shard_file).expanduser()
    records: list[dict] = []
    seen: set[str] = set()
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        config = json.loads(line)
        config_id = safe_config_id(config.get("config_id"))
        if config_id in seen:
            raise ValueError(f"{path}:{number}: duplicate config_id {config_id!r}")
        seen.add(config_id)
        records.append({"config_id": config_id, "text": line,
                        "sha256": hashlib.sha256(line.encode()).hexdigest(),
                        "config": config, "line": number})
    if not records:
        raise ValueError(f"{path}: no configs")
    return records


def publish_config(record: dict, config_dir: Path) -> Path:
    """Write the canonical shard line to <config_dir>/config.json atomically."""
    config_dir.mkdir(parents=True, exist_ok=True)
    target = config_dir / "config.json"
    temp = target.with_name(f"config.json.tmp.{os.getpid()}")
    try:
        temp.write_text(record["text"] + "\n")
        os.replace(temp, target)
    finally:
        temp.unlink(missing_ok=True)
    return target


def verify_config_sha(record: dict) -> None:
    """Fail closed when the shard line is not the canonical form it is hashed as."""
    text = canonical(record["config"])
    if text != record["text"]:
        raise ValueError(f"config {record['config_id']} shard line is not canonical")
    if hashlib.sha256(text.encode()).hexdigest() != record["sha256"]:
        raise ValueError(f"config {record['config_id']} sha256 mismatch")


def already_done(config_dir: Path, sha: str) -> bool:
    """Idempotent skip gate: engine marker + published parquet + matching config sha."""
    marker, parquet = config_dir / DONE_NAME, config_dir / PARQUET_NAME
    published = config_dir / SHA_MARKER
    try:
        if not (marker.is_file() and marker.stat().st_size > 0):
            return False
        if not (parquet.is_file() and parquet.stat().st_size > 0):
            return False
        return published.is_file() and published.read_text() == sha
    except OSError:
        return False


def publish_sha(config_dir: Path, sha: str) -> None:
    temp = config_dir / f"{SHA_MARKER}.tmp.{os.getpid()}"
    try:
        temp.write_text(sha)
        os.replace(temp, config_dir / SHA_MARKER)
    finally:
        temp.unlink(missing_ok=True)


def materialize(config: dict):
    """Config JSON -> engine objects, with the CLI's assembly/annotation fence."""
    if config.get("assembly") != ti.ASSEMBLY or config.get("annotation") != ti.ANNOTATION:
        raise ValueError("GRCh38 / GENCODE v50 assembly fence required")
    windows = [atlas.Window(**dict(w, sites=tuple(atlas.Site(**s) for s in w["sites"])))
               for w in config["windows"]]
    grid = atlas.Grid(**dict(config["grid"], lengths=tuple(config["grid"]["lengths"])))
    return windows, grid


def record_failure(atlas_dir: Path, task_id: int, config_id: str, status: int,
                   reason: str) -> None:
    """Explicit failure rows; one bad locus must not strand the rest of the shard.

    The TSV keeps the two columns jobs/j4_scale_array.sbatch wrote (config_id,
    exit-equivalent status); the JSONL carries the reason for forensics.
    """
    atlas_dir.mkdir(parents=True, exist_ok=True)
    with (atlas_dir / f"failures_task_{task_id}.tsv").open("a") as stream:
        stream.write(f"{config_id}\t{status}\n")
    with (atlas_dir / f"failures_task_{task_id}.jsonl").open("a") as stream:
        stream.write(canonical({"config_id": config_id, "status": status,
                                "reason": reason, "task_id": task_id}) + "\n")


def run_config(record: dict, atlas_dir: Path, genome, scorer, batch_size: int) -> dict:
    """Run one config exactly as the per-config CLI invocation would.

    reference_dir is left at the engine default (``<config_dir>/reference_scores``)
    so the done manifest cites the same artifact paths as the CLI path.
    """
    config_dir = atlas_dir / record["config_id"]
    verify_config_sha(record)
    publish_config(record, config_dir)
    windows, grid = materialize(record["config"])
    return atlas.run_atlas(windows, grid, genome, scorer, config_dir,
                           batch_size=batch_size, shard=0, shards=1)


def drive(records, atlas_dir: Path, genome, scorer, *, batch_size: int, task_id: int,
          checkpoint_path: Path, log=print) -> dict:
    """Process every config in the shard with the already-loaded scorer."""
    state = {"task_id": task_id, "batch_size": batch_size, "configs": len(records),
             "next": 0, "done": 0, "skipped": 0, "failed": 0, "failures": []}

    def flush(_signum=None):
        ck.save_checkpoint(state, checkpoint_path)

    flush()
    with ck.checkpoint_on_signal(flush, requeue=False):
        for index, record in enumerate(records):
            config_dir = atlas_dir / record["config_id"]
            if already_done(config_dir, record["sha256"]):
                state["skipped"] += 1
                state["next"] = index + 1
                # Flush on the skip path too: an all-skip resume must still leave
                # a checkpoint describing what this task actually accounted for.
                flush()
                continue
            start = time.perf_counter()
            try:
                done = run_config(record, atlas_dir, genome, scorer, batch_size)
            except Exception as exc:
                reason = f"{type(exc).__name__}: {exc}"
                record_failure(atlas_dir, task_id, record["config_id"], 1, reason)
                state["failed"] += 1
                state["failures"].append({"config_id": record["config_id"], "reason": reason})
                print(f"ERROR: config {record['config_id']} raised {reason}", file=sys.stderr)
            else:
                # The CLI exits nonzero on failed rows; the shell runner counted
                # that as a shard failure and withheld the sha marker. Match it.
                if done["failures"] > 0:
                    reason = f"complete_with_failures rows={done['failures']}"
                    record_failure(atlas_dir, task_id, record["config_id"], 1, reason)
                    state["failed"] += 1
                    state["failures"].append({"config_id": record["config_id"], "reason": reason})
                    print(f"ERROR: config {record['config_id']} {reason}", file=sys.stderr)
                else:
                    publish_sha(config_dir, record["sha256"])
                    state["done"] += 1
                log(json.dumps({k: done[k] for k in ("status", "expected_rows", "failures", "run_id")},
                               sort_keys=True))
            state["next"] = index + 1
            state["last_seconds"] = time.perf_counter() - start
            flush()
    return state


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--shard-dir", type=Path, default=DEFAULT_SHARD_DIR)
    source.add_argument("--shard-file", type=Path,
                        help="explicit shard jsonl; overrides --shard-dir/--task-id")
    parser.add_argument("--task-id", type=int, default=None,
                        help="array task id; defaults to SLURM_ARRAY_TASK_ID")
    parser.add_argument("--atlas-dir", type=Path, default=DEFAULT_ATLAS_DIR)
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument("--fasta", default=None)
    args = parser.parse_args(argv)

    task_id = args.task_id
    if task_id is None:
        raw = os.environ.get("SLURM_ARRAY_TASK_ID")
        if raw is None:
            parser.error("--task-id is required outside a Slurm array (no SLURM_ARRAY_TASK_ID)")
        task_id = int(raw)
    shard_file = (args.shard_file.expanduser() if args.shard_file
                  else shard_file_for(args.shard_dir, task_id))
    atlas_dir = args.atlas_dir.expanduser()
    fasta = args.fasta if args.fasta is not None else ti.DEFAULT_FASTA

    records = read_plan(shard_file)
    atlas_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_path = atlas_dir / ".j4_driver" / f"task_{task_id}.checkpoint.json"
    print(f"PLAN shard={shard_file.name} configs={len(records)} batch_size={args.batch_size}")

    # Publish readiness only AFTER handlers exist and BEFORE any output work; the
    # shell checks this marker after wait, avoiding a readiness/TERM race.
    with ck.checkpoint_on_signal(lambda _sig: None, requeue=False):
        if ready_file := os.environ.get("J4_READY_FILE"):
            Path(ready_file).resolve().touch(exist_ok=False)
        genome, reference_identity = atlas.grch38_reference(fasta)
        load_start = time.perf_counter()
        scorer = atlas.PangolinScorer()
        model_load_seconds = time.perf_counter() - load_start
        scorer.provenance = dict(scorer.provenance, reference=reference_identity)
        print(f"LOAD model_load_seconds={model_load_seconds:.1f} (once per task)")
        state = drive(records, atlas_dir, genome, scorer, batch_size=args.batch_size,
                      task_id=task_id, checkpoint_path=checkpoint_path)
    print(f"SHARD {task_id} done={state['done']} skipped={state['skipped']} "
          f"failed={state['failed']}")
    return int(state["failed"] > 0)


if __name__ == "__main__":
    raise SystemExit(main())
