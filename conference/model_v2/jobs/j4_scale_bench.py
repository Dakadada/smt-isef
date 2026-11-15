"""J4 scale-out GPU throughput bench driver (real Pangolin; no stubs).

Runs real configs from the J4 config shards through
conference.model_v2.gpu_perturbation_atlas.run_atlas at several --batch-size
values and reports inferences/s and peak GPU memory per sweep, plus the
one-off model-load cost (which the array runner pays once per engine
invocation). Measured numbers choose the array shard size and --time.

Each sweep gets its own output dir AND its own reference_scores dir, so no
sweep is credited with another sweep's cached work and no sweep short-circuits
on an existing .done marker.

Inference accounting: one Pangolin ensemble call per N-masked footprint plus one
per window reference score => inferences = footprints + windows.
Masked rows keep evidence_class sequence_ablation_proxy (set by the engine).
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import time

from conference.model_v2 import gpu_perturbation_atlas as atlas
from conference.model_v2 import j4_build_configs as jb
from conference.model_v2 import trans_intervention as ti
from feng.gpu import checkpoint as ck


def load_configs(shard_dir: Path, count: int, skip_curated: bool = True) -> list[dict]:
    configs: list[dict] = []
    for path in sorted(shard_dir.glob("shard_*.jsonl")):
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            config = json.loads(line)
            if skip_curated and config.get("config_id") == jb.CURATED_ID:
                continue  # stride-1 fine grid is not representative of the atlas grid
            configs.append(config)
            if len(configs) >= count:
                return configs
    if not configs:
        raise ValueError(f"no configs under {shard_dir}")
    return configs


def materialize(config: dict):
    windows = [atlas.Window(**dict(w, sites=tuple(atlas.Site(**s) for s in w["sites"])))
               for w in config["windows"]]
    grid = atlas.Grid(**dict(config["grid"], lengths=tuple(config["grid"]["lengths"])))
    return windows, grid


def peak_memory() -> dict:
    try:
        import torch
        if not torch.cuda.is_available():
            return {"available": False}
        return {"available": True,
                "peak_allocated_gib": torch.cuda.max_memory_allocated() / 2 ** 30,
                "peak_reserved_gib": torch.cuda.max_memory_reserved() / 2 ** 30}
    except Exception as exc:  # pragma: no cover - bench only
        return {"available": False, "error": repr(exc)}


def reset_peak_memory() -> None:
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.reset_peak_memory_stats()
    except Exception:  # pragma: no cover - bench only
        pass


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--shard-dir", type=Path,
                        default=Path("~/orcd/scratch/feng/model_v2/j4_configs/shards"))
    parser.add_argument("--out-dir", type=Path,
                        default=Path("~/orcd/scratch/feng/model_v2/j4_scale_bench"))
    parser.add_argument("--configs", type=int, default=50)
    parser.add_argument("--batch-sizes", default="4,16,32")
    parser.add_argument("--budget-seconds", type=float, default=300.0,
                        help="per-sweep wall budget; the sweep stops early and reports partials")
    parser.add_argument("--fasta", default=ti.DEFAULT_FASTA)
    args = parser.parse_args(argv)

    shard_dir = args.shard_dir.expanduser()
    out_root = args.out_dir.expanduser()
    out_root.mkdir(parents=True, exist_ok=True)
    batch_sizes = [int(n) for n in args.batch_sizes.split(",") if n.strip()]
    configs = load_configs(shard_dir, args.configs)
    planned = [materialize(c) for c in configs]
    footprints_per_config = [
        jb.footprint_count(g.start, g.end, g.lengths, g.stride) for _, g in planned]

    with ck.checkpoint_on_signal(lambda _sig: None, requeue=False):
        if ready := os.environ.get("J4_READY_FILE"):
            Path(ready).resolve().touch(exist_ok=False)
        genome, reference_identity = atlas.grch38_reference(args.fasta)
        load_start = time.perf_counter()
        scorer = atlas.PangolinScorer()
        model_load_seconds = time.perf_counter() - load_start
        scorer.provenance = dict(scorer.provenance, reference=reference_identity)

        sweeps = []
        for batch_size in batch_sizes:
            sweep_dir = out_root / f"bs{batch_size}"
            reset_peak_memory()
            done_configs = footprints = inferences = 0
            start = time.perf_counter()
            for (windows, grid), count, config in zip(planned, footprints_per_config, configs):
                atlas.run_atlas(windows, grid, genome, scorer,
                                sweep_dir / config["config_id"], batch_size=batch_size,
                                reference_dir=sweep_dir / "reference_scores")
                done_configs += 1
                footprints += count
                inferences += count + len(windows)
                if time.perf_counter() - start >= args.budget_seconds:
                    break
            elapsed = time.perf_counter() - start
            sweeps.append({
                "batch_size": batch_size, "configs": done_configs,
                "footprints": footprints, "inferences": inferences,
                "wall_seconds": elapsed,
                "inferences_per_second": inferences / elapsed if elapsed else None,
                "footprints_per_second": footprints / elapsed if elapsed else None,
                "seconds_per_config": elapsed / done_configs if done_configs else None,
                "gpu_memory": peak_memory(),
                "budget_exhausted": done_configs < len(planned),
            })

        report = {
            "shard_dir": str(shard_dir), "configs_requested": args.configs,
            "configs_loaded": len(configs),
            "footprints_per_config": sorted(set(footprints_per_config)),
            "model_load_seconds": model_load_seconds,
            "provenance": scorer.provenance, "sweeps": sweeps,
        }
        ck.save_checkpoint(report, out_root / "j4_scale_bench.json")
        print(json.dumps(report, sort_keys=True, indent=2, default=str))
        for sweep in sweeps:
            print(f"BENCH batch_size={sweep['batch_size']} configs={sweep['configs']} "
                  f"inferences={sweep['inferences']} wall_s={sweep['wall_seconds']:.1f} "
                  f"inferences_per_s={sweep['inferences_per_second']:.2f} "
                  f"peak_gpu={sweep['gpu_memory'].get('peak_allocated_gib')} "
                  f"reserved_gpu={sweep['gpu_memory'].get('peak_reserved_gib')}")
        print(f"BENCH model_load_seconds={model_load_seconds:.1f} "
              f"(paid once per engine invocation by the array runner)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
