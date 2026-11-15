#!/usr/bin/env python
"""J4 completion gate: every emitted config must have a valid .done under the atlas.

Reads j4_configs_manifest.tsv (the builder's own ledger), takes every row with
status=emitted, and checks the config's output tree under <atlas-dir>:

  <config_id>/shard_0.done                 engine done manifest, nonempty
  <config_id>/ablation_maps/shard_0.parquet published map, nonempty
  <config_id>/.j4_config_sha256            == the manifest config_sha256

The sha marker is published by the runner only after a clean run, so a config
that failed, or that was built from a different config text, shows up as a gap.
Exit 1 and report the gaps; exit 0 only when the atlas is complete.

Default mode is deliberately cheap on the shared filesystem: stat + one tiny
read per config, no parquet or manifest parsing. ``--deep`` additionally loads
each done manifest and checks its integrity hash, published row count and
status (that is 20k manifest reads; run it once, not in a loop).

  python conference/model_v2/jobs/j4_scale_gate.py \
      --manifest ~/orcd/scratch/feng/model_v2/j4_configs/j4_configs_manifest.tsv \
      --atlas-dir ~/orcd/scratch/feng/model_v2/j4_atlas_v1
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from conference.model_v2 import gpu_perturbation_atlas as atlas
from conference.model_v2 import j4_build_configs as jb
from feng.gpu import checkpoint as ck

DEFAULT_MANIFEST = Path("~/orcd/scratch/feng/model_v2/j4_configs/j4_configs_manifest.tsv")
DEFAULT_ATLAS_DIR = Path("~/orcd/scratch/feng/model_v2/j4_atlas_v1")
DONE_NAME = "shard_0.done"
PARQUET_NAME = "ablation_maps/shard_0.parquet"
SHA_MARKER = ".j4_config_sha256"


def read_manifest(path) -> list[dict]:
    """Emitted rows of j4_configs_manifest.tsv, fail-closed on shape drift."""
    columns = list(jb.MANIFEST_COLUMNS)
    rows: list[dict] = []
    header_seen = False
    for number, line in enumerate(Path(path).expanduser().read_text().splitlines(), start=1):
        if line.startswith("#"):
            fields = line.lstrip("#").split("\t")
            if fields[0] == columns[0]:
                if fields != columns:
                    raise ValueError(f"{path}:{number}: unexpected manifest columns")
                header_seen = True
            continue
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) != len(columns):
            raise ValueError(f"{path}:{number}: expected {len(columns)} fields, got {len(fields)}")
        row = dict(zip(columns, fields))
        if row["status"] == "emitted":
            rows.append(row)
    if not header_seen:
        raise ValueError(f"{path}: no column header")
    if not rows:
        raise ValueError(f"{path}: no emitted configs")
    return rows


def check_config(atlas_dir: Path, row: dict, deep: bool = False) -> str | None:
    """Return None when complete, else a short gap reason."""
    directory = atlas_dir / row["config_id"]
    if not directory.is_dir():
        return "missing_config_dir"
    marker = directory / DONE_NAME
    if not marker.is_file() or marker.stat().st_size == 0:
        return "missing_done"
    parquet = directory / PARQUET_NAME
    if not parquet.is_file() or parquet.stat().st_size == 0:
        return "missing_parquet"
    published = directory / SHA_MARKER
    if not published.is_file():
        return "missing_config_sha"
    if published.read_text().strip() != row["config_sha256"]:
        return "config_sha_mismatch"
    if deep:
        try:
            done = atlas.load_manifest(marker)
        except (ValueError, json.JSONDecodeError) as exc:
            return f"done_manifest_invalid:{type(exc).__name__}"
        if done is None:
            return "missing_done"
        if done.get("status") != "complete":
            return f"status_{done.get('status')}"
        expected = int(row["rows"])
        if done.get("expected_rows") != expected or done["ablation"]["rows"] != expected:
            return "row_count_mismatch"
        if done.get("failures"):
            return "failed_rows"
    return None


def gate(manifest=DEFAULT_MANIFEST, atlas_dir=DEFAULT_ATLAS_DIR, *, deep: bool = False) -> dict:
    atlas_dir = Path(atlas_dir).expanduser()
    rows = read_manifest(manifest)
    gaps = []
    reasons: dict[str, int] = {}
    for row in rows:
        reason = check_config(atlas_dir, row, deep=deep)
        if reason is not None:
            gaps.append({"config_id": row["config_id"], "shard": row["shard"], "reason": reason})
            reasons[reason] = reasons.get(reason, 0) + 1
    return {"manifest": str(Path(manifest).expanduser()), "atlas_dir": str(atlas_dir),
            "deep": deep, "emitted": len(rows), "complete": len(rows) - len(gaps),
            "gaps": len(gaps), "gap_reasons": reasons, "gap_rows": gaps}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--atlas-dir", type=Path, default=DEFAULT_ATLAS_DIR)
    parser.add_argument("--deep", action="store_true",
                        help="also verify done-manifest integrity, status and row counts")
    parser.add_argument("--report", type=Path, default=None, help="write the full report as JSON")
    parser.add_argument("--list", type=int, default=20, help="gap rows to print (0 = all)")
    args = parser.parse_args(argv)

    report = gate(args.manifest, args.atlas_dir, deep=args.deep)
    if args.report:
        ck.save_checkpoint(report, args.report.expanduser())
    print(json.dumps({k: report[k] for k in
                      ("emitted", "complete", "gaps", "gap_reasons", "deep")}, sort_keys=True))
    shown = report["gap_rows"] if args.list == 0 else report["gap_rows"][:args.list]
    for gap in shown:
        print(f"GAP shard={gap['shard']} {gap['config_id']} {gap['reason']}", file=sys.stderr)
    if report["gaps"] > len(shown):
        print(f"GAP ... {report['gaps'] - len(shown)} more", file=sys.stderr)
    return int(report["gaps"] > 0)


if __name__ == "__main__":
    raise SystemExit(main())
