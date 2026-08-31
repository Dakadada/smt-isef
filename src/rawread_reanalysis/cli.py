import argparse
from pathlib import Path

from . import compute, fetch, harmonize, quantify, resolve
from .common import file_lock, read_tsv, validate_config, write_tsv


FAILURE_FIELDS = ("run_accession", "stage", "error")


def _update_failure(path, accession, stage=None, error=None):
    path = Path(path)
    with file_lock(path.with_name(f".{path.name}.lock")):
        failures = {row["run_accession"]: row for row in read_tsv(path)} if path.exists() else {}
        if stage is None:
            failures.pop(accession, None)
        else:
            failures[accession] = {"run_accession": accession, "stage": stage,
                                   "error": str(error).replace("\t", " ").replace("\r", " ").replace("\n", " ")}
        if failures:
            write_tsv(path, [failures[key] for key in sorted(failures)], FAILURE_FIELDS)
        elif path.exists():
            path.unlink()


def mark_stage_a_completion(cfg, run_dir):
    run_dir = Path(run_dir); run_dir.mkdir(parents=True, exist_ok=True)
    failures = Path(cfg["paths"]["work_dir"]) / "fetch_failures.tsv"
    partial = run_dir / "STAGE_A_PARTIAL"; complete = run_dir / "STAGE_A_COMPLETE"
    if failures.exists() and read_tsv(failures):
        complete.unlink(missing_ok=True); partial.touch(); return partial
    partial.unlink(missing_ok=True); complete.touch(); return complete


def main(argv=None):
    parser = argparse.ArgumentParser(description="Uniform raw-read splicing reanalysis")
    parser.add_argument("stage", choices=["resolve", "fetch", "quantify", "process", "compute", "harmonize", "mark-completion", "all"])
    parser.add_argument("--config", required=True); parser.add_argument("--metadata-dir"); parser.add_argument("--manifest"); parser.add_argument("--run-dir")
    args = parser.parse_args(argv)
    cfg=validate_config(args.config)
    if args.stage in ("resolve", "all"): args.manifest=str(resolve.run(args.config,args.metadata_dir))
    if args.stage in ("process", "all"):
        # Sequential accession processing bounds live FASTQ/alignment storage.
        manifest=Path(args.manifest); rows=read_tsv(manifest)
        failures=Path(cfg["paths"]["work_dir"])/"fetch_failures.tsv"
        for row in rows:
            if row["endpoint_eligible"]!="true": continue
            shard=Path(cfg["paths"]["work_dir"])/"manifests"/f".{row['run_accession']}.tsv"
            write_tsv(shard,[row],rows[0].keys())
            stage="fetch"
            try:
                fetch.run(args.config,str(shard)); stage="quantify"; quantify.run(args.config,str(shard))
            except Exception as exc:
                _update_failure(failures,row["run_accession"],stage,exc)
                continue
            _update_failure(failures,row["run_accession"])
        quantify.rebuild_combined(cfg,manifest)
    else:
        if args.stage == "fetch": fetch.run(args.config,args.manifest)
        if args.stage == "quantify": quantify.run(args.config,args.manifest)
    if args.stage in ("compute", "all"): compute.run(args.config)
    if args.stage in ("harmonize", "all"): harmonize.run(args.config)
    if args.stage == "mark-completion":
        if not args.run_dir: parser.error("mark-completion requires --run-dir")
        mark_stage_a_completion(cfg,args.run_dir)


if __name__ == "__main__": main()
