import os
import shutil
import subprocess
import time
import urllib.parse
import urllib.request
from pathlib import Path

from .common import Checkpoints, file_lock, md5, read_tsv, require_space, validate_config


def _accession_lock(scratch, accession):
    """Serialize fetch and quantification work for one accession."""
    lock_path = Path(scratch) / "fastq" / ".locks" / f"{accession}.lock"
    return file_lock(lock_path)


def _urls(row):
    urls = row["fastq_ftp"].split(";"); sums = row["fastq_md5"].split(";")
    sizes = row.get("fastq_bytes", "").split(";")
    if len(urls) != len(sums):
        raise ValueError(f"{row['run_accession']}: FASTQ URL/MD5 cardinality mismatch")
    return list(zip(urls, sums, sizes + [""] * (len(urls) - len(sizes))))


def run(config_path, manifest=None):
    cfg = validate_config(config_path); paths = cfg["paths"]
    manifest = manifest or str(Path(paths["work_dir"]) / "manifests/all_runs.frozen.tsv")
    cp = Checkpoints(paths["checkpoints"], "fetch", config_path,cfg); qcp = Checkpoints(paths["checkpoints"], "quantify", config_path,cfg)
    scratch = Path(paths["scratch_dir"])
    for row in read_tsv(manifest):
        if row["endpoint_eligible"] != "true": continue
        with _accession_lock(scratch, row["run_accession"]):
            _run_accession(cfg, manifest, cp, qcp, scratch, row)


def _run_accession(cfg, manifest, cp, qcp, scratch, row):
    paths = cfg["paths"]
    run_dir = scratch / "fastq" / row["run_accession"]
    outputs = [run_dir / Path(url).name for url, _, _ in _urls(row)]
    if cp.valid(row["run_accession"], outputs,[manifest],row):
        # Reconcile debris from an older interrupted writer even when the final
        # archives and checkpoint are already valid.
        for debris in run_dir.glob("*.partial*"):
            debris.unlink()
        return
    count_output = Path(paths["work_dir"]) / "counts" / (row["run_accession"] + ".tsv")
    # FASTQs intentionally disappear after a committed quantification; do not fetch them again on resume.
    if qcp.valid(row["run_accession"], [count_output],[manifest],row):
        shutil.rmtree(run_dir,ignore_errors=True); shutil.rmtree(scratch/"alignment"/row["run_accession"],ignore_errors=True); return
    expected_bytes = sum(int(n or 0) for _, _, n in _urls(row))
    require_space(scratch,max(int(expected_bytes*cfg["fetch"]["space_multiplier"]),cfg["fetch"]["minimum_free_bytes"]),cfg["fetch"]["aggregate_scratch_cap_bytes"])
    run_dir.mkdir(parents=True, exist_ok=True)
    for (url, expected, _), dest in zip(_urls(row), outputs):
        if dest.exists() and md5(dest) == expected: continue
        partial = dest.with_suffix(dest.suffix + ".partial")
        normalized = url if "://" in url else "https://" + url
        if normalized.startswith("file://"):
            shutil.copy2(urllib.request.url2pathname(urllib.parse.urlparse(normalized).path), partial)
        else:
            # Preserve .partial and resume byte ranges after preemption.
            for attempt in range(8):
                curl_process = subprocess.run(["curl", "--fail", "--location", "--retry", "8",
                                               "--retry-delay", "10", "--retry-max-time", "600", "--connect-timeout", "30",
                                               "--speed-limit", "1000", "--speed-time", "120", "--continue-at", "-",
                                               "--output", str(partial), normalized], check=False)
                if curl_process.returncode == 0: break
                if attempt == 7:
                    raise IOError(f"{row['run_accession']}: curl failed for {normalized} with exit code {curl_process.returncode}")
                time.sleep(min(10 * 2 ** attempt, 120))
        if md5(partial) != expected:
            partial.rename(partial.with_suffix(partial.suffix + ".corrupt")); raise IOError(f"MD5 mismatch: {dest}")
        os.replace(partial, dest)
    # Debris is useful for resuming failures, but never after all final archives validate.
    for debris in run_dir.glob("*.partial*"):
        debris.unlink()
    cp.complete(row["run_accession"],outputs,inputs=[manifest],resolved_metadata=row,archive_md5=[x[1] for x in _urls(row)])
