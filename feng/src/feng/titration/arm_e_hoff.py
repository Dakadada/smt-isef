"""Checkpointed Pangolin H-off sweep over an Arm E site manifest."""

from __future__ import annotations

import argparse
from itertools import islice
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import subprocess

from feng import source_auth
from feng.gpu.checkpoint import checkpoint_on_signal, load_checkpoint, save_checkpoint
from feng.titration.arm_e_scan import read_site_manifest, validate_burden_artifact
from feng.titration.gpu_scores import (
    DEFAULT_GENOME_FASTA,
    DEFAULT_GPU_PYTHON,
    PANGOLIN_CONTEXT,
    PANGOLIN_MASKING_METHOD,
    PANGOLIN_WEIGHT_FILENAMES,
    pangolin_off_target_burden,
)


HOFF_ALGORITHM_VERSION = "armE-pangolin-hoff-v4"
HOFF_ARTIFACT_SCHEMA_VERSION = 2


class HoffIncompleteError(RuntimeError):
    """A diagnostic artifact was preserved, but the shard did not complete."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _digest(payload: dict) -> str:
    return hashlib.sha256(json.dumps(
        payload, sort_keys=True, separators=(",", ":")
    ).encode()).hexdigest()


def _source_digest_bindings() -> dict:
    """Re-authenticate and return the complete committed source identity."""
    identity = source_auth.verify_source_tree_from_env()
    return {
        "pangolin_distance": DEFAULT_DISTANCE,
        "pangolin_context": PANGOLIN_CONTEXT,
        **identity,
    }


def _execution_provenance() -> dict:
    """Identify this execution without changing reproducible input identity."""
    run_id = os.environ.get("FENG_ARME_RUN_ID", "")
    array_job_id = os.environ.get("SLURM_ARRAY_JOB_ID", "")
    array_task_id = os.environ.get("SLURM_ARRAY_TASK_ID", "")
    launcher_sha = os.environ.get("FENG_ARME_LAUNCHER_BLOB_SHA256", "")
    if not re.fullmatch(r"[0-9]{8}T[0-9]{6}Z", run_id):
        raise ValueError("FENG_ARME_RUN_ID must be a UTC basic timestamp")
    if not re.fullmatch(r"[0-9]+", array_job_id):
        raise ValueError("SLURM_ARRAY_JOB_ID must be numeric")
    if not re.fullmatch(r"[0-9]+", array_task_id):
        raise ValueError("SLURM_ARRAY_TASK_ID must be numeric")
    if not re.fullmatch(r"[0-9a-f]{64}", launcher_sha):
        raise ValueError("FENG_ARME_LAUNCHER_BLOB_SHA256 must be a lowercase SHA-256")
    return {
        "run_id": run_id,
        "slurm_array_job_id": array_job_id,
        "slurm_array_task_id": array_task_id,
        "launcher_blob_sha256": launcher_sha,
    }


def _input_digest(inputs: dict) -> str:
    return _digest(inputs)


def _authenticated_genome_inputs(genome_fasta: Path, genome_pin: Path | None) -> dict:
    if not genome_fasta.is_file() or genome_fasta.stat().st_size == 0:
        raise FileNotFoundError(f"genome FASTA missing or empty: {genome_fasta}")
    actual_genome_sha = _sha256(genome_fasta)
    if genome_pin is None:
        genome_source_sha = None
        genome_pin_sha = None
    else:
        pin_payload = json.loads(genome_pin.read_text())
        genome_source_sha = pin_payload.get("source_genome_fasta_sha256")
        if not isinstance(genome_source_sha, str) or len(genome_source_sha) != 64:
            raise ValueError("genome pin lacks source_genome_fasta_sha256")
        if actual_genome_sha != genome_source_sha:
            raise ValueError(
                "genome FASTA SHA-256 mismatch: "
                f"expected {genome_source_sha}, got {actual_genome_sha}"
            )
        genome_pin_sha = _sha256(genome_pin)
    return {
        "genome_fasta": str(genome_fasta.resolve()),
        "genome_fai_sha256": _sha256(Path(str(genome_fasta) + ".fai"))
        if Path(str(genome_fasta) + ".fai").is_file() else None,
        "genome_fasta_sha256": actual_genome_sha,
        "genome_source_sha256": genome_source_sha,
        "genome_pin_sha256": genome_pin_sha,
    }


def _pangolin_provenance(gpu_python: Path) -> dict:
    """Return installed Pangolin code, version, and weight identities."""
    gpu_python = Path(gpu_python)
    if not gpu_python.is_file():
        raise FileNotFoundError(f"GPU env python missing: {gpu_python}")
    package_identity = source_auth.probe_pangolin_package_identity(gpu_python)
    pinned_package_sha = source_auth.PINNED_PANGOLIN_PACKAGE_SHA256
    if package_identity["pangolin_package_sha256"] != pinned_package_sha:
        raise source_auth.SourceAuthenticationError(
            "live Pangolin package identity disagrees with committed pin"
        )
    script = """
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import sys

names = json.loads(sys.argv[1])
spec = importlib.util.find_spec("pangolin")
if spec is None or not spec.submodule_search_locations:
    raise RuntimeError("Pangolin package location unavailable")
model_dir = Path(next(iter(spec.submodule_search_locations))) / "models"
hashes = {}
for name in names:
    path = model_dir / name
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Pangolin weight missing or empty: {path}")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    hashes[name] = digest.hexdigest()
print(json.dumps({
    "pangolin_package_version": importlib.metadata.version("pangolin"),
    "pangolin_weight_sha256": hashes,
}, sort_keys=True))
"""
    try:
        completed = subprocess.run(
            [str(gpu_python), "-c", script, json.dumps(PANGOLIN_WEIGHT_FILENAMES)],
            check=False, capture_output=True, text=True, timeout=180,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError(f"Pangolin provenance probe failed: {exc}") from exc
    if completed.returncode != 0:
        detail = completed.stderr.strip().splitlines()[-1] if completed.stderr.strip() else "no stderr"
        raise RuntimeError(f"Pangolin provenance probe failed: {detail}")
    try:
        provenance = json.loads(completed.stdout.strip().splitlines()[-1])
    except (IndexError, json.JSONDecodeError) as exc:
        raise RuntimeError("Pangolin provenance probe returned malformed output") from exc
    version = provenance.get("pangolin_package_version")
    hashes = provenance.get("pangolin_weight_sha256")
    expected_names = set(PANGOLIN_WEIGHT_FILENAMES)
    if type(version) is not str or not version.strip():
        raise ValueError("Pangolin package version missing from provenance")
    if type(hashes) is not dict or set(hashes) != expected_names:
        raise ValueError("Pangolin weight provenance is incomplete")
    if any(
        type(value) is not str or len(value) != 64
        or any(char not in "0123456789abcdef" for char in value)
        for value in hashes.values()
    ):
        raise ValueError("Pangolin weight SHA-256 is malformed")
    return {
        **provenance,
        **package_identity,
        "pangolin_package_pinned_sha256": pinned_package_sha,
    }


def _hoff_inputs(*, manifest_sha: str, scan_burden: Path,
                 scan_payload: dict, genome_fasta: Path,
                 genome_pin: Path | None, gpu_python: Path,
                 chunk_size: int) -> dict:
    source_inputs = _source_digest_bindings()
    return {
        "algorithm_version": HOFF_ALGORITHM_VERSION,
        "manifest_sha256": manifest_sha,
        "scan_burden_sha256": _sha256(scan_burden),
        "length": scan_payload["length"],
        "aso_5to3": scan_payload["aso_5to3"],
        **_authenticated_genome_inputs(genome_fasta, genome_pin),
        **_pangolin_provenance(gpu_python),
        **source_inputs,
        "gpu_python": str(gpu_python.resolve()),
        "masking_method": PANGOLIN_MASKING_METHOD,
        "chunk_size": chunk_size,
    }


def validate_output(path: Path, manifest: Path, *,
                    genome_fasta: Path = DEFAULT_GENOME_FASTA,
                    gpu_python: Path = DEFAULT_GPU_PYTHON,
                    genome_pin: Path | None = None,
                    scan_burden: Path | None = None,
                    chunk_size: int = 32) -> dict:
    """Authenticate a successful H-off artifact before creating ``.done``."""
    if not path.is_file() or not manifest.is_file():
        raise FileNotFoundError("H-off output or site manifest missing")
    payload = json.loads(path.read_text())
    if payload.get("schema_version") != HOFF_ARTIFACT_SCHEMA_VERSION:
        raise ValueError("H-off artifact schema mismatch")
    if payload.get("algorithm_version") != HOFF_ALGORITHM_VERSION:
        raise ValueError("H-off artifact algorithm version mismatch")
    if payload.get("state") != "candidate":
        raise ValueError("H-off artifact must remain in candidate state")
    if payload.get("terminal") is not True:
        raise ValueError("H-off artifact is not terminal")
    digest = payload.get("input_digest")
    if not isinstance(digest, str) or len(digest) != 64 or any(
        char not in "0123456789abcdef" for char in digest
    ):
        raise ValueError("H-off artifact input/config digest missing or malformed")
    if payload.get("manifest_sha256") != _sha256(manifest):
        raise ValueError("H-off manifest SHA-256 mismatch")
    if scan_burden is None:
        raise ValueError("scan burden is required to authenticate H-off input")
    scan_payload = validate_burden_artifact(scan_burden, manifest)
    expected_inputs = _hoff_inputs(
        manifest_sha=_sha256(manifest), scan_burden=scan_burden,
        scan_payload=scan_payload, genome_fasta=genome_fasta,
        genome_pin=genome_pin, gpu_python=gpu_python,
        chunk_size=chunk_size,
    )
    expected_digest = _input_digest(expected_inputs)
    if payload["input_digest"] != expected_digest:
        raise ValueError("completed H-off input/config digest mismatch")
    if payload.get("inputs") != expected_inputs:
        raise ValueError("completed H-off input provenance mismatch")
    if payload.get("execution_provenance") != _execution_provenance():
        raise ValueError("completed H-off execution provenance mismatch")
    try:
        rows = list(read_site_manifest(manifest))
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("H-off manifest schema invalid") from exc
    incomplete = payload.get("pangolin_off_target_incomplete")
    if type(incomplete) is not bool:
        raise ValueError("H-off incomplete flag must be a boolean")
    n_sites = payload.get("n_sites")
    if type(n_sites) is not int:
        raise ValueError("H-off site count must be an exact integer")
    if n_sites != len(rows):
        raise ValueError("H-off site count does not match canonical manifest")
    # An empty canonical manifest has two authenticated meanings: a completed
    # scan with zero hits, or a deliberately uninformative short-rung scan.
    if not rows:
        expected_incomplete = scan_payload.get("complete") is False
        if incomplete is not expected_incomplete:
            raise ValueError("empty H-off state disagrees with authenticated scan completeness")
        if expected_incomplete and payload.get("reason") != scan_payload.get("reason"):
            raise ValueError("empty H-off reason disagrees with authenticated scan")
    # Runtime/fetch/orientation failures are terminal=False and rejected above.
    if incomplete:
        if not isinstance(payload.get("reason"), str) or not payload["reason"]:
            raise ValueError("H-off artifact is incomplete")
        if rows or payload.get("pangolin_off_target_burden") is not None:
            raise ValueError("H-off empty-manifest diagnostic is inconsistent")
    else:
        burden = payload.get("pangolin_off_target_burden")
        if type(burden) not in {int, float} or not math.isfinite(burden) or burden < 0:
            raise ValueError("H-off burden is invalid")
        if not rows and burden != 0:
            raise ValueError("authenticated complete-zero scan must have zero H-off burden")
    return payload


def run_sweep(*, manifest: Path, out: Path, checkpoint: Path,
              genome_fasta: Path = DEFAULT_GENOME_FASTA,
              gpu_python: Path = DEFAULT_GPU_PYTHON, genome_pin: Path | None = None,
              scan_burden: Path | None = None,
              chunk_size: int = 32) -> dict:
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    if not manifest.is_file() or manifest.stat().st_size == 0:
        raise FileNotFoundError(f"site manifest missing or empty: {manifest}")
    if scan_burden is None:
        raise ValueError("scan burden is required to authenticate H-off input")
    manifest_sha = _sha256(manifest)
    scan_payload = validate_burden_artifact(scan_burden, manifest)
    inputs = _hoff_inputs(
        manifest_sha=manifest_sha, scan_burden=scan_burden,
        scan_payload=scan_payload, genome_fasta=genome_fasta,
        genome_pin=genome_pin, gpu_python=gpu_python,
        chunk_size=chunk_size,
    )
    input_digest = _input_digest(inputs)
    execution_provenance = _execution_provenance()
    first_iter = read_site_manifest(manifest)
    first = next(first_iter, None)
    if first is None:
        match = re.search(r"sites_L(\d+)\.tsv$", manifest.name)
        length = int(match.group(1)) if match else None
        scan_complete_zero = bool(
            scan_payload is not None
            and scan_payload.get("complete") is True
            and scan_payload.get("off_target_hits") == 0
        )
        if scan_payload is not None:
            length = scan_payload["length"]
        aso = scan_payload.get("aso_5to3") if scan_payload else None
        empty_reason = (
            "completed scan found zero off-target sites"
            if scan_complete_zero else (
                scan_payload.get("reason") if scan_payload is not None
                else "empty site manifest without authenticated scan reason"
            )
        )
        payload = {"schema_version": HOFF_ARTIFACT_SCHEMA_VERSION,
                   "algorithm_version": HOFF_ALGORITHM_VERSION,
                   "length": length, "aso_5to3": aso,
                   "pangolin_off_target_burden": 0.0 if scan_complete_zero else None,
                   "pangolin_off_target_incomplete": not scan_complete_zero,
                   "reason": empty_reason,
                   "n_sites": 0,
                   "manifest_sha256": manifest_sha, "terminal": True,
                   "input_digest": input_digest,
                   "inputs": inputs,
                   "execution_provenance": execution_provenance,
                   "masking_method": PANGOLIN_MASKING_METHOD,
                   "state": "candidate"}
        _atomic_json(payload, out)
        return payload
    length, aso = first["length"], first["aso_5to3"]
    if out.is_file():
        existing = json.loads(out.read_text())
        # A failed diagnostic is intentionally preserved for audit, but it is
        # not a completion cache and must never prevent a checkpointed retry.
        if existing.get("terminal") is True:
            prior = validate_output(
                out, manifest, genome_fasta=genome_fasta,
                gpu_python=gpu_python, genome_pin=genome_pin,
                scan_burden=scan_burden, chunk_size=chunk_size,
            )
            if prior.get("input_digest") != input_digest:
                raise ValueError("completed H-off input/config digest mismatch")
            return prior
    loaded_state = load_checkpoint(checkpoint)
    if loaded_state is None:
        state = {}
        start = 0
        burden = 0.0
    else:
        if type(loaded_state) is not dict:
            raise ValueError("checkpoint state must be an object; refusing stale resume")
        state = loaded_state
        if state.get("input_digest") != input_digest:
            raise ValueError("checkpoint input/config digest mismatch; refusing stale resume")
        prior_burden = state.get("pangolin_off_target_burden")
        if (type(prior_burden) is not float or not math.isfinite(prior_burden)
                or prior_burden < 0):
            raise ValueError(
                "checkpoint pangolin_off_target_burden missing or invalid; refusing stale resume"
            )
        if "site_index" not in state or type(state["site_index"]) is not int:
            raise ValueError(
                "checkpoint site_index missing or invalid; refusing stale resume"
            )
        start = state["site_index"]
        burden = prior_burden
    if start < 0 or not math.isfinite(burden) or burden < 0:
        raise ValueError("checkpoint site_index outside manifest")
    committed = [{"site_index": start, "pangolin_off_target_burden": burden,
                  "input_digest": input_digest, "inputs": inputs}]

    def checkpoint_committed(_signum=None):
        save_checkpoint(dict(committed[0]), checkpoint)

    incomplete = None
    with checkpoint_on_signal(checkpoint_committed, signals=(signal.SIGTERM,), requeue=False):
        rows = read_site_manifest(manifest)
        for _ in range(start):
            if next(rows, None) is None:
                raise ValueError("checkpoint site_index outside manifest")
        while True:
            batch = list(islice(rows, chunk_size))
            if not batch:
                break
            if any(row["length"] != length or row["aso_5to3"] != aso for row in batch):
                raise ValueError("site manifest must contain exactly one length and ASO")
            result = pangolin_off_target_burden(aso, batch,
                                                genome_fasta=genome_fasta,
                                                gpu_python=gpu_python)
            if result["pangolin_off_target_burden"] is None:
                incomplete = result["reason"]
                break
            batch_burden = float(result["pangolin_off_target_burden"])
            if not math.isfinite(batch_burden) or batch_burden < 0:
                incomplete = "worker returned invalid burden"
                break
            next_snapshot = {
                "site_index": committed[0]["site_index"] + len(batch),
                "pangolin_off_target_burden": committed[0]["pangolin_off_target_burden"] + batch_burden,
                "input_digest": input_digest,
                "inputs": inputs,
            }
            save_checkpoint(next_snapshot, checkpoint)
            committed[0] = next_snapshot
    if incomplete is None:
        payload = {"schema_version": HOFF_ARTIFACT_SCHEMA_VERSION,
                   "algorithm_version": HOFF_ALGORITHM_VERSION,
                   "length": length, "aso_5to3": aso,
                   "pangolin_off_target_burden": committed[0]["pangolin_off_target_burden"],
                   "pangolin_off_target_incomplete": False, "reason": None,
                   "n_sites": committed[0]["site_index"], "terminal": True,
                   "manifest_sha256": manifest_sha, "input_digest": input_digest,
                   "inputs": inputs,
                   "execution_provenance": execution_provenance,
                   "masking_method": PANGOLIN_MASKING_METHOD,
                   "state": "candidate"}
    else:
        payload = {"schema_version": HOFF_ARTIFACT_SCHEMA_VERSION,
                   "algorithm_version": HOFF_ALGORITHM_VERSION,
                   "length": length, "aso_5to3": aso,
                   "pangolin_off_target_burden": None,
                   "pangolin_off_target_incomplete": True,
                   "reason": incomplete or "sweep did not reach terminal site",
                   "n_sites_completed": committed[0]["site_index"],
                   "manifest_sha256": manifest_sha, "input_digest": input_digest,
                   "inputs": inputs,
                   "execution_provenance": execution_provenance,
                   "masking_method": PANGOLIN_MASKING_METHOD,
                   "terminal": False, "state": "candidate"}
    _atomic_json(payload, out)
    if incomplete is not None:
        raise HoffIncompleteError(payload["reason"])
    return payload


def _atomic_json(payload: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".partial.{os.getpid()}")
    save_checkpoint(payload, tmp)
    os.replace(tmp, path)


def _main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--genome-fasta", type=Path, default=DEFAULT_GENOME_FASTA)
    parser.add_argument("--gpu-python", type=Path, default=DEFAULT_GPU_PYTHON)
    parser.add_argument("--genome-pin", type=Path)
    parser.add_argument("--scan-burden", type=Path)
    parser.add_argument("--chunk-size", type=int, default=32)
    parser.add_argument("--validate-output", action="store_true")
    args = parser.parse_args()
    if args.validate_output:
        validate_output(
            args.out, args.manifest, genome_fasta=args.genome_fasta,
            gpu_python=args.gpu_python, genome_pin=args.genome_pin,
            scan_burden=args.scan_burden, chunk_size=args.chunk_size,
        )
        return 0
    run_sweep(manifest=args.manifest, out=args.out, checkpoint=args.checkpoint,
              genome_fasta=args.genome_fasta, gpu_python=args.gpu_python,
              genome_pin=args.genome_pin, scan_burden=args.scan_burden,
              chunk_size=args.chunk_size)
    return 0


if __name__ == "__main__":
    raise SystemExit(_main())
