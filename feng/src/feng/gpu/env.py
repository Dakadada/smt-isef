"""Resolve and validate the feng GPU compute environment.

Pure resolution/inspection: this module emits ``module load`` lines and reports
which GPU libraries are importable in the (separately created) GPU conda env. It
NEVER runs a module, installs a package, or imports torch/openmm itself — every
function here is safe to call from the CPU ``feng-env`` and degrades gracefully
when the GPU env is absent.

Confirmed LMOD modules on this cluster (2026-09-03, read-only ``module avail``):
  cuda:      12.9.1, 13.0.1, 13.1.0
  cudnn:     9.8.0.87-cuda12
  gromacs:   2021.6, 2024.1, 2026.2
  miniforge: 23.11.0-0, 24.3.0-0, 25.11.0-0
No AlphaFold/Rosetta module exists.
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

# --- confirmed-available module versions (allowlists; reject anything else so a
# typo can't silently produce a `module load` line that fails only on the node).
CUDA_VERSIONS = frozenset({"12.9.1", "13.0.1", "13.1.0"})
CUDNN_VERSIONS = frozenset({"9.8.0.87-cuda12"})
GROMACS_VERSIONS = frozenset({"2021.6", "2024.1", "2026.2"})
MINIFORGE_VERSIONS = frozenset({"23.11.0-0", "24.3.0-0", "25.11.0-0"})

# --- pinned known-good default combo. cudnn 9.8 is built for cuda 12, so the
# default cuda is the 12.x line (12.9.1), NOT the 13.x modules — pairing cudnn9.8
# with cuda13 would be an unsupported mix. gromacs 2024.1 is the stable current
# release (2026.2 is newest but 2024.1 is the conservative default).
DEFAULT_CUDA = "12.9.1"
DEFAULT_CUDNN = "9.8.0.87-cuda12"
DEFAULT_GROMACS = "2024.1"
DEFAULT_MINIFORGE = "25.11.0-0"

# --- the GPU conda env path (created by setup_gpu_env.sh, not by this module).
GPU_ENV_PREFIX = "~/orcd/pool/feng-gpu-env"

# --- GPU libraries the deep-splice / MD scorers will need.
GPU_PACKAGES = ("torch", "spliceai", "pangolin", "openmm")


def resolve_modules(
    *,
    cuda: str = DEFAULT_CUDA,
    cudnn: str = DEFAULT_CUDNN,
    gromacs: Optional[str] = DEFAULT_GROMACS,
    miniforge: str = DEFAULT_MINIFORGE,
) -> list[str]:
    """Return the ``module load`` lines for a validated module combo.

    Every requested version must be in its confirmed-available allowlist; an
    unknown version raises ``ValueError`` rather than emitting a line that would
    fail only once the job is on a GPU node. ``gromacs=None`` omits the MD
    module (deep-splice-only jobs don't need it).
    """
    _require("cuda", cuda, CUDA_VERSIONS)
    _require("cudnn", cudnn, CUDNN_VERSIONS)
    _require("miniforge", miniforge, MINIFORGE_VERSIONS)
    _require_cuda_cudnn_compatible(cuda, cudnn)
    lines = [
        f"module load cuda/{cuda}",
        f"module load cudnn/{cudnn}",
        f"module load miniforge/{miniforge}",
    ]
    if gromacs is not None:
        _require("gromacs", gromacs, GROMACS_VERSIONS)
        # gromacs after cuda so its GPU build resolves against the loaded toolkit
        lines.insert(2, f"module load gromacs/{gromacs}")
    return lines


def _require(name: str, value: str, allowed: frozenset[str]) -> None:
    if value not in allowed:
        raise ValueError(
            f"unknown {name} module version {value!r}; "
            f"available: {sorted(allowed)}"
        )


def _cudnn_cuda_major(cudnn: str) -> Optional[str]:
    """The CUDA major a cudnn module is built for, from its ``-cuda<major>`` tag.

    e.g. ``9.8.0.87-cuda12`` -> ``"12"``. Returns None if the module name carries
    no such tag (then no cross-constraint can be enforced).
    """
    marker = "-cuda"
    idx = cudnn.rfind(marker)
    if idx == -1:
        return None
    tag = cudnn[idx + len(marker):]
    # keep only the leading digits (e.g. "12" from "12" or "12-foo")
    digits = ""
    for ch in tag:
        if ch.isdigit():
            digits += ch
        else:
            break
    return digits or None


def _require_cuda_cudnn_compatible(cuda: str, cudnn: str) -> None:
    """A ``-cuda<major>`` cudnn requires a cuda line of the same major.

    Both versions are individually in-allowlist by the time this runs; this
    catches an in-allowlist but mutually incompatible pairing (e.g. cuda 13.x
    with the cuda12-built cudnn) that would otherwise fail only on the node.
    Keyed generically off the cuda major, not a hardcoded pair.
    """
    needed_major = _cudnn_cuda_major(cudnn)
    if needed_major is None:
        return
    cuda_major = cuda.split(".", 1)[0]
    if cuda_major != needed_major:
        raise ValueError(
            f"incompatible module combo: cudnn {cudnn!r} is built for CUDA "
            f"{needed_major}.x but cuda {cuda!r} is a {cuda_major}.x line; "
            f"pick a cuda {needed_major}.x module (available: {sorted(CUDA_VERSIONS)})"
        )


def gpu_env_python(prefix: str = GPU_ENV_PREFIX) -> Path:
    """Expected python interpreter path inside the GPU conda env.

    Path-only; does not check existence (see :func:`check_gpu_env`).
    """
    return Path(os.path.expanduser(prefix)) / "bin" / "python"


def check_gpu_env(
    prefix: str = GPU_ENV_PREFIX,
    *,
    packages: tuple[str, ...] = GPU_PACKAGES,
    timeout_s: float = 120.0,
) -> dict:
    """Structured report on the GPU env: does it exist, which libs import.

    ``timeout_s`` is generous (torch's cold import is ~15-25s and slower under
    node CPU load); a probe that times out is reported False with a ``timeout:``
    error string -- inconclusive, not proof of absence -- so callers can tell a
    slow import from a genuinely missing package. This report is ADVISORY: GPU
    scorers do their real import inside their own worker and treat a genuine
    ImportError there as the authoritative "unavailable", never this probe.

    NEVER raises, NEVER installs. If the env python is absent, every package is
    reported missing. If present, each package is probed by a short subprocess
    ``import`` in that interpreter (isolated from this CPU process). Used to
    gate GPU scorers before they run — a missing env yields a clean
    "unavailable" report, not a crash.
    """
    py = gpu_env_python(prefix)
    report = {
        "prefix": os.path.expanduser(prefix),
        "python": str(py),
        "env_exists": py.exists(),
        "packages": {pkg: False for pkg in packages},
        "available": False,
        "errors": {},
    }
    if not py.exists():
        return report
    for pkg in packages:
        try:
            ok, err = _probe_import(py, pkg, timeout_s)
        except Exception as exc:  # never let a probe break the gating contract
            ok, err = False, f"{type(exc).__name__}: {exc}"
        report["packages"][pkg] = ok
        if err:
            report["errors"][pkg] = err
    report["available"] = all(report["packages"].values())
    return report


def _probe_import(py: Path, module: str, timeout_s: float) -> tuple[bool, str]:
    """Try ``python -c 'import <module>'`` in the GPU env; never raise."""
    try:
        proc = subprocess.run(
            [str(py), "-c", f"import {module}"],
            capture_output=True,
            text=True,
            timeout=timeout_s,
        )
    except (subprocess.TimeoutExpired, OSError) as exc:  # pragma: no cover - env-dependent
        return False, f"{type(exc).__name__}: {exc}"
    if proc.returncode == 0:
        return True, ""
    # Last non-empty line of stderr/stdout, or the exit code if the process said
    # nothing (or only whitespace) — must never IndexError on empty output.
    lines = (proc.stderr or proc.stdout or "").strip().splitlines()
    return False, lines[-1] if lines else f"exit {proc.returncode}"


def slurm_gpu_context() -> dict:
    """Read the current SLURM/GPU context from environment variables (logging).

    Returns a plain dict; empty/None fields when not running under SLURM.
    ``n_gpus`` counts the entries in ``CUDA_VISIBLE_DEVICES``.
    """
    cuda_visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    n_gpus = 0
    if cuda_visible:
        n_gpus = len([d for d in cuda_visible.split(",") if d.strip() != ""])
    return {
        "partition": os.environ.get("SLURM_JOB_PARTITION"),
        "job_id": os.environ.get("SLURM_JOB_ID"),
        "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
        "array_task_id": os.environ.get("SLURM_ARRAY_TASK_ID"),
        "cuda_visible_devices": cuda_visible,
        "n_gpus": n_gpus,
        "gres": os.environ.get("SLURM_JOB_GRES"),
        "under_slurm": os.environ.get("SLURM_JOB_ID") is not None,
    }


def module_cmd_available() -> bool:
    """Whether an LMOD ``module`` shim / ``lmod`` is reachable (advisory)."""
    return bool(shutil.which("lmod") or os.environ.get("LMOD_CMD"))
