"""Shared GPU-enablement infrastructure for the feng arms.

This package holds *infrastructure only* — module/env resolution and
preemption-safe checkpoint helpers. It never installs software, never launches
a GPU job, and never imports a heavy GPU library at module load. The actual
GPU compute env lives at ``~/orcd/pool/feng-gpu-env`` and is created by the
documented (not auto-run) ``feng/scripts/setup_gpu_env.sh``.

Reusable by all arms (Arm E MD / deep-splice sweeps first).
"""

from feng.gpu.env import (
    DEFAULT_CUDA,
    DEFAULT_CUDNN,
    DEFAULT_GROMACS,
    GPU_ENV_PREFIX,
    check_gpu_env,
    gpu_env_python,
    resolve_modules,
    slurm_gpu_context,
)
from feng.gpu.checkpoint import (
    checkpoint_on_signal,
    load_checkpoint,
    requeue_self,
    save_checkpoint,
    should_checkpoint,
)

__all__ = [
    "DEFAULT_CUDA",
    "DEFAULT_CUDNN",
    "DEFAULT_GROMACS",
    "GPU_ENV_PREFIX",
    "check_gpu_env",
    "gpu_env_python",
    "resolve_modules",
    "slurm_gpu_context",
    "checkpoint_on_signal",
    "load_checkpoint",
    "requeue_self",
    "save_checkpoint",
    "should_checkpoint",
]
