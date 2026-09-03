"""Tests for the shared feng GPU infrastructure (feng.gpu).

Infra-only: these exercise module resolution, graceful env inspection, atomic
checkpointing, the cadence gate, and SLURM-context parsing. Nothing here runs a
module, installs a package, or imports torch/openmm.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import textwrap
import time
import types
from pathlib import Path

import pytest

from feng.gpu import env as env_mod
from feng.gpu import checkpoint as ckpt_mod
from feng.gpu.env import (
    DEFAULT_CUDA,
    DEFAULT_CUDNN,
    DEFAULT_GROMACS,
    check_gpu_env,
    gpu_env_python,
    resolve_modules,
    slurm_gpu_context,
)
from feng.gpu.checkpoint import (
    load_checkpoint,
    requeue_self,
    save_checkpoint,
    should_checkpoint,
)


# --------------------------------------------------------------------------- #
# env.resolve_modules
# --------------------------------------------------------------------------- #
def test_resolve_modules_defaults():
    lines = resolve_modules()
    assert f"module load cuda/{DEFAULT_CUDA}" in lines
    assert f"module load cudnn/{DEFAULT_CUDNN}" in lines
    assert f"module load gromacs/{DEFAULT_GROMACS}" in lines
    assert any(l.startswith("module load miniforge/") for l in lines)


def test_resolve_modules_default_cuda_is_cuda12_line():
    # cudnn 9.8 is a cuda-12 build; the default cuda must be the 12.x line.
    assert DEFAULT_CUDA.startswith("12.")


def test_resolve_modules_gromacs_none_omits_md_module():
    lines = resolve_modules(gromacs=None)
    assert not any("gromacs" in l for l in lines)
    # cuda/cudnn/miniforge still present
    assert any("cuda/" in l for l in lines)


def test_resolve_modules_rejects_unknown_cuda():
    with pytest.raises(ValueError, match="cuda"):
        resolve_modules(cuda="99.9.9")


def test_resolve_modules_rejects_unknown_gromacs():
    with pytest.raises(ValueError, match="gromacs"):
        resolve_modules(gromacs="1999.9")


def test_resolve_modules_rejects_incompatible_cuda_cudnn(monkeypatch):
    # MED-2: cuda 13.x with the cuda12-built cudnn is in-allowlist per module but
    # a mutually incompatible combo — must raise, not emit a node-only failure.
    monkeypatch.setattr(env_mod, "CUDA_VERSIONS",
                        frozenset(env_mod.CUDA_VERSIONS | {"13.1.0"}))
    with pytest.raises(ValueError, match="incompatible|CUDA 12"):
        resolve_modules(cuda="13.1.0", cudnn="9.8.0.87-cuda12")


def test_resolve_modules_default_combo_is_compatible():
    # the pinned default (cuda 12.9 + cuda12 cudnn) must pass the cross-check
    lines = resolve_modules()
    assert any("cuda/12." in l for l in lines)


def test_cudnn_cuda_major_parses_tag():
    assert env_mod._cudnn_cuda_major("9.8.0.87-cuda12") == "12"
    assert env_mod._cudnn_cuda_major("9.9.0-cuda13") == "13"
    assert env_mod._cudnn_cuda_major("no-tag-here") is None


# --------------------------------------------------------------------------- #
# env.check_gpu_env  (graceful degradation)
# --------------------------------------------------------------------------- #
def test_check_gpu_env_absent_reports_all_missing(tmp_path):
    missing_prefix = str(tmp_path / "no-such-env")
    report = check_gpu_env(missing_prefix)
    assert report["env_exists"] is False
    assert report["available"] is False
    assert set(report["packages"]) == {"torch", "spliceai", "pangolin", "openmm"}
    assert all(v is False for v in report["packages"].values())
    # never raised, structured report returned
    assert report["python"].endswith("bin/python")


def test_check_gpu_env_never_installs_or_raises(tmp_path):
    # Even pointed at a bogus path, it must return cleanly.
    report = check_gpu_env(str(tmp_path / "ghost"), timeout_s=1.0)
    assert isinstance(report, dict)
    assert report["available"] is False


def test_probe_import_no_indexerror_on_whitespace_nonzero(monkeypatch):
    # MED-3: a probe with whitespace-only stderr + nonzero exit must not
    # IndexError — it should fall back to the exit-code message.
    fake = types.SimpleNamespace(returncode=1, stderr="   \n  ", stdout="")
    monkeypatch.setattr(env_mod.subprocess, "run", lambda *a, **k: fake)
    ok, msg = env_mod._probe_import(Path("/nonexistent/python"), "torch", 1.0)
    assert ok is False
    assert msg == "exit 1"


def test_check_gpu_env_never_raises_even_if_probe_raises(tmp_path, monkeypatch):
    # MED-3: even if a probe blows up, check_gpu_env must return a report.
    real_py = tmp_path / "bin" / "python"
    real_py.parent.mkdir(parents=True)
    real_py.write_text("#!/bin/sh\n")

    def boom(*a, **k):
        raise RuntimeError("probe exploded")

    monkeypatch.setattr(env_mod, "_probe_import", boom)
    report = check_gpu_env(str(tmp_path))  # env "exists" so it reaches the probe loop
    assert report["env_exists"] is True
    assert report["available"] is False
    assert all(v is False for v in report["packages"].values())


def test_gpu_env_python_path_shape():
    py = gpu_env_python("~/orcd/pool/feng-gpu-env")
    assert py.name == "python"
    assert py.parent.name == "bin"
    assert "feng-gpu-env" in str(py)
    assert "~" not in str(py)  # expanded


# --------------------------------------------------------------------------- #
# env.slurm_gpu_context
# --------------------------------------------------------------------------- #
def test_slurm_context_parses_env(monkeypatch):
    monkeypatch.setenv("SLURM_JOB_PARTITION", "mit_preemptable")
    monkeypatch.setenv("SLURM_JOB_ID", "12345")
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0,1")
    ctx = slurm_gpu_context()
    assert ctx["partition"] == "mit_preemptable"
    assert ctx["job_id"] == "12345"
    assert ctx["n_gpus"] == 2
    assert ctx["under_slurm"] is True


def test_slurm_context_no_slurm(monkeypatch):
    for var in ("SLURM_JOB_PARTITION", "SLURM_JOB_ID", "CUDA_VISIBLE_DEVICES",
                "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID", "SLURM_JOB_GRES"):
        monkeypatch.delenv(var, raising=False)
    ctx = slurm_gpu_context()
    assert ctx["under_slurm"] is False
    assert ctx["n_gpus"] == 0
    assert ctx["partition"] is None


def test_slurm_context_empty_cuda_visible(monkeypatch):
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "")
    assert slurm_gpu_context()["n_gpus"] == 0


# --------------------------------------------------------------------------- #
# checkpoint save/load round-trip + atomicity
# --------------------------------------------------------------------------- #
def test_checkpoint_round_trip(tmp_path):
    path = tmp_path / "ckpt.json"
    state = {"done_shards": [1, 3, 5], "last_index": 5, "ts": 1000.0}
    save_checkpoint(state, path)
    assert load_checkpoint(path) == state


def test_load_checkpoint_absent_returns_none(tmp_path):
    assert load_checkpoint(tmp_path / "nope.json") is None


def test_checkpoint_atomic_no_partial_on_failure(tmp_path, monkeypatch):
    path = tmp_path / "ckpt.json"
    save_checkpoint({"v": 1}, path)  # good baseline
    assert load_checkpoint(path) == {"v": 1}

    # Force json.dump to blow up mid-write; the target must be untouched and no
    # stray .tmp file left behind.
    def boom(*a, **k):
        raise RuntimeError("disk full")

    monkeypatch.setattr(ckpt_mod.json, "dump", boom)
    with pytest.raises(RuntimeError):
        save_checkpoint({"v": 2}, path)
    assert load_checkpoint(path) == {"v": 1}  # old checkpoint survived
    leftovers = list(tmp_path.glob("ckpt.json.tmp*"))
    assert leftovers == []


def test_checkpoint_creates_parent_dirs(tmp_path):
    path = tmp_path / "nested" / "deep" / "ckpt.json"
    save_checkpoint({"x": 1}, path)
    assert load_checkpoint(path) == {"x": 1}


# --------------------------------------------------------------------------- #
# should_checkpoint cadence gate
# --------------------------------------------------------------------------- #
def test_should_checkpoint_first_time_true():
    assert should_checkpoint(None, 600) is True


def test_should_checkpoint_before_interval_false():
    assert should_checkpoint(1000.0, 600, now=1300.0) is False


def test_should_checkpoint_after_interval_true():
    assert should_checkpoint(1000.0, 600, now=1600.0) is True
    assert should_checkpoint(1000.0, 600, now=1601.0) is True


def test_should_checkpoint_rejects_nonpositive_interval():
    with pytest.raises(ValueError):
        should_checkpoint(1000.0, 0)


# --------------------------------------------------------------------------- #
# requeue_self
# --------------------------------------------------------------------------- #
def test_requeue_self_no_job_id_returns_false(monkeypatch):
    monkeypatch.delenv("SLURM_JOB_ID", raising=False)
    assert requeue_self() is False


def test_requeue_self_dry_run_with_job(monkeypatch):
    monkeypatch.setenv("SLURM_JOB_ID", "999")
    assert requeue_self(dry_run=True) is True


def test_requeue_self_dry_run_explicit_job():
    assert requeue_self(job_id="777", dry_run=True) is True


# --------------------------------------------------------------------------- #
# checkpoint_on_signal restores prior handlers
# --------------------------------------------------------------------------- #
def test_checkpoint_on_signal_restores_handlers():
    original = signal.getsignal(signal.SIGUSR1)
    calls = []
    # exit_after/requeue off so the handler is inert if it ever fired.
    with ckpt_mod.checkpoint_on_signal(
        lambda s: calls.append(s),
        signals=(signal.SIGUSR1,),
        requeue=False,
        exit_after=False,
    ):
        assert signal.getsignal(signal.SIGUSR1) is not original
    assert signal.getsignal(signal.SIGUSR1) is original


def test_checkpoint_on_signal_invokes_callback(monkeypatch):
    calls = []
    with ckpt_mod.checkpoint_on_signal(
        lambda s: calls.append(s),
        signals=(signal.SIGUSR1,),
        requeue=False,
        exit_after=False,
    ):
        os.kill(os.getpid(), signal.SIGUSR1)
    assert calls == [signal.SIGUSR1]


def test_sigterm_flushes_final_checkpoint_before_exit(tmp_path):
    """MED-1: a REAL forwarded SIGTERM must flush a final checkpoint then exit.

    Mirrors the gpu_preemptable.sbatch wiring: the payload installs
    checkpoint_on_signal(requeue=False, exit_after=True) and, on SIGTERM
    (forwarded by the batch-shell trap), writes its final checkpoint before the
    process dies. Spawns a real subprocess so the signal path is genuinely
    exercised, not called in-process.
    """
    ckpt = tmp_path / "shard.ckpt.json"
    ready = tmp_path / "ready"
    # src root = .../feng/src  (checkpoint.py is feng/src/feng/gpu/checkpoint.py)
    src_root = Path(ckpt_mod.__file__).resolve().parents[3]
    script = textwrap.dedent(f"""
        import signal, time
        from pathlib import Path
        from feng.gpu.checkpoint import checkpoint_on_signal, save_checkpoint

        CKPT = Path({str(ckpt)!r})
        READY = Path({str(ready)!r})

        def final_flush(signum):
            save_checkpoint({{"final": True, "signum": signum, "last_index": 42}}, CKPT)

        with checkpoint_on_signal(final_flush, signals=(signal.SIGTERM,),
                                  requeue=False, exit_after=True):
            READY.write_text("up")
            for _ in range(200):   # ~20s max; parent kills well before
                time.sleep(0.1)
    """)
    env = dict(os.environ)
    env["PYTHONPATH"] = str(src_root) + os.pathsep + env.get("PYTHONPATH", "")
    env.pop("SLURM_JOB_ID", None)  # requeue_self is a no-op off SLURM anyway
    proc = subprocess.Popen([sys.executable, "-c", script], env=env)
    try:
        for _ in range(100):  # wait up to 10s for handler-installed readiness
            if ready.exists():
                break
            time.sleep(0.1)
        assert ready.exists(), "child never signalled ready"
        proc.send_signal(signal.SIGTERM)
        proc.wait(timeout=10)
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait(timeout=5)
    assert proc.returncode == 0, f"child exited {proc.returncode}, expected clean 0"
    assert ckpt.exists(), "final checkpoint was not flushed on SIGTERM"
    state = json.loads(ckpt.read_text())
    assert state["final"] is True and state["last_index"] == 42
