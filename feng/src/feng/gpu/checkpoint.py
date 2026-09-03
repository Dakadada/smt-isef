"""Preemption-safe checkpoint / requeue helpers for feng GPU jobs.

``mit_preemptable`` (the only partition with real GPU headroom for this user)
can evict a job at any time; SLURM sends ``SIGTERM`` ~60 s ahead when the sbatch
declares ``--signal=B:SIGTERM@60``. These helpers let a long GPU sweep persist
its per-shard progress on a 5-15 min cadence and, on the eviction warning, flush
one final checkpoint and ask SLURM to requeue THIS job so it resumes from the
last checkpoint instead of restarting from scratch.

Serialization is JSON, deliberately: checkpoint state here is shard-progress
bookkeeping (which shards are ``.done``, last-index, timestamps) — plain data.
JSON is human-inspectable, diffable, and carries no arbitrary-code / version
fragility that ``pickle`` would; anything not JSON-serializable does not belong
in a checkpoint. Writes are atomic (tmp in the same dir + ``os.replace``), so a
crash mid-write never corrupts or truncates the existing checkpoint.

No RNG; deterministic given its inputs.
"""

from __future__ import annotations

import json
import os
import signal
import subprocess
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Callable, Iterator, Optional


def save_checkpoint(state: dict, path: os.PathLike | str) -> Path:
    """Atomically write ``state`` (JSON) to ``path``.

    Writes to a sibling temp file then ``os.replace``s it over the target, so a
    reader never sees a partial file and an existing checkpoint survives a
    failure during serialization. Returns the checkpoint path.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    try:
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(state, fh, sort_keys=True, indent=2)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    finally:
        # Never leave a stray tmp file behind if json.dump raised mid-write.
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:  # pragma: no cover - best-effort cleanup
                pass
    return path


def load_checkpoint(path: os.PathLike | str) -> Optional[dict]:
    """Load a checkpoint, or ``None`` if it does not exist.

    A present-but-corrupt checkpoint raises ``json.JSONDecodeError`` (a real
    problem the caller must see) rather than being silently treated as absent.
    """
    path = Path(path)
    if not path.exists():
        return None
    with open(path, "r", encoding="utf-8") as fh:
        return json.load(fh)


def should_checkpoint(
    last_checkpoint_ts: Optional[float],
    interval_s: float,
    *,
    now: Optional[float] = None,
) -> bool:
    """Cadence gate: True if ``interval_s`` has elapsed since the last flush.

    ``last_checkpoint_ts=None`` (never checkpointed) returns True so the first
    checkpoint fires promptly. ``now`` is injectable for deterministic tests.
    """
    if interval_s <= 0:
        raise ValueError(f"interval_s must be positive, got {interval_s}")
    if last_checkpoint_ts is None:
        return True
    current = time.time() if now is None else now
    return (current - last_checkpoint_ts) >= interval_s


def requeue_self(*, job_id: Optional[str] = None, dry_run: bool = False) -> bool:
    """Ask SLURM to requeue the current job so it resumes after preemption.

    Resolves the job id from ``SLURM_JOB_ID`` when not given. Returns True if a
    requeue was issued (or would be, under ``dry_run``), False when there is no
    job id (e.g. running outside SLURM) or the ``scontrol`` call fails. Never
    raises — a failed requeue must not mask the checkpoint that preceded it.
    """
    jid = job_id or os.environ.get("SLURM_JOB_ID")
    if not jid:
        return False
    if dry_run:
        return True
    try:
        subprocess.run(["scontrol", "requeue", jid], check=True,
                       capture_output=True, text=True)
        return True
    except (subprocess.CalledProcessError, OSError):
        return False


@contextmanager
def checkpoint_on_signal(
    on_signal: Callable[[int], Any],
    *,
    signals: tuple[int, ...] = (signal.SIGTERM, signal.SIGUSR1),
    requeue: bool = True,
    exit_after: bool = True,
) -> Iterator[None]:
    """Install handlers that flush a checkpoint on the preemption warning.

    On any of ``signals`` (default SIGTERM + SIGUSR1, matching
    ``--signal=B:SIGTERM@60``), ``on_signal(signum)`` is called — the caller's
    final-flush callback — then, if ``requeue``, :func:`requeue_self` is invoked
    and (if ``exit_after``) the process exits 0 so SLURM's requeue takes effect
    cleanly. Prior handlers are restored on exit. Must be entered from the main
    thread (Python signal constraint).

    Under the ``gpu_preemptable.sbatch`` template, ``B:SIGTERM`` reaches only the
    batch shell, whose trap forwards SIGTERM to this payload and owns the SINGLE
    ``scontrol requeue``. In that wiring pass ``requeue=False`` here (flush-only)
    so requeue does not double-fire; ``exit_after=True`` still lets the payload
    exit cleanly once its final checkpoint is written.
    """
    previous: dict[int, Any] = {}

    def _handler(signum, _frame):
        try:
            on_signal(signum)
        finally:
            if requeue:
                requeue_self()
            if exit_after:
                # 0 + a pending requeue = clean resume; nonzero could be read as
                # a real failure by afterok dependents.
                os._exit(0)

    try:
        for sig in signals:
            previous[sig] = signal.signal(sig, _handler)
        yield
    finally:
        for sig, prev in previous.items():
            try:
                signal.signal(sig, prev)
            except (OSError, ValueError, TypeError):  # pragma: no cover
                pass
