"""Verify the inert command-line surface for the clean-room Feng arm."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


FENG_ROOT = Path(__file__).resolve().parents[1]


def _run_cli(*arguments: str) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(FENG_ROOT / "src")
    return subprocess.run(
        [sys.executable, "-m", "feng", *arguments],
        cwd=FENG_ROOT.parent,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )


def test_module_help_exits_successfully_and_prints_usage() -> None:
    result = _run_cli("--help")
    assert result.returncode == 0
    assert "usage: feng" in result.stdout
    assert "Clean-room Feng arm foundation" in result.stdout


def test_no_arguments_only_prints_help_and_exits_successfully() -> None:
    result = _run_cli()
    assert result.returncode == 0
    assert "usage: feng" in result.stdout
    assert result.stderr == ""

