"""Guard sealed terminal results from pre-freeze arm implementation code."""

from __future__ import annotations

import ast
from pathlib import Path
import subprocess

import pytest


FENG_ROOT = Path(__file__).resolve().parents[1]
GUARDED_ROOTS = tuple(FENG_ROOT / "src" / "feng" / name for name in ("field", "rbp", "gate", "bench"))
def _freeze_tag_exists() -> bool:
    result = subprocess.run(
        ["git", "tag", "--list", "feng-freeze", "feng-freeze-*"],
        cwd=FENG_ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return bool(result.stdout.strip())


def _terminal_read_violations(path: Path) -> list[int]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    docstring_node = None
    if tree.body and isinstance(tree.body[0], ast.Expr):
        docstring_node = tree.body[0].value
    lines: list[int] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            if any("terminal" in alias.name.casefold() for alias in node.names):
                lines.append(node.lineno)
        elif isinstance(node, ast.ImportFrom):
            names = [node.module or "", *(alias.name for alias in node.names)]
            if any("terminal" in name.casefold() for name in names):
                lines.append(node.lineno)
        elif (
            isinstance(node, ast.Constant)
            and isinstance(node.value, str)
            and "terminal" in node.value.casefold()
            and node is not docstring_node
        ):
            lines.append(node.lineno)
    return lines


@pytest.mark.skipif(_freeze_tag_exists(), reason="terminal results unlocked by freeze tag")
@pytest.mark.parametrize(
    "source_file",
    sorted(path for root in GUARDED_ROOTS for path in root.rglob("*.py")),
)
def test_guarded_code_cannot_read_terminal_before_freeze(source_file: Path) -> None:
    lines = _terminal_read_violations(source_file)
    assert not lines, f"pre-freeze terminal access in {source_file.relative_to(FENG_ROOT)} at {lines}"
