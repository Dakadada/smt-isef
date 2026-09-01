"""Guard the clean-room Feng arm against legacy code dependencies."""

from pathlib import Path

import pytest


SOURCE_ROOT = Path(__file__).resolve().parents[1] / "src" / "feng"
FORBIDDEN_REFERENCES = (
    "src/",
    "pre/src/",
    "design_engine",
    "offtarget",
    "physics_adapter",
    "design_panel",
)


@pytest.mark.parametrize("source_file", sorted(SOURCE_ROOT.rglob("*.py")))
def test_implementation_has_no_legacy_reference(source_file: Path) -> None:
    text = source_file.read_text(encoding="utf-8")
    hits = [token for token in FORBIDDEN_REFERENCES if token in text]
    assert not hits, f"{source_file.relative_to(SOURCE_ROOT)} references {hits}"

