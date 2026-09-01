"""Test-path setup for the clean-room Feng arm source checkout."""

from pathlib import Path
import sys


FENG_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FENG_ROOT / "src"))

