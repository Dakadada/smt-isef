"""Inert command-line foundation for the clean-room Feng arm."""

from __future__ import annotations

import argparse
from collections.abc import Sequence


PROGRAM_DESCRIPTION = (
    "Clean-room Feng arm foundation for bounded candidate benchmarks; "
    "it makes no clinical safety, efficacy, or patient-level claims."
)


def build_parser() -> argparse.ArgumentParser:
    """Build the foundation parser without registering pipeline operations."""
    return argparse.ArgumentParser(prog="feng", description=PROGRAM_DESCRIPTION)


def main(argv: Sequence[str] | None = None) -> int:
    """Display the inert CLI surface and exit successfully."""
    parser = build_parser()
    parser.parse_args(argv)
    parser.print_help()
    return 0

