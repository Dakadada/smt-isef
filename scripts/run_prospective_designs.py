#!/usr/bin/env python3
"""Deposit the committed-context SMN2 prospective design set."""
import argparse
from datetime import datetime
from pathlib import Path

from design_set import build_design_set, deposit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timestamp", required=True, help="timezone-aware ISO-8601 timestamp")
    parser.add_argument("--out", type=Path, default=Path("artifacts/prospective-designs"))
    args = parser.parse_args()
    timestamp = datetime.fromisoformat(args.timestamp.replace("Z", "+00:00"))
    report = build_design_set(clock=lambda: timestamp)
    for path in deposit(report, args.out):
        print(path)


if __name__ == "__main__":
    main()
