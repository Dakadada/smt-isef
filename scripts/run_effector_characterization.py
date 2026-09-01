#!/usr/bin/env python3
"""Deposit the authenticated SMN2 engine-internal score characterization."""
import argparse
from datetime import datetime
from pathlib import Path

from effector_characterization import build_characterization, deposit


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timestamp", required=True,
                        help="timezone-aware ISO-8601 timestamp")
    parser.add_argument("--out", type=Path,
                        default=Path("artifacts/effector-characterization"))
    args = parser.parse_args()
    timestamp = datetime.fromisoformat(args.timestamp.replace("Z", "+00:00"))
    report = build_characterization(clock=lambda: timestamp)
    for path in deposit(report, args.out):
        print(path)


if __name__ == "__main__":
    main()
