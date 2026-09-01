#!/usr/bin/env python3
"""Deposit the checked-in real-context prospective panel."""
import argparse
from datetime import datetime
from pathlib import Path

from design_panel import deposit, run_panel


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--timestamp", required=True, help="timezone-aware ISO-8601 timestamp")
    parser.add_argument("--out", type=Path, default=Path("artifacts/design-panels"))
    args = parser.parse_args()
    timestamp = datetime.fromisoformat(args.timestamp.replace("Z", "+00:00"))
    report = run_panel(clock=lambda: timestamp)
    if not any(item.target_id == "SMN2_exon7_inclusion" for item in report.designs):
        raise RuntimeError("checked-in panel emitted no real SMN2 design")
    for path in deposit(report, args.out): print(path)


if __name__ == "__main__": main()
