#!/usr/bin/env python3
"""Generate the complete Arm A descriptive atlas from the frozen CSV inputs."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


FENG_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FENG_ROOT / "src"))

from feng.atlas.arm_a import OUTPUT_NAMES, build_arm_a_atlas, lint_scientific_outputs


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the descriptive SMN2 field atlas")
    parser.add_argument(
        "--guides", type=Path,
        default=Path("/home/dashii/orcd/scratch/feng/phase0/splicerush_dcas13d_field.csv"),
    )
    parser.add_argument(
        "--asos", type=Path,
        default=FENG_ROOT / "artifacts" / "phase0_arm_a_overlap" / "harmonized_aso_positions.csv",
    )
    parser.add_argument("--output", type=Path, default=FENG_ROOT / "artifacts" / "arm_a_descriptive")
    parser.add_argument("--lint-only", action="store_true")
    args = parser.parse_args()
    if args.lint_only:
        missing = [name for name in OUTPUT_NAMES if not (args.output / name).is_file()]
        violations = missing or lint_scientific_outputs(args.output)
        if violations:
            print("FAIL: " + "; ".join(violations))
            return 1
        print("PASS: forbidden-string scan found no violations outside the fixed claim-boundary section")
        return 0
    data = build_arm_a_atlas(args.guides, args.asos, args.output)
    print(f"Wrote {len(OUTPUT_NAMES)} outputs; {len(data.guides)} guides, {len(data.windows)} ASO records, {len(data.blocks)} block-study rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
