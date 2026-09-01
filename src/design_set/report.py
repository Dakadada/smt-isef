"""Honesty-gated renderers and immutable timestamped design-set deposits."""
from __future__ import annotations

import json
from pathlib import Path

from design_engine import PROSPECTIVE_BANNER, PROSPECTIVE_STATUS

from .models import ProspectiveDesignSet


def _validate_render(report: ProspectiveDesignSet) -> None:
    report.validate_for_rendering()
    if report.banner != PROSPECTIVE_BANNER or report.prediction_status != PROSPECTIVE_STATUS:
        raise ValueError("renderer refuses a non-prospective design set")


def render_json(report: ProspectiveDesignSet) -> str:
    _validate_render(report)
    return json.dumps(report.to_dict(), indent=2, sort_keys=True, ensure_ascii=False) + "\n"


def render_text(report: ProspectiveDesignSet) -> str:
    _validate_render(report)
    reference = report.reference
    lines = [report.banner, f"Status: {report.prediction_status}", report.claim_boundary, "",
             f"Frozen design-set timestamp: {report.generated_at}",
             f"Target: {report.target_id} ({report.locus})",
             f"Coordinates: {report.coordinate_convention}",
             f"Off-target disclosure: {report.provenance['off_target_disclosure']}",
             f"Reference: {reference.label} | {reference.coordinate_status}"]
    if reference.candidate_id is not None:
        lines.append(f"Reference candidate: {reference.candidate_id} | "
                     f"{reference.effector_class} | {reference.chemistry} | "
                     f"start={reference.start} end={reference.end}")
    else:
        lines.append(f"Reference limitation: {reference.limitation}")
    lines.extend(["", "PREDICTED-SUPERIOR TO REFERENCE:",
                  "candidate_id | effector_class | chemistry | start | end | "
                  "directional_confidence | off_target_burden | span_robustness | rank | pareto"])
    superior = [row for row in report.candidates if row.predicted_superior_to_reference]
    if not superior:
        lines.append("(none)")
    for row in superior:
        lines.append(f"{row.candidate_id} | {row.effector_class} | {row.chemistry} | "
                     f"{row.start} | {row.end} | {row.directional_confidence:.6f} | "
                     f"{row.off_target_burden:.6f} | {row.span_robustness:.6f} | "
                     f"{row.rank} | {'yes' if row.pareto else 'no'}")
    lines.extend(["", "ALL ENGINE-PREDICTED MOLECULES:",
                  "candidate_id | effector_class | chemistry | start | end | target_sequence | "
                  "effector_sequence | directional_confidence | off_target_burden | "
                  "span_robustness | rank | pareto | reference"])
    for row in report.candidates:
        rank = "—" if row.rank is None else str(row.rank)
        lines.append(f"{row.candidate_id} | {row.effector_class} | {row.chemistry} | "
                     f"{row.start} | {row.end} | {row.target_sequence} | {row.effector_sequence} | "
                     f"{row.directional_confidence:.6f} | {row.off_target_burden:.6f} | "
                     f"{row.span_robustness:.6f} | {rank} | "
                     f"{'yes' if row.pareto else 'no'} | {'yes' if row.is_reference else 'no'}")
    lines.extend(["", "Provenance:", json.dumps(report.to_dict()["provenance"],
                                                  sort_keys=True, ensure_ascii=False)])
    return "\n".join(lines) + "\n"


def deposit(report: ProspectiveDesignSet, directory: Path) -> tuple[Path, Path]:
    """Write timestamped JSON/text exactly once; refuse either-name collisions."""
    _validate_render(report)
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    stamp = report.generated_at.replace(":", "").replace("+00:00", "Z")
    json_path = directory / f"prospective-designs-{stamp}.json"
    text_path = directory / f"prospective-designs-{stamp}.txt"
    if json_path.exists() or text_path.exists():
        raise FileExistsError(f"refusing to overwrite prospective design-set deposit {stamp}")
    with json_path.open("x", encoding="utf-8") as handle:
        handle.write(render_json(report))
    try:
        with text_path.open("x", encoding="utf-8") as handle:
            handle.write(render_text(report))
    except Exception:
        json_path.unlink()
        raise
    return json_path, text_path
