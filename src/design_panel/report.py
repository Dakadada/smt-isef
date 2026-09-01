"""Stable JSON/text renderers and no-overwrite timestamped deposits."""
from __future__ import annotations

import json
from pathlib import Path

from .models import PanelReport


# Independently defined at each output trust boundary; do not derive from caller data.
RENDER_HONESTY_BANNER = (
    "PROSPECTIVE PREDICTION — NOT EXPERIMENTALLY VALIDATED. "
    "Candidates require independent experimental and safety validation."
)
RENDER_PROSPECTIVE_STATUS = "prospective-not-experimentally-validated"
DEPOSIT_HONESTY_BANNER = (
    "PROSPECTIVE PREDICTION — NOT EXPERIMENTALLY VALIDATED. "
    "Candidates require independent experimental and safety validation."
)
DEPOSIT_PROSPECTIVE_STATUS = "prospective-not-experimentally-validated"


def _validate_render(report: PanelReport) -> None:
    report.validate_for_rendering()
    if report.banner != RENDER_HONESTY_BANNER or report.prediction_status != RENDER_PROSPECTIVE_STATUS:
        raise ValueError("renderer refuses a non-prospective panel")


def render_json(report: PanelReport) -> str:
    _validate_render(report)
    return json.dumps(report.to_dict(), indent=2, sort_keys=True, ensure_ascii=False) + "\n"


def render_text(report: PanelReport) -> str:
    _validate_render(report)
    lines = [report.banner, "", f"Frozen panel timestamp: {report.generated_at}",
             "SMN2 is the only locus intended for eventual real benchmark support; all other targets are exploratory.",
             "", "Designability map (prospective predictions; not a validation ranking)",
             "target | goal | state | verified | best candidates | best rank | confidence | best-pick off-target | abstention | Pareto | status"]
    for row in report.designability_map:
        confidence = "—" if row.best_directional_confidence is None else f"{row.best_directional_confidence:.2f}"
        burden = "—" if row.best_pick_off_target_burden is None else f"{row.best_pick_off_target_burden:.3f}"
        abstention = "—" if row.abstention_rate is None else f"{row.abstention_rate:.1%}"
        rank = "—" if row.best_candidate_rank is None else str(row.best_candidate_rank)
        candidates = ",".join(row.best_candidate_ids) or "—"
        lines.append(f"{row.target_id} | {row.therapeutic_goal.value} | {row.state.value} | "
                     f"{'yes' if row.coordinates_verified else 'no'} | {candidates} | "
                     f"{rank} | {confidence} | {burden} | {abstention} | "
                     f"{row.pareto_candidate_count} | {row.status}")
        lines.append(f"  coordinate provenance: {row.coordinate_provenance}")
    lines.extend(["", f"Generated designs: {len(report.designs)}",
                  f"Pending coordinate verification: {len(report.pending_targets)}", "", "Provenance:",
                  json.dumps(report.to_dict()["provenance"], sort_keys=True, ensure_ascii=False)])
    return "\n".join(lines) + "\n"


def deposit(report: PanelReport, directory: Path) -> tuple[Path, Path]:
    """Create immutable-by-convention artifacts; refuse to overwrite either file."""
    report.validate_for_rendering()
    if report.banner != DEPOSIT_HONESTY_BANNER or report.prediction_status != DEPOSIT_PROSPECTIVE_STATUS:
        raise ValueError("deposit refuses a non-prospective panel")
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    stamp = report.generated_at.replace(":", "").replace("+00:00", "Z")
    json_path = directory / f"design-panel-{stamp}.json"
    text_path = directory / f"design-panel-{stamp}.txt"
    # Preflight prevents the common partial deposit when either name is occupied.
    if json_path.exists() or text_path.exists():
        raise FileExistsError(f"refusing to overwrite frozen panel deposit {stamp}")
    with json_path.open("x", encoding="utf-8") as handle:
        handle.write(render_json(report))
    try:
        with text_path.open("x", encoding="utf-8") as handle:
            handle.write(render_text(report))
    except Exception:
        json_path.unlink()
        raise
    return json_path, text_path
