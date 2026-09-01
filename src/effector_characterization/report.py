"""Render and immutably deposit engine-internal characterizations."""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

from .models import EffectorScoreCharacterization


def render_json(report: EffectorScoreCharacterization) -> str:
    report.validate_for_rendering()
    return json.dumps(report.to_dict(), indent=2, sort_keys=True, ensure_ascii=False) + "\n"


def render_text(report: EffectorScoreCharacterization) -> str:
    report.validate_for_rendering()
    lines = [report.banner, report.claim_boundary, f"P1 gate: {report.p1_status}", "",
             f"Frozen characterization timestamp: {report.generated_at}",
             f"Target: {report.target_id} ({report.locus}); committed context: "
             f"{report.context_length} nt",
             f"Score: {report.score_name}", f"Formula: {report.score_formula}",
             f"Class bases: {dict(report.class_bases)}",
             f"ASO chemistry deltas: {dict(report.aso_chemistry_deltas)}",
             f"gRNA shadow: {report.grna_shadow_nt} nt on each side, clipped to context; "
             "this changes the physics footprint, not the manufacturability scalar.", "",
             f"Tiling limitation: {report.provenance['tiling_disclosure']}", "",
             "CLASS DISTRIBUTIONS:",
             "effector_class | count | min | max | mean | availability | limitation"]
    for row in report.class_summaries:
        values = ["—" if value is None else f"{value:.6f}"
                  for value in (row.minimum, row.maximum, row.mean)]
        lines.append(f"{row.effector_class} | {row.candidate_count} | "
                     f"{' | '.join(values)} | {row.availability} | {row.limitation or '—'}")
    lines.extend(["", "CONFIGURED CHEMISTRY/LENGTH DISTRIBUTIONS:",
                  "effector_class | chemistry | length | count | observed_scores | min | max | mean"])
    for row in report.distributions:
        values = ["—" if value is None else f"{value:.6f}"
                  for value in (row.minimum, row.maximum, row.mean)]
        lines.append(f"{row.effector_class} | {row.chemistry} | {row.length} | "
                     f"{row.candidate_count} | {list(row.observed_scores)} | {' | '.join(values)}")
    lines.extend(["", "ACTUAL ENGINE CANDIDATES:",
                  "candidate_id | class | chemistry | start | end | length | base | "
                  "chemistry_delta | length_penalty | engine_score | gRNA_shadow"])
    for row in report.candidates:
        shadow = ("—" if row.grna_shadow_start is None
                  else f"[{row.grna_shadow_start},{row.grna_shadow_end})")
        lines.append(f"{row.candidate_id} | {row.effector_class} | {row.chemistry} | "
                     f"{row.start} | {row.end} | {row.length} | {row.class_base:.6f} | "
                     f"{row.chemistry_delta:.6f} | {row.length_penalty:.6f} | "
                     f"{row.engine_score:.6f} | {shadow}")
    lines.extend(["", "Provenance:", json.dumps(report.to_dict()["provenance"],
                                                  sort_keys=True, ensure_ascii=False)])
    return "\n".join(lines) + "\n"


def deposit(report: EffectorScoreCharacterization, directory: Path) -> tuple[Path, Path]:
    report.validate_for_rendering()
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.fromisoformat(report.generated_at).astimezone(timezone.utc)
    stamp = timestamp.strftime("%Y-%m-%dT%H%M%SZ")
    json_path = directory / f"effector-score-characterization-{stamp}.json"
    text_path = directory / f"effector-score-characterization-{stamp}.txt"
    if json_path.exists() or text_path.exists():
        raise FileExistsError(f"refusing to overwrite characterization deposit {stamp}")
    with json_path.open("x", encoding="utf-8") as handle:
        handle.write(render_json(report))
    try:
        with text_path.open("x", encoding="utf-8") as handle:
            handle.write(render_text(report))
    except Exception:
        json_path.unlink()
        raise
    return json_path, text_path
