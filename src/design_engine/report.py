"""Immutable-friendly JSON and concise human-readable report rendering."""
from __future__ import annotations

import json

from .models import DesignReport, PROSPECTIVE_BANNER, PROSPECTIVE_STATUS


def render_json(report: DesignReport) -> str:
    report.validate_for_rendering()
    payload = report.to_dict()
    # Constants are emitted independently of the report's fields.
    payload["banner"] = PROSPECTIVE_BANNER
    payload["prediction_status"] = PROSPECTIVE_STATUS
    return json.dumps(payload, indent=2, sort_keys=True) + "\n"


def render_text(report: DesignReport) -> str:
    report.validate_for_rendering()
    lines = [PROSPECTIVE_BANNER, f"Status: {PROSPECTIVE_STATUS}", "",
             f"Target: {report.target_id} ({report.locus})",
             f"Therapeutic goal: exon {report.therapeutic_polarity}",
             f"Candidates: {len(report.candidates)}; Pareto equivalence tier: {len(report.pareto_candidate_ids)}",
             "", "RANKED GOAL-SUPPORTING EQUIVALENCE TIERS:"]
    ranked = [row for row in report.candidates if row.rank is not None]
    if not ranked:
        lines.append("(none; no rank-1 candidate exists)")
    for row in ranked:
        claim = "direction-only; magnitude structurally absent" if row.magnitude_claim is None else "AUDITED MAGNITUDE CLAIM"
        lines.append(f"tier {row.rank:2d} {row.candidate.candidate_id} | {row.goal_alignment} | "
                     f"confidence={row.directional_confidence:.2f} | off-target={row.off_target_burden:.3f} | "
                     f"pareto-tier={'yes' if row.pareto else 'no'} | {claim}")
    lines.extend(["", "UNRANKED DIAGNOSTICS (ABSTAINED OR GOAL-OPPOSING):"])
    diagnostics = [row for row in report.candidates if row.rank is None]
    if not diagnostics:
        lines.append("(none)")
    for row in diagnostics:
        lines.append(f"UNRANKED {row.candidate.candidate_id} | {row.goal_alignment} | "
                     f"confidence={row.directional_confidence:.2f} | off-target={row.off_target_burden:.3f} | "
                     f"{'ABSTAIN' if row.abstained else 'OPPOSES GOAL'}")
    lines.extend(["", "Provenance:", json.dumps(report.to_dict()["provenance"], sort_keys=True)])
    return "\n".join(lines) + "\n"
