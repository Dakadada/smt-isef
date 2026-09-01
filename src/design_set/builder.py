"""Build the SMN2 design set from the authenticated committed context."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from typing import Callable, Iterable

from design_engine import (PROSPECTIVE_BANNER, PROSPECTIVE_STATUS, DesignReport,
                           ScoredCandidate, design_for)
from design_panel.models import RealTranscriptContext
from design_panel.panel import committed_real_contexts

from .models import (CLAIM_BOUNDARY, EMPTY_OFF_TARGET_DISCLOSURE,
                     DesignSetCandidate, ProspectiveDesignSet, ReferencePlacement)


TARGET_ID = "SMN2_exon7_inclusion"
REFERENCE_LENGTH = 18
REFERENCE_EFFECTOR_CLASS = "ASO"
REFERENCE_CHEMISTRY = "2MOE"
REFERENCE_TARGET_SEQUENCE = "CCAGCATTATGAAAGTGA"
REFERENCE_EFFECTOR_SEQUENCE = "TCACTTTCATAATGCTGG"
REFERENCE_GENOMIC_START = 70076584
REFERENCE_GENOMIC_END = 70076601


def _axes(row: ScoredCandidate) -> tuple[float, float, float]:
    return (row.directional_confidence, -row.off_target_burden, row.span_robustness)


def prediction_axis_dominates(a: ScoredCandidate, b: ScoredCandidate) -> bool:
    """Strict Pareto dominance on confidence, burden, and span robustness only."""
    left, right = _axes(a), _axes(b)
    return all(x >= y for x, y in zip(left, right)) and any(x > y for x, y in zip(left, right))


def prediction_axis_tiers(candidates: Iterable[ScoredCandidate]) -> dict[str, int]:
    """Return deterministic three-axis Pareto tiers for goal-supporting predictions."""
    remaining = [row for row in candidates
                 if not row.abstained and row.goal_alignment == "supports_goal"]
    identifiers = [row.candidate.candidate_id for row in remaining]
    if len(identifiers) != len(set(identifiers)):
        raise ValueError("candidate IDs must be unique before design-set Pareto analysis")
    tiers: dict[str, int] = {}
    tier = 1
    while remaining:
        members = [row for row in remaining
                   if not any(prediction_axis_dominates(other, row)
                              for other in remaining if other is not row)]
        for row in members:
            tiers[row.candidate.candidate_id] = tier
        member_ids = {id(row) for row in members}
        remaining = [row for row in remaining if id(row) not in member_ids]
        tier += 1
    return tiers


def superior_candidate_ids(candidates: Iterable[ScoredCandidate],
                           reference: ScoredCandidate | None) -> tuple[str, ...]:
    if reference is None:
        return ()
    return tuple(sorted(row.candidate.candidate_id for row in candidates
                        if not row.abstained and row.goal_alignment == "supports_goal"
                        and prediction_axis_dominates(row, reference)))


def _reference_placement(report: DesignReport,
                         context: RealTranscriptContext) -> tuple[ReferencePlacement,
                                                                  ScoredCandidate | None]:
    label = "nusinersen-class ISS-N1 steric ASO reference"
    limitation = (
        "The committed context cannot verify the expected ISS-N1 18-mer sequence and "
        "GRCh38 genomic-to-transcript mapping; only the reference position within the "
        "committed window is asserted. Coordinates and predicted-superiority comparisons "
        "are withheld."
    )
    definition = context.definition
    expected_genomic = tuple(range(REFERENCE_GENOMIC_START, REFERENCE_GENOMIC_END + 1))
    mapping = dict(context.genomic_to_transcript_mapping)
    transcript_positions = tuple(mapping.get(position) for position in expected_genomic)
    mapping_resolves = (None not in transcript_positions
                        and transcript_positions == tuple(range(
                            transcript_positions[0],
                            transcript_positions[0] + REFERENCE_LENGTH)))
    start = transcript_positions[0] if mapping_resolves else None
    end = start + REFERENCE_LENGTH if start is not None else None
    context_matches = (
        definition.target_id == TARGET_ID and definition.locus == "SMN2"
        and context.assembly == "GRCh38" and context.chrom == "chr5"
        and context.strand == "+" and context.transcript_version == "NM_017411.4"
        and start is not None and end is not None
        and end <= len(definition.context.sequence)
        and definition.context.sequence[start:end].upper().replace("U", "T")
            == REFERENCE_TARGET_SEQUENCE)
    if not context_matches:
        return ReferencePlacement(label, REFERENCE_EFFECTOR_CLASS, REFERENCE_CHEMISTRY,
                                  REFERENCE_LENGTH, "reference-position-within-committed-window",
                                  None, None, None, limitation), None
    matches = [row for row in report.candidates
               if row.candidate.effector_class == REFERENCE_EFFECTOR_CLASS
               and row.candidate.chemistry == REFERENCE_CHEMISTRY
               and row.candidate.start == start and row.candidate.end == end
               and row.candidate.target_sequence == REFERENCE_TARGET_SEQUENCE
               and row.candidate.effector_sequence == REFERENCE_EFFECTOR_SEQUENCE]
    if len(matches) != 1:
        return ReferencePlacement(label, REFERENCE_EFFECTOR_CLASS, REFERENCE_CHEMISTRY,
                                  REFERENCE_LENGTH, "reference-position-within-committed-window",
                                  None, None, None, limitation), None
    row = matches[0]
    return ReferencePlacement(label, REFERENCE_EFFECTOR_CLASS, REFERENCE_CHEMISTRY,
                              REFERENCE_LENGTH, "resolved-full-18mer", row.candidate.start,
                              row.candidate.end, row.candidate.candidate_id, None), row


def _sha256_json(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def build_design_set(*, clock: Callable[[], datetime] | None = None) -> ProspectiveDesignSet:
    timestamp = (clock or (lambda: datetime.now(timezone.utc)))()
    if timestamp.tzinfo is None:
        raise ValueError("design-set clock must return a timezone-aware datetime")
    timestamp = timestamp.astimezone(timezone.utc)
    contexts = committed_real_contexts()
    context = contexts.get(TARGET_ID)
    if context is None:
        raise ValueError("frozen committed SMN2 real context is required")
    definition = context.definition
    if definition.target_id != TARGET_ID or definition.locus != "SMN2":
        raise ValueError("committed context identity mismatch")
    report = design_for(definition, clock=lambda: timestamp)
    report.validate_for_rendering()
    placement, reference = _reference_placement(report, context)
    tiers = prediction_axis_tiers(report.candidates)
    superior = superior_candidate_ids(report.candidates, reference)
    superior_set = set(superior)
    off_target_corpus_size = len(definition.offtarget_transcripts)
    rows = tuple(DesignSetCandidate(
        PROSPECTIVE_BANNER, PROSPECTIVE_STATUS, scored.candidate.candidate_id,
        scored.candidate.effector_class, scored.candidate.chemistry,
        scored.candidate.start, scored.candidate.end, scored.candidate.target_sequence,
        scored.candidate.effector_sequence, scored.directional_confidence,
        scored.off_target_burden, scored.span_robustness,
        tiers.get(scored.candidate.candidate_id), tiers.get(scored.candidate.candidate_id) == 1,
        scored.candidate.candidate_id in superior_set,
        scored.candidate.candidate_id == placement.candidate_id)
        for scored in report.candidates)
    provenance = {
        "prediction_status": PROSPECTIVE_STATUS,
        "engine_api": "design_engine.engine.design_for",
        "engine_report_sha256": _sha256_json(report.to_dict()),
        "engine_run_manifest_sha256": report.provenance["run_manifest_sha256"],
        "committed_context_source": "design_panel.panel.committed_real_contexts",
        "binding_sequence_sha256": context.binding_sequence_sha256,
        "mapping_sha256": context.mapping_sha256,
        "sequence_provenance": context.sequence_provenance,
        "reference_basis": (
            f"Verified committed {context.assembly} {context.chrom}:"
            f"{REFERENCE_GENOMIC_START}-{REFERENCE_GENOMIC_END} maps to zero-based "
            f"half-open [{placement.start},{placement.end}) and matches the expected "
            "ISS-N1 18-mer sequence."
            if placement.candidate_id is not None else placement.limitation
        ),
        "off_target_corpus_size": off_target_corpus_size,
        "off_target_burden_evaluated": off_target_corpus_size > 0,
        "off_target_burden_discriminating": (
            off_target_corpus_size > 0
            and len({row.off_target_burden for row in report.candidates}) > 1
        ),
        "off_target_disclosure": (
            EMPTY_OFF_TARGET_DISCLOSURE if off_target_corpus_size == 0
            else "Off-target burden was evaluated against the committed transcript corpus."
        ),
        "superiority_axes": [
            "directional_confidence: higher",
            "off_target_burden: lower",
            "span_robustness: higher",
        ],
        "superiority_rule": "weakly dominates reference on all axes and strictly beats it on at least one",
    }
    return ProspectiveDesignSet("prospective-design-set/1", PROSPECTIVE_BANNER,
                                PROSPECTIVE_STATUS, CLAIM_BOUNDARY, timestamp.isoformat(),
                                TARGET_ID, "SMN2", "zero-based half-open within committed context",
                                placement, rows, superior, provenance)
