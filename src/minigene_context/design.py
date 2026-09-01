"""Reporter-specific design boundary and serialization."""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import json
from pathlib import Path
from typing import Mapping

from design_engine import (DesignEngine, DesignReport, DeterministicStubFolder,
                           GeneratedCandidate, GenerationConfig, generate_candidates)
from design_engine.targets import DEFAULT_EVENT_CONFIG
from offtarget import EnumeratorConfig, enumerate_offtargets
from physics_adapter import ModelConfig

from .models import CLAIM_BOUNDARY, CompressedReporterContext, ContextDomain, _plain


def _touches_join(start: int, end: int, join: int) -> bool:
    # Endpoint equality is intentionally rejected: "touches" is broader than
    # only crossing a join.
    return start <= join <= end


def _validate_candidate_joins(candidates, joins: tuple[int, ...], grna_shadow: int) -> None:
    for value in candidates:
        candidate: GeneratedCandidate = value.candidate if hasattr(value, "candidate") else value
        for join in joins:
            if _touches_join(candidate.start, candidate.end, join):
                raise ValueError(
                    f"candidate footprint touches artificial reporter junction: {candidate.candidate_id}")
            if candidate.effector_class == "dCas13_gRNA" and _touches_join(
                    candidate.start - grna_shadow, candidate.end + grna_shadow, join):
                raise ValueError(
                    f"dCas13 shadow touches artificial reporter junction: {candidate.candidate_id}")


def _assert_approved_engine(engine: DesignEngine,
        _engine_class=DesignEngine, _folder_class=DeterministicStubFolder,
        _generation=GenerationConfig(), _adapter=ModelConfig(),
        _offtarget=EnumeratorConfig(max_edits=4, max_indel_bases=1),
        _enumerator=enumerate_offtargets, _event_config=DEFAULT_EVENT_CONFIG) -> None:
    """Pin the exact deterministic stack accepted by reporter serialization."""
    if type(engine) is not _engine_class or type(engine.folder) is not _folder_class:
        raise ValueError("reporter design requires the approved deterministic engine/backend")
    if (engine.generation_config != _generation or engine.adapter_config != _adapter
            or engine.offtarget_config != _offtarget or engine.enumerator is not _enumerator
            or Path(engine.event_config).resolve() != Path(_event_config).resolve()
            or engine.magnitude_evidence or engine.magnitude_claims):
        raise ValueError("reporter design engine configuration is not approved")


def _rerun_approved(target, generated_at: str,
        _engine_class=DesignEngine, _folder_class=DeterministicStubFolder,
        _generation=GenerationConfig(), _adapter=ModelConfig(),
        _offtarget=EnumeratorConfig(max_edits=4, max_indel_bases=1),
        _enumerator=enumerate_offtargets, _event_config=DEFAULT_EVENT_CONFIG) -> DesignReport:
    try:
        timestamp = datetime.fromisoformat(generated_at)
    except (TypeError, ValueError) as error:
        raise ValueError("reporter design timestamp is invalid") from error
    if timestamp.tzinfo is None:
        raise ValueError("reporter design timestamp must be timezone-aware")
    engine = _engine_class(folder=_folder_class(), generation_config=_generation,
        adapter_config=_adapter, offtarget_config=_offtarget, enumerator=_enumerator,
        event_config=_event_config, clock=lambda: timestamp)
    return engine.design_for(target)


@dataclass(frozen=True)
class ReporterDesignReport:
    """Design output whose reporter domain and limits cannot be omitted."""

    schema_version: str
    domain_label: ContextDomain
    claim_boundary: str
    compression_receipt: CompressedReporterContext
    design_report: DesignReport

    def __post_init__(self) -> None:
        self.validate_for_rendering()

    def validate_for_rendering(self, _rerun=_rerun_approved,
                               _join_validator=_validate_candidate_joins) -> None:
        self.compression_receipt.validate_for_use()
        self.design_report.validate_for_rendering()
        if self.schema_version != "minigene-reporter-design-report/1":
            raise ValueError("unsupported minigene reporter design-report schema")
        if self.domain_label is not ContextDomain.MINIGENE_REPORTER:
            raise ValueError("design report domain must remain minigene_reporter")
        if self.claim_boundary != CLAIM_BOUNDARY:
            raise ValueError("design report claim boundary changed")
        target = self.compression_receipt.target
        if (self.design_report.target_id != target.target_id
                or self.design_report.locus != target.locus
                or self.design_report.objective != target.objective):
            raise ValueError("design report identity does not match its compression receipt")
        manifest = self.design_report.provenance["run_manifest"]
        if _plain(manifest.get("target_definition")) != _plain(target):
            raise ValueError("design report target definition differs from compression receipt")
        config = manifest.get("generation_config", {})
        if not isinstance(config, Mapping) or "grna_shadow" not in config:
            raise ValueError("design report lacks the candidate shadow configuration")
        _join_validator(self.design_report.candidates,
            target.context.artificial_junction_offsets, int(config["grna_shadow"]))
        expected = _rerun(target, self.design_report.provenance.get("generated_at", ""))
        if _plain(expected) != _plain(self.design_report):
            raise ValueError("design report is not the deterministic approved-engine output")

    def to_dict(self) -> dict:
        self.validate_for_rendering()
        return _plain(self)


def design_compressed_context(
        compressed: CompressedReporterContext,
        engine: DesignEngine | None = None, *,
        _engine_class=DesignEngine, _engine_validator=_assert_approved_engine,
        _generator=generate_candidates,
        _join_validator=_validate_candidate_joins) -> ReporterDesignReport:
    """Run the engine only after reporter-junction safety checks pass."""
    if not isinstance(compressed, CompressedReporterContext):
        raise TypeError("compressed must be a CompressedReporterContext")
    compressed.validate_for_use()
    selected_engine = engine or _engine_class()
    _engine_validator(selected_engine)
    generated = _generator(compressed.target, selected_engine.generation_config)
    _join_validator(generated,
        compressed.target.context.artificial_junction_offsets,
        selected_engine.generation_config.grna_shadow)
    report = selected_engine.design_for(compressed.target)
    return ReporterDesignReport("minigene-reporter-design-report/1",
        ContextDomain.MINIGENE_REPORTER, CLAIM_BOUNDARY, compressed, report)


def render_json(report: ReporterDesignReport) -> str:
    return json.dumps(report.to_dict(), sort_keys=True, indent=2) + "\n"


def render_text(report: ReporterDesignReport) -> str:
    report.validate_for_rendering()
    context = report.compression_receipt.target.context
    return "\n".join((
        f"DOMAIN: {report.domain_label.value}",
        f"CLAIM BOUNDARY: {report.claim_boundary}",
        f"SOURCE SEQUENCE SHA-256: {context.source_sequence_sha256}",
        f"SOURCE PROVENANCE: {context.source_sequence_provenance}",
        f"COMPRESSION MANIFEST SHA-256: {context.compression_manifest_sha256}",
        f"ARTIFICIAL JUNCTION OFFSETS: {list(context.artificial_junction_offsets)}",
        f"CANDIDATES: {len(report.design_report.candidates)}",
        "",
    ))
