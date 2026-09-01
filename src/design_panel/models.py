"""Immutable public types and trust boundaries for cross-disease panels."""
from __future__ import annotations

from dataclasses import asdict, dataclass, fields, is_dataclass
from enum import Enum
from types import MappingProxyType
from typing import Any, Mapping
import hashlib
import json
from datetime import datetime, timezone

from design_engine import (PROSPECTIVE_BANNER, PROSPECTIVE_STATUS, DesignReport,
                           ObjectiveDirection, TargetDefinition)

PENDING_COORDINATES_STATUS = "pending coordinate verification"
PENDING_CONTEXT_STATUS = "pending real transcript context"
PENDING_STATUS = PENDING_COORDINATES_STATUS
PANEL_PREDICTION_STATUS = PROSPECTIVE_STATUS
PANEL_BANNER = PROSPECTIVE_BANNER


def _freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType({str(k): _freeze(v) for k, v in value.items()})
    if isinstance(value, (list, tuple)): return tuple(_freeze(v) for v in value)
    if isinstance(value, set): return frozenset(_freeze(v) for v in value)
    return value


def _plain(value: Any) -> Any:
    if is_dataclass(value): return {f.name: _plain(getattr(value, f.name)) for f in fields(value)}
    if isinstance(value, Mapping): return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list, set, frozenset)): return [_plain(v) for v in value]
    if isinstance(value, Enum): return value.value
    return value


def _manifest_hash(value: Any) -> str:
    def default(item: Any) -> Any:
        if isinstance(item, Enum): return item.value
        if is_dataclass(item): return asdict(item)
        return repr(item)
    raw = json.dumps(_plain(value), default=default, sort_keys=True,
                     separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(raw.encode()).hexdigest()


@dataclass(frozen=True)
class PanelTarget:
    target_id: str
    locus: str
    display_name: str
    coordinate_status: str
    coordinate_provenance: str
    event_id: str | None
    assembly: str | None = None
    chrom: str | None = None
    strand: str | None = None
    transcript_version: str | None = None
    provenance_digest: str | None = None
    sequence_sha256: str | None = None
    mapping_sha256: str | None = None

    def __post_init__(self) -> None:
        if self.coordinate_status not in {"verified", "to-verify"}:
            raise ValueError(f"invalid coordinate_status for {self.target_id}")
        if not self.coordinate_provenance: raise ValueError("coordinate_provenance is required")
        if self.coordinate_status == "verified":
            if not all((self.event_id, self.assembly, self.chrom, self.strand,
                        self.transcript_version, self.provenance_digest)):
                raise ValueError(f"verified target {self.target_id} lacks registry binding")
        elif self.event_id is not None or "REQUIRED-INPUT" not in self.coordinate_provenance:
            raise ValueError(f"unverified target {self.target_id} must be REQUIRED-INPUT without event_id")

    @property
    def coordinates_verified(self) -> bool: return self.coordinate_status == "verified"


@dataclass(frozen=True)
class RealTranscriptContext:
    """Operator-supplied, sequence-bound context from a networked session."""
    definition: TargetDefinition
    assembly: str
    chrom: str
    strand: str
    transcript_version: str
    sequence_provenance: str
    genomic_to_transcript_mapping: tuple[tuple[int, int], ...]
    binding_sequence_sha256: str
    mapping_sha256: str

    def __post_init__(self) -> None:
        if not isinstance(self.definition, TargetDefinition): raise TypeError("real context requires TargetDefinition")
        if self.definition.context.transcript_id != self.transcript_version:
            raise ValueError("transcript ID must be transcript-version exact")
        if not self.sequence_provenance.strip() or not self.genomic_to_transcript_mapping:
            raise ValueError("real sequence provenance and genomic mapping are required")
        if any(len(x) != 2 or x[0] < 0 or x[1] < 0 for x in self.genomic_to_transcript_mapping):
            raise ValueError("invalid genomic-to-transcript mapping")
        if not self.chrom.startswith("chr") or self.strand not in {"+", "-"}:
            raise ValueError("real context requires chromosome and strand")
        transcript = tuple(x[1] for x in self.genomic_to_transcript_mapping)
        genomic = tuple(x[0] for x in self.genomic_to_transcript_mapping)
        expected_genomic_step = 1 if self.strand == "+" else -1
        if (transcript != tuple(range(len(self.definition.context.sequence)))
                or any(b - a != expected_genomic_step for a, b in zip(genomic, genomic[1:]))):
            raise ValueError("genomic mapping must exactly and directionally cover the transcript sequence")
        for value in (self.binding_sequence_sha256, self.mapping_sha256):
            if len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
                raise ValueError("real context hashes must be SHA-256 digests")


@dataclass(frozen=True)
class PendingTarget:
    target_id: str
    locus: str
    display_name: str
    status: str
    coordinate_provenance: str


class DesignabilityState(str, Enum):
    PENDING = "PENDING"
    NON_DESIGNABLE = "NON_DESIGNABLE"
    DESIGNABLE = "DESIGNABLE"


@dataclass(frozen=True)
class DesignabilityRow:
    target_id: str
    locus: str
    therapeutic_goal: ObjectiveDirection
    coordinates_verified: bool
    state: DesignabilityState
    status: str
    best_candidate_ids: tuple[str, ...]
    best_candidate_rank: int | None
    best_directional_confidence: float | None
    best_pick_off_target_burden: float | None
    abstention_rate: float | None
    pareto_candidate_count: int
    coordinate_provenance: str


@dataclass(frozen=True)
class PanelReport:
    schema_version: str
    banner: str
    prediction_status: str
    generated_at: str
    designs: tuple[DesignReport, ...]
    pending_targets: tuple[PendingTarget, ...]
    designability_map: tuple[DesignabilityRow, ...]
    provenance: Mapping[str, Any]
    panel_manifest_sha256: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "provenance", _freeze(self.provenance))
        self.validate_for_rendering()

    def validate_for_rendering(self) -> None:
        # Lazy import avoids models -> panel -> models at module-import time.  These
        # authenticated, frozen files are the trust anchor; provenance is not.
        from . import panel as frozen_panel
        frozen_targets = frozen_panel.load_panel(frozen_panel.DEFAULT_PANEL_CONFIG,
                                                  frozen_panel.DEFAULT_COORDINATE_REGISTRY)
        frozen_ordered_targets = tuple({**asdict(target), "expected_goal":
            frozen_panel.EXPECTED_GOALS.get(target.event_id or target.target_id,
                                            ObjectiveDirection.INCREASE_INCLUSION).value}
            for target in frozen_targets)
        if self.banner != PANEL_BANNER or self.prediction_status != PANEL_PREDICTION_STATUS:
            raise ValueError("panel honesty banner/status must remain exactly prospective")
        if self.provenance.get("prediction_status") != PANEL_PREDICTION_STATUS:
            raise ValueError("panel provenance must carry prospective status")
        value = self.panel_manifest_sha256
        if len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
            raise ValueError("panel manifest requires a SHA-256 digest")
        if self.provenance.get("panel_manifest_sha256") != value:
            raise ValueError("stored panel manifest digest mismatch")
        if _manifest_hash(self.provenance.get("panel_manifest")) != value:
            raise ValueError("panel manifest content does not match its digest")
        for design in self.designs: design.validate_for_rendering()
        manifest = self.provenance["panel_manifest"]
        if (manifest.get("panel_version") != frozen_panel.PANEL_VERSION
                or manifest.get("engine_api") != frozen_panel.APPROVED_ENGINE_API
                or not isinstance(manifest.get("options"), Mapping)
                or dict(manifest["options"]) != {}
                or manifest.get("coordinate_registry_sha256")
                    != frozen_panel.EXPECTED_REGISTRY_SHA256):
            raise ValueError("panel manifest metadata does not match frozen policy")
        manifest_ordered_targets = manifest.get("ordered_targets")
        if (not isinstance(manifest_ordered_targets, (tuple, list))
                or tuple(_plain(item) for item in manifest_ordered_targets) != frozen_ordered_targets):
            raise ValueError("panel manifest ordered targets do not match frozen panel")
        try:
            parsed_time = datetime.fromisoformat(self.generated_at)
        except (TypeError, ValueError) as exc:
            raise ValueError("panel generated_at must be canonical UTC ISO-8601") from exc
        if (parsed_time.tzinfo is None or parsed_time.utcoffset() != timezone.utc.utcoffset(parsed_time)
                or parsed_time.isoformat() != self.generated_at
                or manifest.get("generated_at") != self.generated_at):
            raise ValueError("panel timestamps must be identical canonical UTC ISO-8601")
        if any(design.provenance.get("generated_at") != self.generated_at for design in self.designs):
            raise ValueError("nested report timestamps must equal panel generated_at")
        report_entries = tuple({"target_id": design.target_id,
                                "run_manifest_sha256": design.provenance["run_manifest_sha256"]}
                               for design in self.designs)
        if tuple(dict(x) for x in manifest.get("report_manifests", ())) != report_entries:
            raise ValueError("panel designs do not match report manifest list")

        ordered_targets = frozen_ordered_targets
        target_ids = tuple(item.get("target_id") for item in ordered_targets)
        row_ids = tuple(row.target_id for row in self.designability_map)
        if (not target_ids or len(target_ids) != len(set(target_ids))
                or row_ids != target_ids):
            raise ValueError("designability map does not exactly cover ordered panel targets")
        design_ids = tuple(design.target_id for design in self.designs)
        pending_ids = tuple(item.target_id for item in self.pending_targets)
        if len(design_ids) != len(set(design_ids)) or len(pending_ids) != len(set(pending_ids)):
            raise ValueError("panel output target IDs must be unique")
        if set(design_ids) & set(pending_ids) or set(design_ids) | set(pending_ids) != set(target_ids):
            raise ValueError("every panel target must be exactly designed or pending")
        if design_ids != tuple(target_id for target_id in target_ids if target_id in set(design_ids)):
            raise ValueError("panel designs must follow frozen target order")
        rows = {row.target_id: row for row in self.designability_map}
        if any(rows[target_id].state is DesignabilityState.PENDING for target_id in design_ids):
            raise ValueError("designed targets cannot be marked pending")
        if any(rows[target_id].state is not DesignabilityState.PENDING for target_id in pending_ids):
            raise ValueError("pending targets cannot be marked designable")
        if any(not rows[target_id].coordinates_verified for target_id in design_ids):
            raise ValueError("unverified targets cannot carry designs")
        frozen_by_id = {target.target_id: target for target in frozen_targets}
        if any(not frozen_by_id[target_id].sequence_sha256 for target_id in design_ids):
            raise ValueError("targets without a frozen committed sequence digest cannot carry designs")
        # The independently validated design reports, not caller-held manifest
        # provenance, are the authority for which targets were designed.
        design_by_id = {item.target_id: item for item in self.designs}
        for target_id, report in design_by_id.items():
            target = frozen_by_id[target_id]
            run_manifest = report.provenance["run_manifest"]
            definition = run_manifest.get("target_definition", {})
            if not isinstance(definition, Mapping):
                definition = {}
            context = definition.get("context", {})
            sequence = context.get("sequence") if isinstance(context, Mapping) else None
            if (not isinstance(sequence, str)
                    or frozen_panel._sequence_sha256(sequence) != target.sequence_sha256):
                raise ValueError(
                    f"design report sequence does not match frozen committed digest for {target_id}")
            backend = run_manifest.get("backend", {})
            if (not isinstance(backend, Mapping)
                    or backend.get("identity_sha256")
                        != frozen_panel.APPROVED_BACKEND_IDENTITY_SHA256
                    or report.provenance.get("folding_backend_identity_sha256")
                        != frozen_panel.APPROVED_BACKEND_IDENTITY_SHA256):
                raise ValueError("design report lacks approved backend receipt")
            objective = definition.get("objective", {})
            event = definition.get("event", {})
            expected = frozen_panel.EXPECTED_GOALS[target.event_id or target_id]
            if (report.target_id != target_id or report.locus != target.locus
                    or report.objective.event_id != target.event_id
                    or report.objective.locus != target.locus
                    or report.objective.direction is not expected
                    or definition.get("target_id") != target_id
                    or definition.get("locus") != target.locus
                    or not isinstance(objective, Mapping)
                    or objective.get("event_id") != target.event_id
                    or objective.get("locus") != target.locus
                    or objective.get("direction") != expected.value
                    or not isinstance(event, Mapping)
                    or event.get("event_id") != target.event_id):
                raise ValueError(f"design report identity/polarity mismatch for {target_id}")
        # Authenticity anchor: re-derive each design from the FROZEN committed context
        # (frozen_panel.committed_real_contexts, never the caller-held report) with the
        # approved engine (frozen_panel.APPROVED_DESIGN_ENGINE).  The engine is
        # deterministic (default DeterministicStubFolder + frozen configs), so a genuine
        # report reproduces byte-for-byte while a hand-built report with fabricated or
        # mutated candidates, scores, ranks, or Pareto set cannot.  The committed sequence
        # being public in the repo is not enough to forge a design: the caller would also
        # have to reproduce the exact deterministic engine output, which IS the real run.
        committed_contexts = frozen_panel.committed_real_contexts()
        for target_id, report in design_by_id.items():
            context = committed_contexts.get(target_id)
            if context is None:
                raise ValueError(f"no committed real context to re-derive design for {target_id}")
            if context.binding_sequence_sha256 != frozen_by_id[target_id].sequence_sha256:
                raise ValueError(f"committed context does not match frozen digest for {target_id}")
            expected = frozen_panel.APPROVED_DESIGN_ENGINE(context.definition, clock=lambda: parsed_time)
            if _manifest_hash(expected) != _manifest_hash(report):
                raise ValueError(f"design report is not the deterministic engine output for {target_id}")
        # Accepted downgrade-only residual: frozen files do not record whether an
        # operator supplied a runtime RealTranscriptContext, so a withheld context
        # is indistinguishable here from legitimate no-context.  That can only leave
        # a real design PENDING; it cannot create the validated engine/backend receipt
        # required for a designable claim, which is the safe direction for this harness.
        expected_pending = []
        expected_rows = []
        for raw in ordered_targets:
            target_id = raw["target_id"]
            if target_id in design_by_id:
                report = design_by_id[target_id]
                abstention = (sum(row.abstained for row in report.candidates) / len(report.candidates)
                              if report.candidates else 1.0)
                supported = [row for row in report.candidates if not row.abstained
                             and row.goal_alignment == "supports_goal" and row.rank is not None]
                reasons = []
                if not supported: reasons.append("no goal-supporting candidates")
                if abstention > .50: reasons.append("abstention threshold exceeded")
                eligible = [row for row in supported if row.off_target_burden <= 1.0]
                if supported and not eligible: reasons.append("off-target burden threshold exceeded")
                designable = not reasons and bool(eligible)
                rank = min((row.rank for row in eligible), default=None)
                best = [row for row in eligible if row.rank == rank]
                expected_rows.append(DesignabilityRow(target_id, raw["locus"], report.objective.direction,
                    True, DesignabilityState.DESIGNABLE if designable else DesignabilityState.NON_DESIGNABLE,
                    "prospective designs meet declared thresholds" if designable else "; ".join(reasons),
                    tuple(sorted(x.candidate.candidate_id for x in best)) if designable else (),
                    rank if designable else None,
                    max((x.directional_confidence for x in best), default=None) if designable else None,
                    min((x.off_target_burden for x in best), default=None) if designable else None,
                    round(abstention, 6), len(report.pareto_candidate_ids), raw["coordinate_provenance"]))
            else:
                if raw["coordinate_status"] != "verified": status = PENDING_COORDINATES_STATUS
                elif not raw.get("sequence_sha256"): status = "pending real transcript context (no committed digest)"
                else:
                    # Downgrade-only: this cannot fabricate a validated/designable claim.
                    # Runtime context-availability is not recorded in frozen files, so
                    # pending-context is not render-time falsifiable.
                    status = PENDING_CONTEXT_STATUS
                pending = PendingTarget(target_id, raw["locus"], raw["display_name"], status,
                                        raw["coordinate_provenance"])
                expected_pending.append(pending)
                goal = ObjectiveDirection(raw.get("expected_goal", "inclusion"))
                expected_rows.append(DesignabilityRow(target_id, raw["locus"], goal,
                    raw["coordinate_status"] == "verified", DesignabilityState.PENDING, status,
                    (), None, None, None, None, 0, raw["coordinate_provenance"]))
        if tuple(expected_pending) != self.pending_targets:
            raise ValueError("pending target records do not match ordered panel targets")
        if tuple(expected_rows) != self.designability_map:
            raise ValueError("designability rows do not match deterministic panel outputs")

    def to_dict(self) -> dict[str, Any]:
        self.validate_for_rendering()
        return _plain(self)
