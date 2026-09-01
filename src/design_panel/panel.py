"""Fail-closed registry resolution and deterministic prospective aggregation."""
from __future__ import annotations

import csv
import io
from dataclasses import asdict, fields, is_dataclass, replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any, Callable, Mapping

from design_engine import DesignReport, ObjectiveDirection, design_for
from .models import (PANEL_BANNER, PANEL_PREDICTION_STATUS, PENDING_CONTEXT_STATUS,
                     PENDING_COORDINATES_STATUS, DesignabilityRow, DesignabilityState,
                     PanelReport, PanelTarget, PendingTarget, RealTranscriptContext)

DEFAULT_PANEL_CONFIG = Path(__file__).with_name("target_panel.json")
DEFAULT_COORDINATE_REGISTRY = Path(__file__).resolve().parents[2] / "config" / "event_junctions.GRCh38.tsv"
PANEL_VERSION = "prospective-design-panel/2"
MAX_ABSTENTION_RATE = 0.50
MAX_OFFTARGET_BURDEN = 1.0
EXPECTED_GOALS = {"SMN2_exon7_inclusion": ObjectiveDirection.INCREASE_INCLUSION,
                  "DMD_exon53_inclusion": ObjectiveDirection.DECREASE_INCLUSION}
EXPECTED_REGISTRY_SHA256 = "499ae7dd7eee6715eeadf8c2a7a2a68f6aa0de1a9f679afbb8bf5098d98f0f51"
EXPECTED_PANEL_SHA256 = "0de5a15d1725f0743fe238b008ecceacce3a2c38c9a33ba6bf6d3cfb672d1a9d"
APPROVED_ENGINE_API = "design_engine.engine.design_for"
APPROVED_BACKEND_IDENTITY_SHA256 = "013a655b7a3e6e3708153410e202dc2b40933eb180698b9916b089c9a6ac4077"
PENDING_UNCOMMITTED_CONTEXT_STATUS = "pending real transcript context (no committed digest)"

# Frozen seams for render-time re-derivation.  These name the approved engine and the
# authenticated committed-context source that PanelReport.validate_for_rendering re-runs
# to prove a designed report IS the deterministic engine output, not a hand-built forgery.
# They mirror the DEFAULT_*/EXPECTED_* frozen anchors: production uses the real values;
# the test suite patches them (alongside the registry/panel/backend anchors) so genuine
# test panels re-derive while fabricated candidates do not.
APPROVED_DESIGN_ENGINE = design_for


def committed_real_contexts() -> Mapping[str, RealTranscriptContext]:
    """Authenticated operator contexts used ONLY for render-time re-derivation.

    Sourced from the frozen committed files (never from the caller-held report), so the
    re-run's target definition — event intervals, off-target inputs, window — cannot be
    steered by a forged report.  Lazy import: real_contexts -> panel -> models."""
    from .real_contexts import default_real_contexts
    return default_real_contexts()


def _sha256_bytes(raw: bytes) -> str: return hashlib.sha256(raw).hexdigest()
def _sha256(path: Path) -> str: return _sha256_bytes(path.read_bytes())
def _sequence_sha256(sequence: str) -> str: return hashlib.sha256(sequence.encode("utf-8")).hexdigest()


def _hash(value: Any) -> str:
    def plain(item: Any) -> Any:
        if is_dataclass(item): return {field.name: plain(getattr(item, field.name)) for field in fields(item)}
        if isinstance(item, Mapping): return {str(key): plain(val) for key, val in item.items()}
        if isinstance(item, (tuple, list, set, frozenset)): return [plain(val) for val in item]
        if hasattr(item, "value"): return item.value
        return item
    def default(item: Any) -> Any:
        if isinstance(item, Path): return str(item)
        if callable(item): return f"{item.__module__}.{getattr(item, '__qualname__', item.__class__.__qualname__)}"
        return repr(item)
    raw = json.dumps(plain(value), default=default, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(raw.encode()).hexdigest()


def registry_provenance_digest(row: Mapping[str, str], transcript_version: str) -> str:
    return _hash({"assembly": "GRCh38", "row": dict(row), "transcript_version": transcript_version})


def bind_real_transcript_context(*, definition, assembly: str, transcript_version: str,
                                 chrom: str, strand: str,
                                 sequence_provenance: str,
                                 genomic_to_transcript_mapping: tuple[tuple[int, int], ...]) -> RealTranscriptContext:
    """Bind operator-retrieved sequence/mapping bytes to a typed production input."""
    return RealTranscriptContext(definition, assembly, chrom, strand, transcript_version, sequence_provenance,
        genomic_to_transcript_mapping, _sequence_sha256(definition.context.sequence),
        _hash(genomic_to_transcript_mapping))


def load_panel(path: Path = DEFAULT_PANEL_CONFIG,
               registry_path: Path = DEFAULT_COORDINATE_REGISTRY, *,
               _test_only_allow_unfrozen_config: bool = False) -> tuple[PanelTarget, ...]:
    registry_path = Path(registry_path).resolve()
    frozen_registry = DEFAULT_COORDINATE_REGISTRY.resolve()
    registry_bytes = registry_path.read_bytes()
    registry_hash = _sha256_bytes(registry_bytes)
    if registry_path != frozen_registry:
        raise ValueError("coordinate registry must be the frozen project registry")
    if registry_hash != EXPECTED_REGISTRY_SHA256:
        raise ValueError("frozen coordinate registry digest mismatch")
    panel_path = Path(path).resolve()
    panel_bytes = panel_path.read_bytes()
    if not _test_only_allow_unfrozen_config:
        if panel_path != DEFAULT_PANEL_CONFIG.resolve():
            raise ValueError("target panel must be the frozen project panel")
        if _sha256_bytes(panel_bytes) != EXPECTED_PANEL_SHA256:
            raise ValueError("frozen target panel digest mismatch")
    config = json.loads(panel_bytes.decode("utf-8"))
    if config.get("schema_version") != "design-target-panel/2" or config.get("assembly") != "GRCh38":
        raise ValueError("unsupported target panel schema/assembly")
    with io.StringIO(registry_bytes.decode("utf-8"), newline="") as handle:
        registry = {row["event"]: row for row in csv.DictReader(handle, delimiter="\t")}
    targets = tuple(PanelTarget(**item) for item in config.get("targets", ()))
    ids = [item.target_id for item in targets]
    if not targets or len(ids) != len(set(ids)): raise ValueError("target panel must contain unique targets")
    if not _test_only_allow_unfrozen_config:
        # panel_bytes are the already authenticated bytes of DEFAULT_PANEL_CONFIG.
        committed = json.loads(panel_bytes.decode("utf-8"))["targets"]
        if config.get("targets") != committed:
            raise ValueError("target panel does not match exact committed targets and statuses")
    resolved_targets = []
    for target in targets:
        if not target.coordinates_verified:
            resolved_targets.append(target)
            continue
        row = registry.get(target.event_id or "")
        expected_provenance = f"event_junctions.GRCh38.tsv:sha256:{registry_hash}:{target.event_id}"
        if (row is None or target.target_id != target.event_id or target.assembly != "GRCh38"
                or row["chrom"] != target.chrom or row["strand"] != target.strand
                or row.get("transcript_version") != target.transcript_version
                or target.event_id not in EXPECTED_GOALS
                or target.provenance_digest != registry_provenance_digest(row, target.transcript_version or "")
                or target.coordinate_provenance != expected_provenance):
            raise ValueError(f"target {target.target_id} is not exactly resolved by frozen GRCh38 registry")
        resolved_targets.append(replace(target, sequence_sha256=row.get("sequence_sha256") or None))
    return tuple(resolved_targets)


def _pending(target: PanelTarget, status: str) -> tuple[PendingTarget, DesignabilityRow]:
    goal = EXPECTED_GOALS.get(target.event_id or target.target_id, ObjectiveDirection.INCREASE_INCLUSION)
    return (PendingTarget(target.target_id, target.locus, target.display_name, status, target.coordinate_provenance),
            DesignabilityRow(target.target_id, target.locus, goal, target.coordinates_verified,
                DesignabilityState.PENDING, status, (), None, None, None, None, 0,
                target.coordinate_provenance))


def _summary(target: PanelTarget, report: DesignReport) -> DesignabilityRow:
    abstention_rate = (sum(row.abstained for row in report.candidates) / len(report.candidates)) if report.candidates else 1.0
    supported = [row for row in report.candidates if not row.abstained
                 and row.goal_alignment == "supports_goal" and row.rank is not None]
    reasons = []
    if not supported: reasons.append("no goal-supporting candidates")
    if abstention_rate > MAX_ABSTENTION_RATE: reasons.append("abstention threshold exceeded")
    eligible = [row for row in supported if row.off_target_burden <= MAX_OFFTARGET_BURDEN]
    if supported and not eligible: reasons.append("off-target burden threshold exceeded")
    designable = not reasons and bool(eligible)
    best_rank = min((row.rank for row in eligible), default=None)
    best_rows = [row for row in eligible if row.rank == best_rank]
    tied = tuple(sorted(row.candidate.candidate_id for row in best_rows))
    return DesignabilityRow(target.target_id, target.locus, report.objective.direction, True,
        DesignabilityState.DESIGNABLE if designable else DesignabilityState.NON_DESIGNABLE,
        "prospective designs meet declared thresholds" if designable else "; ".join(reasons),
        tied if designable else (), best_rank if designable else None,
        max((x.directional_confidence for x in best_rows), default=None) if designable else None,
        min((x.off_target_burden for x in best_rows), default=None) if designable else None,
        round(abstention_rate, 6), len(report.pareto_candidate_ids), target.coordinate_provenance)


def run_panel(*, config_path: Path = DEFAULT_PANEL_CONFIG,
              registry_path: Path = DEFAULT_COORDINATE_REGISTRY,
              target_contexts: Mapping[str, RealTranscriptContext] | None = None,
              clock: Callable[[], datetime] | None = None,
              design_api: Callable[..., DesignReport] = design_for,
              design_options: dict | None = None,
              _test_only_allow_unapproved_engine: bool = False,
              _test_only_allow_unfrozen_config: bool = False) -> PanelReport:
    timestamp = (clock or (lambda: datetime.now(timezone.utc)))()
    if timestamp.tzinfo is None: raise ValueError("panel clock must return a timezone-aware datetime")
    timestamp = timestamp.astimezone(timezone.utc)
    targets = load_panel(config_path, registry_path,
                         _test_only_allow_unfrozen_config=_test_only_allow_unfrozen_config)
    if target_contexts is None:
        # Import lazily to avoid a module cycle: real_contexts uses the public binder.
        from .real_contexts import default_real_contexts
        contexts = default_real_contexts()
    else:
        contexts = dict(target_contexts)
    unknown = set(contexts) - {x.target_id for x in targets}
    if unknown: raise ValueError(f"contexts supplied for unknown panel targets: {sorted(unknown)}")
    reports: list[DesignReport] = []; pending = []; summaries = []
    options = dict(design_options or {})
    engine_api = f"{design_api.__module__}.{getattr(design_api, '__qualname__', design_api.__class__.__qualname__)}"
    if not _test_only_allow_unapproved_engine and (design_api is not design_for or options):
        raise ValueError("production panel requires the approved immutable engine manifest")
    if "clock" in options: raise ValueError("panel owns the frozen engine clock")
    options["clock"] = lambda: timestamp
    for target in targets:
        if not target.coordinates_verified:
            item, row = _pending(target, PENDING_COORDINATES_STATUS); pending.append(item); summaries.append(row); continue
        if not target.sequence_sha256:
            item, row = _pending(target, PENDING_UNCOMMITTED_CONTEXT_STATUS)
            pending.append(item); summaries.append(row); continue
        bound = contexts.get(target.target_id)
        if bound is None:
            item, row = _pending(target, PENDING_CONTEXT_STATUS); pending.append(item); summaries.append(row); continue
        definition = bound.definition
        if (_sequence_sha256(definition.context.sequence) != bound.binding_sequence_sha256
                or _hash(bound.genomic_to_transcript_mapping) != bound.mapping_sha256):
            raise ValueError(f"real context hash mismatch for {target.target_id}")
        if bound.binding_sequence_sha256 != target.sequence_sha256:
            raise ValueError(f"real context sequence does not match committed registry digest for {target.target_id}")
        if (bound.assembly != target.assembly or bound.chrom != target.chrom or bound.strand != target.strand
                or bound.transcript_version != target.transcript_version
                or definition.target_id != target.target_id or definition.objective.event_id != target.event_id
                or definition.locus != target.locus or definition.objective.locus != target.locus):
            raise ValueError(f"real context identity mismatch for {target.target_id}")
        expected = EXPECTED_GOALS[target.event_id]
        if definition.objective.direction is not expected: raise ValueError(f"therapeutic polarity mismatch for {target.target_id}")
        report = design_api(definition, **options)
        report.validate_for_rendering()
        backend = report.provenance["run_manifest"].get("backend", {})
        if (not _test_only_allow_unapproved_engine
                and (backend.get("identity_sha256") != APPROVED_BACKEND_IDENTITY_SHA256
                     or report.provenance.get("folding_backend_identity_sha256") != APPROVED_BACKEND_IDENTITY_SHA256)):
            raise ValueError("design report lacks approved backend receipt")
        if bound.mapping_sha256 != target.mapping_sha256:
            raise ValueError(f"real context mapping does not match committed digest for {target.target_id}")
        reported_definition = report.provenance["run_manifest"].get("target_definition")
        if _hash(reported_definition) != _hash(definition):
            raise ValueError(f"design report is not bound to supplied real context for {target.target_id}")
        if (report.target_id != target.target_id or report.objective.event_id != target.event_id
                or report.locus != target.locus or report.objective.locus != target.locus
                or report.objective.direction is not expected):
            raise ValueError(f"design report identity/polarity mismatch for {target.target_id}")
        if report.provenance.get("generated_at") != timestamp.isoformat():
            raise ValueError(f"design report ignored frozen panel clock for {target.target_id}")
        reports.append(report); summaries.append(_summary(target, report))
    manifest = {"panel_version": PANEL_VERSION,
        "engine_api": engine_api,
        "options": {k: v for k, v in options.items() if k != "clock"},
        "ordered_targets": [{**asdict(x), "expected_goal":
                             EXPECTED_GOALS.get(x.event_id or x.target_id,
                                                ObjectiveDirection.INCREASE_INCLUSION).value}
                            for x in targets],
        "bound_target_definitions": [{
            "target_id": report.target_id,
            "definition": contexts[report.target_id].definition,
            "sequence_provenance": contexts[report.target_id].sequence_provenance,
            "binding_sequence_sha256": contexts[report.target_id].binding_sequence_sha256,
            "genomic_to_transcript_mapping": contexts[report.target_id].genomic_to_transcript_mapping,
            "mapping_sha256": contexts[report.target_id].mapping_sha256,
            "chrom": contexts[report.target_id].chrom,
            "strand": contexts[report.target_id].strand,
        } for report in reports],
        # load_panel authenticated the exact bytes against this committed digest.
        "coordinate_registry_sha256": EXPECTED_REGISTRY_SHA256,
        "report_manifests": [{"target_id": x.target_id,
                              "run_manifest_sha256": x.provenance["run_manifest_sha256"]} for x in reports],
        "generated_at": timestamp.isoformat()}
    manifest_hash = _hash(manifest)
    provenance = {"panel_manifest": manifest, "panel_manifest_sha256": manifest_hash,
                  "prediction_status": PANEL_PREDICTION_STATUS,
                  "coordinate_policy": "fail closed: registry exact and real sequence bound"}
    return PanelReport("prospective-design-panel-report/2", PANEL_BANNER, PANEL_PREDICTION_STATUS,
                       timestamp.isoformat(), tuple(reports), tuple(pending), tuple(summaries),
                       provenance, manifest_hash)
