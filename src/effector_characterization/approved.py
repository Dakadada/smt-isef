"""Derive score structure from the authenticated SMN2 engine run."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import inspect
import json
from pathlib import Path
from typing import Any

from design_engine import (PROSPECTIVE_BANNER, PROSPECTIVE_STATUS,
                           GeneratedCandidate, design_for)
from design_engine.engine import _manufacturability
from design_panel.panel import (_hash as panel_hash, _sequence_sha256,
                                committed_real_contexts, load_panel)

from .policy import CLAIM_BOUNDARY, P1_STATUS


TARGET_ID = "SMN2_exon7_inclusion"
CLASS_BASES = {"ASO": 0.82, "U7": 0.67, "dCas13_gRNA": 0.58}
ASO_CHEMISTRY_DELTAS = {"PMO": 0.08, "2MOE": 0.06, "OMe": 0.04, "LNA": -0.02}
LENGTH_PENALTY_THRESHOLD = 20
LENGTH_PENALTY_PER_NT = 0.01
SCORE_FORMULA = (
    "round(clamp(class_base + ASO_chemistry_delta - "
    "max(0, length-20)*0.01, 0, 1), 6)"
)


def _sha256_json(value: Any) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _stats(scores: list[float]) -> tuple[tuple[float, ...], float | None,
                                          float | None, float | None]:
    observed = tuple(sorted(scores))
    if not observed:
        return observed, None, None, None
    return observed, min(observed), max(observed), round(sum(observed) / len(observed), 6)


def _verify_disclosed_score_rule(*, _score_fn=_manufacturability,
                                 _candidate_type=GeneratedCandidate) -> str:
    """Exercise the implementation, including branches absent from the 21-nt tiling."""
    probes = []
    chemistries = {
        "ASO": tuple(ASO_CHEMISTRY_DELTAS),
        "U7": ("U7 snRNA cassette",),
        "dCas13_gRNA": ("RNA-guided RNP",),
    }
    for effector_class, class_base in CLASS_BASES.items():
        for chemistry in chemistries[effector_class]:
            for length in (18, 20, 21, 24):
                candidate = _candidate_type(f"rule-probe:{effector_class}:{chemistry}:{length}",
                    effector_class, chemistry, 0, length, "A" * length, "T" * length)
                delta = ASO_CHEMISTRY_DELTAS.get(chemistry, 0.0)
                penalty = max(0, length - LENGTH_PENALTY_THRESHOLD) * LENGTH_PENALTY_PER_NT
                expected = round(max(0.0, min(1.0, class_base + delta - penalty)), 6)
                actual = _score_fn(candidate)
                if actual != expected:
                    raise ValueError("engine implementation failed disclosed score-rule probe")
                probes.append({"effector_class": effector_class, "chemistry": chemistry,
                               "length": length, "expected": expected, "actual": actual})
    return _sha256_json(probes)


def _authenticated_context(*, _panel_loader=load_panel,
                           _context_loader=committed_real_contexts,
                           _sequence_hash=_sequence_sha256,
                           _mapping_hash=panel_hash):
    targets = {row.target_id: row for row in _panel_loader()}
    target = targets.get(TARGET_ID)
    context = _context_loader().get(TARGET_ID)
    if target is None or context is None or not target.coordinates_verified:
        raise ValueError("frozen committed SMN2 target and real context are required")
    definition = context.definition
    if (target.sequence_sha256 is None or target.mapping_sha256 is None
            or context.binding_sequence_sha256 != target.sequence_sha256
            or context.mapping_sha256 != target.mapping_sha256
            or _sequence_hash(definition.context.sequence) != context.binding_sequence_sha256
            or _mapping_hash(context.genomic_to_transcript_mapping) != context.mapping_sha256
            or definition.target_id != TARGET_ID or definition.locus != "SMN2"
            or context.assembly != target.assembly or context.chrom != target.chrom
            or context.strand != target.strand
            or context.transcript_version != target.transcript_version):
        raise ValueError("committed SMN2 context failed registry/hash authentication")
    return target, context


def derive_characterization(timestamp: datetime, *,
                            _context_fn=_authenticated_context,
                            _probe_fn=_verify_disclosed_score_rule,
                            _design_fn=design_for,
                            _score_fn=_manufacturability) -> dict[str, Any]:
    """Return canonical constructor fields; does not trust caller-held summaries."""
    if timestamp.tzinfo is None:
        raise ValueError("characterization clock must return a timezone-aware datetime")
    timestamp = timestamp.astimezone(timezone.utc)
    target, context = _context_fn()
    score_rule_probe_sha256 = _probe_fn()
    definition = context.definition
    report = _design_fn(definition, clock=lambda: timestamp)
    report.validate_for_rendering()
    config = report.provenance["run_manifest"]["generation_config"]
    grna_shadow = int(config["grna_shadow"])

    candidate_rows: list[dict[str, Any]] = []
    for scored in report.candidates:
        candidate = scored.candidate
        length = len(candidate.effector_sequence)
        base = CLASS_BASES[candidate.effector_class]
        delta = ASO_CHEMISTRY_DELTAS.get(candidate.chemistry, 0.0)
        penalty = round(max(0, length - LENGTH_PENALTY_THRESHOLD)
                        * LENGTH_PENALTY_PER_NT, 6)
        expected = round(max(0.0, min(1.0, base + delta - penalty)), 6)
        # Independent disclosed-component reconstruction must agree with the actual
        # deterministic engine output before any aggregate is emitted.
        if scored.manufacturability != expected or _score_fn(candidate) != expected:
            raise ValueError("engine score no longer matches disclosed score structure")
        shadow_start = shadow_end = None
        if candidate.effector_class == "dCas13_gRNA":
            shadow_start = max(0, candidate.start - grna_shadow)
            shadow_end = min(len(definition.context.sequence), candidate.end + grna_shadow)
        candidate_rows.append({"candidate_id": candidate.candidate_id,
            "effector_class": candidate.effector_class, "chemistry": candidate.chemistry,
            "start": candidate.start, "end": candidate.end, "length": length,
            "class_base": base, "chemistry_delta": delta, "length_penalty": penalty,
            "engine_score": scored.manufacturability,
            "grna_shadow_start": shadow_start, "grna_shadow_end": shadow_end})
    candidate_rows.sort(key=lambda row: (row["start"], row["end"], row["effector_class"],
                                         row["chemistry"], row["candidate_id"]))

    configured = {
        "ASO": tuple(int(x) for x in config["aso_lengths"]),
        "U7": tuple(int(x) for x in config["u7_lengths"]),
        "dCas13_gRNA": tuple(int(x) for x in config["grna_lengths"]),
    }
    chemistries = {
        "ASO": tuple(str(x) for x in config["aso_chemistries"]),
        "U7": ("U7 snRNA cassette",),
        "dCas13_gRNA": ("RNA-guided RNP",),
    }
    distributions: list[dict[str, Any]] = []
    for effector_class in ("ASO", "U7", "dCas13_gRNA"):
        for chemistry in chemistries[effector_class]:
            for length in configured[effector_class]:
                scores = [row["engine_score"] for row in candidate_rows
                          if row["effector_class"] == effector_class
                          and row["chemistry"] == chemistry and row["length"] == length]
                observed, minimum, maximum, mean = _stats(scores)
                distributions.append({"effector_class": effector_class,
                    "chemistry": chemistry, "length": length,
                    "candidate_count": len(scores), "observed_scores": observed,
                    "minimum": minimum, "maximum": maximum, "mean": mean})

    class_summaries: list[dict[str, Any]] = []
    for effector_class in ("ASO", "U7", "dCas13_gRNA"):
        scores = [row["engine_score"] for row in candidate_rows
                  if row["effector_class"] == effector_class]
        observed, minimum, maximum, mean = _stats(scores)
        limitation = None
        if not scores:
            lengths = ", ".join(str(value) for value in configured[effector_class])
            limitation = (
                f"No {effector_class} candidates were generated: committed context length "
                f"{len(definition.context.sequence)} nt is shorter than configured lengths "
                f"{lengths} nt. The configured {grna_shadow}-nt gRNA shadow therefore was "
                "not instantiated and no observed class distribution exists."
            )
        class_summaries.append({"effector_class": effector_class,
            "configured_lengths": configured[effector_class], "candidate_count": len(scores),
            "observed_scores": observed, "minimum": minimum, "maximum": maximum,
            "mean": mean,
            "availability": "observed-in-engine-output" if scores else "not-generated",
            "limitation": limitation})

    engine_path = Path(inspect.getsourcefile(_manufacturability) or "")
    provenance = {
        "prediction_status": PROSPECTIVE_STATUS,
        "engine_api": "design_engine.engine.design_for",
        "engine_report_sha256": _sha256_json(report.to_dict()),
        "engine_run_manifest_sha256": report.provenance["run_manifest_sha256"],
        "manufacturability_function_source_sha256": hashlib.sha256(
            inspect.getsource(_manufacturability).encode("utf-8")).hexdigest(),
        "score_rule_probe_sha256": score_rule_probe_sha256,
        "design_engine_module_sha256": _sha256_file(engine_path),
        "generation_config": dict(config),
        "committed_context_source": "design_panel.panel.committed_real_contexts",
        "coordinate_registry_sha256": report.provenance.get("event_config_sha256"),
        "panel_coordinate_provenance": target.coordinate_provenance,
        "binding_sequence_sha256": context.binding_sequence_sha256,
        "mapping_sha256": context.mapping_sha256,
        "sequence_provenance": context.sequence_provenance,
        "candidate_count": len(candidate_rows),
        "tiling_start_positions": sorted({row["start"] for row in candidate_rows}),
        "tiling_disclosure": (
            "The committed 21-nt window yields only start 0 for each generated 18/20-nt "
            "length. Counts above 2 within ASO reflect chemistry variants, not independent "
            "positions or biological replicates."
        ),
        "p1_activation_criteria": (
            "repeated interventions in at least 2 effector classes at at least 2 loci"
        ),
        "p1_activation_met": False,
    }
    manifest = {
        "generated_at": timestamp.isoformat(),
        "target_id": TARGET_ID,
        "candidate_rows": candidate_rows,
        "distributions": distributions,
        "class_summaries": class_summaries,
        "engine_run_manifest_sha256": report.provenance["run_manifest_sha256"],
        "binding_sequence_sha256": context.binding_sequence_sha256,
        "mapping_sha256": context.mapping_sha256,
    }
    provenance["characterization_manifest_sha256"] = _sha256_json(manifest)
    provenance["approved_derivation_source_sha256"] = hashlib.sha256(
        inspect.getsource(derive_characterization).encode("utf-8")).hexdigest()
    provenance["approved_derivation_module_sha256"] = _sha256_file(Path(__file__))
    return {
        "schema_version": "engine-effector-score-characterization/1",
        "banner": PROSPECTIVE_BANNER,
        "prediction_status": PROSPECTIVE_STATUS,
        "claim_boundary": CLAIM_BOUNDARY,
        "p1_status": P1_STATUS,
        "generated_at": timestamp.isoformat(),
        "target_id": TARGET_ID,
        "locus": "SMN2",
        "context_length": len(definition.context.sequence),
        "score_name": "manufacturability (engine heuristic)",
        "score_formula": SCORE_FORMULA,
        "class_bases": CLASS_BASES,
        "aso_chemistry_deltas": ASO_CHEMISTRY_DELTAS,
        "length_penalty_threshold": LENGTH_PENALTY_THRESHOLD,
        "length_penalty_per_nt": LENGTH_PENALTY_PER_NT,
        "grna_shadow_nt": grna_shadow,
        "candidates": tuple(candidate_rows),
        "distributions": tuple(distributions),
        "class_summaries": tuple(class_summaries),
        "provenance": provenance,
    }
