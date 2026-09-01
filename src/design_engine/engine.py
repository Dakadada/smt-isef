"""Candidate generation, adapter-backed scoring, Pareto tiers, and reporting."""
from __future__ import annotations

from dataclasses import asdict, is_dataclass, replace
from datetime import datetime, timezone
from enum import Enum
import hashlib
import inspect
import json
from pathlib import Path
from typing import Callable, Mapping

from offtarget import EnumeratorConfig, enumerate_offtargets
from physics_adapter import (CoordinateConvention, DuplexOrientation,
    EffectorPhysicsAdapter, EffectorSpec, IntermolecularComplexSpec,
    IntermolecularPair, Interval, ModelConfig, MagnitudeEvidence,
    prediction_input_hash, protocol_hash)
from offtarget.alignment import reverse_complement

from .backends import DeterministicStubFolder
from .models import (DesignReport, GeneratedCandidate, GenerationConfig,
                     MagnitudeClaim, PROSPECTIVE_BANNER, PROSPECTIVE_STATUS,
                     ScoredCandidate, TargetDefinition)
from .targets import DEFAULT_EVENT_CONFIG, offline_target


def _canonical(value):
    if is_dataclass(value):
        return {key: _canonical(item) for key, item in asdict(value).items()}
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Mapping):
        return {str(key): _canonical(item) for key, item in sorted(value.items(), key=lambda x: str(x[0]))}
    if isinstance(value, (tuple, list, set, frozenset)):
        items = [_canonical(item) for item in value]
        return sorted(items, key=lambda x: json.dumps(x, sort_keys=True)) if isinstance(value, (set, frozenset)) else items
    if isinstance(value, Path):
        return str(value)
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return repr(value)


def _hash(value) -> str:
    payload = json.dumps(_canonical(value), sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode()).hexdigest()


def _file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _implementation_identity(obj) -> dict[str, object]:
    target = obj if inspect.isfunction(obj) or inspect.isclass(obj) else type(obj)
    module = inspect.getmodule(target)
    source_file = inspect.getsourcefile(target)
    try:
        source = inspect.getsource(target)
    except (OSError, TypeError):
        source = repr(target)
    identity: dict[str, object] = {
        "qualified_name": f"{target.__module__}.{target.__qualname__}",
        "source_sha256": hashlib.sha256(source.encode()).hexdigest(),
    }
    if source_file and Path(source_file).is_file():
        identity["module_artifact_sha256"] = _file_hash(Path(source_file))
    if not inspect.isfunction(obj) and not inspect.isclass(obj):
        identity["instance_config"] = _canonical(getattr(obj, "__dict__", {}))
    identity["identity_sha256"] = _hash(identity)
    return identity


def _package_artifacts(obj) -> dict[str, str]:
    module = inspect.getmodule(obj if inspect.isfunction(obj) or inspect.isclass(obj) else type(obj))
    module_file = Path(module.__file__).resolve() if module and module.__file__ else None
    if module_file is None:
        return {}
    package_dir = module_file.parent
    return {path.name: _file_hash(path) for path in sorted(package_dir.glob("*.py"))}


def _effector_spec(item: GeneratedCandidate, transcript_id: str, context_length: int,
                   shadow_size: int) -> EffectorSpec:
    target = Interval(item.start, item.end)
    pairs = tuple(IntermolecularPair(item.start + i, len(item.effector_sequence) - 1 - i)
                  for i in range(len(item.effector_sequence)))
    complex_spec = IntermolecularComplexSpec(
        target, CoordinateConvention.ZERO_BASED_HALF_OPEN, item.effector_sequence,
        transcript_id, item.candidate_id, DuplexOrientation.ANTIPARALLEL,
        item.start, item.chemistry, pairs)
    shadow = () if item.effector_class != "dCas13_gRNA" else (
        Interval(max(0, item.start - shadow_size), min(context_length, item.end + shadow_size)),)
    return EffectorSpec(item.candidate_id, target, complex_spec,
                        item.effector_class == "dCas13_gRNA", shadow)


def generate_candidates(target: TargetDefinition,
                        config: GenerationConfig = GenerationConfig()) -> tuple[GeneratedCandidate, ...]:
    candidates: list[GeneratedCandidate] = []
    classes = (("ASO", config.aso_lengths, config.aso_chemistries),
               ("dCas13_gRNA", config.grna_lengths, ("RNA-guided RNP",)),
               ("U7", config.u7_lengths, ("U7 snRNA cassette",)))
    for effector_class, lengths, chemistries in classes:
        for length in sorted(set(lengths)):
            for start in range(target.design_start, target.design_end - length + 1, config.stride):
                sequence = target.context.sequence[start:start + length].upper().replace("U", "T")
                for chemistry in chemistries:
                    identifier = f"{target.target_id}:{effector_class}:{chemistry}:{start}-{start+length}"
                    candidates.append(GeneratedCandidate(identifier, effector_class, chemistry,
                        start, start + length, sequence, reverse_complement(sequence)))
    candidates.sort(key=lambda x: (x.start, x.end, x.effector_class, x.chemistry, x.candidate_id))
    if len(candidates) <= config.max_candidates:
        return tuple(candidates)
    if config.max_candidates == 1:
        return (candidates[0],)
    indices = [i * (len(candidates) - 1) // (config.max_candidates - 1)
               for i in range(config.max_candidates)]
    return tuple(candidates[i] for i in indices)


def _manufacturability(item: GeneratedCandidate) -> float:
    base = {"ASO": 0.82, "U7": 0.67, "dCas13_gRNA": 0.58}[item.effector_class]
    chemistry = {"PMO": .08, "2MOE": .06, "OMe": .04, "LNA": -.02}.get(item.chemistry, 0)
    return round(max(0.0, min(1.0, base + chemistry - max(0, len(item.effector_sequence) - 20) * .01)), 6)


def _candidate_enumerator_config(item: GeneratedCandidate, base: EnumeratorConfig) -> EnumeratorConfig:
    chemistry = item.chemistry if item.effector_class == "ASO" else "OMe"
    return replace(base, chemistry=chemistry)


def _offtarget(item: GeneratedCandidate, target: TargetDefinition,
               base_config: EnumeratorConfig, enumerator) -> tuple[float, int]:
    if not target.offtarget_transcripts:
        return 0.0, 0
    hits = enumerator(item.effector_sequence, target.offtarget_transcripts,
                      _candidate_enumerator_config(item, base_config), target.offtarget_evidence)
    burden = sum(hit.recall_score * hit.expression_localization_weight for hit in hits)
    return round(burden, 6), len(hits)


def _vector(row: ScoredCandidate) -> tuple[float, float, float, float]:
    return (row.directional_confidence, -row.off_target_burden,
            row.manufacturability, row.span_robustness)


def _dominates(a: ScoredCandidate, b: ScoredCandidate) -> bool:
    return all(x >= y for x, y in zip(_vector(a), _vector(b))) and any(
        x > y for x, y in zip(_vector(a), _vector(b)))


def pareto_front(candidates: tuple[ScoredCandidate, ...]) -> tuple[str, ...]:
    ids = [x.candidate.candidate_id for x in candidates]
    if len(ids) != len(set(ids)):
        raise ValueError("candidate IDs must be unique before Pareto analysis")
    eligible = tuple(x for x in candidates if not x.abstained and x.goal_alignment == "supports_goal")
    # Sorting is serialization-only. Every member has the same semantic tier.
    return tuple(sorted(x.candidate.candidate_id for x in eligible
                        if not any(_dominates(other, x) for other in eligible if other is not x)))


def _dominance_tiers(candidates: tuple[ScoredCandidate, ...]) -> dict[str, int]:
    remaining = list(candidates)
    tiers: dict[str, int] = {}
    tier = 1
    while remaining:
        members = [row for row in remaining
                   if not any(_dominates(other, row) for other in remaining if other is not row)]
        for row in members:
            tiers[row.candidate.candidate_id] = tier
        member_ids = {id(row) for row in members}
        remaining = [row for row in remaining if id(row) not in member_ids]
        tier += 1
    return tiers


def pareto_tiers(candidates: tuple[ScoredCandidate, ...]) -> dict[str, int]:
    """Return nondominated equivalence tiers for eligible candidates only."""
    ids = [x.candidate.candidate_id for x in candidates]
    if len(ids) != len(set(ids)):
        raise ValueError("candidate IDs must be unique before Pareto analysis")
    eligible = tuple(x for x in candidates if not x.abstained and x.goal_alignment == "supports_goal")
    return _dominance_tiers(eligible)


class DesignEngine:
    def __init__(self, folder=None, generation_config: GenerationConfig | None = None,
                 *, clock: Callable[[], datetime] | None = None,
                 event_config: Path = DEFAULT_EVENT_CONFIG,
                 adapter_config: ModelConfig | None = None,
                 offtarget_config: EnumeratorConfig | None = None,
                 enumerator=enumerate_offtargets,
                 magnitude_evidence: Mapping[str, MagnitudeEvidence] | None = None,
                 magnitude_claims: Mapping[str, MagnitudeClaim] | None = None):
        self.folder = folder or DeterministicStubFolder()
        self.generation_config = generation_config or GenerationConfig()
        self.clock = clock or (lambda: datetime.now(timezone.utc))
        self.event_config = Path(event_config)
        self.adapter_config = adapter_config or ModelConfig()
        self.offtarget_config = offtarget_config or EnumeratorConfig(max_edits=4, max_indel_bases=1)
        self.enumerator = enumerator
        self.magnitude_evidence = dict(magnitude_evidence or {})
        self.magnitude_claims = dict(magnitude_claims or {})

    def design_for(self, target: str | TargetDefinition) -> DesignReport:
        definition = offline_target(target, self.event_config) if isinstance(target, str) else target
        generated = generate_candidates(definition, self.generation_config)
        if len({row.candidate_id for row in generated}) != len(generated):
            raise ValueError("generated candidate IDs must be unique")
        adapter = EffectorPhysicsAdapter(self.folder, self.adapter_config)
        provisional: list[ScoredCandidate] = []
        directions: dict[tuple[str, str, int, int], int | None] = {}
        for item in generated:
            spec = _effector_spec(item, definition.context.transcript_id,
                                  len(definition.context.sequence), self.generation_config.grna_shadow)
            evidence = self.magnitude_evidence.get(item.candidate_id)
            claim = self.magnitude_claims.get(item.candidate_id)
            if claim is not None and evidence is None:
                raise ValueError("a magnitude claim requires matching MagnitudeEvidence")
            if claim is not None and (claim.evidence_provenance != evidence.provenance
                                      or claim.input_hash != evidence.input_hash
                                      or claim.protocol_hash != evidence.protocol_hash):
                raise ValueError("magnitude claim and MagnitudeEvidence provenance do not match")
            prediction = adapter.predict(definition.context, spec, definition.event,
                                         magnitude_evidence=evidence)
            direction = prediction.primary_direction
            desired = definition.objective.direction.sign
            alignment = ("abstained" if direction is None else
                         "supports_goal" if direction == desired else "opposes_goal")
            confidence = 0.0 if direction is None else (0.9 if prediction.directional_concordance else 0.6)
            burden, hit_count = _offtarget(item, definition, self.offtarget_config, self.enumerator)
            uncertainty = list(prediction.warnings)
            if not prediction.directional_concordance and direction is not None:
                uncertainty.append("primary/secondary model directions are not concordant")
            uncertainty.append("off-target enumeration is recall-oriented; hit effects have unknown sign")
            key = (item.effector_class, item.chemistry, item.start, item.end)
            directions[key] = direction
            input_digest = prediction_input_hash(definition.context, spec, definition.event)
            protocol_digest = protocol_hash(adapter.config)
            provisional.append(ScoredCandidate(
                item, {1: "increases_inclusion", -1: "decreases_inclusion", None: "abstain"}[direction],
                definition.objective.direction, alignment, confidence, claim,
                input_digest, protocol_digest, direction is None, tuple(uncertainty),
                burden, hit_count, "unknown", _manufacturability(item), 0.0))
        scored: list[ScoredCandidate] = []
        for row in provisional:
            item = row.candidate
            direction = directions[(item.effector_class, item.chemistry, item.start, item.end)]
            neighbors = [value for (cls, chem, start, end), value in directions.items()
                         if cls == item.effector_class and chem == item.chemistry
                         and abs(start - item.start) <= self.generation_config.stride
                         and (start, end) != (item.start, item.end)]
            robustness = 0.0 if direction is None else round(
                sum(value == direction for value in neighbors) / max(1, len(neighbors)), 6)
            scored.append(replace(row, span_robustness=robustness))

        eligible = tuple(row for row in scored if not row.abstained and row.goal_alignment == "supports_goal")
        tiers = pareto_tiers(tuple(scored))
        front = set(pareto_front(tuple(scored)))
        scored = [replace(row, pareto=row.candidate.candidate_id in front,
                          pareto_tier=tiers.get(row.candidate.candidate_id),
                          rank=tiers.get(row.candidate.candidate_id)) for row in scored]
        # Candidate ID only stabilizes serialization within an equivalence tier; it never breaks the tier.
        scored.sort(key=lambda row: (row.rank is None, row.rank if row.rank is not None else 0,
                                     row.goal_alignment, row.candidate.candidate_id))

        backend_identity = _implementation_identity(self.folder)
        adapter_identity = _implementation_identity(EffectorPhysicsAdapter)
        adapter_identity["package_artifacts"] = _package_artifacts(EffectorPhysicsAdapter)
        adapter_identity["protocol_sha256"] = protocol_hash(self.adapter_config)
        adapter_identity["identity_sha256"] = _hash(adapter_identity)
        enumerator_identity = _implementation_identity(self.enumerator)
        enumerator_identity["package_artifacts"] = _package_artifacts(self.enumerator)
        enumerator_identity["identity_sha256"] = _hash(enumerator_identity)
        manifest = {
            "target_definition": definition,
            "generated_candidates": generated,
            "scored_candidates": scored,
            "pareto_candidate_ids": tuple(sorted(front)),
            "generation_config": self.generation_config,
            "adapter": adapter_identity,
            "offtarget": {
                "config": self.offtarget_config,
                "corpus_sha256": _hash(definition.offtarget_transcripts),
                "evidence_sha256": _hash(definition.offtarget_evidence),
            },
            "backend": backend_identity,
            "enumerator": {**enumerator_identity, "config": self.offtarget_config},
            "magnitude_evidence": self.magnitude_evidence,
            "magnitude_claims": self.magnitude_claims,
            "event_config_sha256": _file_hash(self.event_config),
        }
        provenance = {
            "run_manifest_sha256": _hash(manifest),
            "run_manifest": _canonical(manifest),
            "physics_adapter_identity_sha256": adapter_identity["identity_sha256"],
            "folding_backend_identity_sha256": backend_identity["identity_sha256"],
            "offtarget_enumerator_identity_sha256": enumerator_identity["identity_sha256"],
            "event_config_sha256": _file_hash(self.event_config),
            "generated_at": self.clock().isoformat(),
            "prediction_status": PROSPECTIVE_STATUS,
        }
        return DesignReport("design-report/2", PROSPECTIVE_BANNER, PROSPECTIVE_STATUS,
                            definition.target_id, definition.locus, definition.objective,
                            tuple(scored), tuple(sorted(front)), provenance)


def design_for(target: str | TargetDefinition, **engine_options) -> DesignReport:
    return DesignEngine(**engine_options).design_for(target)
