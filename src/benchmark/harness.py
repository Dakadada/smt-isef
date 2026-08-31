"""Authenticated decisive cis-vs-TRANS benchmark boundary."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
import hashlib
import json
from pathlib import Path
from types import MappingProxyType
from typing import Mapping

from baselines.evaluation import (
    BenchmarkHarness, EvaluationConfig, _cluster_uncertainty, _metrics,
    _validate_labels, identity_disjoint_component_splits,
)
from baselines.models import CapacityMatchedLinear
from physics_adapter import MagnitudeEvidence, ModelClass

PROTOCOL_SCHEMA = "decisive-cis-vs-trans-benchmark/v2"
PROTOCOL_SEED = 20260828
PROTOCOL_NEUTRAL_THRESHOLD = 5.0
PROTOCOL_SENSITIVITY_THRESHOLDS = (2.0, 10.0)
PROTOCOL_DIRECTION_MARGIN = 0.05
PROTOCOL_MAX_EDF = 8.0
PROTOCOL_CAPACITY_TOLERANCE_EDF = 1.0
PROTOCOL_MIN_CI_COMPONENTS = 8
PROTOCOL_POWERED_MIN_COMPONENTS = 10
PROTOCOL_HASH = hashlib.sha256(Path(__file__).parents[2].joinpath("docs/benchmark-protocol.md").read_bytes()).hexdigest()
OCCLUSION_TRANSFORM = "model-derived-sequence-occlusion-not-native-treatment-output"
OCCLUSION_TRANSFORM_HASH = hashlib.sha256(OCCLUSION_TRANSFORM.encode()).hexdigest()
FROZEN_CIS_BACKENDS = MappingProxyType({
    "spliceai": ("1.3.1", hashlib.sha256(b"spliceai-frozen-weights-v1").hexdigest()),
    "pangolin": ("1.0.2", hashlib.sha256(b"pangolin-frozen-weights-v1").hexdigest()),
})


def _plain(value):
    if isinstance(value, Enum): return value.value
    if isinstance(value, Mapping): return {str(k): _plain(v) for k, v in sorted(value.items(), key=lambda x: str(x[0]))}
    if isinstance(value, (tuple, list, set, frozenset)): return [_plain(v) for v in value]
    if hasattr(value, "__dataclass_fields__"): return {k: _plain(v) for k, v in asdict(value).items()}
    return value


def _hash(value):
    return hashlib.sha256(json.dumps(_plain(value), sort_keys=True, separators=(",", ":"), default=repr).encode()).hexdigest()


@dataclass(frozen=True)
class BenchmarkConfig:
    seed: int = PROTOCOL_SEED
    neutral_threshold: float = PROTOCOL_NEUTRAL_THRESHOLD
    sensitivity_thresholds: tuple[float, ...] = PROTOCOL_SENSITIVITY_THRESHOLDS
    direction_margin: float = PROTOCOL_DIRECTION_MARGIN
    max_effective_df: float = PROTOCOL_MAX_EDF
    capacity_tolerance_df: float = PROTOCOL_CAPACITY_TOLERANCE_EDF
    min_ci_families: int = PROTOCOL_MIN_CI_COMPONENTS
    powered_min_families: int = PROTOCOL_POWERED_MIN_COMPONENTS
    directional_power_validated: bool = False

    def __post_init__(self):
        expected = (PROTOCOL_SEED, PROTOCOL_NEUTRAL_THRESHOLD, PROTOCOL_SENSITIVITY_THRESHOLDS,
                    PROTOCOL_DIRECTION_MARGIN, PROTOCOL_MAX_EDF, PROTOCOL_CAPACITY_TOLERANCE_EDF,
                    PROTOCOL_MIN_CI_COMPONENTS, PROTOCOL_POWERED_MIN_COMPONENTS, False)
        actual = (self.seed, self.neutral_threshold, tuple(self.sensitivity_thresholds),
                  self.direction_margin, self.max_effective_df, self.capacity_tolerance_df,
                  self.min_ci_families, self.powered_min_families, self.directional_power_validated)
        if self.direction_margin <= 0: raise ValueError("protocol direction margin must be positive")
        if actual != expected: raise ValueError("decision constants are protocol-locked and non-lowerable")

    @property
    def digest(self): return _hash({"protocol": PROTOCOL_HASH, "config": self})

    @classmethod
    def from_path(cls, path):
        values = json.loads(Path(path).read_text())
        if values.pop("protocol_locked", None) is not True:
            raise ValueError("benchmark config must declare the protocol lock")
        values["sensitivity_thresholds"] = tuple(values.get("sensitivity_thresholds", PROTOCOL_SENSITIVITY_THRESHOLDS))
        return cls(**values)


class FeatureOrigin(str, Enum):
    NON_LABEL_INPUT = "non-label-input/v1"
    FROZEN_CIS_BACKEND = "frozen-cis-backend/v1"
    LABEL_DERIVED = "label-derived/v1"


@dataclass(frozen=True)
class BackendProvenance:
    backend: str
    version: str
    weight_hash: str
    occlusion_transform_hash: str

    def verify(self):
        if self.backend not in FROZEN_CIS_BACKENDS: raise ValueError("unverified cis backend")
        if (self.version, self.weight_hash) != FROZEN_CIS_BACKENDS[self.backend]:
            raise ValueError("unverified cis backend version or weights")
        if self.occlusion_transform_hash != OCCLUSION_TRANSFORM_HASH:
            raise ValueError("unverified occlusion transform")


@dataclass(frozen=True)
class FeatureValue:
    value: float
    origin: FeatureOrigin
    source_fields: tuple[str, ...]
    backend: BackendProvenance | None = None

    def verify(self):
        if isinstance(self.value, bool) or not isinstance(self.value, (int, float)): raise ValueError("baseline feature must be numeric")
        if not isinstance(self.origin, FeatureOrigin): raise TypeError("feature origin must be a typed FeatureOrigin")
        if not self.source_fields: raise ValueError("feature provenance requires source fields")
        if self.origin is FeatureOrigin.LABEL_DERIVED: raise ValueError("label-derived baseline feature rejected structurally")
        if self.origin is FeatureOrigin.FROZEN_CIS_BACKEND:
            if self.backend is None: raise ValueError("cis feature lacks verified backend provenance")
            self.backend.verify()
        elif self.backend is not None: raise ValueError("backend proof only applies to frozen cis features")


@dataclass(frozen=True)
class ThetaEvidence:
    evidence: MagnitudeEvidence
    case_input_hash: str
    protocol_hash: str = PROTOCOL_HASH

    def verify(self, case_hash):
        if self.protocol_hash != PROTOCOL_HASH or self.evidence.protocol_hash != PROTOCOL_HASH:
            raise ValueError("theta evidence is not bound to the frozen protocol")
        if self.case_input_hash != case_hash or self.evidence.input_hash != case_hash:
            raise ValueError("theta evidence is not hash-bound to physics inputs")


@dataclass(frozen=True)
class PhysicsCase:
    fixture_id: str
    context: object | None = None
    effector: object | None = None
    event: object | None = None
    theta_evidence: ThetaEvidence | None = None

    def input_hash(self):
        return _hash({"fixture_id": self.fixture_id, "context": repr(self.context),
                      "effector": repr(self.effector), "event": repr(self.event), "protocol": PROTOCOL_HASH})

    def __getitem__(self, key):
        raise TypeError(f"outcome-blind PhysicsCase does not expose {key!r}")


@dataclass(frozen=True)
class PhysicsPrediction:
    direction: int | None
    magnitude: float | None = None
    theta_evidence: ThetaEvidence | None = None
    model_class: str = ModelClass.EQUILIBRIUM_FULL_TRANSCRIPT.value


def _measure_edf(rows, ridge):
    if not rows or not rows[0] or any(len(row) != len(rows[0]) for row in rows): raise ValueError("physics EDF requires rectangular nonempty design")
    features = [{f"x{i}": value for i, value in enumerate(row)} for row in rows]
    return CapacityMatchedLinear(PROTOCOL_MAX_EDF, ridge).fit(features, [0.0] * len(rows)).effective_degrees_of_freedom


@dataclass(frozen=True)
class PhysicsCapacityEvidence:
    design_rows: tuple[tuple[float, ...], ...]
    ridge: float
    claimed_effective_df: float
    input_hash: str
    protocol_hash: str = PROTOCOL_HASH

    @classmethod
    def create(cls, design_rows, ridge=1.0):
        rows = tuple(tuple(map(float, row)) for row in design_rows); ridge = float(ridge)
        digest = _hash({"design": rows, "ridge": ridge, "protocol": PROTOCOL_HASH})
        return cls(rows, ridge, _measure_edf(rows, ridge), digest)

    def verified_edf(self):
        if self.protocol_hash != PROTOCOL_HASH or self.input_hash != _hash({"design": self.design_rows, "ridge": self.ridge, "protocol": PROTOCOL_HASH}):
            raise ValueError("physics EDF evidence hash/protocol mismatch")
        measured = _measure_edf(self.design_rows, self.ridge)
        if abs(measured - self.claimed_effective_df) > 1e-10: raise ValueError("physics EDF under-report or measure mismatch")
        return measured


@dataclass(frozen=True)
class RegistryEntry:
    intervention: str; family: str; study: str; locus: str; samples: tuple[str, ...]; chemistry: str


@dataclass(frozen=True)
class FrozenFamilyRegistry:
    entries: tuple[RegistryEntry, ...]
    manifest_hash: str

    @classmethod
    def freeze(cls, entries):
        ordered = tuple(sorted(entries, key=lambda x: x.intervention))
        return cls(ordered, _hash({"registry": ordered, "protocol": PROTOCOL_HASH}))

    def verify(self):
        if self.manifest_hash != _hash({"registry": self.entries, "protocol": PROTOCOL_HASH}): raise ValueError("family registry manifest hash mismatch")
        if len({x.intervention for x in self.entries}) != len(self.entries): raise ValueError("family registry duplicate intervention")
        return {x.intervention: x for x in self.entries}


@dataclass(frozen=True)
class BenchmarkInputs:
    labels: Mapping[str, Mapping]
    features: Mapping[str, Mapping[str, FeatureValue]]
    physics_cases: Mapping[str, PhysicsCase]
    registry: FrozenFamilyRegistry
    digest: str

    @classmethod
    def freeze(cls, labels, features, physics_cases, registry):
        labels_f = MappingProxyType({k: MappingProxyType(dict(v)) for k, v in labels.items()})
        features_f = MappingProxyType({k: MappingProxyType(dict(v)) for k, v in features.items()})
        cases_f = MappingProxyType(dict(physics_cases))
        digest = _hash({"labels": labels, "features": features, "cases": physics_cases,
                        "registry": registry.manifest_hash, "protocol": PROTOCOL_HASH})
        return cls(labels_f, features_f, cases_f, registry, digest)

    def verify(self):
        digest = _hash({"labels": self.labels, "features": self.features, "cases": self.physics_cases,
                        "registry": self.registry.manifest_hash, "protocol": PROTOCOL_HASH})
        if digest != self.digest: raise ValueError("benchmark input hash mismatch")


class AdapterPhysicsRunner:
    def __init__(self, adapter, capacity): self.adapter, self.capacity = adapter, capacity
    def predict(self, case):
        result = self.adapter.predict(case.context, case.effector, case.event,
            magnitude_evidence=case.theta_evidence.evidence if case.theta_evidence else None)
        return PhysicsPrediction(result.primary_direction, theta_evidence=case.theta_evidence if result.magnitude_claim_allowed else None)


class DeterministicPhysicsStub:
    def __init__(self, predictions, capacity): self.predictions, self.capacity = MappingProxyType(dict(predictions)), capacity
    def predict(self, case): return PhysicsPrediction(self.predictions.get(case.fixture_id))


class TierState(str, Enum):
    INSUFFICIENT = "INSUFFICIENT"
    N2_REPRODUCIBILITY = "N=2"
    N3_FEASIBILITY = "N=3-FEASIBILITY"
    N4_STRESS_TEST = "N=4-FEASIBILITY"
    ABOVE_N4_DESCRIPTIVE = "ABOVE-N=4-DESCRIPTIVE"
    POWERED_CONFIRMATORY = "POWERED-CONFIRMATORY"


def _tier(n, loci, interventions):
    if n < 2 or interventions < 2: return TierState.INSUFFICIENT
    if n == 2: return TierState.N2_REPRODUCIBILITY
    if n == 3: return TierState.N3_FEASIBILITY if loci >= 2 and interventions >= 3 else TierState.INSUFFICIENT
    if n == 4: return TierState.N4_STRESS_TEST if loci >= 2 and interventions >= 5 else TierState.INSUFFICIENT
    return TierState.ABOVE_N4_DESCRIPTIVE


def _freeze_tree(value):
    if isinstance(value, dict): return MappingProxyType({k: _freeze_tree(v) for k, v in value.items()})
    if isinstance(value, list): return tuple(_freeze_tree(v) for v in value)
    return value


def _thaw(value):
    if isinstance(value, Mapping): return {k: _thaw(v) for k, v in value.items()}
    if isinstance(value, tuple): return [_thaw(v) for v in value]
    if isinstance(value, Enum): return value.value
    return value


@dataclass(frozen=True)
class BenchmarkResult:
    tier: TierState; claim_class: str; verdict: str; observed_margin: float | None
    coverage: float; component_count: int; ci_present: bool; payload: Mapping
    input_hash: str; config_hash: str; model_hash: str; manifest_hash: str

    def __getitem__(self, key): return self.to_dict()[key]
    def to_dict(self):
        data = _thaw(self.payload)
        data.update({"schema": PROTOCOL_SCHEMA, "tier": self.tier.value, "claim_class": self.claim_class,
                     "hashes": {"input": self.input_hash, "config": self.config_hash,
                     "model": self.model_hash, "manifest": self.manifest_hash, "protocol": PROTOCOL_HASH}})
        data["gate"]["verdict"] = self.verdict
        return data


def _derived_decision(tier, margin, coverage, count, ci_present):
    win = coverage == 1 and margin is not None and margin >= PROTOCOL_DIRECTION_MARGIN
    powered = tier is TierState.POWERED_CONFIRMATORY
    if powered and (count < PROTOCOL_POWERED_MIN_COMPONENTS or not ci_present): raise ValueError("contradictory powered tier/count/CI state")
    claim = "powered-generalization" if powered else "feasibility/descriptive"
    verdict = "pass" if powered and win else ("feasibility-only" if win and tier is not TierState.INSUFFICIENT else "fail")
    return claim, verdict


def validate_result(result):
    claim, verdict = _derived_decision(result.tier, result.observed_margin, result.coverage, result.component_count, result.ci_present)
    if (claim, verdict) != (result.claim_class, result.verdict): raise ValueError("contradictory tier/counts/CI/claim/verdict result")
    payload = _thaw(result.payload)
    if payload["counts"]["identity_components"] != result.component_count: raise ValueError("result count invariant failed")
    if payload["gate"]["observed_margin"] != result.observed_margin: raise ValueError("result margin invariant failed")
    if payload["physics"]["coverage"] != result.coverage: raise ValueError("result coverage invariant failed")


class DecisiveBenchmark:
    def __init__(self, config=None): self.config = config or BenchmarkConfig()

    def run(self, inputs, physics_runner):
        if not isinstance(inputs, BenchmarkInputs): raise TypeError("decisive benchmark requires frozen typed BenchmarkInputs")
        inputs.verify(); registry = inputs.registry.verify()
        labels = {k: dict(v) for k, v in inputs.labels.items()}; _validate_labels(labels); keys = set(labels)
        if set(inputs.features) != keys or set(inputs.physics_cases) != keys or set(registry) != keys: raise ValueError("input and registry coverage differ")
        for key, row in labels.items():
            entry = registry[key]
            authenticated = (entry.family, entry.study, entry.locus, entry.samples, entry.intervention, entry.chemistry)
            claimed = (row["family"], row["study"], row["locus"], tuple(row["sample"]), row["intervention"], row["chemistry"])
            if authenticated != claimed: raise ValueError("forged family/study/sample identity rejected")
        numeric = {}
        for key, row in inputs.features.items():
            numeric[key] = {}
            for name, feature in row.items():
                if not isinstance(feature, FeatureValue): raise TypeError("baseline features require typed provenance")
                feature.verify(); numeric[key][name] = float(feature.value)
        if not isinstance(getattr(physics_runner, "capacity", None), PhysicsCapacityEvidence): raise ValueError("physics requires authenticated EDF evidence")
        physics_edf = physics_runner.capacity.verified_edf()
        records = [labels[k] for k in labels]; splits = identity_disjoint_component_splits(records)
        if len(splits) < 2: raise ValueError("at least two identity-disjoint connected components required")
        # All weighting, nested selection, CIs, and counts use the same
        # authenticated connected-component identity. Registry family names
        # remain provenance only and can never inflate independent replication.
        for component_id, (_, test, _) in enumerate(splits):
            for row in test:
                labels[row["intervention"]]["family"] = f"identity-component-{component_id}"
        records = [labels[k] for k in labels]
        splits = identity_disjoint_component_splits(records)
        evaluation = EvaluationConfig(neutral_threshold=PROTOCOL_NEUTRAL_THRESHOLD,
            sensitivity_thresholds=PROTOCOL_SENSITIVITY_THRESHOLDS, min_ci_families=PROTOCOL_MIN_CI_COMPONENTS,
            max_effective_df=PROTOCOL_MAX_EDF, seed=PROTOCOL_SEED)
        baselines = BenchmarkHarness(evaluation).run(labels, numeric, physics_effective_degrees_of_freedom=physics_edf)
        physics_rows, folds = [], []
        for fold_id, (_, test, _) in enumerate(splits):
            fold_ids = []
            for row in test:
                case = inputs.physics_cases[row["intervention"]]
                if not isinstance(case, PhysicsCase): raise TypeError("physics cases must be typed and outcome-blind")
                if bool(row["magnitude_eligible"]) != bool(case.theta_evidence): raise ValueError("magnitude eligibility requires typed theta evidence")
                if case.theta_evidence: case.theta_evidence.verify(case.input_hash())
                try: prediction = physics_runner.predict(case)
                except (KeyError, TypeError, AttributeError) as exc: raise ValueError("physics runner violated outcome-blind boundary") from exc
                if not isinstance(prediction, PhysicsPrediction): raise TypeError("physics runner must return typed PhysicsPrediction")
                if prediction.direction not in (-1, 0, 1, None): raise ValueError("physics direction invalid")
                if prediction.magnitude is not None and (case.theta_evidence is None or prediction.theta_evidence != case.theta_evidence):
                    raise ValueError("physics magnitude lacks matching typed theta evidence")
                predicted = None if prediction.direction is None else prediction.direction * PROTOCOL_NEUTRAL_THRESHOLD
                physics_rows.append(dict(row, predicted=predicted, predicted_magnitude=prediction.magnitude))
                fold_ids.append(row["intervention"])
            folds.append({"held_out_component": fold_id, "interventions": fold_ids, "identity_disjoint": True})
        covered = [r for r in physics_rows if r["predicted"] is not None]
        direction = _metrics(covered, PROTOCOL_NEUTRAL_THRESHOLD) if covered else {"directional_balanced_accuracy": None, "directional_concordance": None}
        direction.update({"sensitivities": {str(t): _metrics(covered, t) for t in PROTOCOL_SENSITIVITY_THRESHOLDS} if covered else {},
                          "coverage": len(covered)/len(physics_rows), "effective_degrees_of_freedom": physics_edf,
                          "folds": folds, "mae": None, "calibration_slope": None,
                          "magnitude_status": "withheld unless typed-theta metric floors pass"})
        for name, result in baselines.items():
            for fold in result["folds"]:
                difference = abs(fold["baseline_effective_degrees_of_freedom"] - physics_edf)
                if difference > PROTOCOL_CAPACITY_TOLERANCE_EDF + 1e-8: raise ValueError(f"matched-capacity violation for {name}")
            result["mae"] = result["calibration_slope"] = None
        selected_rows, selections = [], []
        for fold_id, (_, _, held) in enumerate(splits):
            training_keys = [r["intervention"] for i, r in enumerate(records) if i not in held]
            test_keys = [r["intervention"] for i, r in enumerate(records) if i in held]
            try:
                inner = BenchmarkHarness(evaluation).run({k: labels[k] for k in training_keys}, {k: numeric[k] for k in training_keys}, physics_effective_degrees_of_freedom=physics_edf)
            except ValueError as error:
                selections.append({"held_out_component": fold_id, "selected_baseline": None, "status": f"not-selectable: {error}"}); continue
            chosen = max(inner, key=lambda n: inner[n]["directional_concordance"])
            predictions = {r["intervention"]: r for r in baselines[chosen]["predictions"]}
            selected_rows.extend(predictions[k] for k in test_keys)
            selections.append({"held_out_component": fold_id, "selected_baseline": chosen, "status": "selected-on-training-components-only"})
        selected = _metrics(selected_rows, PROTOCOL_NEUTRAL_THRESHOLD) if len(selected_rows) == len(labels) else {"directional_balanced_accuracy": None, "directional_concordance": None}
        selected.update({"fold_selections": selections, "coverage": len(selected_rows)/len(labels)})
        best = selected["directional_concordance"]
        margin = None if direction["directional_concordance"] is None or best is None else direction["directional_concordance"] - best
        count = len(splits); loci = len({x.locus for x in registry.values()}); tier = _tier(count, loci, len(labels)); ci = count >= PROTOCOL_MIN_CI_COMPONENTS
        direction["cluster_uncertainty"] = _cluster_uncertainty(covered, PROTOCOL_NEUTRAL_THRESHOLD, PROTOCOL_SEED) if ci else None
        direction["ci_status"] = "component-floor-met" if ci else "withheld-fewer-than-8-independent-components"
        claim, verdict = _derived_decision(tier, margin, direction["coverage"], count, ci)
        gate = {"verdict": verdict, "best_baseline": "outer-fold-selected-best-baseline",
                "required_directional_concordance_margin": PROTOCOL_DIRECTION_MARGIN, "observed_margin": margin,
                "margin_met": direction["coverage"] == 1 and margin is not None and margin >= PROTOCOL_DIRECTION_MARGIN,
                "powered_directional_simulation_validated": False}
        payload = _freeze_tree({"seed": PROTOCOL_SEED, "endpoint": "direction-of-treated-minus-control-delta-PSI",
            "gate": gate, "physics": direction, "baselines": baselines, "outer_fold_selected_baseline": selected,
            "counts": {"identity_components": count, "loci": loci, "interventions": len(labels)},
            "cluster_ci": "available" if ci else "withheld: fewer than 8 independent components",
            "honest_limitations": ["Direction is the sole primary endpoint.", "Confirmatory claims require separately frozen directional power validation."]})
        model_hash = _hash({"capacity": physics_runner.capacity, "predictions": [r["predicted"] for r in physics_rows]})
        result = BenchmarkResult(tier, claim, verdict, margin, direction["coverage"], count, ci, payload,
                                 inputs.digest, self.config.digest, model_hash, inputs.registry.manifest_hash)
        validate_result(result); return result


def _human(data):
    return ("DECISIVE CIS-vs-TRANS BENCHMARK\n\n" + f"Tier: {data['tier']}\nClaim class: {data['claim_class']}\n"
            + f"Gate verdict: {data['gate']['verdict']}\nObserved direction margin: {data['gate']['observed_margin']}\n"
            + f"Independent components: {data['counts']['identity_components']}\n")


def write_report(result, json_path, text_path):
    if not isinstance(result, BenchmarkResult): raise TypeError("report writer accepts only immutable BenchmarkResult")
    validate_result(result)
    claim, verdict = _derived_decision(result.tier, result.observed_margin, result.coverage, result.component_count, result.ci_present)
    if (claim, verdict) != (result.claim_class, result.verdict): raise ValueError("report mutation/contradiction detected")
    data = result.to_dict(); data["gate"]["verdict"] = verdict
    Path(json_path).write_text(json.dumps(data, indent=2, sort_keys=True) + "\n"); Path(text_path).write_text(_human(data))
