import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

from baselines.evaluation import identity_disjoint_component_splits
from benchmark import (BackendProvenance, BenchmarkConfig, BenchmarkInputs,
    DecisiveBenchmark, DeterministicPhysicsStub, FeatureOrigin, FeatureValue,
    FrozenFamilyRegistry, PhysicsCapacityEvidence, PhysicsCase, PhysicsPrediction,
    RegistryEntry, ThetaEvidence, TierState, write_report)
from benchmark.harness import FROZEN_CIS_BACKENDS, OCCLUSION_TRANSFORM_HASH, PROTOCOL_HASH
from physics_adapter import MagnitudeEvidence, OccupancySource


def capacity():
    return PhysicsCapacityEvidence.create(((-1.0,), (1.0,)), ridge=1.0)


def frozen_inputs(mode="perfect", families=4, theta=False):
    labels, features, cases, planted, entries = {}, {}, {}, {}, []
    version, weight_hash = FROZEN_CIS_BACKENDS["spliceai"]
    backend = BackendProvenance("spliceai", version, weight_hash, OCCLUSION_TRANSFORM_HASH)
    for family_index in range(families):
        for within in range(2):
            index = family_index * 2 + within; key = f"int-{index}"; sign = 1 if within else -1
            samples = (f"sample-{family_index}-{within}",)
            labels[key] = {"family": f"family-{family_index}", "study": f"study-{family_index}",
                "locus": f"locus-{family_index}", "sample": samples, "intervention": key,
                "chemistry": f"chemistry-{family_index}", "replicate_clusters": (f"rep-{index}",),
                "event": f"event-{index}", "observed": float(sign * 20), "magnitude_eligible": theta}
            ordinary = lambda value, source: FeatureValue(value, FeatureOrigin.NON_LABEL_INPUT, (source,))
            features[key] = {"affinity_dg": ordinary(0.0, "oligo-sequence"),
                "accessibility": ordinary(0.0, "reference-sequence"),
                "motif_ESE_fixture_mean": ordinary(0.0, "reference-sequence"),
                "splice_distance_donor_include": ordinary(0.0, "event-coordinate"),
                "kmer_k3_AAA": ordinary(0.0, "reference-sequence"),
                "length": ordinary(18.0, "oligo-sequence"), "dose": ordinary(1.0, "protocol-dose"),
                "chemistry_fixture": ordinary(1.0, "chemistry-registry"),
                "spliceai_endpoint_occlusion": FeatureValue(0.0, FeatureOrigin.FROZEN_CIS_BACKEND,
                                                             ("reference-sequence",), backend)}
            case = PhysicsCase(key)
            if theta:
                digest = case.input_hash()
                raw = MagnitudeEvidence("assay:fixture", MagnitudeEvidence.MATCHING_CRITERION,
                    OccupancySource.MATCHED_EXPERIMENTAL, digest, PROTOCOL_HASH)
                case = replace(case, theta_evidence=ThetaEvidence(raw, digest))
            cases[key] = case; planted[key] = sign if mode == "perfect" else -sign
            entries.append(RegistryEntry(key, f"family-{family_index}", f"study-{family_index}",
                                         f"locus-{family_index}", samples, f"chemistry-{family_index}"))
    registry = FrozenFamilyRegistry.freeze(entries)
    return BenchmarkInputs.freeze(labels, features, cases, registry), planted


def run_fixture(mode="perfect", families=4, theta=False):
    inputs, planted = frozen_inputs(mode, families, theta)
    return DecisiveBenchmark().run(inputs, DeterministicPhysicsStub(planted, capacity()))


class DecisiveBenchmarkTest(unittest.TestCase):
    def test_outcome_blindness_rejects_peeking_runner(self):
        inputs, _ = frozen_inputs()
        class Peeker:
            capacity = capacity()
            def predict(self, case): return PhysicsPrediction(1 if case["observed"] > 0 else -1)
        with self.assertRaisesRegex(ValueError, "outcome-blind"):
            DecisiveBenchmark().run(inputs, Peeker())

    def test_zero_margin_and_weakened_constants_rejected(self):
        with self.assertRaisesRegex(ValueError, "positive"): BenchmarkConfig(direction_margin=0)
        with self.assertRaisesRegex(ValueError, "protocol-locked"): BenchmarkConfig(neutral_threshold=0)
        config_path = Path(__file__).parents[1]/"src/benchmark/default_config.json"
        self.assertEqual(BenchmarkConfig.from_path(config_path), BenchmarkConfig())

    def test_connected_component_never_splits(self):
        inputs, _ = frozen_inputs(); rows = [dict(x) for x in inputs.labels.values()]
        rows[1]["sample"] = rows[0]["sample"]
        for _, _, held in identity_disjoint_component_splits(rows): self.assertEqual(0 in held, 1 in held)

    def test_physics_edf_underreport_rejected(self):
        inputs, planted = frozen_inputs(); proof = capacity(); forged = replace(proof, claimed_effective_df=proof.claimed_effective_df-.5)
        with self.assertRaisesRegex(ValueError, "under-report"):
            DecisiveBenchmark().run(inputs, DeterministicPhysicsStub(planted, forged))

    def test_blacklist_bypass_label_leakage_rejected_structurally(self):
        inputs, planted = frozen_inputs(); labels = {k: dict(v) for k,v in inputs.labels.items()}; features = {k: dict(v) for k,v in inputs.features.items()}
        features["int-0"]["innocent_x7"] = FeatureValue(20.0, FeatureOrigin.LABEL_DERIVED, ("endpoint-table-column-9",))
        forged = BenchmarkInputs.freeze(labels, features, dict(inputs.physics_cases), inputs.registry)
        with self.assertRaisesRegex(ValueError, "structurally"):
            DecisiveBenchmark().run(forged, DeterministicPhysicsStub(planted, capacity()))

    def test_unverified_occlusion_backend_rejected(self):
        inputs, planted = frozen_inputs(); labels = {k: dict(v) for k,v in inputs.labels.items()}; features = {k: dict(v) for k,v in inputs.features.items()}
        bad = BackendProvenance("spliceai", "latest", "0"*64, OCCLUSION_TRANSFORM_HASH)
        features["int-0"]["spliceai_endpoint_occlusion"] = FeatureValue(0.0, FeatureOrigin.FROZEN_CIS_BACKEND, ("reference",), bad)
        forged = BenchmarkInputs.freeze(labels, features, dict(inputs.physics_cases), inputs.registry)
        with self.assertRaisesRegex(ValueError, "unverified"):
            DecisiveBenchmark().run(forged, DeterministicPhysicsStub(planted, capacity()))

    def test_typed_theta_evidence_end_to_end(self):
        inputs, planted = frozen_inputs(theta=True)
        class MagnitudeRunner:
            capacity = capacity()
            def predict(self, case):
                return PhysicsPrediction(planted[case.fixture_id], 12.0, case.theta_evidence)
        result = DecisiveBenchmark().run(inputs, MagnitudeRunner())
        self.assertIn("typed-theta", result["physics"]["magnitude_status"])

    def test_forged_family_identity_rejected(self):
        inputs, planted = frozen_inputs(); labels = {k: dict(v) for k,v in inputs.labels.items()}; labels["int-0"]["family"] = "family-forged-9"
        forged = BenchmarkInputs.freeze(labels, {k:dict(v) for k,v in inputs.features.items()}, dict(inputs.physics_cases), inputs.registry)
        with self.assertRaisesRegex(ValueError, "forged family"):
            DecisiveBenchmark().run(forged, DeterministicPhysicsStub(planted, capacity()))

    def test_contradictory_tier_rejected(self):
        result = run_fixture(); contradictory = replace(result, tier=TierState.POWERED_CONFIRMATORY, claim_class="powered-generalization", verdict="pass")
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, "contradictory powered"):
                write_report(contradictory, Path(directory)/"x.json", Path(directory)/"x.txt")

    def test_report_mutation_rejected_and_hashes_emitted(self):
        result = run_fixture(); mutated = replace(result, verdict="pass")
        with tempfile.TemporaryDirectory() as directory:
            jp, tp = Path(directory)/"report.json", Path(directory)/"report.txt"
            with self.assertRaisesRegex(ValueError, "contradictory"): write_report(mutated, jp, tp)
            write_report(result, jp, tp); report = json.loads(jp.read_text())
            self.assertEqual(set(report["hashes"]), {"input","config","model","manifest","protocol"})

    def test_win_stays_feasibility_at_current_floor(self):
        result = run_fixture(); self.assertEqual(result.tier, TierState.N4_STRESS_TEST)
        self.assertEqual(result.verdict, "feasibility-only"); self.assertEqual(result.claim_class, "feasibility/descriptive")

    def test_wrong_physics_cannot_fabricate_win(self):
        self.assertEqual(run_fixture("wrong").verdict, "fail")


if __name__ == "__main__": unittest.main()
