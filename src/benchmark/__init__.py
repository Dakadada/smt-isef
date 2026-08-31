"""Decisive, leakage-resistant cis-vs-trans benchmark orchestration."""

from .harness import (AdapterPhysicsRunner, BackendProvenance, BenchmarkConfig,
                      BenchmarkInputs, BenchmarkResult, DecisiveBenchmark,
                      DeterministicPhysicsStub, FeatureOrigin, FeatureValue,
                      FrozenFamilyRegistry, PhysicsCapacityEvidence, PhysicsCase,
                      PhysicsPrediction, RegistryEntry, ThetaEvidence, TierState,
                      write_report)

__all__ = ["AdapterPhysicsRunner", "BackendProvenance", "BenchmarkConfig",
           "BenchmarkInputs", "BenchmarkResult", "DecisiveBenchmark",
           "DeterministicPhysicsStub", "FeatureOrigin", "FeatureValue",
           "FrozenFamilyRegistry", "PhysicsCapacityEvidence", "PhysicsCase",
           "PhysicsPrediction", "RegistryEntry", "ThetaEvidence", "TierState",
           "write_report"]
