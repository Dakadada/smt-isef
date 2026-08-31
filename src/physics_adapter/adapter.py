"""Stable entry point for benchmarking and prospective design."""
from __future__ import annotations

import hashlib

from .folding import Folder
from .models import (Cotranscriptional, EclipProbingAblation,
                     EquilibriumFullTranscript, ModelConfig)
from .types import (EffectorSpec, EventDefinition, MagnitudeEvidence, ModelClass,
                    Prediction, TranscriptContext)


def _digest(parts: tuple[object, ...]) -> str:
    return hashlib.sha256(repr(parts).encode()).hexdigest()


def prediction_input_hash(context: TranscriptContext, effector: EffectorSpec,
                          event: EventDefinition) -> str:
    return _digest((context.transcript_id, context.sequence.upper(), tuple(sorted(context.probing.items())),
                    tuple(sorted(context.eclip_repressive_occupancy.items())), effector, event))


def protocol_hash(config: ModelConfig) -> str:
    return _digest(("effector-physics/v1", config.theta_grid, config.elongation_rates,
                    config.cotranscriptional_max_span))


class EffectorPhysicsAdapter:
    def __init__(self, folder: Folder, config: ModelConfig | None = None):
        self.config = config or ModelConfig()
        self.models = (EquilibriumFullTranscript(folder, self.config),
                       Cotranscriptional(folder, self.config),
                       EclipProbingAblation(folder, self.config))

    def predict(self, context: TranscriptContext, effector: EffectorSpec,
                event: EventDefinition, *, magnitude_evidence: MagnitudeEvidence | None = None) -> Prediction:
        evidence_valid = magnitude_evidence is not None
        if magnitude_evidence is not None:
            if magnitude_evidence.input_hash != prediction_input_hash(context, effector, event):
                raise ValueError("magnitude evidence does not match prediction inputs")
            if magnitude_evidence.protocol_hash != protocol_hash(self.config):
                raise ValueError("magnitude evidence does not match prediction protocol")
        results = {m.model_class: m.evaluate(context, effector, event) for m in self.models}
        primary = results[ModelClass.EQUILIBRIUM_FULL_TRANSCRIPT].robust_direction
        # Diagnostic only: no vote changes the primary result. Ablation is displayed,
        # but concordance of independent classes compares primary and secondary only.
        diagnostic = {k: v.robust_direction for k, v in results.items()}
        independent = [diagnostic[ModelClass.EQUILIBRIUM_FULL_TRANSCRIPT],
                       diagnostic[ModelClass.COTRANSCRIPTIONAL]]
        concordant = None not in independent and len(set(independent)) == 1
        warnings = []
        if not evidence_valid:
            warnings.append("theta is grid-bounded only: direction is primary; magnitude/calibration claims are prohibited")
        if primary is None:
            warnings.append("primary direction is not robust over the frozen theta grid")
        return Prediction(results, primary, concordant, diagnostic,
                          evidence_valid, self.config.theta_grid, tuple(warnings))
