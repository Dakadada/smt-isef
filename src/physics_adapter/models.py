"""Three frozen model classes and the occupancy-weighted keystone estimand."""
from __future__ import annotations

from dataclasses import dataclass

from .folding import EnsembleRequest, Folder, verified_probability
from .types import (ClassResult, EffectorSpec, EventDefinition, ModelClass,
                    Probability, ProductiveState, TranscriptContext)

FROZEN_THETA_GRID = (0.0, 0.25, 0.5, 0.75, 1.0)


@dataclass(frozen=True)
class ModelConfig:
    theta_grid: tuple[float, ...] = FROZEN_THETA_GRID
    elongation_rates: tuple[tuple[float, float], ...] = ((20.0, .25), (40.0, .5), (60.0, .25))
    cotranscriptional_max_span: int = 200

    def __post_init__(self) -> None:
        if self.theta_grid != FROZEN_THETA_GRID:
            raise ValueError("theta grid is frozen by the benchmark protocol")
        if abs(sum(w for _, w in self.elongation_rates) - 1) > 1e-12:
            raise ValueError("elongation distribution weights must sum to one")
        if any(rate <= 0 or weight < 0 for rate, weight in self.elongation_rates):
            raise ValueError("elongation rates must be positive and weights non-negative")

    @classmethod
    def from_mapping(cls, config: dict[str, object]) -> "ModelConfig":
        """Build from JSON/YAML-decoded data, rejecting undeclared options."""
        allowed = {"theta_grid", "elongation_rates", "cotranscriptional_max_span"}
        unknown = set(config) - allowed
        if unknown:
            raise ValueError(f"unknown model configuration: {sorted(unknown)}")
        values = dict(config)
        if "theta_grid" in values:
            values["theta_grid"] = tuple(float(x) for x in values["theta_grid"])
        if "elongation_rates" in values:
            values["elongation_rates"] = tuple((float(x[0]), float(x[1])) for x in values["elongation_rates"])
        return cls(**values)


def _mix(a: Probability, b: Probability, theta: float) -> Probability:
    return Probability((1-theta)*a.estimate + theta*b.estimate,
                       (1-theta)*a.lower + theta*b.lower,
                       (1-theta)*a.upper + theta*b.upper)


def _sign(x: float, tolerance: float = 1e-12) -> int:
    return 1 if x > tolerance else (-1 if x < -tolerance else 0)


def robust_interval_direction(delta_intervals: dict[float, tuple[float, float]]) -> int | None:
    """Require every occupied-theta interval to exclude zero with one common sign."""
    signs = {1 if low > 0 else (-1 if high < 0 else 0)
             for theta, (low, high) in delta_intervals.items() if theta > 0}
    return next(iter(signs)) if len(signs) == 1 and 0 not in signs else None


class PhysicsModel:
    model_class: ModelClass
    role: str

    def __init__(self, folder: Folder, config: ModelConfig):
        self.folder, self.config = folder, config

    def evaluate(self, context: TranscriptContext, effector: EffectorSpec,
                 event: EventDefinition) -> ClassResult:
        n = len(context.sequence)
        all_intervals = event.accessible_elements + event.silencers + (effector.target,) + effector.shadow
        if any(iv.end > n for iv in all_intervals):
            raise ValueError("event/effector interval lies outside transcript")
        unbound = ProductiveState("effector_unbound", event.accessible_elements,
                                  event.silencers, (), ())
        unavailable = (effector.target,) + (effector.shadow if effector.bulky_rnp else ())
        bound = ProductiveState("aso_duplexed", event.accessible_elements,
                                event.silencers, (effector.target,), unavailable)
        # A splice contact inside occupied/shadowed sequence is explicitly impossible.
        occluded = any(a.overlaps(x) for a in event.accessible_elements for x in unavailable)
        pu = verified_probability(self.folder, self._request(context, effector, unbound))
        pb = Probability(0, 0, 0) if occluded else verified_probability(self.folder, self._request(context, effector, bound))
        probs = {t: _mix(pu, pb, t) for t in self.config.theta_grid}
        deltas = {t: p.estimate - pu.estimate for t, p in probs.items()}
        delta_intervals = {t: (p.lower - pu.upper, p.upper - pu.lower)
                           for t, p in probs.items()}
        directions = {t: (1 if low > 0 else (-1 if high < 0 else 0))
                      for t, (low, high) in delta_intervals.items()}
        # theta=0 is necessarily no-effect; robustness concerns every occupied point.
        robust = robust_interval_direction(delta_intervals)
        return ClassResult(self.model_class, self.role, unbound, bound, pu, pb,
                           probs, deltas, directions, robust)

    def _request(self, context: TranscriptContext, effector: EffectorSpec,
                 state: ProductiveState) -> EnsembleRequest:
        return EnsembleRequest(context, self.model_class, state.required_accessible,
            state.required_silencer_free, effector.bound_complex if state.target_paired_to_aso else None,
            state.unavailable_to_factors)


class EquilibriumFullTranscript(PhysicsModel):
    model_class = ModelClass.EQUILIBRIUM_FULL_TRANSCRIPT
    role = "PRIMARY"


class Cotranscriptional(PhysicsModel):
    model_class = ModelClass.COTRANSCRIPTIONAL
    role = "SECONDARY"

    def _request(self, context: TranscriptContext, effector: EffectorSpec,
                 state: ProductiveState) -> EnsembleRequest:
        r = super()._request(context, effector, state)
        return EnsembleRequest(**{**r.__dict__, "max_bp_span": self.config.cotranscriptional_max_span,
                                  "elongation_rates": self.config.elongation_rates})


class EclipProbingAblation(PhysicsModel):
    model_class = ModelClass.ECLIP_PROBING_ABLATION
    role = "ABLATION_FEATURE_VARIANT"

    def _request(self, context: TranscriptContext, effector: EffectorSpec,
                 state: ProductiveState) -> EnsembleRequest:
        r = super()._request(context, effector, state)
        return EnsembleRequest(**{**r.__dict__, "feature_variant": True})
