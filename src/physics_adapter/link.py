"""Capacity contract for the cross-fitted PSI link and matched baselines."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math
from typing import Mapping, Sequence


class CapacityError(ValueError): pass


class Baseline(str, Enum):
    UNTREATED_PROPENSITY = "untreated_propensity"
    OCCLUSION_IMPORTANCE = "propensity_plus_occlusion_importance"
    MEASURABLE_FEATURES = "measurable_features"
    MATCHED_STACK = "matched_capacity_stack"


@dataclass(frozen=True)
class LinkCapacity:
    interior_knots: tuple[float, ...] = (.25, .5, .75)
    max_group_levels: int = 12
    min_observations_per_group: int = 2
    min_observations_per_degree_of_freedom: int = 4
    min_group_label_fraction: float = .2
    ridge_precision: float = 4.0

    def __post_init__(self) -> None:
        if len(self.interior_knots) > 3 or tuple(sorted(self.interior_knots)) != self.interior_knots:
            raise CapacityError("monotone PSI link permits at most three ordered knots")
        if self.ridge_precision <= 0:
            raise CapacityError("group intercepts require a zero-centering penalty")
        if self.min_observations_per_degree_of_freedom < 2:
            raise CapacityError("parameter-to-sample ratio is too permissive")
        if not 0 < self.min_group_label_fraction < .5:
            raise CapacityError("group label fraction must be between zero and one half")

    @property
    def global_degrees_of_freedom(self) -> int:
        return 2 + len(self.interior_knots)


@dataclass(frozen=True)
class CapacityConstrainedLink:
    capacity: LinkCapacity = LinkCapacity()
    forbidden_fields: tuple[str, ...] = ("locus_id", "oligo_id", "study_by_chemistry", "label")

    def validate_training_design(self, rows: Sequence[Mapping[str, object]]) -> None:
        if not rows:
            raise CapacityError("training rows cannot be empty")
        present = set().union(*(r.keys() for r in rows))
        bad = present.intersection(self.forbidden_fields)
        if bad:
            raise CapacityError(f"forbidden lookup-table fields: {sorted(bad)}")
        required = {"p", "psi", "study_family", "chemistry", "fold"}
        missing = required - present
        if missing:
            raise CapacityError(f"missing training-design fields: {sorted(missing)}")
        counts_by_field: dict[str, dict[object, int]] = {}
        for field in ("study_family", "chemistry"):
            counts: dict[object, int] = {}
            for row in rows:
                value = row.get(field)
                counts[value] = counts.get(value, 0) + 1
            if len(counts) > self.capacity.max_group_levels:
                raise CapacityError(f"too many {field} levels for frozen capacity")
            if min(counts.values()) < self.capacity.min_observations_per_group:
                raise CapacityError(f"{field} cannot identify individual rows/labels")
            counts_by_field[field] = counts

        studies = set(counts_by_field["study_family"])
        chemistries = set(counts_by_field["chemistry"])
        study_chems = {s: {r["chemistry"] for r in rows if r["study_family"] == s}
                       for s in studies}
        chem_studies = {c: {r["study_family"] for r in rows if r["chemistry"] == c}
                        for c in chemistries}
        if any(len(values) < 2 for values in study_chems.values()) or any(
                len(values) < 2 for values in chem_studies.values()):
            raise CapacityError("study x chemistry design must be crossed, not sparse identities")
        seen_studies, seen_chems = set(), set()
        frontier_studies = {next(iter(studies))}
        while frontier_studies:
            seen_studies |= frontier_studies
            new_chems = set().union(*(study_chems[s] for s in frontier_studies)) - seen_chems
            seen_chems |= new_chems
            frontier_studies = set().union(*(chem_studies[c] for c in new_chems)) - seen_studies if new_chems else set()
        if seen_studies != studies or seen_chems != chemistries:
            raise CapacityError("study x chemistry design must form one connected component")

        for field, levels in (("study_family", studies), ("chemistry", chemistries)):
            for level in levels:
                labels = [float(r["psi"]) for r in rows if r[field] == level]
                prevalence = sum(labels) / len(labels)
                floor = self.capacity.min_group_label_fraction
                if prevalence < floor or prevalence > 1 - floor:
                    raise CapacityError(f"{field} group is label-pure or near-label-pure")

        folds = {r["fold"] for r in rows}
        if len(folds) < 2:
            raise CapacityError("group-disjoint cross-fitting requires at least two folds")
        for study_name in studies:
            if len({r["fold"] for r in rows if r["study_family"] == study_name}) != 1:
                raise CapacityError("study groups must be disjoint across cross-fitting folds")
        for chemistry in chemistries:
            if len({r["fold"] for r in rows if r["chemistry"] == chemistry}) < 2:
                raise CapacityError("each chemistry must be represented across cross-fitting folds")
        for fold in folds:
            labels = {float(r["psi"]) >= .5 for r in rows if r["fold"] == fold}
            if len(labels) < 2:
                raise CapacityError("each cross-fitting fold must contain both label classes")

        fitted_df = (self.capacity.global_degrees_of_freedom + len(studies) - 1
                     + len(chemistries) - 1)
        minimum = fitted_df * self.capacity.min_observations_per_degree_of_freedom
        if len(rows) < minimum:
            raise CapacityError(f"insufficient observations per fitted degree of freedom: need {minimum}")

    def baseline_capacity(self, baseline: Baseline) -> Mapping[str, object]:
        """The required matched budget; feature content differs, capacity does not."""
        return {"baseline": baseline.value,
                "global_df": self.capacity.global_degrees_of_freedom,
                "max_group_levels": self.capacity.max_group_levels,
                "min_observations_per_degree_of_freedom": self.capacity.min_observations_per_degree_of_freedom,
                "ridge_precision": self.capacity.ridge_precision}

    def fit(self, rows: Sequence[Mapping[str, object]], *, iterations: int = 2000,
            learning_rate: float = .03) -> "FittedLink":
        """Fit a deterministic monotone logistic spline with ridge group effects.

        Required fields are ``p`` in [0,1], binary/fractional ``psi`` in [0,1],
        ``study_family`` and ``chemistry``. Identity-like groups fail validation.
        Positive coefficients are parameterised with softplus, guaranteeing the
        fitted map is monotone rather than repairing it after outcome access.
        """
        self.validate_training_design(rows)
        studies = sorted({str(r["study_family"]) for r in rows})
        chemistries = sorted({str(r["chemistry"]) for r in rows})
        data = []
        for r in rows:
            x, y = float(r["p"]), float(r["psi"])
            if not (0 <= x <= 1 and 0 <= y <= 1):
                raise ValueError("p and psi must lie in [0,1]")
            data.append((x, y, str(r["study_family"]), str(r["chemistry"])))
        intercept, raw_slope = 0.0, 0.0
        raw_hinges = [0.0] * len(self.capacity.interior_knots)
        study = {x: 0.0 for x in studies}
        chemistry = {x: 0.0 for x in chemistries}
        n = float(len(data))
        for _ in range(iterations):
            gi = grs = 0.0
            grh = [0.0] * len(raw_hinges)
            gs = {x: 0.0 for x in studies}; gc = {x: 0.0 for x in chemistries}
            slope = _softplus(raw_slope); hinges = [_softplus(x) for x in raw_hinges]
            for x, y, s, c in data:
                z = intercept + slope*x + sum(a*max(0.0, x-k) for a,k in zip(hinges, self.capacity.interior_knots)) + study[s] + chemistry[c]
                error = _sigmoid(z) - y
                gi += error; grs += error*x*_sigmoid(raw_slope)
                for j, knot in enumerate(self.capacity.interior_knots):
                    grh[j] += error*max(0.0, x-knot)*_sigmoid(raw_hinges[j])
                gs[s] += error; gc[c] += error
            intercept -= learning_rate*gi/n; raw_slope -= learning_rate*grs/n
            raw_hinges = [v-learning_rate*g/n for v,g in zip(raw_hinges,grh)]
            for key in study:
                study[key] -= learning_rate*(gs[key]/n + self.capacity.ridge_precision*study[key]/n)
            for key in chemistry:
                chemistry[key] -= learning_rate*(gc[key]/n + self.capacity.ridge_precision*chemistry[key]/n)
            # Zero centering makes the global intercept identifiable.
            sm = sum(study.values())/len(study); cm = sum(chemistry.values())/len(chemistry)
            study = {k: v-sm for k,v in study.items()}; chemistry = {k: v-cm for k,v in chemistry.items()}
        return FittedLink(intercept, _softplus(raw_slope),
                          tuple(_softplus(x) for x in raw_hinges),
                          self.capacity.interior_knots, study, chemistry)


def _sigmoid(x: float) -> float:
    if x >= 0:
        return 1/(1+math.exp(-x))
    e = math.exp(x); return e/(1+e)


def _softplus(x: float) -> float:
    return math.log1p(math.exp(-abs(x))) + max(x, 0.0)


@dataclass(frozen=True)
class FittedLink:
    intercept: float
    slope: float
    hinge_coefficients: tuple[float, ...]
    knots: tuple[float, ...]
    study_intercepts: Mapping[str, float]
    chemistry_intercepts: Mapping[str, float]

    def predict(self, p: float, study_family: str | None = None,
                chemistry: str | None = None) -> float:
        if not 0 <= p <= 1:
            raise ValueError("p must lie in [0,1]")
        z = self.intercept + self.slope*p
        z += sum(a*max(0.0, p-k) for a,k in zip(self.hinge_coefficients, self.knots))
        z += self.study_intercepts.get(study_family, 0.0)
        z += self.chemistry_intercepts.get(chemistry, 0.0)
        return _sigmoid(z)
