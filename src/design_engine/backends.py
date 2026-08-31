"""Explicit deterministic offline backend; never a biological predictor."""
from __future__ import annotations

from physics_adapter import FoldingResult, Probability, expected_receipt


class DeterministicStubFolder:
    """Fixture backend based only on declared event/occupancy constraints.

    It exists for end-to-end tests and demonstrations without weights. Reports
    using it remain prospective and identify this backend in provenance.
    """

    backend_version = "design-engine-deterministic-stub/1"

    def productive_probability(self, request):
        base = 0.55
        if request.feature_variant:
            base -= 0.03
        if request.max_bp_span is not None:
            base -= 0.02
        bound = request.bound_complex is not None
        value = base
        if bound:
            target = request.bound_complex.transcript_coordinates
            accessible_overlap = any(target.overlaps(x) for x in request.required_accessible)
            silencer_overlap = any(target.overlaps(x) for x in request.required_silencer_free)
            value += (0.22 if silencer_overlap else 0) - (0.22 if accessible_overlap else 0)
        width = 0.025
        probability = Probability(value, value - width, value + width)
        return FoldingResult(probability, expected_receipt(request, self.backend_version))

