"""A2: the A_5ss -> PSI link function. DELIBERATELY UNFIT.

This module declares the *shape* of the mapping from computed 5' splice site
accessibility to measured exon 7 inclusion (PSI) and nothing else. It ships
with NO default parameter values and it will refuse to run without explicit
ones.

Why: the parameters of this link are exactly where a placeholder constant would
silently become a "prediction". Until real experimental PSI data exists for the
actual construct the wet lab runs (docs/lab-request.md, Ask 1 and Ask 2), any
number here would be invented. A model that reports delta_A honestly is useful;
a model that converts delta_A into a fabricated PSI is not.

When data arrives, fit :class:`LinkParams` by maximum likelihood against
measured PSI and store the fit -- with its provenance, n, and CI -- in a
version-controlled file. Do not hardcode it here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


class UnfitLinkError(NotImplementedError):
    """Raised when the link function is called without explicit parameters."""


@dataclass(frozen=True)
class LinkParams:
    """Parameters of the logistic link.

    psi = psi_min + (psi_max - psi_min) / (1 + exp(-slope * (A_5ss - a_half)))

    slope:
        Steepness, dimensionless (A_5ss is a probability in [0, 1]).
    a_half:
        A_5ss at the half-maximal point, dimensionless probability.
    psi_min, psi_max:
        Inclusion floor and ceiling in percent, 0..100. Both are properties of
        the construct and cell line, not of the effector, and must be measured.

    Every field is required. There are no defaults, by design.
    """

    slope: float
    a_half: float
    psi_min: float
    psi_max: float
    provenance: str  # what data this was fit to; empty string is rejected

    def __post_init__(self) -> None:
        if not self.provenance.strip():
            raise UnfitLinkError(
                "LinkParams.provenance must name the experimental dataset this "
                "fit came from; unattributed parameters are not permitted."
            )
        if not 0.0 <= self.a_half <= 1.0:
            raise ValueError("a_half is an A_5ss probability and must lie in [0, 1]")
        if not 0.0 <= self.psi_min < self.psi_max <= 100.0:
            raise ValueError("require 0 <= psi_min < psi_max <= 100 (percent)")


def psi_from_a5ss(a_5ss: float, params: Optional[LinkParams] = None) -> float:
    """Map A_5ss to PSI (percent).

    Raises
    ------
    UnfitLinkError
        Always, unless ``params`` is supplied explicitly. There is no fitted
        default and there will not be one until measured PSI exists.
    """
    if params is None:
        raise UnfitLinkError(
            "The A_5ss -> PSI link is unfit. No calibration constants are shipped "
            "with this project. Supply LinkParams fit to real measured PSI, or "
            "report delta_A directly."
        )
    if not 0.0 <= a_5ss <= 1.0:
        raise ValueError("a_5ss must be a probability in [0, 1]")
    import math

    span = params.psi_max - params.psi_min
    return params.psi_min + span / (1.0 + math.exp(-params.slope * (a_5ss - params.a_half)))


def delta_psi(
    a_baseline: float, a_blocked: float, params: Optional[LinkParams] = None
) -> float:
    """PSI(blocked) - PSI(baseline), percent. Same unfit refusal applies."""
    return psi_from_a5ss(a_blocked, params) - psi_from_a5ss(a_baseline, params)
