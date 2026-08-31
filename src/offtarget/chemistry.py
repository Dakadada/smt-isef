"""Transparent chemistry-specific scoring parameters.

These are coarse enumeration priors, not a substitute for an experimentally
validated nearest-neighbor model. Unknown chemistries are explicitly flagged.
"""

from dataclasses import dataclass
from typing import Dict, Optional


@dataclass(frozen=True)
class ChemistryParameters:
    name: str
    match_delta_g: float
    mismatch_penalty: float
    indel_penalty: float
    known: bool
    heuristic: bool = False


PARAMETERS: Dict[str, ChemistryParameters] = {
    "2MOE": ChemistryParameters("2MOE", -1.75, 1.35, 2.10, True),
    "OME": ChemistryParameters("OMe", -1.60, 1.30, 2.00, True),
    "PMO": ChemistryParameters("PMO", -1.45, 1.20, 1.90, True),
    "LNA": ChemistryParameters("LNA", -2.10, 1.55, 2.35, True),
}

HEURISTIC_FALLBACK = ChemistryParameters("unknown", -1.50, 1.50, 2.25, False, True)


def parameters_for(name: str, allow_heuristic: bool = False) -> Optional[ChemistryParameters]:
    key = name.upper().replace("'", "").replace("-", "")
    aliases = {"2MOE": "2MOE", "MOE": "2MOE", "2OME": "OME", "OME": "OME"}
    known = PARAMETERS.get(aliases.get(key, key))
    if known is not None:
        return known
    if not allow_heuristic:
        return None
    return ChemistryParameters(name, HEURISTIC_FALLBACK.match_delta_g,
                               HEURISTIC_FALLBACK.mismatch_penalty,
                               HEURISTIC_FALLBACK.indel_penalty, False, True)


def estimate_delta_g(length: int, mismatches: int, indels: int,
                     params: ChemistryParameters) -> float:
    paired = max(0, length - mismatches - indels)
    return (paired * params.match_delta_g + mismatches * params.mismatch_penalty
            + indels * params.indel_penalty)
