"""Allowed evidence states for every reported result in the clean-room Feng arm."""

from enum import Enum


class ResultState(str, Enum):
    """The exhaustive evidence-state vocabulary permitted by the arm."""

    CANDIDATE = "candidate"
    PAPER_SUPPORTED = "paper_supported"
    CLEAN_REANALYSIS_REPRODUCED = "clean_reanalysis_reproduced"


ALLOWED_RESULT_STATES = frozenset(state.value for state in ResultState)

