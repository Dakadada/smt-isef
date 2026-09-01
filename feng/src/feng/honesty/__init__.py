"""Claim-boundary enforcement for every result in the clean-room Feng arm."""

from .claims import Claim
from .linter import ClaimLintError, LintViolation, lint_claim
from .states import ALLOWED_RESULT_STATES, ResultState

__all__ = [
    "ALLOWED_RESULT_STATES",
    "Claim",
    "ClaimLintError",
    "LintViolation",
    "ResultState",
    "lint_claim",
]

