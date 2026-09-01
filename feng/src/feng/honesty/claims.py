"""Validated claim record for result reporting in the clean-room Feng arm."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .linter import LintViolation, lint_claim
from .states import ResultState


@dataclass(frozen=True, slots=True, init=False)
class Claim:
    """A claim coupled to its evidence state and declared output fields."""

    text: str
    state: ResultState
    field_names: tuple[str, ...] = ()
    lint_violations: tuple[LintViolation, ...] = ()

    def __init__(
        self,
        text: str,
        state: ResultState | str,
        field_names: Iterable[str] = (),
        *,
        override_reason: str | None = None,
        override_log: str | Path | None = None,
    ) -> None:
        """Construct a claim only after validating its state, text, and schema."""
        parsed_state = ResultState(state)
        fields = tuple(field_names)
        violations = lint_claim(
            text,
            fields,
            override_reason=override_reason,
            override_log=override_log,
        )
        object.__setattr__(self, "text", text)
        object.__setattr__(self, "state", parsed_state)
        object.__setattr__(self, "field_names", fields)
        object.__setattr__(self, "lint_violations", violations)

    @classmethod
    def create(
        cls,
        text: str,
        state: ResultState | str,
        field_names: Iterable[str] = (),
        *,
        override_reason: str | None = None,
        override_log: str | Path | None = None,
    ) -> tuple[Claim, tuple[LintViolation, ...]]:
        """Construct a claim only after its text and schema pass the linter."""
        claim = cls(
            text,
            state,
            field_names,
            override_reason=override_reason,
            override_log=override_log,
        )
        return claim, claim.lint_violations
