"""Audited claim-language linter for the clean-room Feng arm."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
from typing import Iterable


FORBIDDEN_STRINGS = (
    "first",
    "genome-wide safety map",
    "clinical response",
    "beats nusinersen",
    "generalizable",
)
FORBIDDEN_FIELD_NAMES = frozenset(
    {
        "clinical_response",
        "patient_responder",
        "safe",
        "efficacious",
        "superior_to_nusinersen",
    }
)


@dataclass(frozen=True, slots=True)
class LintViolation:
    """One forbidden phrase or schema field found by the arm linter."""

    kind: str
    value: str


class ClaimLintError(ValueError):
    """Raised when forbidden claim content lacks a logged override."""

    def __init__(self, violations: Iterable[LintViolation]) -> None:
        self.violations = tuple(violations)
        detail = ", ".join(f"{item.kind}:{item.value}" for item in self.violations)
        super().__init__(f"claim boundary violation: {detail}")


def _normalized_text(text: str) -> str:
    return re.sub(r"\s+", " ", text.casefold()).strip()


def find_violations(text: str, field_names: Iterable[str] = ()) -> tuple[LintViolation, ...]:
    """Return every forbidden phrase and exact forbidden field name."""
    normalized = _normalized_text(text)
    violations = [
        LintViolation("string", phrase)
        for phrase in FORBIDDEN_STRINGS
        if phrase in normalized
    ]
    violations.extend(
        LintViolation("field", field)
        for field in sorted({name.casefold() for name in field_names})
        if field in FORBIDDEN_FIELD_NAMES
    )
    return tuple(violations)


def lint_claim(
    text: str,
    field_names: Iterable[str] = (),
    *,
    override_reason: str | None = None,
    override_log: str | Path | None = None,
) -> tuple[LintViolation, ...]:
    """Validate claim content, requiring a durable log for any explicit override."""
    violations = find_violations(text, field_names)
    if not violations:
        return ()
    if not override_reason or not override_reason.strip() or override_log is None:
        raise ClaimLintError(violations)

    record = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "reason": override_reason.strip(),
        "claim_sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
        "violations": [asdict(item) for item in violations],
    }
    log_path = Path(override_log)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")
    return violations

