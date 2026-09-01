"""Public construction wrapper around the captured approved derivation."""
from __future__ import annotations

from datetime import datetime, timezone

from .approved import derive_characterization as _APPROVED_DERIVATION
from .models import (CandidateScore, EffectorClassSummary,
                     EffectorScoreCharacterization, ScoreDistribution)


def derive_characterization(timestamp: datetime) -> dict:
    """Public inspection seam; render validation does not resolve through this name."""
    return _APPROVED_DERIVATION(timestamp)


def build_characterization(*, clock=None) -> EffectorScoreCharacterization:
    timestamp = (clock or (lambda: datetime.now(timezone.utc)))()
    raw = _APPROVED_DERIVATION(timestamp)
    raw["candidates"] = tuple(CandidateScore(**row) for row in raw["candidates"])
    raw["distributions"] = tuple(ScoreDistribution(**row) for row in raw["distributions"])
    raw["class_summaries"] = tuple(EffectorClassSummary(**row)
                                    for row in raw["class_summaries"])
    return EffectorScoreCharacterization(**raw)
