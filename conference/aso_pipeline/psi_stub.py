"""Item 2 plug-in boundary; intentionally performs no PSI analysis."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from typing import Any, NoReturn


EXPECTED_PSI_COLUMNS = ("aso_id", "length", "replicate", "PSI_percent")


class PsiDataUnavailable(RuntimeError):
    """Raised because wet-lab PSI is not yet available for this candidate work."""


def fit_dg_psi(
    rows: Iterable[Mapping[str, Any]], dg_by_aso: Mapping[str, float]
) -> NoReturn:
    """Validate the future replicate-level schema, then refuse to fit.

    Expected row schema: ``aso_id: str``, ``length: int``, ``replicate: str``,
    and ``PSI_percent: float``. ``dg_by_aso`` will map each ASO ID to one frozen
    candidate ΔG value. A later, separately authorized task may replace this
    refusal with a pre-specified regression; this task never reads or invents
    PSI.
    """
    materialized = list(rows)
    if not materialized:
        raise ValueError("PSI rows are empty")
    for index, row in enumerate(materialized):
        missing = [column for column in EXPECTED_PSI_COLUMNS if column not in row]
        if missing:
            raise ValueError(f"PSI row {index} missing required columns: {missing}")
        psi = float(row["PSI_percent"])
        if not 0.0 <= psi <= 100.0:
            raise ValueError(f"PSI row {index} is outside [0, 100]")
        if str(row["aso_id"]) not in dg_by_aso:
            raise ValueError(f"PSI row {index} has no frozen ΔG value")
    raise PsiDataUnavailable(
        "wet-lab PSI data are not available; Item 2 is an interface stub only"
    )

