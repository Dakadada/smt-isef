"""Arm D Phase-9 matched-stratum permutation test -- GUARDED SHUT until freeze.

The locked benchmark's confirmatory readout is a one-sided, matched-stratum label
permutation over 10^6 shuffles with a frozen seed, comparing the baseline score
(``min 18-nt Hamming distance``) against the test score (accessibility-aware
IntaRNA interaction energy). The effect size is::

    delta_r = mean over positives of ( percentile_test - percentile_baseline )

where each percentile is the positive's rank within its own matched stratum.

FREEZE DISCIPLINE. Computing delta_r or its permutation p-value inspects the
benchmark OUTCOME and must not happen before preregistration is written and the
git freeze tag is cut (Phase 8), after which the test runs ONCE (Phase 9). This
module therefore refuses to run unless an explicit unlock is present:

    FENG_PHASE9_UNLOCK=1  (environment)   -- or --
    a sentinel file whose path is passed to ``assert_phase9_unlocked``.

Every outcome entry point calls :func:`assert_phase9_unlocked` first. The pure
building block :func:`percentile_within_stratum` is safe to unit-test on
synthetic data; the aggregate / permutation functions are never run pre-freeze
and are intentionally left inert.
"""

from __future__ import annotations

from pathlib import Path
import os

import numpy as np
import pandas as pd

from feng.bench.strata import FreezeGuardError

PHASE9_UNLOCK_ENV = "FENG_PHASE9_UNLOCK"
DEFAULT_N_PERMUTATIONS = 1_000_000
DEFAULT_SEED = 20260901


def phase9_is_unlocked(sentinel: Path | None = None) -> bool:
    """True iff the Phase-9 unlock is present (env flag or sentinel file)."""
    if os.environ.get(PHASE9_UNLOCK_ENV) == "1":
        return True
    return sentinel is not None and Path(sentinel).is_file()


def assert_phase9_unlocked(sentinel: Path | None = None) -> None:
    """Refuse to run any outcome computation before the Phase-9 unlock."""
    if not phase9_is_unlocked(sentinel):
        raise FreezeGuardError(
            "Arm D permutation test is sealed until Phase 9. Set "
            f"{PHASE9_UNLOCK_ENV}=1 (or pass an existing sentinel file) only "
            "after the preregistration freeze tag is cut. Refusing to compute "
            "delta_r or a permutation p-value pre-freeze."
        )


def percentile_within_stratum(
    positive_value: float, background_values: np.ndarray, *, lower_is_better: bool
) -> float:
    """Percentile rank of a positive against its matched background.

    Returns the fraction of the (positive + background) set that the positive
    ranks at or better than, in [0, 1], where "better" means lower when
    ``lower_is_better`` (Hamming distance, IntaRNA energy) and higher otherwise.
    A pure ranking primitive -- safe to test on synthetic data. It does not, by
    itself, read any real score.
    """
    background = np.asarray(background_values, dtype=float)
    if background.size == 0:
        raise ValueError("cannot rank a positive against an empty stratum")
    total = background.size + 1
    if lower_is_better:
        rank_better = int(np.sum(background > positive_value))
    else:
        rank_better = int(np.sum(background < positive_value))
    ties = int(np.sum(background == positive_value))
    # Mid-rank convention for ties, counting the positive itself.
    return (rank_better + (ties + 1) / 2) / total


def run_phase9_test(
    scored_strata: pd.DataFrame,
    *,
    baseline_col: str = "baseline_min_hamming",
    test_col: str = "intarna_min_energy",
    n_permutations: int = DEFAULT_N_PERMUTATIONS,
    seed: int = DEFAULT_SEED,
    sentinel: Path | None = None,
) -> dict:  # pragma: no cover - Phase 9 only, never run pre-freeze
    """Run the one-shot matched-stratum permutation test.

    SEALED until Phase 9. Calls :func:`assert_phase9_unlocked` before touching
    any score, then would compute delta_r and its one-sided permutation p-value
    over ``n_permutations`` shuffles with ``seed``. The body is intentionally
    left unimplemented so it cannot run pre-freeze even if the guard were
    bypassed.
    """
    assert_phase9_unlocked(sentinel)
    raise NotImplementedError(
        "Phase 9 execution: implement against the frozen preregistered score "
        "columns after the freeze tag is cut."
    )
