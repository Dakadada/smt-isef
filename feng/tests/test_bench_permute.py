"""Phase-9 permutation-test seal + pure ranking-primitive tests."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from feng.bench.permute import (
    PHASE9_UNLOCK_ENV,
    assert_phase9_unlocked,
    percentile_within_stratum,
    phase9_is_unlocked,
    run_phase9_test,
)
from feng.bench.strata import FreezeGuardError


# --- pure ranking primitive (safe to test) -----------------------------------


def test_percentile_lower_is_better() -> None:
    bg = np.array([5.0, 6.0, 7.0, 8.0])
    # positive better (lower) than all background -> high percentile
    assert percentile_within_stratum(1.0, bg, lower_is_better=True) == pytest.approx(0.9)
    # positive worse (higher) than all -> low percentile
    assert percentile_within_stratum(9.0, bg, lower_is_better=True) == pytest.approx(0.1)


def test_percentile_higher_is_better_symmetry() -> None:
    bg = np.array([5.0, 6.0, 7.0, 8.0])
    assert percentile_within_stratum(9.0, bg, lower_is_better=False) == pytest.approx(0.9)


def test_percentile_ties_use_midrank() -> None:
    bg = np.array([1.0, 1.0, 3.0])
    # positive ties two backgrounds; mid-rank keeps it centred
    val = percentile_within_stratum(1.0, bg, lower_is_better=True)
    assert 0.0 < val < 1.0


def test_percentile_empty_stratum_raises() -> None:
    with pytest.raises(ValueError):
        percentile_within_stratum(1.0, np.array([]), lower_is_better=True)


# --- freeze seal on the outcome computation ----------------------------------


def test_assert_phase9_unlocked_raises_prefreeze(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(PHASE9_UNLOCK_ENV, raising=False)
    assert not phase9_is_unlocked()
    with pytest.raises(FreezeGuardError):
        assert_phase9_unlocked()


def test_phase9_unlock_via_env(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(PHASE9_UNLOCK_ENV, "1")
    assert phase9_is_unlocked()
    assert_phase9_unlocked()  # does not raise once unlocked


def test_phase9_unlock_via_sentinel(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(PHASE9_UNLOCK_ENV, raising=False)
    sentinel = tmp_path / "unlock.flag"
    assert not phase9_is_unlocked(sentinel)  # absent
    sentinel.write_text("ok")
    assert phase9_is_unlocked(sentinel)


def test_run_phase9_test_locked_raises_freezeguard(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(PHASE9_UNLOCK_ENV, raising=False)
    scored = pd.DataFrame({"positive_event_id": ["POS1"], "baseline_min_hamming": [3], "intarna_min_energy": [-5.0]})
    with pytest.raises(FreezeGuardError):
        run_phase9_test(scored)
