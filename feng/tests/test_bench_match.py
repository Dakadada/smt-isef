"""Matching + kill-gate + control-only tests for the Arm D benchmark strata."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from feng.bench.match import (
    COVARIATES,
    MatchConfig,
    build_strata,
    check_kill_gates,
    map_positives_to_universe,
)
from feng.bench.strata import FreezeGuardError


REPO_ROOT = Path(__file__).resolve().parents[2]
UNIVERSE = REPO_ROOT / "feng" / "artifacts" / "bench" / "detectable_universe_v50.parquet"
POSITIVES = REPO_ROOT / "feng" / "artifacts" / "bench" / "positives_v50.tsv"


def _universe_row(exon_key: str, chrom: str, start0: int, end: int, strand: str, **cov) -> dict:
    row = {
        "exon_key": exon_key,
        "chrom": chrom,
        "exon_start_0base": start0,
        "exon_end": end,
        "strand": strand,
        "detectable": cov.pop("detectable", True),
    }
    # default covariate values, overridable
    defaults = {
        "baseline_psi": 0.5,
        "mean_total_junction_reads": 50.0,
        "control_mean_cpm": 10.0,
        "exon_length": end - start0,
        "gc_fraction": 0.5,
    }
    defaults.update(cov)
    row.update(defaults)
    return row


def _synthetic_universe(n_close: int, n_far: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Universe with one positive, ``n_close`` matchable + ``n_far`` unmatchable backgrounds."""
    rows = []
    # positive at chr1:100-200:+ -> exon_key chr1:99-200:+ (1-based 100..200)
    pos_key = "chr1:99-200:+"
    rows.append(_universe_row(pos_key, "chr1", 99, 200, "+"))
    # close backgrounds: identical covariates to the positive -> standardised dist 0
    for i in range(n_close):
        s0 = 1000 + i * 500
        rows.append(_universe_row(f"chr1:{s0}-{s0 + 101}:+", "chr1", s0, s0 + 101, "+"))
    # far backgrounds: extreme covariates -> outside the caliper
    for i in range(n_far):
        s0 = 500000 + i * 500
        rows.append(
            _universe_row(
                f"chr1:{s0}-{s0 + 300}:+", "chr1", s0, s0 + 300, "+",
                baseline_psi=0.99, mean_total_junction_reads=5000.0,
                control_mean_cpm=5000.0, exon_length=300, gc_fraction=0.9,
            )
        )
    universe = pd.DataFrame(rows)
    positives = pd.DataFrame(
        [{"event_id": "POS1", "chrom": "chr1", "start": 100, "end": 200, "strand": "+"}]
    )
    return universe, positives


def test_map_positives_to_universe_flags_presence() -> None:
    universe, positives = _synthetic_universe(5, 0)
    mapping = map_positives_to_universe(universe, positives)
    assert mapping.loc[0, "exon_key"] == "chr1:99-200:+"
    assert bool(mapping.loc[0, "in_universe"])
    assert bool(mapping.loc[0, "detectable"])


def test_map_positive_absent_from_universe() -> None:
    universe, positives = _synthetic_universe(5, 0)
    positives.loc[0, "end"] = 250  # shift so the exon_key no longer exists
    mapping = map_positives_to_universe(universe, positives)
    assert not bool(mapping.loc[0, "in_universe"])
    assert not bool(mapping.loc[0, "detectable"])


def test_build_strata_excludes_positive_and_selects_close_backgrounds() -> None:
    universe, positives = _synthetic_universe(n_close=60, n_far=40)
    strata, summary = build_strata(universe, positives, MatchConfig())
    # positive appears once as a flagged anchor, never as a background row
    backgrounds = strata.loc[~strata["is_positive"]]
    assert "chr1:99-200:+" not in set(backgrounds["exon_key"])
    anchors = strata.loc[strata["is_positive"]]
    assert set(anchors["exon_key"]) == {"chr1:99-200:+"}
    # only the close group is eligible; selection capped at n_background
    row = summary.iloc[0]
    assert int(row["eligible_background"]) == 60
    assert int(row["n_selected"]) == 60  # fewer than n_background=100 available
    assert bool(row["kill_gate_pass"])
    # every selected background is from the close group (exon_length 101), never
    # the far group (exon_length 300) -- eligible_background == 60 already proves
    # only the close group cleared the caliper.
    assert (strata["exon_length"] == 101).all()


def test_kill_gate_fails_below_floor() -> None:
    universe, positives = _synthetic_universe(n_close=10, n_far=40)
    strata, summary = build_strata(universe, positives, MatchConfig())
    row = summary.iloc[0]
    assert int(row["eligible_background"]) == 10
    assert not bool(row["kill_gate_pass"])
    assert not check_kill_gates(summary)
    assert "10 eligible" in row["kill_reason"]


def test_undetectable_positive_is_killed_with_no_stratum() -> None:
    # n_far>0 so the background pool has real covariate spread; the positive
    # itself is flipped undetectable below.
    universe, positives = _synthetic_universe(n_close=60, n_far=10)
    universe.loc[universe["exon_key"] == "chr1:99-200:+", "detectable"] = False
    strata, summary = build_strata(universe, positives, MatchConfig())
    row = summary.iloc[0]
    assert not bool(row["kill_gate_pass"])
    assert int(row["n_selected"]) == 0
    assert strata.empty


def test_contrast_field_in_universe_is_refused() -> None:
    universe, positives = _synthetic_universe(n_close=60, n_far=0)
    universe["IncLevelDifference"] = 0.0
    with pytest.raises(FreezeGuardError):
        build_strata(universe, positives, MatchConfig())


def test_selection_is_deterministic() -> None:
    universe, positives = _synthetic_universe(n_close=60, n_far=40)
    a, _ = build_strata(universe, positives, MatchConfig(seed=1))
    b, _ = build_strata(universe, positives, MatchConfig(seed=999))
    # selection does not depend on the seed (deterministic nearest + tiebreak)
    assert list(a["exon_key"]) == list(b["exon_key"])


@pytest.mark.skipif(not UNIVERSE.is_file(), reason="committed universe parquet required")
def test_all_committed_positives_pass_kill_gate() -> None:
    universe = pd.read_parquet(UNIVERSE)
    positives = pd.read_csv(POSITIVES, sep="\t")
    strata, summary = build_strata(universe, positives, MatchConfig())
    assert len(summary) == 8
    assert bool((summary["in_universe"] & summary["detectable"]).all())
    assert bool((summary["eligible_background"] >= 50).all())
    assert bool((summary["n_selected"] == 100).all())
    assert check_kill_gates(summary)
    # each covariate is present in the strata output
    for cov in COVARIATES:
        assert cov in strata.columns
