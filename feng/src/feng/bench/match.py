"""Matched-stratum construction for the Arm D locked nusinersen benchmark.

For each experimentally supported nusinersen off-target exon (a *positive*), this
module draws a stratum of background internal cassette exons from the
control-only detectable universe, matched on five baseline covariates measured in
the GSE307431 control samples only:

* ``baseline_psi``              -- baseline inclusion level
* ``mean_total_junction_reads`` -- junction depth
* ``control_mean_cpm``          -- host-gene expression
* ``exon_length``               -- exon length
* ``gc_fraction``               -- GC content

Every covariate is a control-derived quantity that already lives in the
committed universe parquet. This module NEVER reads a treatment (F18OMe /
F18MOE / F20PMO) alignment, a treatment-vs-control contrast, or the sealed
ΔPSI. It also refuses to build if the universe carries any forbidden contrast
column.

FREEZE BOUNDARY. This module builds the *matching* only -- the strata and their
covariates. It does not compute any score, rank, percentile, ΔR, or permutation
p-value. Those live in :mod:`feng.bench.scores` (frozen sequence scores) and
:mod:`feng.bench.permute` (the Phase 9 permutation test, guarded shut). Nothing
here inspects a benchmark outcome.

Allowed evidence states for anything reported here are ``candidate`` and
``paper_supported`` only.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Mapping
import argparse
import json

import numpy as np
import pandas as pd

from feng.bench.strata import (
    FreezeGuardError,
    assert_no_contrast_fields,
    positive_interval_to_exon_key,
)

# --- covariate contract ------------------------------------------------------

# Transform applied to each covariate BEFORE standardisation. The three
# heavy-tailed count/size covariates are matched on the log scale so that a
# handful of very deep or very long exons cannot dominate the Euclidean
# distance; the two already-bounded fractions are matched linearly. These
# choices are fixed here so the specification is auditable and stable.
LINEAR = "linear"
LOG1P = "log1p"

COVARIATE_TRANSFORMS: dict[str, str] = {
    "baseline_psi": LINEAR,
    "mean_total_junction_reads": LOG1P,
    "control_mean_cpm": LOG1P,
    "exon_length": LOG1P,
    "gc_fraction": LINEAR,
}
COVARIATES: tuple[str, ...] = tuple(COVARIATE_TRANSFORMS)

# Exon-identity columns carried from the universe into every strata row so the
# strata parquet is directly consumable by feng.bench.scores (which needs
# chrom/exon_start_0base/exon_end/strand to window each exon). Without these the
# match->scores seam breaks with a missing-column error.
COORD_COLUMNS: tuple[str, ...] = ("chrom", "exon_start_0base", "exon_end", "strand")

# The default matched-stratum size and the caliper defining an eligible match.
# A background exon is *eligible* for a positive only if, on every covariate, it
# lies within ``caliper_sd`` standardised units of that positive. The kill-gate
# is expressed against the eligible count, so the caliper is what gives the
# >=50-background floor its meaning.
DEFAULT_N_BACKGROUND = 100
DEFAULT_CALIPER_SD = 1.0
DEFAULT_MIN_BACKGROUND = 50
# Frozen seed. Selection below is fully deterministic (nearest by distance with
# an ``exon_key`` tiebreak), so the seed does not affect the strata; it is
# recorded and carried forward for the Phase 9 permutation test, whose null
# shuffling is the only stochastic step in the benchmark.
DEFAULT_SEED = 20260901


def _apply_transform(values: pd.Series, transform: str) -> pd.Series:
    if transform == LINEAR:
        return values.astype(float)
    if transform == LOG1P:
        if (values < 0).any():
            raise ValueError(f"log1p transform requires non-negative values, got negatives in {values.name!r}")
        return np.log1p(values.astype(float))
    raise ValueError(f"unknown covariate transform {transform!r}")


@dataclass(frozen=True, slots=True)
class MatchConfig:
    """Frozen matching specification for the locked benchmark."""

    n_background: int = DEFAULT_N_BACKGROUND
    caliper_sd: float = DEFAULT_CALIPER_SD
    min_background: int = DEFAULT_MIN_BACKGROUND
    seed: int = DEFAULT_SEED
    covariates: tuple[str, ...] = COVARIATES
    transforms: Mapping[str, str] = field(default_factory=lambda: dict(COVARIATE_TRANSFORMS))


def map_positives_to_universe(
    universe: pd.DataFrame, positives: pd.DataFrame
) -> pd.DataFrame:
    """Return per-positive mapping status against the detectable universe.

    Positives carry 1-based inclusive coordinates; the conversion to the
    universe ``exon_key`` frame goes through
    :func:`strata.positive_interval_to_exon_key` so the join cannot silently
    shift by one base. The kill condition "any positive unmappable to the locked
    universe" is checked here (``in_universe`` / ``detectable``); it is not an
    outcome, only a coordinate-join audit.
    """
    required_pos = {"event_id", "chrom", "start", "end", "strand"}
    missing = required_pos.difference(positives.columns)
    if missing:
        raise ValueError(f"positives table lacks columns: {sorted(missing)}")
    if "exon_key" not in universe.columns:
        raise ValueError("universe lacks an exon_key column")

    detect_by_key = (
        universe.set_index("exon_key")["detectable"].to_dict()
        if "detectable" in universe.columns
        else {}
    )
    key_set = set(universe["exon_key"])
    rows = []
    for row in positives.itertuples(index=False):
        key = positive_interval_to_exon_key(row.chrom, int(row.start), int(row.end), row.strand)
        rows.append(
            {
                "event_id": row.event_id,
                "exon_key": key,
                "in_universe": key in key_set,
                "detectable": bool(detect_by_key.get(key, False)),
            }
        )
    return pd.DataFrame(rows)


def _standardise(
    pool: pd.DataFrame, config: MatchConfig
) -> tuple[np.ndarray, dict[str, float], dict[str, float]]:
    """Transform + standardise the covariate matrix over the candidate pool.

    Returns the standardised matrix (rows aligned to ``pool``) plus the
    per-covariate centre and scale used, so the same transform can be applied to
    each positive's covariate vector.
    """
    centres: dict[str, float] = {}
    scales: dict[str, float] = {}
    columns: list[np.ndarray] = []
    for name in config.covariates:
        if name not in pool.columns:
            raise ValueError(f"pool lacks covariate column {name!r}")
        transformed = _apply_transform(pool[name], config.transforms[name])
        if transformed.isna().any():
            raise FreezeGuardError(
                f"covariate {name!r} has null values in the candidate pool; "
                "matching refuses to guess a covariate"
            )
        centre = float(transformed.mean())
        scale = float(transformed.std(ddof=0))
        if scale == 0.0 or not np.isfinite(scale):
            raise FreezeGuardError(
                f"covariate {name!r} has zero or non-finite spread in the pool; "
                "cannot standardise"
            )
        centres[name] = centre
        scales[name] = scale
        columns.append(((transformed - centre) / scale).to_numpy())
    return np.column_stack(columns), centres, scales


def _standardise_vector(
    values: Mapping[str, float],
    config: MatchConfig,
    centres: Mapping[str, float],
    scales: Mapping[str, float],
) -> np.ndarray:
    out = []
    for name in config.covariates:
        raw = float(values[name])
        if config.transforms[name] == LOG1P:
            if raw < 0:
                raise ValueError(f"log1p transform requires non-negative {name!r}, got {raw}")
            raw = float(np.log1p(raw))
        out.append((raw - centres[name]) / scales[name])
    return np.asarray(out, dtype=float)


def build_strata(
    universe: pd.DataFrame,
    positives: pd.DataFrame,
    config: MatchConfig | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build one matched background stratum per mappable positive.

    Returns ``(strata, summary)``:

    * ``strata`` -- long format. Each stratum contributes one ``is_positive``
      anchor row (the positive exon itself) plus one row per matched background
      exon (``is_positive`` False). Every row carries the exon-identity columns
      (``COORD_COLUMNS``) and the five covariates, so the strata parquet is at
      once the scoring input (feng.bench.scores windows each exon) and the
      ranking input (feng.bench.permute ranks each positive within its matched
      backgrounds). Including positive anchor rows is freeze-safe: the scores
      are pure sequence functions and the ranking stays sealed in permute.py.
    * ``summary`` -- one row per positive with ``in_universe`` / ``detectable``,
      the eligible-background count (within caliper), the number actually
      selected, and the per-positive kill-gate verdict.

    Freeze guards: aborts if the universe carries any contrast column, and
    always removes every positive exon_key from the background pool so a
    positive can never match itself or another positive.
    """
    config = config or MatchConfig()
    assert_no_contrast_fields(universe)

    mapping = map_positives_to_universe(universe, positives)
    positive_keys = set(mapping["exon_key"])

    # Candidate background pool: detectable exons, excluding every positive.
    if "detectable" in universe.columns:
        pool = universe.loc[universe["detectable"]].copy()
    else:
        pool = universe.copy()
    pool = pool.loc[~pool["exon_key"].isin(positive_keys)].copy()
    if pool.empty:
        raise FreezeGuardError("background pool is empty after removing positives")
    pool = pool.sort_values("exon_key", kind="mergesort").reset_index(drop=True)

    matrix, centres, scales = _standardise(pool, config)
    pool_keys = pool["exon_key"].to_numpy()

    # A universe row keyed by exon_key so each positive's covariate vector is
    # read from its own row (never inferred).
    universe_by_key = universe.set_index("exon_key")

    strata_rows: list[pd.DataFrame] = []
    summary_rows: list[dict] = []
    for map_row in mapping.itertuples(index=False):
        event_id = map_row.event_id
        key = map_row.exon_key
        base = {
            "event_id": event_id,
            "positive_exon_key": key,
            "in_universe": bool(map_row.in_universe),
            "detectable": bool(map_row.detectable),
        }
        if not map_row.in_universe or not map_row.detectable:
            # Unmappable / undetectable positive: no stratum can be built. This
            # is a kill condition, surfaced explicitly rather than skipped.
            summary_rows.append(
                {
                    **base,
                    "eligible_background": 0,
                    "n_selected": 0,
                    "kill_gate_pass": False,
                    "kill_reason": "positive not in detectable universe",
                }
            )
            continue

        # Read the positive's own coordinate + covariate columns from its own
        # universe row (never inferred).
        pos_row = universe_by_key.loc[key, [*COORD_COLUMNS, *config.covariates]]
        if isinstance(pos_row, pd.DataFrame):  # duplicate key guard
            raise FreezeGuardError(f"positive {event_id} maps to duplicate universe rows for {key}")
        pos_vec = _standardise_vector(
            pos_row[list(config.covariates)].to_dict(), config, centres, scales
        )

        # Per-covariate absolute standardised deviation; eligibility requires
        # every covariate within the caliper.
        deviations = np.abs(matrix - pos_vec)
        eligible_mask = np.all(deviations <= config.caliper_sd, axis=1)
        eligible_idx = np.flatnonzero(eligible_mask)
        eligible_count = int(eligible_idx.size)

        # Euclidean distance in standardised space, restricted to eligible
        # candidates; deterministic nearest-N with an exon_key tiebreak.
        distances = np.sqrt(np.sum((matrix[eligible_idx] - pos_vec) ** 2, axis=1))
        order = np.lexsort((pool_keys[eligible_idx], distances))
        take = eligible_idx[order][: config.n_background]
        take_dist = distances[order][: config.n_background]

        # Background rows: exon_key + coordinates + covariates, plus the anchor
        # columns. The universe carries its own rMATS ``event_id``; the positive
        # anchor is named ``positive_event_id`` to avoid clobbering it.
        block = pool.iloc[take][["exon_key", *COORD_COLUMNS, *config.covariates]].copy()
        block.insert(0, "match_distance", take_dist)
        block.insert(0, "is_positive", False)

        # Positive anchor row for this stratum: same schema, is_positive True,
        # distance 0. This is what the frozen score sweep scores and what the
        # sealed ranker ranks against the backgrounds.
        anchor = pd.DataFrame(
            [{
                "match_distance": 0.0,
                "is_positive": True,
                "exon_key": key,
                **{c: pos_row[c] for c in COORD_COLUMNS},
                **{c: pos_row[c] for c in config.covariates},
            }]
        )

        stratum = pd.concat([anchor, block], ignore_index=True)
        stratum.insert(0, "positive_exon_key", key)
        stratum.insert(0, "positive_event_id", event_id)
        strata_rows.append(stratum)

        summary_rows.append(
            {
                **base,
                "eligible_background": eligible_count,
                "n_selected": int(take.size),
                "kill_gate_pass": eligible_count >= config.min_background,
                "kill_reason": ""
                if eligible_count >= config.min_background
                else f"only {eligible_count} eligible backgrounds (< {config.min_background})",
            }
        )

    strata = (
        pd.concat(strata_rows, ignore_index=True)
        if strata_rows
        else pd.DataFrame(
            columns=[
                "positive_event_id", "positive_exon_key", "is_positive",
                "match_distance", "exon_key", *COORD_COLUMNS, *config.covariates,
            ]
        )
    )
    summary = pd.DataFrame(summary_rows)
    strata.attrs["seed"] = config.seed
    strata.attrs["caliper_sd"] = config.caliper_sd
    strata.attrs["n_background"] = config.n_background
    return strata, summary


def check_kill_gates(summary: pd.DataFrame, min_background: int = DEFAULT_MIN_BACKGROUND) -> bool:
    """Return True iff every positive cleared the >=min_background eligibility floor."""
    if summary.empty:
        return False
    return bool((summary["eligible_background"] >= min_background).all())


# --- CLI ---------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build Arm D matched background strata (control-only).")
    parser.add_argument("--universe", type=Path, required=True, help="detectable universe parquet")
    parser.add_argument("--positives", type=Path, required=True, help="positives_v50.tsv")
    parser.add_argument("--out-strata", type=Path, required=True, help="output strata parquet")
    parser.add_argument("--out-summary", type=Path, required=True, help="output per-positive summary tsv")
    parser.add_argument("--n-background", type=int, default=DEFAULT_N_BACKGROUND)
    parser.add_argument("--caliper-sd", type=float, default=DEFAULT_CALIPER_SD)
    parser.add_argument("--min-background", type=int, default=DEFAULT_MIN_BACKGROUND)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args(argv)

    universe = pd.read_parquet(args.universe)
    positives = pd.read_csv(args.positives, sep="\t")
    config = MatchConfig(
        n_background=args.n_background,
        caliper_sd=args.caliper_sd,
        min_background=args.min_background,
        seed=args.seed,
    )
    strata, summary = build_strata(universe, positives, config)

    args.out_strata.parent.mkdir(parents=True, exist_ok=True)
    strata.to_parquet(args.out_strata, index=False)
    summary.to_csv(args.out_summary, sep="\t", index=False)

    all_pass = check_kill_gates(summary, args.min_background)
    print(json.dumps({"positives": len(summary), "all_kill_gates_pass": all_pass}))
    print(summary.to_string(index=False))
    if not all_pass:
        print("KILL-GATE FAIL: one or more positives below the eligible-background floor", flush=True)
        return 1
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
