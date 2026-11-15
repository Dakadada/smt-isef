"""Generate Item 6 thermodynamic null artifacts and alignment query manifest."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
import statistics
from typing import Iterable

from . import CANDIDATE
from .item1 import ASO_ORDER, PRIMARY_SALT_MOLAR, _write_tsv
from .sequences import ASOS, derive_target
from .shuffle import dinucleotide_shuffles
from .thermo import score_pair


RNG_SEED = 20260926
NULL_DRAWS_PER_ASO = 100


def _linear_quantile(sorted_values: list[float], probability: float) -> float:
    position = (len(sorted_values) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return sorted_values[lower]
    fraction = position - lower
    return sorted_values[lower] * (1 - fraction) + sorted_values[upper] * fraction


def summarize_null(
    real_value: float,
    values: Iterable[float],
    *,
    lower_is_stronger: bool,
) -> dict[str, float | int]:
    numeric = sorted(float(value) for value in values if value is not None)
    if not numeric:
        raise ValueError("null distribution has no numeric values")
    below = sum(value < real_value for value in numeric)
    equal = sum(value == real_value for value in numeric)
    if lower_is_stronger:
        extreme = sum(value <= real_value for value in numeric)
    else:
        extreme = sum(value >= real_value for value in numeric)
    return {
        "null_n_numeric": len(numeric),
        "null_mean": statistics.fmean(numeric),
        "null_sd_population": statistics.pstdev(numeric),
        "null_min": numeric[0],
        "null_q05": _linear_quantile(numeric, 0.05),
        "null_q25": _linear_quantile(numeric, 0.25),
        "null_median": _linear_quantile(numeric, 0.50),
        "null_q75": _linear_quantile(numeric, 0.75),
        "null_q95": _linear_quantile(numeric, 0.95),
        "null_max": numeric[-1],
        "real_minus_null_mean": real_value - statistics.fmean(numeric),
        "real_percentile_midrank": 100.0 * (below + 0.5 * equal) / len(numeric),
        "one_sided_enrichment_p": (extreme + 1) / (len(numeric) + 1),
    }


def censored_no_duplex_lower_tail_p(
    real_value: float, numeric_values: Iterable[float], no_duplex_count: int
) -> float:
    """Empirical lower-tail p with no-duplex draws censored as weaker.

    No energy is assigned to those draws. For a negative real interaction they
    cannot be at least as strong in the lower (more-negative) tail, but they do
    remain valid null draws in the denominator.
    """
    numeric = [float(value) for value in numeric_values]
    if no_duplex_count < 0:
        raise ValueError("no_duplex_count must be nonnegative")
    denominator = len(numeric) + no_duplex_count + 1
    return (sum(value <= real_value for value in numeric) + 1) / denominator


def run_item6_thermo(
    output_dir: Path,
    *,
    intarna: Path,
    rnaduplex: Path,
    seed: int = RNG_SEED,
) -> tuple[list[dict], list[dict]]:
    raw_rows: list[dict] = []
    summary_rows: list[dict] = []
    query_rows: list[dict] = []

    for aso_id in ASO_ORDER:
        aso = ASOS[aso_id]
        target = derive_target(aso.rna_sequence)
        draws = dinucleotide_shuffles(
            aso.rna_sequence, n=NULL_DRAWS_PER_ASO, seed=seed
        )
        real = score_pair(
            aso.rna_sequence,
            target.sequence,
            intarna=intarna,
            rnaduplex=rnaduplex,
            salt_molar=PRIMARY_SALT_MOLAR,
        )
        score_cache = {}
        for sequence in sorted(set(draws)):
            score_cache[sequence] = score_pair(
                sequence,
                target.sequence,
                intarna=intarna,
                rnaduplex=rnaduplex,
                salt_molar=PRIMARY_SALT_MOLAR,
            )
        for draw_id, sequence in enumerate(draws, 1):
            score = score_cache[sequence]
            raw_rows.append(
                {
                    "state": CANDIDATE,
                    "aso_id": aso_id,
                    "draw_id": draw_id,
                    "rng_seed": seed,
                    "shuffled_aso_rna_5to3": sequence,
                    "target_rna_5to3": target.sequence,
                    "intarna_dg_kcal_mol": (
                        "" if score.intarna_no_duplex else score.intarna_dg_kcal_mol
                    ),
                    "intarna_no_duplex": score.intarna_no_duplex,
                    "nn_reference_dg_kcal_mol": score.nn_raw_reference_dg_kcal_mol,
                    "nn_salt_correction_kcal_mol": score.nn_salt_correction_kcal_mol,
                    "nn_tbe_dg_kcal_mol": score.nn_dg_kcal_mol,
                    "nn_tbe_salt_molar": PRIMARY_SALT_MOLAR,
                }
            )
            query_rows.append(
                {
                    "state": CANDIDATE,
                    "query_id": f"{aso_id}__shuffle_{draw_id:03d}",
                    "aso_id": aso_id,
                    "draw_id": draw_id,
                    "is_real": False,
                    "sequence_rna_5to3": sequence,
                }
            )
        query_rows.append(
            {
                "state": CANDIDATE,
                "query_id": f"{aso_id}__real",
                "aso_id": aso_id,
                "draw_id": "",
                "is_real": True,
                "sequence_rna_5to3": aso.rna_sequence,
            }
        )

        intarna_null = [
            score_cache[sequence].intarna_dg_kcal_mol
            for sequence in draws
            if not score_cache[sequence].intarna_no_duplex
        ]
        for model, real_value, null_values, no_duplex_count in (
            (
                "IntaRNA-default",
                real.intarna_dg_kcal_mol,
                intarna_null,
                len(draws) - len(intarna_null),
            ),
            (
                "Turner2004-NN-TBE",
                real.nn_dg_kcal_mol,
                [score_cache[sequence].nn_dg_kcal_mol for sequence in draws],
                0,
            ),
        ):
            summary = summarize_null(real_value, null_values, lower_is_stronger=True)
            censored_p = (
                censored_no_duplex_lower_tail_p(
                    real_value, null_values, no_duplex_count
                )
                if model == "IntaRNA-default"
                else summary["one_sided_enrichment_p"]
            )
            summary_rows.append(
                {
                    "state": CANDIDATE,
                    "aso_id": aso_id,
                    "length_nt": aso.length,
                    "model": model,
                    "real_dg_kcal_mol": real_value,
                    "null_n_requested": NULL_DRAWS_PER_ASO,
                    "null_n_no_duplex": no_duplex_count,
                    "null_n_unique_sequences": len(set(draws)),
                    **summary,
                    "one_sided_enrichment_p_all_draws_censored": censored_p,
                }
            )

    output_dir = Path(output_dir)
    _write_tsv(output_dir / "item6_thermo_null_raw.tsv", raw_rows)
    _write_tsv(output_dir / "item6_thermo_null_summary.tsv", summary_rows)
    _write_tsv(output_dir / "item6_alignment_queries.tsv", query_rows)
    return raw_rows, summary_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--intarna", type=Path, required=True)
    parser.add_argument("--rnaduplex", type=Path, required=True)
    args = parser.parse_args()
    run_item6_thermo(
        args.output_dir, intarna=args.intarna, rnaduplex=args.rnaduplex
    )


if __name__ == "__main__":
    main()
