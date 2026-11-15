"""Generate Item 1 candidate artifacts from the locked ASO set."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess

from . import CANDIDATE
from .sequences import ASOS, ISS_N1_CONTEXT, derive_target
from .thermo import effective_monovalent_from_mixed_salt, rnaduplex_score, score_pair


# Operational buffer (user-specified 2026-09-30): 50 mM monovalent + 1.5 mM Mg2+.
# The divalent is folded into an effective monovalent via the von Ahsen (2001)
# Na+-equivalence, then consumed by the existing SantaLucia monovalent correction.
PRIMARY_MONOVALENT_MOLAR = 0.050
PRIMARY_DIVALENT_MOLAR = 0.0015
PRIMARY_SALT_MOLAR = effective_monovalent_from_mixed_salt(
    PRIMARY_MONOVALENT_MOLAR, PRIMARY_DIVALENT_MOLAR
)
SALT_GRID_MOLAR = (0.025, 0.050, 0.089, 0.100, 0.150)
ASO_ORDER = ("F18", "F14", "F10", "3UP8")


def _write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def _version(executable: Path) -> str:
    result = subprocess.run(
        [str(executable), "--version"], capture_output=True, text=True, check=True
    )
    return (result.stdout or result.stderr).splitlines()[0].strip()


def run_item1(
    output_dir: Path,
    *,
    intarna: Path,
    rnaduplex: Path,
    parameter_file: Path | str | None = None,
    chemistry: str = "RNA:RNA",
) -> list[dict]:
    raw_rows: list[dict] = []
    sensitivity_rows: list[dict] = []
    for aso_id in ASO_ORDER:
        aso = ASOS[aso_id]
        target = derive_target(aso.rna_sequence, ISS_N1_CONTEXT)
        primary = score_pair(
            aso.rna_sequence,
            target.sequence,
            intarna=intarna,
            rnaduplex=rnaduplex,
            salt_molar=PRIMARY_SALT_MOLAR,
            parameter_file=parameter_file,
            chemistry=chemistry,
        )
        raw_rows.append(
            {
                "state": CANDIDATE,
                "aso_id": aso_id,
                "chemistry": chemistry,
                "length_nt": aso.length,
                "aso_rna_5to3": aso.rna_sequence,
                "target_rna_5to3": target.sequence,
                "target_start_1based": target.start_1based,
                "target_end_1based": target.end_1based,
                "target_requires_human_confirmation": target.requires_human_confirmation,
                "intarna_default_dg_kcal_mol": primary.intarna_dg_kcal_mol,
                "intarna_no_duplex": primary.intarna_no_duplex,
                "nn_reference_dg_kcal_mol": primary.nn_raw_reference_dg_kcal_mol,
                "nn_salt_correction_kcal_mol": primary.nn_salt_correction_kcal_mol,
                "nn_buffer_dg_kcal_mol": primary.nn_dg_kcal_mol,
                "nn_buffer_monovalent_molar": PRIMARY_MONOVALENT_MOLAR,
                "nn_buffer_divalent_molar": PRIMARY_DIVALENT_MOLAR,
                "nn_buffer_effective_monovalent_molar": PRIMARY_SALT_MOLAR,
                "temperature_c": primary.temperature_c,
            }
        )
        for salt in SALT_GRID_MOLAR:
            nn = rnaduplex_score(
                aso.rna_sequence,
                target.sequence,
                rnaduplex=rnaduplex,
                salt_molar=salt,
                parameter_file=parameter_file,
                chemistry=chemistry,
            )
            sensitivity_rows.append(
                {
                    "state": CANDIDATE,
                    "aso_id": aso_id,
                    "chemistry": chemistry,
                    "length_nt": aso.length,
                    "salt_molar": salt,
                    "nn_reference_dg_kcal_mol": nn.raw_reference_dg_kcal_mol,
                    "nn_salt_correction_kcal_mol": nn.salt_correction_kcal_mol,
                    "nn_dg_kcal_mol": nn.energy_kcal_mol,
                    "nn_structure": nn.structure,
                    "nn_paired_bases": nn.paired_bases,
                    "temperature_c": nn.temperature_c,
                }
            )

    intarna_ranked = sorted(
        raw_rows, key=lambda row: float(row["intarna_default_dg_kcal_mol"])
    )
    nn_ranked = sorted(raw_rows, key=lambda row: float(row["nn_buffer_dg_kcal_mol"]))
    intarna_rank = {row["aso_id"]: rank for rank, row in enumerate(intarna_ranked, 1)}
    nn_rank = {row["aso_id"]: rank for rank, row in enumerate(nn_ranked, 1)}
    rank_rows = [
        {
            "state": CANDIDATE,
            "aso_id": aso_id,
            "intarna_rank_most_negative_first": intarna_rank[aso_id],
            "nn_buffer_rank_most_negative_first": nn_rank[aso_id],
            "ranks_identical": intarna_rank[aso_id] == nn_rank[aso_id],
        }
        for aso_id in ASO_ORDER
    ]

    output_dir = Path(output_dir)
    _write_tsv(output_dir / "item1_dg_raw.tsv", raw_rows)
    _write_tsv(output_dir / "item1_salt_sensitivity.tsv", sensitivity_rows)
    _write_tsv(output_dir / "item1_dg_ranks.tsv", rank_rows)
    provenance = {
        "state": CANDIDATE,
        "chemistry": chemistry,
        "temperature_c": 37.0,
        "intarna_salt": "default; no explicit salt",
        "nn_parameter_set": "ViennaRNA default Turner 2004",
        "nn_reference_salt_molar": 1.021,
        "salt_correction": (
            "SantaLucia entropy form: dS=0.368*(paired_bases-1)*"
            "ln([M+]/1.021 M); dG=-T*dS"
        ),
        "primary_buffer": "50 mM monovalent + 1.5 mM Mg2+ (user-specified 2026-09-30)",
        "primary_monovalent_molar": PRIMARY_MONOVALENT_MOLAR,
        "primary_divalent_molar": PRIMARY_DIVALENT_MOLAR,
        "divalent_equivalence": (
            "von Ahsen 2001 Na+-equivalent: [Na+]_eq(mM)=[mono](mM)+120*"
            "sqrt([Mg2+](mM)); DNA-derived empirical rule applied to RNA:RNA "
            "-> heuristic salt-adjusted dG, not a validated RNA measurement"
        ),
        "primary_effective_monovalent_molar": PRIMARY_SALT_MOLAR,
        "salt_grid_molar": list(SALT_GRID_MOLAR),
        "intarna_version": _version(Path(intarna)),
        "rnaduplex_version": _version(Path(rnaduplex)),
        "context_rna_5to3": ISS_N1_CONTEXT,
        "ranks_identical_all_asos": all(row["ranks_identical"] for row in rank_rows),
    }
    (output_dir / "item1_provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )
    return raw_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--intarna", type=Path, required=True)
    parser.add_argument("--rnaduplex", type=Path, required=True)
    parser.add_argument("--parameter-file", type=Path)
    parser.add_argument("--chemistry", choices=("RNA:RNA", "DNA:RNA"), default="RNA:RNA")
    args = parser.parse_args()
    run_item1(
        args.output_dir,
        intarna=args.intarna,
        rnaduplex=args.rnaduplex,
        parameter_file=args.parameter_file,
        chemistry=args.chemistry,
    )


if __name__ == "__main__":
    main()
