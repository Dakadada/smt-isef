"""Build the Arm E (StericRuler) length-titration artifact package.

This module renders the per-mer therapeutic-window panel from
:func:`feng.titration.window.compute_window` into a self-describing artifact
bundle: a claim-bounded README, the per-mer score/window table (tsv + parquet),
a length-axis figure, and a provenance record. It mirrors the honesty framing of
the Arm A atlas (``feng.atlas.arm_a``): every value is tagged ``candidate``, a
fixed claim-boundary section states what Arm E does and does not claim, and a
forbidden-language lint runs over the emitted text before the bundle is accepted.

Freeze-safe: consumes only the candidate sequence-derived window rows. It reads
no treatment alignment, contrast, rank, or sealed ΔPSI. The preregistration lock
is reported as an unsealed **DRAFT** -- cutting the real seal is a separate gated
step (see :mod:`feng.titration.window`).
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
from typing import Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from feng.honesty.linter import find_violations  # noqa: E402
from feng.honesty.states import ResultState  # noqa: E402
from feng.titration.window import (  # noqa: E402
    AXIS_DEFINITIONS,
    FORMULA_VERSION,
    normalize_window_rows,
    prereg_hash,
    prereg_payload,
)

plt.rcParams["svg.fonttype"] = "none"  # keep figure text machine-auditable

_CANDIDATE = ResultState.CANDIDATE.value

# --- output package ----------------------------------------------------------
REPORT_NAME = "arm_e_titration_report.md"
TABLE_TSV = "arm_e_per_mer_table.tsv"
TABLE_PARQUET = "arm_e_per_mer_table.parquet"
FIG_PDF = "arm_e_length_axis.pdf"
FIG_SVG = "arm_e_length_axis.svg"
PROVENANCE = "arm_e_provenance.json"
OUTPUT_NAMES = (REPORT_NAME, TABLE_TSV, TABLE_PARQUET, FIG_PDF, FIG_SVG, PROVENANCE)

# Fixed table schema (ladder + evaluate + window keys), stable column order.
TABLE_COLUMNS = (
    "length", "chrom", "start", "end", "strand", "target_5to3", "aso_5to3",
    "intarna_dg", "intarna_no_duplex",
    "rnaduplex_dg", "rnacofold_dg", "rnaup_dg",
    "rnaplfold_accessibility", "rnafold_mfe",
    "off_target_hits", "off_target_burden",
    "potency", "potency_incomplete",
    "therapeutic_window", "window_incomplete",
    "state",
)

DEFAULT_HYPOTHESIS = (
    "a shorter mer improves the therapeutic window: the predicted "
    "potency-over-off-target-burden ratio favours a length below the "
    "18-nt nusinersen anchor"
)

CLAIM_BOUNDARY_HEADING = "## What Arm E does not claim"
CLAIM_BOUNDARY_TEXT = """## What Arm E does not claim

- Arm E does not claim that any predicted therapeutic-window value is a measured potency, a measured specificity, or a wet-lab outcome.
- Arm E does not claim that the sequence-derived CPU scores (duplex energy, accessibility, self-structure, off-target count) establish mechanism, causality, or efficacy.
- Arm E does not claim a clinical response, a dosing recommendation, or that any rung is safe or effective in cells, animals, or patients.
- Arm E does not claim that a shorter or longer oligo is superior to nusinersen; the ladder is scored on its own pre-declared axes only.
- Arm E does not claim the preregistration is sealed: the reported lock is an unsealed DRAFT until the GPU potency axes are included and the seal is cut with explicit sign-off.
- Arm E does not claim the directional prediction is confirmed; confirmation requires the tandem wet-lab length series, read only after the seal."""


def _row_frame(window_rows: Sequence[Mapping]) -> pd.DataFrame:
    """Ordered per-mer table with a fixed schema, sorted by length."""
    records = []
    for row in window_rows:
        rec = {col: row.get(col, None) for col in TABLE_COLUMNS}
        rec["state"] = row.get("state", _CANDIDATE)
        records.append(rec)
    frame = pd.DataFrame.from_records(records, columns=list(TABLE_COLUMNS))
    return frame.sort_values("length").reset_index(drop=True)


def _draw_length_axis(frame: pd.DataFrame, pdf_path: Path, svg_path: Path) -> list[str]:
    """Three stacked length-axis panels: potency, therapeutic window, burden."""
    fig, axes = plt.subplots(
        3, 1, figsize=(9, 9), sharex=True, gridspec_kw={"hspace": 0.18}
    )
    lengths = frame["length"].to_numpy()

    def _plot(axis, column, ylabel, color, title):
        sub = frame.loc[frame[column].notna(), ["length", column]]
        if len(sub):
            axis.plot(sub["length"], sub[column], marker="o", color=color, lw=1.4)
        else:
            axis.text(
                0.5, 0.5, "axis not yet run (off-target scan / GPU pending)",
                ha="center", va="center", transform=axis.transAxes,
                fontsize=9, color="#888888",
            )
        axis.axvline(18, color="#333333", ls=":", lw=0.9)
        axis.set_ylabel(ylabel)
        axis.set_title(title, loc="left", fontweight="bold", fontsize=10)
        return title

    titles = [
        _plot(axes[0], "potency", "potency\n(-mean on-target ΔG)", "#08519c",
              "A  On-target potency vs length [candidate]"),
        _plot(axes[1], "therapeutic_window", "therapeutic window\n(potency / burden)",
              "#238b45", "B  Therapeutic window vs length [candidate]"),
        _plot(axes[2], "off_target_burden", "off-target burden\n(near-complement hits)",
              "#d7301f", "C  Off-target burden vs length [candidate]"),
    ]
    axes[2].set_xlabel("oligo length (nt); dotted line = 18-nt nusinersen anchor")
    axes[0].set_xlim(min(lengths) - 1, max(lengths) + 1)
    caption = (
        "Sequence-derived candidate scores over the symmetric nusinersen length "
        "ladder. The dotted line marks the 18-nt drug anchor. Values are model "
        "outputs, not measurements. [candidate]"
    )
    fig.text(0.01, 0.005, caption, ha="left", va="bottom", fontsize=8)
    fig.suptitle("Arm E StericRuler length titration", x=0.01, ha="left",
                 fontsize=14, fontweight="bold")
    fig.subplots_adjust(left=0.14, right=0.98, top=0.93, bottom=0.09)
    fig.savefig(pdf_path)
    fig.savefig(svg_path)
    plt.close(fig)
    return titles + [caption, "Arm E StericRuler length titration"]


def _markdown_table(frame: pd.DataFrame) -> str:
    columns = list(frame.columns)
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for row in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(v).replace("|", "\\|") for v in row) + " |")
    return "\n".join(lines)


def _summary_columns(frame: pd.DataFrame) -> pd.DataFrame:
    keep = ["length", "aso_5to3", "intarna_dg", "intarna_no_duplex", "rnaup_dg",
            "potency", "off_target_burden", "therapeutic_window", "window_incomplete"]
    return frame[[c for c in keep if c in frame.columns]]


def write_report(
    frame: pd.DataFrame,
    *,
    hypothesis: str,
    prereg_hash_hex: str,
    output_path: Path,
) -> None:
    """Render the claim-bounded Arm E README."""
    versions = {
        name: importlib.metadata.version(name)
        for name in ("numpy", "pandas", "matplotlib")
    }
    n_mers = len(frame)
    n_window = int(frame["therapeutic_window"].notna().sum())
    axis_lines = "\n".join(
        f"- [candidate] **{key}**: {value}"
        for key, value in AXIS_DEFINITIONS.items()
        if key != "formula_version"
    )
    report = f"""# Arm E StericRuler length-titration report

> Claim boundary: predictive, pre-registration-eligible, and hypothesis-driven. Every value in this package is a sequence-derived model output tagged `candidate`. Nothing here is a measurement or a wet-lab outcome.

## What Arm E claims

- [candidate] Over a symmetric truncation/extension ladder of the true nusinersen 18-mer (even lengths, centred on the ISS-N1 target `chr5:70076584-70076601`, GRCh38 v50), Arm E emits one pre-declared **therapeutic window** value per length: predicted on-target potency divided by predicted off-target burden.
- [candidate] The predicted window-vs-length curve and the directional hypothesis are captured in a preregistration payload whose SHA-256 is recorded here, to be sealed and compared against a tandem wet-lab length series.
- [candidate] Directional hypothesis under test: {hypothesis}.

## Methods

- [candidate] Ladder frame: each even length `L` takes the `L`-nt window centred on the 18-nt target, trimming (`L<18`) or extending (`L>18`) one nucleotide per end; the ASO is the reverse complement of that genomic window. The `L=18` rung reproduces nusinersen exactly (fail-closed in-code assertion).
- [candidate] On-target axes: IntaRNA duplex ΔG and RNAup accessibility-corrected ΔG form the potency backbone. IntaRNA's `no_duplex` (0.0) and sub-seed `None` readings are never averaged as real binding. RNAduplex/RNAcofold (a perfect reverse-complement duplex by construction) and RNAplfold accessibility (already folded into RNAup) are carried as secondary/diagnostic columns and are excluded from the backbone to avoid double-counting.
- [candidate] Off-target axis: near-complement hit count from the seed-and-extend scan. A missing scan is reported as `None`, never as zero burden.
- [candidate] Endpoint definitions (all candidate, pre-registration-eligible):

{axis_lines}

- [candidate] Coverage: {n_mers} rungs scored; {n_window} have a computable therapeutic window (the rest await the off-target scan and/or GPU potency axes).
- [candidate] Preregistration formula version `{FORMULA_VERSION}`; payload SHA-256 `{prereg_hash_hex}`. The lock is an unsealed DRAFT.
- [candidate] Software: Python {platform.python_version()}, NumPy {versions['numpy']}, pandas {versions['pandas']}, matplotlib {versions['matplotlib']}.
- [candidate] File-access audit: this generator consumes only the in-memory candidate window rows. It reads no sealed material and reuses no earlier-arm code.

## Per-mer summary [candidate]

{_markdown_table(_summary_columns(frame))}

## Figure caption

- [candidate] `arm_e_length_axis`: on-target potency, therapeutic window, and off-target burden over the symmetric nusinersen length ladder. The dotted line marks the 18-nt drug anchor. Panels whose axis has not been run are annotated as pending; values are model outputs, not measurements.

## Limitations

- [candidate] CPU-only potency (duplex energy + accessibility) rises with length by construction; the directional test lives in the full therapeutic window (÷ off-target burden) and the GPU potency axes, not in raw binding energy.
- [candidate] RNAduplex/RNAcofold score a perfect reverse-complement duplex and carry little signal beyond composition; they are diagnostic only.
- [candidate] The off-target burden counts near-complement sequence sites and does not model expression, chromatin, or delivery; it is a sequence-specificity proxy.

{CLAIM_BOUNDARY_TEXT}
"""
    output_path.write_text(report, encoding="utf-8")


def lint_scientific_outputs(output_dir: Path, figure_text: Iterable[str] = ()) -> list[str]:
    """Return claim-boundary and forbidden-language violations, or an empty list."""
    output_dir = Path(output_dir)
    violations: list[str] = []

    md = (output_dir / REPORT_NAME).read_text(encoding="utf-8")
    if CLAIM_BOUNDARY_TEXT not in md:
        violations.append(f"{REPORT_NAME}: fixed claim-boundary text differs")
    if CLAIM_BOUNDARY_HEADING not in md:
        violations.append(f"{REPORT_NAME}: missing 'does not claim' section")

    # The fixed claim-boundary box is the one exempt region: it must name the
    # things Arm E does NOT claim (clinical response, safe, superior...), so it
    # is removed before the forbidden-language scan. Everything else -- before or
    # after the box -- is still scanned.
    md_scanned = md.replace(CLAIM_BOUNDARY_TEXT, "")
    texts: list[tuple[str, str]] = [(REPORT_NAME, md_scanned)]
    tsv_path = output_dir / TABLE_TSV
    if tsv_path.exists():
        header = tsv_path.read_text(encoding="utf-8").splitlines()[:1]
        texts.append((TABLE_TSV, "\n".join(header)))
    texts.append(("figure text", "\n".join(figure_text)))

    # Forbidden-string scan over the emitted prose.
    field_names = list(TABLE_COLUMNS)
    for label, text in texts:
        for item in find_violations(text, field_names if label == REPORT_NAME else ()):
            violations.append(f"{label}: {item.kind}:{item.value}")
    return violations


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_arm_e_report(
    window_rows: Sequence[Mapping],
    *,
    out_dir: str | Path,
    hypothesis: str = DEFAULT_HYPOTHESIS,
    input_provenance: Mapping | None = None,
) -> pd.DataFrame:
    """Generate the six-file Arm E length-titration artifact package.

    ``window_rows`` are the output of :func:`feng.titration.window.compute_window`.
    Writes README + per-mer table (tsv/parquet) + length-axis figure (pdf/svg) +
    provenance under ``out_dir``, runs the honesty lint, and returns the table.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stale = [p for p in out_dir.iterdir() if p.name not in OUTPUT_NAMES]
    if stale:
        raise ValueError(
            f"output directory contains non-package files: {[p.name for p in stale]}"
        )

    normalized_rows = normalize_window_rows(window_rows)
    frame = _row_frame(normalized_rows)

    # DRAFT preregistration payload hash (no seal cut here).
    payload = prereg_payload(normalized_rows, hypothesis=hypothesis)
    payload_hash = prereg_hash(payload)

    frame.to_csv(out_dir / TABLE_TSV, sep="\t", index=False, float_format="%.6g")
    frame.to_parquet(out_dir / TABLE_PARQUET, index=False)
    figure_text = _draw_length_axis(frame, out_dir / FIG_PDF, out_dir / FIG_SVG)
    write_report(
        frame, hypothesis=hypothesis, prereg_hash_hex=payload_hash,
        output_path=out_dir / REPORT_NAME,
    )

    provenance = {
        "artifact": "arm_e_titration",
        "formula_version": FORMULA_VERSION,
        "hypothesis": hypothesis,
        "prereg_hash": payload_hash,
        "prereg_locked": False,
        "prereg_note": (
            "unsealed DRAFT; the real locked=true seal is a gated step pending "
            "the GPU potency axes and explicit sign-off"
        ),
        "directional_prediction": payload["directional_prediction"],
        "directional_prediction_axis": payload["directional_prediction_axis"],
        "curve_complete": payload["curve_complete"],
        "n_rungs": int(len(frame)),
        "n_with_window": int(frame["therapeutic_window"].notna().sum()),
        "table_sha256": _sha256_file(out_dir / TABLE_PARQUET),
        "software": {
            "python": platform.python_version(),
            **{
                name: importlib.metadata.version(name)
                for name in ("numpy", "pandas", "matplotlib")
            },
        },
        "inputs": dict(input_provenance) if input_provenance else {},
        "state": _CANDIDATE,
    }
    (out_dir / PROVENANCE).write_text(json.dumps(provenance, indent=2, sort_keys=True))

    violations = lint_scientific_outputs(out_dir, figure_text)
    if violations:
        raise ValueError("forbidden scientific language: " + "; ".join(violations))
    return frame
