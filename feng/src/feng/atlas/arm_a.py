"""Build the descriptive Arm A SMN2 dCas13d/ASO field atlas.

This module deliberately performs only coordinate alignment and descriptive summaries.
The three publication sources remain separate throughout.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import hashlib
import importlib.metadata
import math
from pathlib import Path
import platform
import re
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd

plt.rcParams["svg.fonttype"] = "none"  # Keep figure language machine-auditable.


SPAN_START = 260
SPAN_END = 434
ARTIFICIAL_G = 356
QUANTILE_METHOD = "linear"
STUDIES = ("Hua 2007 PLoS Biology", "Hua 2008 AJHG", "Singh 2010 RNA")
REGIONS = (("I6", 260, 319), ("E7", 320, 374), ("I7", 375, 434))
BLOCKS = (
    ("I6-1", "I6", 260, 274), ("I6-2", "I6", 275, 289),
    ("I6-3", "I6", 290, 304), ("I6-4", "I6", 305, 319),
    ("E7-1", "E7", 320, 334), ("E7-2", "E7", 335, 350),
    ("E7-3", "E7", 351, 374), ("I7-1", "I7", 375, 389),
    ("I7-2", "I7", 390, 404), ("I7-3", "I7", 405, 419),
    ("I7-4", "I7", 420, 434),
)
DIRECTION_MAP = {
    "increase": "inclusion increase",
    "decrease": "inclusion decrease",
    "neutral": "neutral",
    "discordant": "assay-discordant",
    "neutral_or_decrease": "mixed/uncertain",
    "figure_reported": "unresolved, not digitized",
}
DIRECTION_COLUMNS = (
    "inclusion increase", "inclusion decrease", "neutral", "assay-discordant",
    "mixed/uncertain", "unresolved, not digitized",
)
DIRECTION_COLORS = {
    "inclusion increase": "#2b8cbe",
    "inclusion decrease": "#d7301f",
    "neutral": "#777777",
    "assay-discordant": "#7b3294",
    "mixed/uncertain": "#e08214",
    "unresolved, not digitized": "#bdbdbd",
}
DIRECTION_HATCHES = {
    "assay-discordant": "///",
    "mixed/uncertain": "xxx",
    "unresolved, not digitized": "...",
}
FORBIDDEN_SCIENTIFIC_STRINGS = (
    "predict", "prediction", "transfer", "held-out", "precision", "fdr", "p-value"
)
CLAIM_BOUNDARY_HEADING = "## What Arm A no longer claims"
CLAIM_BOUNDARY_TEXT = """## What Arm A no longer claims

- Arm A no longer claims that dCas13d guide activity predicts ASO splice modulation.
- Arm A no longer claims a cross-effector transfer law or a kernel validated on held-out ASOs.
- Arm A no longer claims that a correlation, p-value, precision, or FDR from the historical walks establishes predictive performance.
- Arm A no longer claims that 115 exact ASO windows are 115 independent observations or satisfy the minimum effective-sample gate.
- Arm A no longer claims that Hua 2007, Hua 2008, and Singh 2010 form one poolable numeric response.
- Arm A no longer claims that descriptive same-direction co-localization establishes mechanism, causality, portability beyond this reporter field, or clinical utility."""
SIGN_SOURCE = (
    "Recinos et al. 2024, Nature Communications, Figure 5a and Methods equations 3-5; "
    "https://doi.org/10.1038/s41467-024-47140-8"
)
OUTPUT_NAMES = (
    "arm_a_shared_field_map.pdf", "arm_a_shared_field_map.svg",
    "arm_a_coarse_block_map.pdf", "arm_a_coarse_block_map.svg",
    "arm_a_window_alignment.tsv", "arm_a_block_summary.tsv",
    "arm_a_descriptive_report.md",
)


@dataclass(frozen=True)
class AtlasData:
    """Validated inputs and derived descriptive tables."""

    guides: pd.DataFrame
    asos: pd.DataFrame
    positional: pd.DataFrame
    windows: pd.DataFrame
    blocks: pd.DataFrame


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _quantiles(values: Iterable[float]) -> tuple[float, float, float]:
    arr = np.asarray(list(values), dtype=float)
    return tuple(float(x) for x in np.quantile(arr, (0.25, 0.5, 0.75), method=QUANTILE_METHOD))


def parse_segments(value: str) -> tuple[tuple[int, int], ...]:
    segments = []
    for part in str(value).split(";"):
        match = re.fullmatch(r"(\d+)-(\d+)", part.strip())
        if not match:
            raise ValueError(f"invalid transformed segment: {value!r}")
        start, end = map(int, match.groups())
        if start > end:
            raise ValueError(f"reversed transformed segment: {value!r}")
        segments.append((start, end))
    return tuple(segments)


def _published_interval(value: str) -> tuple[int, int]:
    numbers = [int(x) for x in re.findall(r"[+-]?\d+", str(value))]
    if len(numbers) != 2:
        raise ValueError(f"invalid original coordinate: {value!r}")
    return tuple(numbers)


def expected_s4_geometry(region: str, original: str) -> tuple[int, int, str]:
    """Return corrected S4 hull/segments directly from a paper-local interval."""
    left, right = _published_interval(original)
    if region == "I6":
        start, end = 320 + left, 320 + right
        segments = f"{start}-{end}"
    elif region == "I7":
        start, end = 374 + left, 374 + right
        segments = f"{start}-{end}"
    elif region == "E7":
        def convert(pos: int) -> int:
            return 319 + pos if pos <= 36 else 320 + pos

        start, end = convert(left), convert(right)
        if left <= 36 < right:
            segments = f"{start}-355;357-{end}"
        else:
            segments = f"{start}-{end}"
    else:
        raise ValueError(f"unknown ASO region: {region!r}")
    return start, end, segments


def load_inputs(guide_csv: Path, aso_csv: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    guides = pd.read_csv(guide_csv)
    required_guides = {"OligoID", "minigene", "start_pos", "end_pos", "length", "region", "z_combined"}
    missing = required_guides - set(guides.columns)
    if missing:
        raise ValueError(f"guide input missing columns: {sorted(missing)}")
    guides = guides.loc[
        (guides["minigene"] == "SMN2")
        & (guides["start_pos"] <= SPAN_END)
        & (guides["end_pos"] >= SPAN_START)
    ].copy()
    guides[["start_pos", "end_pos", "length"]] = guides[["start_pos", "end_pos", "length"]].astype(int)
    if guides.empty or guides["z_combined"].isna().any():
        raise ValueError("retained guide field is empty or has missing combined Z values")
    if not set(guides["length"]).issubset({22, 25}):
        raise ValueError(f"unexpected guide lengths: {sorted(set(guides['length']))}")

    asos = pd.read_csv(aso_csv, dtype={"aso_id": str})
    required_asos = {
        "source", "source_locator", "aso_id", "original_coordinate", "region",
        "transformed_minigene_start", "transformed_minigene_end", "transformed_segments",
        "effect_direction", "effect_magnitude", "distinct_window_key",
    }
    missing = required_asos - set(asos.columns)
    if missing:
        raise ValueError(f"ASO input missing columns: {sorted(missing)}")
    if len(asos) != 125:
        raise ValueError(f"expected 125 ASO study records, found {len(asos)}")
    if tuple(asos["source"].drop_duplicates()) != STUDIES:
        raise ValueError("ASO studies or their frozen order differ from the specification")
    unknown = sorted(set(asos["effect_direction"].dropna()) - set(DIRECTION_MAP))
    if unknown or asos["effect_direction"].isna().any():
        raise ValueError(f"unrecognized ASO direction values: {unknown or ['<missing>']}")

    # Load-bearing +9 check: recompute every coordinate from the paper-local frame.
    errors = []
    for idx, row in asos.iterrows():
        expected = expected_s4_geometry(row["region"], row["original_coordinate"])
        observed = (
            int(row["transformed_minigene_start"]),
            int(row["transformed_minigene_end"]),
            str(row["transformed_segments"]),
        )
        if observed != expected:
            errors.append((idx + 2, observed, expected))
    if errors:
        line, observed, expected = errors[0]
        raise ValueError(
            "ASO coordinates are not in the corrected S4 frame; "
            f"first mismatch at CSV line {line}: observed {observed}, expected {expected}. "
            "The pre-correction frame is nine nucleotides lower."
        )
    for value in asos["transformed_segments"]:
        parse_segments(value)
    if asos["distinct_window_key"].nunique() != 115:
        raise ValueError("expected 115 exact S4 windows")
    return guides.reset_index(drop=True), asos.reset_index(drop=True)


def positional_summary(guides: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for position in range(SPAN_START, SPAN_END + 1):
        values = guides.loc[
            (guides["start_pos"] <= position) & (guides["end_pos"] >= position), "z_combined"
        ].to_numpy()
        if not len(values):
            raise ValueError(f"no guide covers S4 position {position}")
        q1, median, q3 = _quantiles(values)
        rows.append({"position": position, "z_q1": q1, "z_median": median, "z_q3": q3, "coverage": len(values)})
    return pd.DataFrame(rows)


def assign_block(start: int, end: int) -> str:
    midpoint = (start + end) / 2
    for index, (block_id, _, block_start, block_end) in enumerate(BLOCKS):
        # Integer coordinates denote base centers. Base-centered edges make adjacent
        # bins exhaustive and unique for half-integer interval midpoints.
        left_edge, right_edge = block_start - 0.5, block_end + 0.5
        if left_edge <= midpoint < right_edge or (index == len(BLOCKS) - 1 and midpoint == right_edge):
            return block_id
    raise ValueError(f"ASO midpoint {midpoint} is outside the fixed display blocks")


def _relation(direction: str, median: float) -> str:
    if direction not in {"increase", "decrease"} or median == 0 or math.isnan(median):
        return "not_evaluable"
    positive = median > 0
    return "same" if (direction == "increase") == positive else "opposite"


def window_alignment(guides: pd.DataFrame, asos: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, aso in asos.iterrows():
        start = int(aso["transformed_minigene_start"])
        end = int(aso["transformed_minigene_end"])
        local = guides.loc[(guides["start_pos"] <= start) & (guides["end_pos"] >= end), "z_combined"]
        if local.empty:
            q1 = median = q3 = minimum = maximum = np.nan
        else:
            q1, median, q3 = _quantiles(local)
            minimum, maximum = float(local.min()), float(local.max())
        rows.append({
            "source": aso["source"], "source_locator": aso["source_locator"], "aso_id": aso["aso_id"],
            "region": aso["region"], "start": start, "end": end,
            "segments": aso["transformed_segments"], "distinct_window_key": aso["distinct_window_key"],
            "block_id": assign_block(start, end), "display_direction": DIRECTION_MAP[aso["effect_direction"]],
            "direction_eligible": str(aso["effect_direction"] in {"increase", "decrease"}).lower(),
            "n_guides_containing_hull": len(local), "local_z_median": median, "local_z_q1": q1,
            "local_z_q3": q3, "local_z_min": minimum, "local_z_max": maximum,
            "local_z_positive_n": int((local > 0).sum()), "local_z_negative_n": int((local < 0).sum()),
            "local_z_zero_n": int((local == 0).sum()),
            "direction_relation": _relation(aso["effect_direction"], median),
            "result_state": "candidate",
        })
    return pd.DataFrame(rows)


def block_summary(positional: pd.DataFrame, windows: pd.DataFrame) -> pd.DataFrame:
    rows = []
    reverse_display = {value: key for key, value in DIRECTION_MAP.items()}
    for block_id, region, start, end in BLOCKS:
        p = positional.loc[positional["position"].between(start, end)]
        zq1, zmedian, zq3 = _quantiles(p["z_median"])
        coverage = p["coverage"].to_numpy()
        for study in STUDIES:
            cell = windows.loc[(windows["block_id"] == block_id) & (windows["source"] == study)]
            counts = Counter(cell["display_direction"])
            eligible = int(counts["inclusion increase"] + counts["inclusion decrease"])
            if zmedian > 0:
                same, opposite = int(counts["inclusion increase"]), int(counts["inclusion decrease"])
            elif zmedian < 0:
                same, opposite = int(counts["inclusion decrease"]), int(counts["inclusion increase"])
            else:
                same = opposite = 0
            row = {
                "block_id": block_id, "region": region, "start": start, "end": end, "study": study,
                "total_aso_records": len(cell), "distinct_exact_windows": cell["distinct_window_key"].nunique(),
            }
            for label in DIRECTION_COLUMNS:
                row[f"{reverse_display[label]}_n"] = counts[label]
            row.update({
                "eligible_direction_n": eligible, "same_n": same, "opposite_n": opposite,
                "same_direction_fraction": same / eligible if eligible else "NA",
                "block_z_positional_median": zmedian, "block_z_positional_q1": zq1,
                "block_z_positional_q3": zq3, "block_z_positional_min": float(p["z_median"].min()),
                "block_z_positional_max": float(p["z_median"].max()), "n_integer_positions": len(p),
                "coverage_median": float(np.median(coverage)), "coverage_min": int(coverage.min()),
                "coverage_max": int(coverage.max()), "result_state": "candidate",
            })
            rows.append(row)
    return pd.DataFrame(rows)


def prepare_atlas(guide_csv: Path, aso_csv: Path) -> AtlasData:
    guides, asos = load_inputs(guide_csv, aso_csv)
    positional = positional_summary(guides)
    windows = window_alignment(guides, asos)
    blocks = block_summary(positional, windows)
    return AtlasData(guides, asos, positional, windows, blocks)


def _style_axis(axis) -> None:
    axis.set_xlim(SPAN_START - 0.5, SPAN_END + 0.5)
    axis.xaxis.set_major_locator(MultipleLocator(5))
    axis.xaxis.set_minor_locator(MultipleLocator(1))
    axis.axvline(319.5, color="#333333", lw=0.8)
    axis.axvline(374.5, color="#333333", lw=0.8)
    axis.axvline(ARTIFICIAL_G, color="#8856a7", lw=0.8, ls=":")


def _draw_landmarks(axis) -> None:
    colors = {"I6": "#deebf7", "E7": "#c7e9c0", "I7": "#fee0d2"}
    for region, start, end in REGIONS:
        axis.add_patch(Rectangle((start - 0.5, 0), end - start + 1, 1, facecolor=colors[region], edgecolor="#555555"))
        axis.text((start + end) / 2, 0.5, f"{region} {start}-{end}", ha="center", va="center", fontsize=9)
    axis.text(ARTIFICIAL_G, 1.08, "artificial exon-7 G (356)", ha="center", va="bottom", fontsize=7, color="#6a3d9a")
    axis.set_ylim(0, 1.35)
    axis.set_yticks([])
    axis.set_title("A  Construct landmarks [paper_supported]", loc="left", fontweight="bold")
    _style_axis(axis)
    axis.tick_params(axis="x", which="both", labelbottom=True)


def _stack_rows(study_asos: pd.DataFrame) -> list[int]:
    # A dedicated row keeps the complete source label legible in both vector formats.
    # It is a conservative interval-stack: overlapping records can never share a row.
    return list(range(len(study_asos)))


def draw_shared_field_map(data: AtlasData, pdf_path: Path, svg_path: Path) -> list[str]:
    study_frames = [data.asos.loc[data.asos["source"] == study].copy() for study in STUDIES]
    lane_counts = [max(_stack_rows(frame)) + 1 for frame in study_frames]
    height_ratios = [1.2, 6, 3, 1.1] + [max(3, n * 0.42) for n in lane_counts]
    fig, axes = plt.subplots(
        7, 1, figsize=(30, sum(height_ratios) * 0.60), sharex=True,
        gridspec_kw={"height_ratios": height_ratios, "hspace": 0.32},
    )
    _draw_landmarks(axes[0])

    guides = data.guides
    vmax = float(np.abs(guides["z_combined"]).max())
    norm = mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    for length, marker in ((22, "o"), (25, "^")):
        subset = guides.loc[guides["length"] == length]
        axes[1].scatter(
            (subset["start_pos"] + subset["end_pos"]) / 2, subset["z_combined"],
            c=subset["z_combined"], cmap="RdBu", norm=norm, marker=marker,
            s=15, alpha=0.72, linewidths=0.2, edgecolors="#222222", rasterized=False,
        )
    axes[1].axhline(0, color="#222222", lw=0.8)
    axes[1].set_ylabel("raw combined Z")
    axes[1].set_title(f"B  All dCas13d guides: {len(guides)} retained [paper_supported]", loc="left", fontweight="bold")
    axes[1].legend(handles=[Line2D([], [], marker="o", ls="", color="#555555", label="22 nt"), Line2D([], [], marker="^", ls="", color="#555555", label="25 nt")], loc="upper right")
    colorbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap="RdBu"), ax=axes[1], pad=0.005, fraction=0.012)
    colorbar.set_label("raw combined Z", fontsize=8)
    _style_axis(axes[1])

    p = data.positional
    axes[2].fill_between(p["position"], p["z_q1"], p["z_q3"], color="#9ecae1", alpha=0.65, label="Q1-Q3")
    axes[2].plot(p["position"], p["z_median"], color="#08519c", lw=1.2, label="median")
    axes[2].axhline(0, color="#222222", lw=0.7)
    axes[2].set_ylabel("raw combined Z")
    axes[2].set_title("C  Per-base guide summary (no smoothing) [candidate]", loc="left", fontweight="bold")
    axes[2].legend(loc="upper right")
    _style_axis(axes[2])

    axes[3].fill_between(p["position"], 0, p["coverage"], step="mid", color="#636363")
    axes[3].set_ylabel("guides")
    axes[3].set_title("Guide interval coverage |Gx| [candidate]", loc="left", fontsize=9)
    _style_axis(axes[3])

    panel_names = ("D", "E", "F")
    for axis, study, frame, panel in zip(axes[4:], STUDIES, study_frames, panel_names):
        lanes = _stack_rows(frame)
        for (_, row), lane in zip(frame.iterrows(), lanes):
            display = DIRECTION_MAP[row["effect_direction"]]
            color = DIRECTION_COLORS[display]
            hatch = DIRECTION_HATCHES.get(display)
            for start, end in parse_segments(row["transformed_segments"]):
                axis.add_patch(Rectangle((start - 0.5, lane - 0.34), end - start + 1, 0.68, facecolor=color, edgecolor="#222222", hatch=hatch, lw=0.45))
            label = f"{row['aso_id']} | {row['original_coordinate']} | {row['source_locator']} | {row['effect_magnitude']}"
            axis.text(SPAN_END + 2, lane, label, va="center", fontsize=5.5, clip_on=False)
        axis.set_ylim(max(lanes) + 0.8, -0.8)
        axis.set_yticks(range(max(lanes) + 1))
        axis.set_yticklabels([f"stack {x + 1}" for x in range(max(lanes) + 1)], fontsize=6)
        axis.set_ylabel("overlap stack")
        axis.set_title(f"{panel}  {study}: {len(frame)} study records [paper_supported]", loc="left", fontweight="bold")
        _style_axis(axis)
    axes[-1].set_xlabel("corrected S4 SMN2 minigene coordinate (inclusive)")

    legend = [Patch(facecolor=DIRECTION_COLORS[label], edgecolor="#222222", hatch=DIRECTION_HATCHES.get(label), label=label) for label in DIRECTION_COLUMNS]
    axes[4].legend(handles=legend, loc="upper left", bbox_to_anchor=(1.0, 1.25), fontsize=7)
    caption = (
        "dCas13d values are numeric source measurements; ASO colors are source-reported qualitative classes. "
        "The studies differ in assay and chemistry. Visual proximity or same-direction color is hypothesis-generating only. "
        "[candidate]"
    )
    fig.text(0.01, 0.003, caption, ha="left", va="bottom", fontsize=8)
    fig.suptitle("Arm A descriptive SMN2 field atlas", x=0.01, ha="left", fontsize=15, fontweight="bold")
    fig.subplots_adjust(left=0.055, right=0.55, top=0.98, bottom=0.025)
    fig.savefig(pdf_path)
    fig.savefig(svg_path)
    plt.close(fig)
    return [axis.get_title() for axis in axes] + [caption, "Arm A descriptive SMN2 field atlas"]


def draw_block_map(data: AtlasData, pdf_path: Path, svg_path: Path) -> list[str]:
    fig, axes = plt.subplots(7, 1, figsize=(18, 11), sharex=True, gridspec_kw={"height_ratios": [1.0, 1.7, 1.7, 1.7, 1, 1, 1], "hspace": 0.16})
    zvalues = []
    for block_id, _, start, end in BLOCKS:
        zvalues.append(float(data.blocks.loc[data.blocks["block_id"] == block_id, "block_z_positional_median"].iloc[0]))
    vmax = max(abs(x) for x in zvalues)
    axes[0].imshow([zvalues], aspect="auto", cmap="RdBu", norm=mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax), extent=(-0.5, 10.5, -0.5, 0.5))
    for col, value in enumerate(zvalues):
        axes[0].text(col, 0, f"{value:.2f}", ha="center", va="center", fontsize=8)
    axes[0].set_yticks([0], ["block median\nof positional medians"])
    axes[0].set_title("Eleven coarse display blocks [candidate]", loc="left", fontweight="bold")

    short = {"Hua 2007 PLoS Biology": "Hua 2007", "Hua 2008 AJHG": "Hua 2008", "Singh 2010 RNA": "Singh 2010"}
    abbreviations = (("increase_n", "inc"), ("decrease_n", "dec"), ("neutral_n", "neu"), ("discordant_n", "disc"), ("neutral_or_decrease_n", "mixed"), ("figure_reported_n", "unres"))
    for axis, study in zip(axes[1:4], STUDIES):
        axis.set_facecolor("#fafafa")
        for col, (block_id, _, _, _) in enumerate(BLOCKS):
            row = data.blocks.loc[(data.blocks["block_id"] == block_id) & (data.blocks["study"] == study)].iloc[0]
            text = "\n".join(f"{label} {int(row[field])}" for field, label in abbreviations)
            axis.text(col, 0, text, ha="center", va="center", fontsize=7)
        axis.set_yticks([0], [f"{short[study]}\ndirection counts"])
        axis.set_ylim(-0.5, 0.5)

    for axis, study in zip(axes[4:], STUDIES):
        axis.set_facecolor("#eef5f9")
        for col, (block_id, _, _, _) in enumerate(BLOCKS):
            row = data.blocks.loc[(data.blocks["block_id"] == block_id) & (data.blocks["study"] == study)].iloc[0]
            if int(row["eligible_direction_n"]) == 0:
                text = "NA (0 eligible)\nopp 0"
            else:
                text = f"{int(row['same_n'])}/{int(row['eligible_direction_n'])}\nopp {int(row['opposite_n'])}"
            axis.text(col, 0, text, ha="center", va="center", fontsize=8)
        axis.set_yticks([0], [f"{short[study]}\nsame fraction"])
        axis.set_ylim(-0.5, 0.5)
    for axis in axes:
        axis.set_xlim(-0.5, 10.5)
        axis.set_xticks(np.arange(11))
        axis.grid(axis="x", color="white", lw=1.5)
        for spine in axis.spines.values():
            spine.set_visible(False)
    labels = [f"{bid}\n{start}-{end}" for bid, _, start, end in BLOCKS]
    axes[-1].set_xticklabels(labels, fontsize=8)
    axes[-1].set_xlabel("fixed coarse display block (inclusive corrected S4 coordinates)")
    caption = "Counts and sign relations are descriptive and study-separated; each block is a display bin, not an independent experiment. [candidate]"
    fig.text(0.01, 0.01, caption, fontsize=8)
    fig.subplots_adjust(left=0.14, right=0.99, top=0.96, bottom=0.10)
    fig.savefig(pdf_path)
    fig.savefig(svg_path)
    plt.close(fig)
    return [axis.get_title() for axis in axes] + [caption]


def _markdown_table(frame: pd.DataFrame) -> str:
    columns = list(frame.columns)
    lines = ["| " + " | ".join(columns) + " |", "| " + " | ".join("---" for _ in columns) + " |"]
    for row in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(str(value).replace("|", "\\|") for value in row) + " |")
    return "\n".join(lines)


def _inventory_tables(data: AtlasData) -> str:
    aso = data.asos.copy()
    aso["display_direction"] = aso["effect_direction"].map(DIRECTION_MAP)
    records = aso.groupby(["source", "region"], observed=True).size().rename("records").reset_index()
    exact = aso.groupby(["source", "region"], observed=True)["distinct_window_key"].nunique().rename("exact_windows").reset_index()
    inv = records.merge(exact, on=["source", "region"])
    all_cells = pd.MultiIndex.from_product([STUDIES, [r[0] for r in REGIONS]], names=["source", "region"]).to_frame(index=False)
    inv = all_cells.merge(inv, how="left").fillna(0)
    inv[["records", "exact_windows"]] = inv[["records", "exact_windows"]].astype(int)

    dirs = aso.groupby(["source", "region", "display_direction"], observed=True).size().unstack(fill_value=0)
    dirs = dirs.reindex(pd.MultiIndex.from_product([STUDIES, [r[0] for r in REGIONS]], names=["source", "region"]), fill_value=0)
    dirs = dirs.reindex(columns=DIRECTION_COLUMNS, fill_value=0).reset_index()
    dirs["non_evaluable"] = dirs[["neutral", "assay-discordant", "mixed/uncertain", "unresolved, not digitized"]].sum(axis=1)
    dirs["study_region_records_denominator"] = dirs[list(DIRECTION_COLUMNS)].sum(axis=1)

    guides = data.guides.groupby(["region", "length"], dropna=False).size().rename("guides").reset_index()
    guides["length"] = guides["length"].astype(int)
    guides["retained_guides_denominator"] = len(data.guides)
    coverage_rows = []
    for region, start, end in REGIONS:
        c = data.positional.loc[data.positional["position"].between(start, end), "coverage"]
        coverage_rows.append({"region": region, "positions": len(c), "coverage_median": float(c.median()), "coverage_min": int(c.min()), "coverage_max": int(c.max())})
    return (
        "### ASO records and exact windows by study and region [paper_supported]\n\n" + _markdown_table(inv) +
        "\n\n### Direction classes and non-evaluable records by study and region [paper_supported]\n\n" + _markdown_table(dirs) +
        "\n\n### Retained guides by source annotation and length [paper_supported]\n\n" + _markdown_table(guides) +
        "\n\n### Per-base guide coverage by construct region [candidate]\n\n" + _markdown_table(pd.DataFrame(coverage_rows))
    )


def write_report(data: AtlasData, guide_csv: Path, aso_csv: Path, output_path: Path) -> None:
    versions = {name: importlib.metadata.version(name) for name in ("numpy", "pandas", "matplotlib")}
    study_rel = data.blocks[["block_id", "study", "eligible_direction_n", "same_n", "opposite_n", "same_direction_fraction"]]
    # Derive the inventory counts from data so the prose cannot silently drift from the artifact.
    n_records = len(data.windows)
    n_exact_windows = int(data.windows["distinct_window_key"].nunique())
    n_guide_covered_windows = int(
        data.windows.loc[data.windows["n_guides_containing_hull"] > 0, "distinct_window_key"].nunique()
    )
    report = f"""# Arm A descriptive SMN2 field atlas

> Claim boundary: descriptive and hypothesis-generating only. Allowed result states in this package are `candidate` and `paper_supported`.

## Manifest and coordinate reconciliation

- [paper_supported] Input ASO inventory: {n_records} records / {n_exact_windows} exact windows / {n_guide_covered_windows} guide-covered exact windows. The exact-window count has no independence interpretation.
- [paper_supported] The source audit found that the earlier ASO artifact was nine nucleotides below the SpliceRUSH S4 field frame because it omitted the S4 5-prime prefix `ATGGCTAGC`. This package uses the corrected input coordinates: I6 maps as `320 + q` (within I6 the offset q is negative, placing I6 positions below the E7 start at 320), natural E7 starts at 320 with the artificial G at 356, and I7 maps as `374 + q`.
- [paper_supported] Sequence-alignment evidence recorded in the source memo: guide 1001 contains the prefix before the Dataset S1 exon-6 sequence; guide 1340 aligns at S1-concatenated 331-355 and S4 340-364; guide 1384 aligns at S1-concatenated 375-399 and S4 384-408. The generator independently recomputes all {n_records} corrected hulls and segments from paper-local coordinates and stops on any mismatch.
- [candidate] Display span: inclusive corrected S4 positions {SPAN_START}-{SPAN_END}; construct regions I6 260-319, E7 320-374, and I7 375-434.
- [paper_supported] dCas13d sign semantics are verified: positive combined Z denotes exon-7 inclusion activation and negative combined Z denotes inhibition. Source: {SIGN_SOURCE}.
- [candidate] Quantiles use NumPy method `{QUANTILE_METHOD}` at every base, window, and coarse display block. There is no smoothing or interpolation.
- [candidate] Input SHA-256: dCas13d CSV `{sha256_file(guide_csv)}`; ASO CSV `{sha256_file(aso_csv)}`.
- [candidate] Input rows: {len(pd.read_csv(guide_csv))} total dCas13d rows, {len(data.guides)} retained SMN2 guides intersecting 260-434, and {len(data.asos)} ASO study records.
- [candidate] Software: Python {platform.python_version()}, NumPy {versions['numpy']}, pandas {versions['pandas']}, matplotlib {versions['matplotlib']}.
- [candidate] File-access audit: the generator reads only the two named CSV inputs. It does not read sealed terminal material or reuse earlier-arm code.
- [candidate] Block-edge assumption: the specification gives integer-inclusive blocks but also requires unique assignment of half-integer hull midpoints. Blocks are therefore interpreted with base-centered, left-closed/right-open edges; a midpoint exactly at 334.5 enters E7-2 rather than E7-1, and likewise at every adjacent boundary. The final right edge is closed.

## Methods

The raw guide scatter uses interval midpoints only for point placement. At each integer coordinate, the central track and quartile ribbon use every guide whose full inclusive interval covers that base. ASO window summaries instead use only guides that fully contain the inclusive ASO hull. Targets crossing the artificial exon-7 G are drawn as their two native segments, while containment uses the conservative hull. ASO directions come only from the fixed categorical input field; magnitude prose is carried as annotation and is not converted to numbers.

Each ASO study record is assigned once to one of eleven fixed coarse display blocks by its hull midpoint. Study lanes and summary rows remain separate. Duplicate coordinates in different papers remain separate records. No effect direction is averaged or voted across papers.

## Inventory

{_inventory_tables(data)}

### Exact-window coverage [candidate]

- {n_guide_covered_windows} / {n_exact_windows} distinct exact windows have at least one fully containing retained guide.

### Block-by-study descriptive sign counts [candidate]

The denominator is the number of `increase` or `decrease` records in that named study and block. `same_direction_fraction` is `same_n / eligible_direction_n`; zero-denominator cells are explicit `NA`.
Two distinct reference medians appear in this package and are not interchangeable: the *window-local containing-guide median* (`local_z_median` in the window-alignment table — the median combined Z over the guides that fully contain an individual ASO hull) governs each window's own `direction_relation`, whereas the *coarse-block positional median* (`block_z_positional_median` — the median of the per-base positional Z medians across a display block) governs the same/opposite counts in the block-by-study table below.
When the block median is exactly zero, eligible ASOs are in the denominator but are neither same nor opposite; the printed fraction is therefore zero and both sign counts are zero.

{_markdown_table(study_rel)}

## Figure captions

- [candidate] `arm_a_shared_field_map`: dCas13d values are numeric source measurements, while ASO colors are source-reported qualitative classes. The studies differ in assay and chemistry. Visual proximity or same-direction color is hypothesis-generating only.
- [candidate] `arm_a_coarse_block_map`: eleven fixed coarse display blocks summarize positional Z medians and six ASO direction classes in separate study rows. Blocks are display bins, not independent experiments.

## Limitations

- [paper_supported] ASO outcomes are categorical transcriptions. `figure_reported` is unresolved, not neutral and not assumed randomly missing; assay-discordant and mixed/uncertain records remain separate.
- [candidate] Study identity is confounded with sampled region and assay context. Empty study-region cells remain visible, and no cross-study effect estimate is made.
- [candidate] Local and block sign relations are coordinate juxtapositions, not mechanistic or causal evidence.

{CLAIM_BOUNDARY_TEXT}
"""
    output_path.write_text(report, encoding="utf-8")


def lint_scientific_outputs(output_dir: Path, figure_text: Iterable[str] = ()) -> list[str]:
    """Return forbidden-language violations outside the one required exception box."""
    violations = []
    texts: list[tuple[str, str]] = []
    for name in OUTPUT_NAMES:
        path = output_dir / name
        if path.suffix in {".md", ".tsv", ".svg"}:
            text = path.read_text(encoding="utf-8")
            if path.suffix == ".md":
                if CLAIM_BOUNDARY_TEXT not in text:
                    violations.append(f"{name}: fixed claim-boundary text differs")
                before, marker, after = text.partition(CLAIM_BOUNDARY_HEADING)
                if not marker:
                    violations.append(f"{name}: missing fixed claim-boundary section")
                text = before
            texts.append((name, text))
    texts.append(("figure text", "\n".join(figure_text)))
    for label, text in texts:
        lowered = text.casefold()
        for token in FORBIDDEN_SCIENTIFIC_STRINGS:
            if token in lowered:
                violations.append(f"{label}: {token}")
    return violations


def build_arm_a_atlas(guide_csv: Path, aso_csv: Path, output_dir: Path) -> AtlasData:
    """Generate the complete seven-file Arm A descriptive output package."""
    guide_csv, aso_csv, output_dir = Path(guide_csv), Path(aso_csv), Path(output_dir)
    data = prepare_atlas(guide_csv, aso_csv)
    output_dir.mkdir(parents=True, exist_ok=True)
    stale = [path for path in output_dir.iterdir() if path.name not in OUTPUT_NAMES]
    if stale:
        raise ValueError(f"output directory contains non-package files: {[p.name for p in stale]}")
    figure_text = []
    figure_text.extend(draw_shared_field_map(data, output_dir / OUTPUT_NAMES[0], output_dir / OUTPUT_NAMES[1]))
    figure_text.extend(draw_block_map(data, output_dir / OUTPUT_NAMES[2], output_dir / OUTPUT_NAMES[3]))
    data.windows.to_csv(output_dir / OUTPUT_NAMES[4], sep="\t", index=False, float_format="%.6g")
    data.blocks.to_csv(output_dir / OUTPUT_NAMES[5], sep="\t", index=False, float_format="%.6g")
    write_report(data, guide_csv, aso_csv, output_dir / OUTPUT_NAMES[6])
    violations = lint_scientific_outputs(output_dir, figure_text)
    if violations:
        raise ValueError("forbidden scientific language: " + "; ".join(violations))
    return data
