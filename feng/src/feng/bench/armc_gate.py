"""Arm C bind->act gate: leave-one-SBO-out (LOSO) cross-validation.

Tests whether an actionability model trained on BINDING features alone (never
the label, never S4's own energy values -- see ``armc_build``'s leakage-
discipline note) can predict, for bound exons held out by SBO, which ones the
paper actually reports as spliced (``actionable``). Compares the model against
two single-feature baselines under the same LOSO protocol:

* ``baseline_edit_distance`` -- rank by ``-edit_distance`` (closer sequence
  match = more likely to act). Orientation fixed a priori, not tuned to the
  outcome.
* ``baseline_duplex_dg`` -- rank by ``-duplex_dg`` (more favourable predicted
  duplex energy = more likely to act). Same fixed-orientation rule.

Model: no ML library is available in feng-env (``sklearn`` is NOT installed,
confirmed by import failure -- a deviation from the task brief, which assumed
it might be present). Rather than silently degrading to a single-feature
ranker, this module implements a small L2-regularized logistic regression by
batch gradient descent (numpy only), deterministic (zero-initialized weights,
no RNG), over the multivariate binding-feature vector.

ROC-AUC and average precision are also implemented from scratch (numpy +
``scipy.stats.rankdata`` for tie handling), matching the standard rank-sum
AUC definition and the sklearn-style step-function AP definition, since
``sklearn.metrics`` is unavailable for the same reason.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata

# SITE / DUPLEX-level binding features, each computed from the SBO's best binding
# window in a target by an identical code path for both labels (see armc_build's
# redesign/leakage notes). Whole-region ``exon_length`` and ``gc_fraction`` were
# dropped: region-positives and exon-negatives differ systematically in
# length/composition (a substrate confound, not binding). ``baseline_psi`` was
# dropped too (unavailable genome-wide, action-adjacent).
FEATURE_COLUMNS = [
    "edit_distance",
    "exact_match_length",
    "duplex_dg",
    "binding_window_gc",
]

# Margin (absolute AUC / AP points) the model must clear over the best
# baseline, pooled across held-out folds, for a GO verdict. Fixed in advance.
GO_MARGIN_AUC = 0.05
GO_MARGIN_AP = 0.05


# --- metrics (numpy/scipy only -- no sklearn in feng-env) -------------------


def roc_auc(y_true: np.ndarray, scores: np.ndarray) -> float | None:
    """Rank-sum (Mann-Whitney U) ROC-AUC. ``None`` if only one class present."""
    y_true = np.asarray(y_true)
    scores = np.asarray(scores, dtype=float)
    n_pos = int((y_true == 1).sum())
    n_neg = int((y_true == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return None
    ranks = rankdata(scores, method="average")
    sum_ranks_pos = ranks[y_true == 1].sum()
    u = sum_ranks_pos - n_pos * (n_pos + 1) / 2.0
    return float(u / (n_pos * n_neg))


def average_precision(y_true: np.ndarray, scores: np.ndarray) -> float | None:
    """Step-function average precision (sklearn convention). ``None`` if no
    positives present.

    Tie handling matters here: the baselines rank by integer ``edit_distance``
    (and by ``duplex_dg``, which is proxy-quantized for many rows), so tied
    scores are common. Precision/recall are therefore evaluated ONLY at the end
    of each equal-score block -- i.e. every sample sharing a score is treated as
    predicted-positive together -- exactly as ``sklearn.metrics.
    average_precision_score`` does. Stepping one sample at a time through a tie
    block (the naive form) makes AP depend on the arbitrary within-tie row order
    and can inflate it substantially (e.g. y=[1,0,1,0], all-tied-in-pairs scores
    -> 0.833 naive vs the correct 0.5), which would silently corrupt the
    GO/NO-GO margin between the model and the tie-heavy baselines.
    """
    y_true = np.asarray(y_true)
    scores = np.asarray(scores, dtype=float)
    n_pos = int((y_true == 1).sum())
    if n_pos == 0:
        return None
    order = np.argsort(-scores, kind="mergesort")
    y_sorted = (y_true[order] == 1).astype(float)
    s_sorted = scores[order]
    tp_cum = np.cumsum(y_sorted)
    n_seen = np.arange(1, len(y_sorted) + 1)
    # Keep only the last index of each equal-score run (a distinct threshold).
    group_end = np.ones(len(s_sorted), dtype=bool)
    group_end[:-1] = s_sorted[:-1] != s_sorted[1:]
    precision = tp_cum[group_end] / n_seen[group_end]
    recall = tp_cum[group_end] / n_pos
    recall_prev = np.concatenate(([0.0], recall[:-1]))
    return float(np.sum((recall - recall_prev) * precision))


# --- model: small L2 logistic regression (numpy only) ------------------------


@dataclass
class LogisticRegressionGD:
    """Deterministic (zero-init, no RNG) batch-gradient-descent logistic
    regression with L2 regularization and z-score standardization."""

    l2: float = 1.0
    lr: float = 0.5
    n_iter: int = 2000
    mean_: np.ndarray = field(default=None, repr=False)
    std_: np.ndarray = field(default=None, repr=False)
    weights_: np.ndarray = field(default=None, repr=False)
    bias_: float = 0.0

    def fit(self, X: np.ndarray, y: np.ndarray) -> "LogisticRegressionGD":
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)
        self.mean_ = X.mean(axis=0)
        self.std_ = X.std(axis=0)
        self.std_[self.std_ == 0] = 1.0
        Xs = (X - self.mean_) / self.std_
        n, d = Xs.shape
        w = np.zeros(d)
        b = 0.0
        for _ in range(self.n_iter):
            z = Xs @ w + b
            p = 1.0 / (1.0 + np.exp(-z))
            grad_w = Xs.T @ (p - y) / n + (self.l2 / n) * w
            grad_b = float(np.mean(p - y))
            w -= self.lr * grad_w
            b -= self.lr * grad_b
        self.weights_ = w
        self.bias_ = b
        return self

    def predict_proba(self, X: np.ndarray) -> np.ndarray:
        X = np.asarray(X, dtype=float)
        Xs = (X - self.mean_) / self.std_
        z = Xs @ self.weights_ + self.bias_
        return 1.0 / (1.0 + np.exp(-z))


# --- LOSO CV ------------------------------------------------------------


@dataclass
class FoldResult:
    group: str
    n_pos: int
    n_neg: int
    trained: bool
    evaluable: bool
    model_auc: float | None = None
    model_ap: float | None = None
    baseline_edit_auc: float | None = None
    baseline_edit_ap: float | None = None
    baseline_dg_auc: float | None = None
    baseline_dg_ap: float | None = None


@dataclass
class GateResult:
    folds: list[FoldResult]
    pooled_y: np.ndarray
    pooled_model_score: np.ndarray
    pooled_baseline_edit_score: np.ndarray
    pooled_baseline_dg_score: np.ndarray
    pooled_model_auc: float | None
    pooled_model_ap: float | None
    pooled_baseline_edit_auc: float | None
    pooled_baseline_edit_ap: float | None
    pooled_baseline_dg_auc: float | None
    pooled_baseline_dg_ap: float | None
    n_groups_total: int
    n_groups_skipped_train: int
    n_groups_skipped_eval: int
    n_rows_total: int
    n_positive_total: int
    verdict: str


def run_loso_gate(
    features: pd.DataFrame,
    *,
    feature_columns: list[str] = FEATURE_COLUMNS,
    l2: float = 1.0,
    lr: float = 0.5,
    n_iter: int = 2000,
) -> GateResult:
    required = {"sbo_group", "actionable", "edit_distance", "duplex_dg", *feature_columns}
    missing = required - set(features.columns)
    if missing:
        raise ValueError(f"features table missing columns: {sorted(missing)}")

    df = features.dropna(subset=feature_columns + ["actionable"]).reset_index(drop=True)
    groups = sorted(df["sbo_group"].unique())

    folds: list[FoldResult] = []
    pooled_y: list[int] = []
    pooled_model: list[float] = []
    pooled_edit: list[float] = []
    pooled_dg: list[float] = []
    n_skipped_train = 0
    n_skipped_eval = 0

    for g in groups:
        test_mask = df["sbo_group"] == g
        train_df = df.loc[~test_mask]
        test_df = df.loc[test_mask]
        n_pos = int((test_df["actionable"] == 1).sum())
        n_neg = int((test_df["actionable"] == 0).sum())
        evaluable = n_pos > 0 and n_neg > 0

        if train_df["actionable"].nunique() < 2:
            n_skipped_train += 1
            folds.append(FoldResult(group=g, n_pos=n_pos, n_neg=n_neg, trained=False, evaluable=False))
            continue
        if not evaluable:
            n_skipped_eval += 1

        model = LogisticRegressionGD(l2=l2, lr=lr, n_iter=n_iter).fit(
            train_df[feature_columns].to_numpy(), train_df["actionable"].to_numpy()
        )
        model_scores = model.predict_proba(test_df[feature_columns].to_numpy())
        edit_scores = -test_df["edit_distance"].to_numpy(dtype=float)
        dg_scores = -test_df["duplex_dg"].to_numpy(dtype=float)
        y_test = test_df["actionable"].to_numpy()

        pooled_y.extend(y_test.tolist())
        pooled_model.extend(model_scores.tolist())
        pooled_edit.extend(edit_scores.tolist())
        pooled_dg.extend(dg_scores.tolist())

        fold = FoldResult(
            group=g, n_pos=n_pos, n_neg=n_neg, trained=True, evaluable=evaluable,
            model_auc=roc_auc(y_test, model_scores) if evaluable else None,
            model_ap=average_precision(y_test, model_scores) if evaluable else None,
            baseline_edit_auc=roc_auc(y_test, edit_scores) if evaluable else None,
            baseline_edit_ap=average_precision(y_test, edit_scores) if evaluable else None,
            baseline_dg_auc=roc_auc(y_test, dg_scores) if evaluable else None,
            baseline_dg_ap=average_precision(y_test, dg_scores) if evaluable else None,
        )
        folds.append(fold)

    pooled_y_arr = np.array(pooled_y)
    pooled_model_arr = np.array(pooled_model)
    pooled_edit_arr = np.array(pooled_edit)
    pooled_dg_arr = np.array(pooled_dg)

    pooled_model_auc = roc_auc(pooled_y_arr, pooled_model_arr)
    pooled_model_ap = average_precision(pooled_y_arr, pooled_model_arr)
    pooled_edit_auc = roc_auc(pooled_y_arr, pooled_edit_arr)
    pooled_edit_ap = average_precision(pooled_y_arr, pooled_edit_arr)
    pooled_dg_auc = roc_auc(pooled_y_arr, pooled_dg_arr)
    pooled_dg_ap = average_precision(pooled_y_arr, pooled_dg_arr)

    verdict = "NO-GO (insufficient held-out coverage)"
    if pooled_model_auc is not None and pooled_model_ap is not None:
        best_baseline_auc = max(v for v in (pooled_edit_auc, pooled_dg_auc) if v is not None)
        best_baseline_ap = max(v for v in (pooled_edit_ap, pooled_dg_ap) if v is not None)
        clears_auc = pooled_model_auc - best_baseline_auc >= GO_MARGIN_AUC
        clears_ap = pooled_model_ap - best_baseline_ap >= GO_MARGIN_AP
        verdict = "GO" if (clears_auc and clears_ap) else "NO-GO"

    return GateResult(
        folds=folds,
        pooled_y=pooled_y_arr,
        pooled_model_score=pooled_model_arr,
        pooled_baseline_edit_score=pooled_edit_arr,
        pooled_baseline_dg_score=pooled_dg_arr,
        pooled_model_auc=pooled_model_auc,
        pooled_model_ap=pooled_model_ap,
        pooled_baseline_edit_auc=pooled_edit_auc,
        pooled_baseline_edit_ap=pooled_edit_ap,
        pooled_baseline_dg_auc=pooled_dg_auc,
        pooled_baseline_dg_ap=pooled_dg_ap,
        n_groups_total=len(groups),
        n_groups_skipped_train=n_skipped_train,
        n_groups_skipped_eval=n_skipped_eval,
        n_rows_total=len(df),
        n_positive_total=int((df["actionable"] == 1).sum()),
        verdict=verdict,
    )


# --- report -------------------------------------------------------------


def _fmt(v: float | None) -> str:
    return "n/a" if v is None else f"{v:.3f}"


def render_gate_result(result: GateResult) -> str:
    pos_rate = result.n_positive_total / result.n_rows_total if result.n_rows_total else float("nan")
    lines = [
        "# Arm C gate result -- bind->act, leave-one-SBO-out CV",
        "",
        f"**Verdict: {result.verdict}**",
        "",
        "Result is `candidate` / `paper_supported` only -- no sealed ground-truth "
        "directory or Phase-9 data was read to produce it. Model = site-level binding "
        "features only (edit_distance, exact_match_length, duplex_dg, binding_window_gc); "
        "no whole-region length/GC, no baseline_psi, no SpliceAI/Pangolin/SPLAM. Positive "
        "class = direct-bound-and-acted (S4 numeric edit_dist<=5, lifted hg19->GRCh38); "
        "negatives = bound-but-not-acted GENCODE v50 pool exons.",
        "",
        "## Coverage",
        "",
        f"- Rows (bound SBO,exon pairs) scored: {result.n_rows_total}",
        f"- Positive rate: {result.n_positive_total}/{result.n_rows_total} = {pos_rate:.3f}",
        f"- SBO groups total: {result.n_groups_total}",
        f"- SBO groups skipped (training fold had 0 pos or 0 neg -- could not fit): "
        f"{result.n_groups_skipped_train}",
        f"- SBO groups with a degenerate held-out set (0 pos or 0 neg in that fold's test "
        f"rows -- excluded from per-fold AUC/AP, still included in pooled metrics): "
        f"{result.n_groups_skipped_eval}",
        "",
        "## Pooled held-out metrics (GO margin: "
        f"AUC >= best baseline + {GO_MARGIN_AUC}, AP >= best baseline + {GO_MARGIN_AP})",
        "",
        "| | ROC-AUC | Average Precision |",
        "|---|---|---|",
        f"| Model (LOSO logistic regression, binding features) | {_fmt(result.pooled_model_auc)} | "
        f"{_fmt(result.pooled_model_ap)} |",
        f"| Baseline: -edit_distance | {_fmt(result.pooled_baseline_edit_auc)} | "
        f"{_fmt(result.pooled_baseline_edit_ap)} |",
        f"| Baseline: -duplex_dg | {_fmt(result.pooled_baseline_dg_auc)} | "
        f"{_fmt(result.pooled_baseline_dg_ap)} |",
        "",
        "## Per-fold (SBO held out)",
        "",
        "| dg_id | n_pos | n_neg | trained | evaluable | model AUC | model AP | "
        "edit AUC | dG AUC |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for f in result.folds:
        lines.append(
            f"| {f.group} | {f.n_pos} | {f.n_neg} | {f.trained} | {f.evaluable} | "
            f"{_fmt(f.model_auc)} | {_fmt(f.model_ap)} | {_fmt(f.baseline_edit_auc)} | "
            f"{_fmt(f.baseline_dg_auc)} |"
        )
    lines.append("")
    return "\n".join(lines)


# --- CLI ------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Arm C LOSO bind->act gate.")
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)

    features = pd.read_parquet(args.features)
    result = run_loso_gate(features)
    report = render_gate_result(result)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(report, encoding="utf-8")
    print(report)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
