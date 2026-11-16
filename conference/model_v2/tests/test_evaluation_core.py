"""Adversarial tests for F3 evaluation_core: leakage, missingness, tiny-N, determinism,
forbidden sources, PU mode. Synthetic rows only - no lab data, no files outside tmp_path."""
from __future__ import annotations

import json
import os
import random
import sys
from dataclasses import replace
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from conference.model_v2 import evaluation_core as ec  # noqa: E402

BASES = "ACGT"
NEG_EVIDENCE = {"assayed_and_powered": True, "detection_limit": 0.05, "basis": "dose_response"}


@pytest.fixture(autouse=True)
def synthetic_chemistry_vocabulary(monkeypatch):
    # Synthetic independent families are explicitly registered by the test caller.
    names = ["chema", "chemy", "chemz"] + [f"chem{i}" for i in range(200)]
    names += [f"c{i}" for i in range(40)] + [f"cw{i}" for i in range(20)]
    names += [f"k{i}_{j}" for i in range(12) for j in range(6)]
    monkeypatch.setattr(ec, "CHEMISTRY_ALIASES", {**ec.CHEMISTRY_ALIASES,
                                               **{name: name for name in names}})


def rand_seq(rng, n=18):
    return "".join(rng.choice(BASES) for _ in range(n))


def row(oid, state="measured_positive", *, seq="ACGTACGTACGTACGTAC", chem="chemA", study="S1",
        unit=None, chrom="chr5", start=1000, end=1017, strand="+", **extra):
    r = {"observation_id": oid, "label_state": state, "assay_id": "assay1", "study_id": study,
         "independent_unit_id": unit or oid, "provenance": {"source": "synthetic:test"},
         "aso_sequence": seq, "chemistry_family": chem, "chrom": chrom, "start": start,
         "end": end, "strand": strand}
    if state == "verified_negative":
        r["negative_evidence"] = dict(NEG_EVIDENCE)
    r.update(extra)
    return r


def independent_rows(n, *, seed=7, states=None, n_studies=None, n_chem=None):
    """n rows that share NO group on any axis unless n_studies / n_chem fold them."""
    rng = random.Random(seed)
    out = []
    for i in range(n):
        out.append(row(
            f"o{i:03d}", (states[i] if states else "measured_positive"), seq=rand_seq(rng),
            chem=f"chem{i % n_chem if n_chem else i}", study=f"S{i % n_studies if n_studies else i}",
            start=1000 + 500 * i, end=1017 + 500 * i))
    return out


def labelled(n_pos, n_neg, n_unl=0, n_unassayed=0, n_censored=0, **kw):
    states = (["measured_positive"] * n_pos + ["verified_negative"] * n_neg
              + ["unlabeled"] * n_unl + ["unassayed"] * n_unassayed + ["censored"] * n_censored)
    return independent_rows(len(states), states=states, **kw)


def assert_no_straddle(plan):
    for fold in plan.folds:
        assert not set(fold.test_idx) & set(fold.train_idx)
        assert not set(fold.purged_idx) & (set(fold.test_idx) | set(fold.train_idx))
        assert len(fold.test_idx) + len(fold.train_idx) + len(fold.purged_idx) == len(plan.observation_ids)
        for axis in plan.axes:
            g = plan.group_ids[axis]
            assert not {g[i] for i in fold.test_idx} & {g[i] for i in fold.train_idx}, (fold.fold_id, axis)


# ------------------------------------------------------------------ label states / validator
def test_label_states_are_exhaustive_and_closed():
    assert set(ec.LABEL_STATES) == {"measured_positive", "verified_negative", "unlabeled",
                                    "unassayed", "censored"}
    with pytest.raises(ec.SchemaError, match="label_state"):
        ec.validate_table([row("a", "negative")])
    with pytest.raises(ec.SchemaError):
        ec.validate_table([row("a", "non_responder")])


@pytest.mark.parametrize("field", ["observation_id", "label_state", "assay_id", "study_id",
                                   "independent_unit_id", "provenance"])
def test_every_required_field_is_enforced(field):
    r = row("a")
    del r[field]
    with pytest.raises(ec.SchemaError):
        ec.validate_table([r])
    r = row("a")
    r[field] = ""
    with pytest.raises(ec.SchemaError):
        ec.validate_table([r])


def test_duplicate_observation_ids_and_bad_provenance_rejected():
    with pytest.raises(ec.SchemaError, match="duplicate"):
        ec.validate_table([row("a"), row("a")])
    with pytest.raises(ec.SchemaError, match="source"):
        ec.validate_table([row("a", provenance={"note": "x"})])
    with pytest.raises(ec.SchemaError, match="evidence_class"):
        ec.validate_table([row("a", evidence_class="measured")])


def test_valid_table_passes_and_all_states_accepted():
    table = ec.validate_table(labelled(1, 1, 1, 1, 1))
    assert [r["label_state"] for r in table.rows] == list(ec.LABEL_STATES)


@pytest.mark.parametrize("evidence", [
    None, {}, {"assayed_and_powered": True}, {"detection_limit": 0.1},
    {"assayed_and_powered": False, "detection_limit": 0.1},
    {"assayed_and_powered": "yes", "detection_limit": 0.1},
    {"assayed_and_powered": 1, "detection_limit": 0.1},
    {"assayed_and_powered": True, "detection_limit": None},
    {"assayed_and_powered": True, "detection_limit": float("nan")},
    {"assayed_and_powered": True, "detection_limit": True},
    {"assayed_and_powered": True, "detection_limit": 0.1, "basis": "absent_from_table"},
    {"assayed_and_powered": True, "detection_limit": 0.1, "basis": "Not_Listed"},
])
def test_verified_negative_needs_explicit_powered_evidence(evidence):
    r = row("n1", "verified_negative")
    r["negative_evidence"] = evidence
    with pytest.raises(ec.NegativeEvidenceError):
        ec.validate_table([r])


def test_negative_from_responder_only_table_is_rejected_even_with_claimed_evidence():
    r = row("n1", "verified_negative", source_table_kind="responder_only")
    with pytest.raises(ec.NegativeEvidenceError, match="responder-only"):
        ec.validate_table([r])
    r = row("n2", "verified_negative", study="holgersen_s4")
    with pytest.raises(ec.NegativeEvidenceError, match="responder-only"):
        ec.validate_table([r], responder_only_studies=["holgersen_s4"])
    # positives and unlabeled rows from the same responder-only study are fine
    ec.validate_table([row("p", study="holgersen_s4"), row("u", "unlabeled", study="holgersen_s4")],
                      responder_only_studies=["holgersen_s4"])


@pytest.mark.parametrize("state", ["unassayed", "unlabeled", "censored"])
@pytest.mark.parametrize("field", ["binary_label", "label", "y"])
def test_non_negative_state_encoded_as_zero_raises(state, field):
    with pytest.raises(ec.LabelStateError, match="never a negative"):
        ec.validate_table([row("a", state, **{field: 0})])


def test_binary_label_must_agree_with_state():
    ec.validate_table([row("p", y=1), row("n", "verified_negative", y=0)])
    with pytest.raises(ec.LabelStateError):
        ec.validate_table([row("u", "unlabeled", y=1)])
    with pytest.raises(ec.LabelStateError):
        ec.validate_table([row("p", y=0)])
    with pytest.raises(ec.SchemaError):
        ec.validate_table([row("p", y=2)])


# ------------------------------------------------------------------ grouping axes
def test_duplicate_aso_under_different_ids_shares_a_sequence_group():
    rows = [row("a", seq="ACGUACGUACGUACGUAC"), row("b", seq="acgt acgt acgt acgt ac"),
            row("c", seq="TTTTTTTTTTGGGGGGGG")]
    g = ec.axis_groups(rows, axes=("sequence",))["sequence"]
    assert g[0] == g[1] != g[2]


def test_near_identical_and_truncated_sequences_cluster():
    base = "ACGTACGTTTGGCCAAGT"
    one_sub = base[:5] + "T" + base[6:]
    truncated = base[2:16]              # 14-mer inside the 18-mer (edit distance 4)
    far = "GGGGGGCCCCCCAAAAAA"
    rows = [row("a", seq=base), row("b", seq=one_sub), row("c", seq=truncated), row("d", seq=far)]
    g = ec.axis_groups(rows, axes=("sequence",))["sequence"]
    assert g[0] == g[1] == g[2] != g[3]
    strict = ec.GroupingConfig(sequence_threshold=0, containment_is_identical=False)
    g = ec.axis_groups(rows, strict, axes=("sequence",))["sequence"]
    assert len(set(g)) == 4
    hamming = ec.GroupingConfig(sequence_metric="hamming", sequence_threshold=1,
                                containment_is_identical=False)
    g = ec.axis_groups(rows, hamming, axes=("sequence",))["sequence"]
    assert g[0] == g[1] and g[2] not in (g[0], g[3])


def test_sequence_clusters_are_single_linkage_transitive():
    a = "AAAAAAAAAAAAAAAAAA"
    b = "AAAAAAAAAAAAAAAACC"   # 2 from a
    c = "AAAAAAAAAAAAAACCCC"   # 2 from b, 4 from a
    g = ec.axis_groups([row("a", seq=a), row("b", seq=b), row("c", seq=c)], axes=("sequence",))
    assert len(set(g["sequence"])) == 1


def test_window_clustering_overlap_padding_strand_and_chrom():
    rows = [row("a", start=100, end=120), row("b", start=120, end=140),   # overlap at 120
            row("c", start=141, end=160),                                  # adjacent, no overlap
            row("d", start=100, end=120, strand="-"),                      # other strand
            row("e", start=100, end=120, chrom="chr1")]                    # other chromosome
    g = ec.axis_groups(rows, axes=("window",))["window"]
    assert g[0] == g[1] and len({g[0], g[2], g[3], g[4]}) == 4
    padded = ec.axis_groups(rows, ec.GroupingConfig(window_padding=1), axes=("window",))["window"]
    assert padded[0] == padded[1] == padded[2] != padded[3]
    unstranded = ec.axis_groups(rows, ec.GroupingConfig(strand_aware=False), axes=("window",))["window"]
    assert unstranded[0] == unstranded[3] != unstranded[4]


def test_window_chain_merges_transitively_and_contained_windows_do_not_reset_span():
    rows = [row("a", start=100, end=300), row("b", start=110, end=120), row("c", start=290, end=310)]
    g = ec.axis_groups(rows, axes=("window",))["window"]
    assert len(set(g)) == 1 and g[0] == "win:chr5:+:100-310"


@pytest.mark.parametrize("bad", [{"start": None}, {"start": 50, "end": 40}, {"strand": "."},
                                 {"chrom": ""}, {"start": 0}, {"start": 1.5}])
def test_window_axis_rejects_missing_or_malformed_coordinates(bad):
    with pytest.raises(ec.SchemaError, match="window axis"):
        ec.axis_groups([row("a", **bad)], axes=("window",))


def test_missing_chemistry_or_bad_sequence_is_not_silently_grouped():
    with pytest.raises(ec.SchemaError):
        ec.axis_groups([row("a", chem=None)], axes=("chemistry",))
    with pytest.raises(ec.SchemaError):
        ec.axis_groups([row("a", seq="ACGX")], axes=("sequence",))


# ------------------------------------------------------------------ leakage traps
def test_trap_duplicate_aso_under_different_ids_never_straddles():
    rows = independent_rows(6)
    rows.append(row("dup", seq=rows[0]["aso_sequence"], chem="chemZ", study="SZ",
                    start=900000, end=900017))
    plan = ec.make_splits(rows, "leave_study_out")
    assert_no_straddle(plan)
    for fold in plan.folds:
        ids = {plan.observation_ids[i] for i in fold.train_idx}
        held = {plan.observation_ids[i] for i in fold.test_idx} | \
               {plan.observation_ids[i] for i in fold.purged_idx}
        assert not ({"o000", "dup"} & ids and {"o000", "dup"} & held)
    # the duplicate is purged (not trained on) when its twin's study is held out
    fold = next(f for f in plan.folds if f.test_group_ids == ("study:S0",))
    assert [plan.observation_ids[i] for i in fold.purged_idx] == ["dup"]


def test_trap_overlapping_windows_across_studies_never_split():
    rows = independent_rows(6)
    rows.append(row("ovl", seq="GGGGGGGGGGCCCCCCCC", chem="chemY", study="SY",
                    start=rows[2]["start"] + 5, end=rows[2]["end"] + 5))
    for mode in ("leave_one_aso_out", "leave_study_out", "region_blocked"):
        plan = ec.make_splits(rows, mode)
        assert_no_straddle(plan)
        for fold in plan.folds:
            train = {plan.observation_ids[i] for i in fold.train_idx}
            assert not ({"o002", "ovl"} & train and {"o002", "ovl"} - train), mode
    plan = ec.make_splits(rows, "region_blocked")
    fold = next(f for f in plan.folds if "o002" in {plan.observation_ids[i] for i in f.test_idx})
    assert {plan.observation_ids[i] for i in fold.test_idx} == {"o002", "ovl"}


def test_trap_chemistry_family_straddle():
    rows = independent_rows(8, n_chem=4)        # 4 families x 2 otherwise-independent ASOs
    plan = ec.make_splits(rows, "leave_one_aso_out")
    assert_no_straddle(plan)
    assert plan.n_independent_groups == 4
    for fold in plan.folds:
        assert len(fold.test_idx) == 1 and len(fold.purged_idx) == 1   # family twin purged
        chem = plan.group_ids["chemistry"]
        assert chem[fold.purged_idx[0]] == chem[fold.test_idx[0]]


def test_relaxing_an_axis_is_explicit_and_recorded():
    rows = independent_rows(6, n_chem=1)        # one chemistry family for everything
    with pytest.raises(ec.SplitRefused, match="independent N too small"):
        ec.make_splits(rows, "leave_one_aso_out")
    plan = ec.make_splits(rows, "leave_one_aso_out", config=ec.GroupingConfig(relax=("chemistry",)))
    assert plan.relaxed_axes == ("chemistry",) and "chemistry" not in plan.axes
    assert len(plan.folds) == 6
    assert_no_straddle(plan)
    with pytest.raises(ec.SchemaError, match="primary axis"):
        ec.make_splits(rows, "leave_one_aso_out", config=ec.GroupingConfig(relax=("sequence",)))


def test_randomised_leakage_fuzz_all_modes():
    rng = random.Random(123)
    for trial in range(12):
        rows = []
        for i in range(40):
            start = rng.randrange(1, 4000)
            rows.append(row(f"r{i}", seq=rand_seq(rng, 16), chem=f"c{rng.randrange(12)}",
                            study=f"s{rng.randrange(12)}", start=start, end=start + rng.randrange(10, 40),
                            strand=rng.choice("+-"), chrom=rng.choice(["chr1", "chr5", "chrX"])))
        for mode in ec.MODES:
            for n_folds in (None, 2):
                try:
                    plan = ec.make_splits(rows, mode, seed=trial, n_folds=n_folds,
                                          min_independent_groups=2)
                except ec.SplitRefused:
                    continue
                assert_no_straddle(plan)
                ec.verify_no_leakage(plan)
                tested = [i for f in plan.folds for i in f.test_idx]
                assert sorted(tested) == list(range(40))      # every row tested exactly once


def test_verify_no_leakage_detects_a_tampered_plan():
    plan = ec.make_splits(independent_rows(5), "leave_one_aso_out")
    f = plan.folds[0]
    bad = ec.Fold(f.fold_id, f.test_idx, f.train_idx + f.test_idx, f.purged_idx,
                  f.test_group_ids, f.held_out_groups)
    tampered = ec.SplitPlan(plan.mode, plan.primary_axis, plan.axes, plan.relaxed_axes, plan.seed,
                            plan.n_independent_groups, plan.observation_ids, plan.group_ids, (bad,))
    with pytest.raises(ec.LeakageError):
        ec.verify_no_leakage(tampered)


# ------------------------------------------------------------------ modes, tiny N, nesting
def test_leave_study_out_holds_out_whole_studies():
    rows = independent_rows(12, n_studies=4)
    plan = ec.make_splits(rows, "leave_study_out")
    assert [f.test_group_ids for f in plan.folds] == [(f"study:S{k}",) for k in range(4)]
    assert all(len(f.test_idx) == 3 and len(f.train_idx) == 9 for f in plan.folds)


@pytest.mark.parametrize("mode", sorted(ec.MODES))
def test_tiny_n_is_refused_with_a_stated_reason(mode):
    with pytest.raises(ec.SplitRefused, match=r"independent N too small - 2 independent group"):
        ec.make_splits(independent_rows(2), mode)
    # three near-identical truncations of one ASO are ONE independent unit, not three
    base = "ACGTACGTTTGGCCAAGT"
    trio = [row("f18", seq=base), row("f14", seq=base[:14], start=1000, end=1013),
            row("f10", seq=base[:10], start=1000, end=1009)]
    with pytest.raises(ec.SplitRefused, match="1 independent group"):
        ec.make_splits(trio, mode)


def test_kfold_refuses_more_folds_than_independent_groups():
    with pytest.raises(ec.SplitRefused, match="n_folds=5"):
        ec.make_splits(independent_rows(4), "leave_one_aso_out", n_folds=5)
    with pytest.raises(ec.SchemaError):
        ec.make_splits(independent_rows(4), "leave_everything_out")
    with pytest.raises(ec.SchemaError, match="make_nested_splits"):
        ec.make_splits(independent_rows(4), "nested")


def test_nested_inner_folds_never_see_outer_test_groups():
    rows = independent_rows(15, n_studies=5)
    rows.append(row("dup", seq=rows[0]["aso_sequence"], chem="chemZ", study="S3",
                    start=700000, end=700017))     # links study S0 and S3 through a duplicate
    nested = ec.make_nested_splits(rows, outer_mode="leave_study_out",
                                   inner_mode="leave_one_aso_out", min_outer_groups=3,
                                   min_inner_groups=3)
    plan = nested.outer
    assert_no_straddle(plan)
    assert len(nested.inner) == len(plan.folds)
    for outer, inner in zip(plan.folds, nested.inner):
        barred_rows = set(outer.test_idx) | set(outer.purged_idx)
        assert inner
        for f in inner:
            seen = set(f.test_idx) | set(f.train_idx) | set(f.purged_idx)
            assert seen == set(outer.train_idx)
            assert not seen & barred_rows
            for axis in plan.axes:
                g = plan.group_ids[axis]
                assert not {g[i] for i in seen} & {g[i] for i in outer.test_idx}
                assert not {g[i] for i in f.test_idx} & {g[i] for i in f.train_idx}


def test_nested_refuses_when_independent_n_is_too_small_and_says_which_level():
    with pytest.raises(ec.SplitRefused, match="independent N too small"):
        ec.make_nested_splits(independent_rows(2), outer_mode="leave_one_aso_out")
    with pytest.raises(ec.SplitRefused, match=r"nested refused: outer fold 0 .*inner"):
        ec.make_nested_splits(independent_rows(4), outer_mode="leave_one_aso_out",
                              min_outer_groups=3, min_inner_groups=4)
    ok = ec.make_nested_splits(independent_rows(4), outer_mode="leave_one_aso_out")
    assert all(len(inner) == 3 for inner in ok.inner)


def test_splits_are_deterministic_and_order_independent_in_group_ids():
    rows = independent_rows(20, n_studies=10)
    a = ec.make_splits(rows, "leave_study_out", seed=11, n_folds=4)
    b = ec.make_splits(rows, "leave_study_out", seed=11, n_folds=4)
    c = ec.make_splits(rows, "leave_study_out", seed=12, n_folds=4)
    assert a == b
    assert [f.test_idx for f in a.folds] != [f.test_idx for f in c.folds]
    shuffled = list(reversed(rows))
    d = ec.make_splits(shuffled, "leave_study_out", seed=11, n_folds=4)
    as_ids = lambda plan: sorted(sorted(plan.observation_ids[i] for i in f.test_idx) for f in plan.folds)
    assert as_ids(a) == as_ids(d)
    n1 = ec.make_nested_splits(rows, outer_mode="leave_study_out", seed=3, outer_n_folds=5, inner_n_folds=2)
    n2 = ec.make_nested_splits(rows, outer_mode="leave_study_out", seed=3, outer_n_folds=5, inner_n_folds=2)
    assert n1 == n2


# ------------------------------------------------------------------ metric guard
POLICY = ec.MetricPolicy(min_verified_negatives=3, min_positives=3)


def scores_for(rows, values):
    return {r["observation_id"]: v for r, v in zip(rows, values)}


@pytest.mark.parametrize("name", ["AUPRC", "auroc", "roc_auc_score", "average_precision", "FDR",
                                  "accuracy", "precision", "recall", "brier", "calibration"])
def test_metrics_blocked_without_verified_negatives(name):
    rows = labelled(5, 0, n_unl=20, n_unassayed=5)
    with pytest.raises(ec.MetricBlocked) as err:
        ec.compute_metric(name, rows, [0.5] * len(rows), policy=POLICY, threshold=0.5)
    assert "0 independent verified_negative" in err.value.reason
    assert "pu_report" in err.value.reason and err.value.metric == name


def test_metric_blocked_when_either_class_is_below_minimum():
    with pytest.raises(ec.MetricBlocked, match="2 independent verified_negative"):
        ec.compute_metric("auroc", labelled(5, 2), [0.1] * 7, policy=POLICY)
    with pytest.raises(ec.MetricBlocked, match="1 independent measured_positive"):
        ec.compute_metric("auroc", labelled(1, 5), [0.1] * 6, policy=POLICY)
    with pytest.raises(ec.MetricBlocked):                      # default minimum is 10 + 10
        ec.compute_metric("auroc", labelled(9, 9), [0.1] * 18)


def test_minimum_counts_independent_units_not_rows():
    rows = labelled(5, 5)
    for r in rows[5:]:
        r["independent_unit_id"] = "one_plate"     # 5 negative rows, a single independent unit
    with pytest.raises(ec.MetricBlocked, match=r"1 independent verified_negative unit\(s\) \(5 rows\)"):
        ec.compute_metric("auprc", rows, list(range(10)), policy=POLICY)


@pytest.mark.parametrize("states", [("verified_negative", "unassayed"), ("unlabeled",),
                                    ("verified_negative", "unlabeled"), ("censored",), ()])
def test_unassayed_or_unlabeled_can_never_be_declared_negative(states):
    rows = labelled(5, 5, n_unl=5, n_unassayed=5)
    with pytest.raises(ec.LabelStateError):
        ec.compute_metric("auroc", rows, list(range(20)), policy=POLICY, negative_states=states)


def test_non_negative_states_are_excluded_not_counted_as_negatives():
    rows = labelled(4, 4, n_unl=6, n_unassayed=3, n_censored=2)
    # positives outscore verified negatives; unlabeled/unassayed/censored get the TOP scores
    values = [10, 11, 12, 13, 1, 2, 3, 4] + [100] * 11
    res = ec.compute_metric("auroc", rows, values, policy=POLICY)
    assert res.value == 1.0 and res.evidence_class == "derived" and res.mode == "supervised"
    assert (res.n_positive_rows, res.n_negative_rows) == (4, 4)
    assert res.excluded_by_state == {"unlabeled": 6, "unassayed": 3, "censored": 2}
    assert ec.compute_metric("auprc", rows, values, policy=POLICY).value == 1.0


def test_metric_values_and_tie_handling():
    rows = labelled(3, 3)
    # pos scores 3,2,2 ; neg scores 2,1,1 -> pairs: 3>all (3), each 2: >1,1 and ties 2 (2.5 each)
    values = [3, 2, 2, 2, 1, 1]
    res = ec.compute_metric("auroc", rows, values, policy=POLICY)
    assert res.value == pytest.approx(8 / 9) and res.n_tied_rows == 3
    # AP: thresholds 3 (P=1, R=1/3), 2 (tp=3, fp=1 -> P=.75, R=1): 1/3*1 + 2/3*.75
    assert ec.compute_metric("auprc", rows, values, policy=POLICY).value == pytest.approx(1 / 3 + 0.5)
    assert ec.compute_metric("auroc", rows, [1] * 6, policy=POLICY).value == 0.5   # all tied
    assert ec.compute_metric("auprc", rows, [1] * 6, policy=POLICY).value == 0.5   # = prevalence
    assert ec.compute_metric("accuracy", rows, values, policy=POLICY, threshold=2).value == pytest.approx(5 / 6)
    assert ec.compute_metric("precision", rows, values, policy=POLICY, threshold=2).value == 0.75
    assert ec.compute_metric("fdr", rows, values, policy=POLICY, threshold=2).value == 0.25
    assert ec.compute_metric("recall", rows, values, policy=POLICY, threshold=3).value == pytest.approx(1 / 3)
    curve = ec.precision_recall_curve(rows, values, policy=POLICY)
    assert curve == [(3.0, 1.0, pytest.approx(1 / 3)), (2.0, 0.75, 1.0), (1.0, 0.5, 1.0)]
    with pytest.raises(ec.SchemaError, match="threshold"):
        ec.compute_metric("fdr", rows, values, policy=POLICY)
    with pytest.raises(ec.MetricBlocked, match="undefined"):
        ec.compute_metric("fdr", rows, values, policy=POLICY, threshold=99)
    with pytest.raises(ec.SchemaError, match="unknown metric"):
        ec.compute_metric("vibes", rows, values, policy=POLICY)


def test_calibration_and_brier():
    rows = labelled(3, 3)
    perfect = [1, 1, 1, 0, 0, 0]
    assert ec.compute_metric("brier", rows, perfect, policy=POLICY).value == 0.0
    assert ec.compute_metric("calibration_error", rows, perfect, policy=POLICY).value == 0.0
    half = [0.5] * 6
    assert ec.compute_metric("brier", rows, half, policy=POLICY).value == 0.25
    assert ec.compute_metric("ece", rows, half, policy=POLICY).value == 0.0
    over = [0.9] * 6                                # one bin: |0.9 - 0.5|
    assert ec.compute_metric("ece", rows, over, policy=POLICY).value == pytest.approx(0.4)
    with pytest.raises(ec.SchemaError, match="probabilities"):
        ec.compute_metric("brier", rows, [2, 1, 1, 0, 0, 0], policy=POLICY)


def test_scores_must_be_finite_complete_and_aligned():
    rows = labelled(3, 3)
    with pytest.raises(ec.SchemaError):
        ec.compute_metric("auroc", rows, [1, 2, 3], policy=POLICY)
    with pytest.raises(ec.SchemaError):
        ec.compute_metric("auroc", rows, [1, 2, 3, 4, 5, float("nan")], policy=POLICY)
    with pytest.raises(ec.SchemaError, match="missing"):
        ec.compute_metric("auroc", rows, {"o000": 1.0}, policy=POLICY)


def test_precision_recall_curve_is_guarded_too():
    with pytest.raises(ec.MetricBlocked):
        ec.precision_recall_curve(labelled(5, 0, n_unl=5), [0.0] * 10, policy=POLICY)


# ------------------------------------------------------------------ PU mode
def test_pu_mode_never_emits_supervised_metrics():
    rows = labelled(4, 0, n_unl=16)
    report = ec.pu_report(rows, list(range(20, 0, -1)), top_k=(4, 10), class_prior=0.1, min_positives=3)
    names = {q.name for q in report.quantities}
    assert names <= set(ec.PU_QUANTITIES)
    banned = set(ec.SUPERVISED_METRICS) | {"auc", "average_precision", "f1", "roc_auc"}
    assert not {ec.metric_key(n) for n in names} & banned
    assert not any(b in n for n in names for b in ("auprc", "auroc", "fdr", "accuracy", "precision"))
    assert report.mode == "positive_unlabeled"
    for q in report.quantities:
        assert q.evidence_class in ec.EVIDENCE_CLASSES
        assert "never counted as negatives" in q.pu_assumption
    for name in ("auprc", "auroc", "fdr", "accuracy", "precision", "recall"):
        with pytest.raises(ec.MetricBlocked, match="positive-unlabeled mode never emits"):
            ec.compute_metric(name, rows, list(range(20)), mode="pu", policy=ec.MetricPolicy(1, 1))


def test_pu_quantities_values():
    rows = labelled(4, 0, n_unl=16)
    report = ec.pu_report(rows, list(range(20, 0, -1)), top_k=(2, 4, 20), min_positives=3)
    assert report.get("pu_recall_at_top_k", 2).value == 0.5
    assert report.get("pu_recall_at_top_k", 4).value == 1.0
    assert report.get("pu_recall_at_top_k", 20).value == 1.0
    assert report.get("pu_positive_vs_unlabeled_rank_probability").value == 1.0
    assert report.get("pu_positive_mean_percentile_rank").value == pytest.approx((19.5 + 18.5 + 17.5 + 16.5) / 80)
    assert report.blocked == {"pu_prior_adjusted_rank_probability": "no class_prior stated"}
    assert (report.n_positive_rows, report.n_universe) == (4, 20)


def test_pu_ties_get_expected_credit_not_optimistic_credit():
    rows = labelled(4, 0, n_unl=16)
    report = ec.pu_report(rows, [1.0] * 20, top_k=(5,), min_positives=3)
    assert report.get("pu_recall_at_top_k", 5).value == pytest.approx(0.25)     # 5/20 by chance
    assert report.get("pu_positive_vs_unlabeled_rank_probability").value == 0.5
    assert report.get("pu_positive_median_percentile_rank").value == 0.5


def test_pu_prior_adjustment_is_tagged_model_assumption_and_states_the_prior():
    rows = labelled(4, 0, n_unl=16)
    values = [20, 19, 5, 4] + list(range(18, 2, -1))
    base = ec.pu_report(rows, values, min_positives=3).get("pu_positive_vs_unlabeled_rank_probability").value
    q = ec.pu_report(rows, values, class_prior=0.2, min_positives=3).get("pu_prior_adjusted_rank_probability")
    assert q.evidence_class == "model_assumption"
    assert "class prior of 0.2" in q.pu_assumption and "SCAR" in q.pu_assumption
    assert q.value == pytest.approx((base - 0.1) / 0.8)
    with pytest.raises(ec.SchemaError):
        ec.pu_report(rows, values, class_prior=1.0, min_positives=3)


def test_pu_universe_excludes_verified_negatives_censored_and_by_default_unassayed():
    rows = labelled(4, 3, n_unl=6, n_unassayed=5, n_censored=2)
    values = list(range(20))
    report = ec.pu_report(rows, values, min_positives=3)
    assert report.n_universe == 10
    assert report.excluded_by_state == {"verified_negative": 3, "unassayed": 5, "censored": 2}
    wider = ec.pu_report(rows, values, min_positives=3, include_unassayed=True)
    assert wider.n_universe == 15 and "unassayed" in wider.universe_states


def test_pu_blocked_with_too_few_independent_positives_and_bad_k():
    rows = labelled(4, 0, n_unl=6)
    for r in rows[:4]:
        r["independent_unit_id"] = "same_animal"
    with pytest.raises(ec.MetricBlocked, match="1 independent measured_positive"):
        ec.pu_report(rows, list(range(10)), min_positives=3)
    with pytest.raises(ec.SchemaError, match="top_k"):
        ec.pu_report(labelled(4, 0, n_unl=6), list(range(10)), top_k=(11,), min_positives=3)
    only_pos = ec.pu_report(labelled(4, 0), [1, 2, 3, 4], min_positives=3)
    assert "pu_positive_vs_unlabeled_rank_probability" in only_pos.blocked


# ------------------------------------------------------------------ clustered uncertainty
def mean_score(values):
    return lambda idx: sum(values[i] for i in idx) / len(idx)


def test_ci_refused_below_minimum_clusters():
    rows = independent_rows(9)
    with pytest.raises(ec.UncertaintyRefused, match="9 independent cluster"):
        ec.cluster_bootstrap_ci(rows, mean_score(list(range(9))), n_boot=50)
    ci = ec.cluster_bootstrap_ci(independent_rows(10), mean_score(list(range(10))), n_boot=200)
    assert ci.n_clusters == 10 and ci.lower <= ci.point <= ci.upper and ci.evidence_class == "derived"


def test_correlated_windows_are_one_cluster_even_with_distinct_unit_ids():
    rows = independent_rows(8)
    for i in range(8):          # 8 more rows overlapping the first window, each its own "unit"
        rows.append(row(f"w{i}", seq=rand_seq(random.Random(900 + i)), chem=f"cw{i}", study=f"sw{i}",
                        start=1002 + i, end=1019 + i))
    assert len(set(r["independent_unit_id"] for r in rows)) == 16
    clusters = ec.effective_clusters(rows)
    assert len(set(clusters)) == 8
    with pytest.raises(ec.UncertaintyRefused, match="8 independent cluster"):
        ec.cluster_bootstrap_ci(rows, mean_score([0.0] * 16), n_boot=20)
    with pytest.raises(ec.UncertaintyRefused, match="sequence.*window"):
        ec.cluster_bootstrap_ci(rows, mean_score([0.0] * 16), n_boot=20, merge_axes=())


def test_shared_unit_ids_and_duplicate_sequences_collapse_clusters():
    rows = independent_rows(12)
    rows[1]["independent_unit_id"] = rows[0]["independent_unit_id"]
    rows[3]["aso_sequence"] = rows[2]["aso_sequence"]
    assert len(set(ec.effective_clusters(rows))) == 10


def test_cluster_bootstrap_is_wider_than_row_bootstrap_for_clustered_data():
    rng = random.Random(5)
    rows, values = [], []
    for c in range(12):
        effect = rng.gauss(0, 1)
        for j in range(6):
            rows.append(row(f"c{c}_{j}", seq=rand_seq(rng), chem=f"k{c}_{j}", study=f"s{c}_{j}",
                            unit=f"cluster{c}", start=1000 * (6 * c + j + 1), end=1000 * (6 * c + j + 1) + 17))
            values.append(effect + rng.gauss(0, 0.05))
    clustered = ec.cluster_bootstrap_ci(rows, mean_score(values), n_boot=400, seed=1)
    for r in rows:
        r["independent_unit_id"] = r["observation_id"]          # (wrongly) call every row independent
    naive = ec.cluster_bootstrap_ci(rows, mean_score(values), n_boot=400, seed=1)
    assert clustered.n_clusters == 12 and naive.n_clusters == 72
    assert (clustered.upper - clustered.lower) > 1.5 * (naive.upper - naive.lower)


def test_bootstrap_is_deterministic_per_seed():
    rows = independent_rows(15)
    stat = mean_score([float(i * i % 7) for i in range(15)])
    a = ec.cluster_bootstrap_ci(rows, stat, n_boot=300, seed=4)
    b = ec.cluster_bootstrap_ci(rows, stat, n_boot=300, seed=4)
    c = ec.cluster_bootstrap_ci(rows, stat, n_boot=300, seed=5)
    assert a == b and (a.lower, a.upper) != (c.lower, c.upper)


def test_metric_ci_respects_guard_and_degenerate_replicates():
    rows = labelled(12, 12)
    values = [i + (18 if i < 12 else 0) for i in range(24)]
    ci = ec.bootstrap_metric_ci("auroc", rows, values, policy=ec.MetricPolicy(10, 10), n_boot=200, seed=2)
    assert 0.5 < ci.lower <= ci.point <= ci.upper <= 1.0 and ci.n_clusters == 24
    with pytest.raises(ec.MetricBlocked):
        ec.bootstrap_metric_ci("auroc", labelled(12, 2), list(range(14)), n_boot=20)
    with pytest.raises(ec.MetricBlocked):
        ec.bootstrap_metric_ci("auprc", labelled(12, 0, n_unl=12), list(range(24)), n_boot=20)
    # 11 positives + 1 negative: passes a (deliberately lax) policy but ~35% of resamples
    # contain no negative at all -> the interval is refused rather than fabricated
    with pytest.raises(ec.UncertaintyRefused, match="degenerate"):
        ec.bootstrap_metric_ci("auroc", labelled(11, 1), list(range(12)),
                               policy=ec.MetricPolicy(1, 1), n_boot=200)


# ------------------------------------------------------------------ forbidden sources
@pytest.mark.parametrize("path", [
    "SMN2_Exon7_7_27+28_Complete_Data.xlsx",
    "/data/lab/SMN2_Exon7_7_27+28_Complete_Data (1).xlsx",
    "exports/smn2_exon7_7_27+28_complete_data_tidy.csv",
    "exports/SMN2_Exon7_7_27%2B28_Complete_Data.parquet",
    "C:\\lab\\SMN2_Exon7_7_27+28_Complete_Data.xlsx",
    "results/psi25/table.csv", "psi25/", "/x/PSI25/y.json", "runs/psi25",
    "feng/terminal/outcomes.parquet", "/home/u/repo/feng/terminal", "../feng/terminal/a/b.csv",
    "feng//terminal/x", "repo/feng/./terminal/x",
])
def test_forbidden_sources_are_refused(path):
    with pytest.raises(ec.ForbiddenSourceError):
        ec.check_source_not_forbidden(path)
    with pytest.raises(ec.ForbiddenSourceError):
        ec.read_table(path)
    with pytest.raises(ec.ForbiddenSourceError):
        ec.validate_table([row("a", provenance={"source": "lab", "file": path})])
    with pytest.raises(ec.ForbiddenSourceError):
        ec.validate_table([row("a", source_path=path)])


@pytest.mark.parametrize("path", ["features.parquet", "feng/src/feng/scan.py", "x/psi250/t.csv",
                                  "feng/terminals/x", "SMN2_Exon7_design_notes.md",
                                  "j9_smn2_v1/ranks.parquet"])
def test_unrelated_paths_are_allowed(path):
    assert ec.check_source_not_forbidden(path) == path
    assert ec.check_source_not_forbidden(Path(path)) == path


def test_symlink_to_forbidden_directory_is_refused(tmp_path):
    target = tmp_path / "psi25"
    target.mkdir()
    (target / "t.jsonl").write_text("")
    link = tmp_path / "innocent_name"
    os.symlink(target, link)
    with pytest.raises(ec.ForbiddenSourceError):
        ec.read_table(link / "t.jsonl")


def test_read_table_round_trip_and_format_guard(tmp_path):
    path = tmp_path / "obs.jsonl"
    path.write_text("\n".join(json.dumps(r) for r in labelled(2, 1)) + "\n")
    assert len(ec.read_table(path)) == 3
    (tmp_path / "obs.csv").write_text("x")
    with pytest.raises(ec.SchemaError, match="unsupported"):
        ec.read_table(tmp_path / "obs.csv")


def test_ablation_refuses_forbidden_source_paths_before_running_anything():
    rows = labelled(6, 6)
    plan = ec.make_splits(rows, "leave_one_aso_out")
    called = []
    with pytest.raises(ec.ForbiddenSourceError):
        ec.run_ablation(rows, {}, "M1", plan, lambda *a: called.append(1),
                        source_paths=["ok.parquet", "lab/SMN2_Exon7_7_27+28_Complete_Data.xlsx"])
    assert not called


# ------------------------------------------------------------------ baselines / ablation
def test_baseline_registry_is_cumulative_and_evidence_tagged():
    assert sorted(ec.BASELINES) == ["M0", "M1", "M2", "M3", "M4", "M5", "M6"]
    assert ec.baseline_feature_columns("M0") == ()
    previous = ()
    for name in ("M1", "M2", "M3", "M4", "M5"):
        cols = ec.baseline_feature_columns(name)
        assert cols[:len(previous)] == previous and len(cols) > len(previous)
        assert all(c.startswith(("m1_", "m2_", "m3_", "m4_", "m5_")) for c in cols)
        previous = cols
    fams = {f.name: f.evidence_class for f in ec.BASELINES["M5"].families}
    assert fams == {"complementarity": "derived", "duplex": "model_assumption",
                    "accessibility": "proxy", "splice_geometry": "derived",
                    "sequence_ablation_proxy": "proxy"}
    with pytest.raises(ec.BaselineUndefined, match="M6"):
        ec.baseline_feature_columns("M6")
    with pytest.raises(ec.SchemaError):
        ec.baseline_feature_columns("M7")


def test_registry_columns_exist_in_the_j9_feature_vocabulary():
    text = (REPO_ROOT / "conference" / "model_v2" / "j9_smn2_prediction.py").read_text()
    for col in ec.baseline_feature_columns("M5"):
        assert f'"{col}"' in text, col


def features_for(rows, rng):
    cols = ec.baseline_feature_columns("M5")
    return {r["observation_id"]: {c: rng.random() for c in cols} for r in rows}


def feature_lineage(features):
    return {"source_paths": ("synthetic/features.json",),
            "feature_sources": {oid: {col: ("synthetic/features.json",) for col in values}
                                for oid, values in features.items()}}


def test_ablation_returns_typed_per_fold_groups_and_blocked_reasons():
    rows = labelled(6, 6)
    plan = ec.make_splits(rows, "leave_one_aso_out", n_folds=3, seed=1)
    feats = features_for(rows, random.Random(0))
    seen = []

    def fit_predict(train, test, columns):
        assert all("label_state" in r for r in train)
        assert all("label_state" not in r for r in test)          # labels never reach prediction
        assert not {r["observation_id"] for r in train} & {r["observation_id"] for r in test}
        assert all(set(r["features"]) == set(columns) for r in train + test)
        seen.append(columns)
        return {r["observation_id"]: sum(v for v in r["features"].values()) for r in test}

    result = ec.run_ablation(rows, feats, "M3", plan, fit_predict, policy=POLICY,
                             **feature_lineage(feats))
    assert isinstance(result, ec.AblationResult) and result.mode == "supervised"
    variants = [c.variant for c in result.pooled]
    assert variants == ["full", "drop:complementarity", "drop:duplex", "drop:accessibility"]
    assert len(result.cells) == 4 * 3 and len(seen) == 12
    drop = next(c for c in result.pooled if c.dropped_family == "duplex")
    assert "m2_duplex_energy_kcal_mol" not in drop.feature_columns
    for cell in result.cells:
        assert isinstance(cell, ec.AblationCell)
        fold = plan.folds[cell.fold_id]
        assert cell.test_group_ids == fold.test_group_ids and cell.test_group_ids
        assert cell.held_out_groups == fold.held_out_groups
        assert cell.n_test == 4 and not cell.metrics            # 4-row folds can't clear 3+3
        assert set(cell.blocked) == {"auroc", "auprc"}
        assert all("minimum 3" in reason for reason in cell.blocked.values())
    for cell in result.pooled:                                   # pooled 6+6 clears the guard
        assert {m.name for m in cell.metrics} == {"auroc", "auprc"} and not cell.blocked
    assert len(result.blocked_reasons()) == 24


def test_ablation_m0_is_a_seeded_permuted_reference_needing_no_model():
    rows = labelled(6, 6)
    plan = ec.make_splits(rows, "leave_one_aso_out")
    a = ec.run_ablation(rows, {}, "M0", plan, policy=POLICY, seed=9)
    b = ec.run_ablation(rows, {}, "M0", plan, policy=POLICY, seed=9)
    c = ec.run_ablation(rows, {}, "M0", plan, policy=POLICY, seed=10)
    assert a == b and a.pooled[0].variant == "permuted_reference"
    assert a.pooled[0].metrics[0].value != c.pooled[0].metrics[0].value
    scores = ec.permuted_reference_scores(["b", "a", "c"], 3)
    assert scores == ec.permuted_reference_scores(["c", "b", "a"], 3)


def test_ablation_pu_mode_never_emits_supervised_metrics_and_records_blocks():
    rows = labelled(6, 0, n_unl=18)
    plan = ec.make_splits(rows, "leave_one_aso_out", n_folds=3)
    feats = features_for(rows, random.Random(1))
    result = ec.run_ablation(rows, feats, "M1", plan,
                             lambda tr, te, cols: {r["observation_id"]: sum(r["features"].values()) for r in te},
                             pu=True, policy=POLICY, top_k=(5, 500), **feature_lineage(feats))
    assert result.mode == "positive_unlabeled"
    for cell in result.cells + result.pooled:
        assert cell.metrics == ()
        assert not set(cell.blocked) & set(ec.SUPERVISED_METRICS)
        assert "pu_recall_at_top_k[500]" in cell.blocked or cell.pu is None
    pooled = result.pooled[0]
    assert pooled.pu is not None and pooled.pu.n_universe == 24
    assert {q.name for q in pooled.pu.quantities} <= set(ec.PU_QUANTITIES)
    # the same data in supervised mode is blocked everywhere, with reasons
    sup = ec.run_ablation(rows, feats, "M1", plan,
                          lambda tr, te, cols: {r["observation_id"]: 0.0 for r in te}, policy=POLICY,
                          **feature_lineage(feats))
    assert all(not c.metrics and set(c.blocked) == {"auroc", "auprc"} for c in sup.cells + sup.pooled)


def test_ablation_input_contract_errors():
    rows = labelled(6, 6)
    plan = ec.make_splits(rows, "leave_one_aso_out")
    feats = features_for(rows, random.Random(2))
    with pytest.raises(ec.BaselineUndefined):
        ec.run_ablation(rows, feats, "M6", plan, lambda *a: {})
    with pytest.raises(ec.SchemaError, match="fit_predict"):
        ec.run_ablation(rows, feats, "M1", plan)
    with pytest.raises(ec.SchemaError, match="features missing"):
        ec.run_ablation(rows, {}, "M1", plan, lambda *a: {})
    with pytest.raises(ec.SchemaError, match="exactly the test rows"):
        ec.run_ablation(rows, feats, "M1", plan, lambda tr, te, cols: {}, **feature_lineage(feats))
    with pytest.raises(ec.SchemaError, match="different table"):
        ec.run_ablation(rows[:-1], feats, "M1", plan, lambda *a: {})
    with pytest.raises(ec.SchemaError, match="unknown metric"):
        ec.run_ablation(rows, feats, "M1", plan, lambda *a: {}, metrics=("p_value",))


def test_register_baseline_defines_m6_explicitly(monkeypatch):
    monkeypatch.setattr(ec, "BASELINES", dict(ec.BASELINES))
    fam = ec.FeatureFamily("offtarget_burden", ("m6_example",), "derived", "M6")
    ec.register_baseline(ec.BaselineSpec("M6", "feature_set", ec.BASELINES["M5"].families + (fam,), "test"))
    assert ec.baseline_feature_columns("M6")[-1] == "m6_example"
    with pytest.raises(ec.SchemaError):
        ec.register_baseline(ec.BaselineSpec("M6", "feature_set", (
            ec.FeatureFamily("x", ("c",), "measured", "M6"),), "bad evidence class"))


# F3 review repros: these must fail on the original implementation.
@pytest.mark.parametrize("metric", ["auprc", "auroc", "fdr", "accuracy"])
def test_r1_public_table_cannot_bypass_negative_evidence(metric):
    rows = labelled(10, 10)
    for r in rows[10:]:
        r.pop("negative_evidence")
    with pytest.raises(ec.NegativeEvidenceError):
        ec.compute_metric(metric, ec.ObservationTable(tuple(rows)), [1.] * 10 + [0.] * 10,
                          threshold=.5)


def test_r1_validated_rows_are_deeply_immutable_and_detached():
    rows = labelled(1, 1)
    table = ec.validate_table(rows)
    rows[1]["negative_evidence"]["assayed_and_powered"] = False
    assert table.rows[1]["negative_evidence"]["assayed_and_powered"] is True
    with pytest.raises(TypeError):
        table.rows[0]["label_state"] = "verified_negative"
    with pytest.raises(TypeError):
        table.rows[1]["negative_evidence"]["detection_limit"] = None


@pytest.mark.parametrize("bad", [0, -1, float("nan"), float("inf"), 1.5, True])
@pytest.mark.parametrize("field", ["min_verified_negatives", "min_positives"])
def test_r2_policy_requires_positive_integer_counts(bad, field):
    with pytest.raises(ec.SchemaError):
        policy = ec.MetricPolicy(**{field: bad})
        ec.compute_metric("auprc", labelled(10, 0), [1.] * 10, policy=policy)


@pytest.mark.parametrize("metadata", [
    {"study_id": "holgersen_s4"},
    {"provenance": {"source": "Holgersen S4"}},
    {"provenance": {"source": "exports/Holgersen_S4.tsv"}},
    {"source_table_kind": "responder-only"},
    {"source_table_kind": "RESPONDERS ONLY"},
    {"negative_evidence": {**NEG_EVIDENCE, "basis": "absence from responder list"}},
    {"negative_evidence": {**NEG_EVIDENCE, "basis": "unrecognized assertion"}},
])
def test_r3_responder_sources_and_absence_evidence_cannot_supply_negatives(metadata):
    rows = labelled(10, 10)
    for r in rows[10:]:
        r.update(metadata)
    with pytest.raises(ec.NegativeEvidenceError):
        ec.compute_metric("auroc", rows, [1.] * 10 + [0.] * 10)


@pytest.mark.parametrize("pu", [False, True])
def test_r4_renamed_assay_copies_do_not_inflate_independent_n(pu):
    rows = [row("p" + str(i)) for i in range(10)] + [
        row("n" + str(i), "unlabeled" if pu else "verified_negative") for i in range(10)]
    with pytest.raises(ec.MetricBlocked, match="1 independent"):
        if pu:
            ec.pu_report(rows, [1.] * 10 + [0.] * 10)
        else:
            ec.compute_metric("auroc", rows, [1.] * 10 + [0.] * 10)


@pytest.mark.parametrize("mode", sorted(ec.MODES))
@pytest.mark.parametrize("n_folds", [None, 2])
def test_r5_independent_unit_is_an_axis_in_every_split(mode, n_folds):
    rows = independent_rows(8)
    rows[1]["independent_unit_id"] = rows[0]["independent_unit_id"]
    nested = ec.make_nested_splits(rows, outer_mode=mode, outer_n_folds=n_folds,
                                   inner_n_folds=2, min_inner_groups=2)
    for outer, inner in zip(nested.outer.folds, nested.inner):
        for fold in (outer,) + inner:
            assert not {rows[i]["independent_unit_id"] for i in fold.train_idx} & {
                rows[i]["independent_unit_id"] for i in fold.test_idx}
        for fold in inner:
            assert not {rows[i]["independent_unit_id"] for i in outer.test_idx} & {
                rows[i]["independent_unit_id"] for i in fold.train_idx + fold.test_idx}


def test_r6_inner_splits_preserve_transitive_outer_groups_after_bridge_purge():
    rows = independent_rows(7)
    for i, seq in enumerate(("A" * 18, "A" * 16 + "CC", "A" * 14 + "CCCC", "G" * 18)):
        rows[i]["aso_sequence"] = seq
    rows[3].update(start=rows[1]["start"], end=rows[1]["end"])
    nested = ec.make_nested_splits(rows, outer_mode="leave_study_out")
    for inner in nested.inner:
        for fold in inner:
            for groups in nested.outer.group_ids.values():
                assert not {groups[i] for i in fold.test_idx} & {groups[i] for i in fold.train_idx}


@pytest.mark.parametrize("alias", ["2'-MOE", "2-MOE", "2′ MOE", "2′-O-methoxyethyl"])
def test_r7_moe_aliases_cannot_straddle(alias):
    rows = independent_rows(6)
    rows[0]["chemistry_family"] = "MOE"
    rows[1]["chemistry_family"] = alias
    plan = ec.make_splits(rows, "leave_study_out")
    assert plan.group_ids["chemistry"][0] == plan.group_ids["chemistry"][1]
    assert 1 not in plan.folds[0].train_idx


@pytest.mark.parametrize("chem", ["MOE/LNA", "modified", "unknown", "PS", "MOE?", "new chemistry"])
def test_r7_ambiguous_or_unknown_chemistry_is_rejected(chem):
    with pytest.raises(ec.SchemaError, match="chemistry"):
        ec.axis_groups([row("a", chem=chem)], axes=("chemistry",))


@pytest.mark.parametrize("field,value", [("aso_sequence", "A" * 18), ("start", 1),
                                        ("independent_unit_id", "new-unit")])
def test_r8_plan_is_bound_to_metadata_not_only_ids(field, value):
    rows = labelled(6, 6)
    plan = ec.make_splits(rows, "leave_study_out", n_folds=3)
    rows[0][field] = value
    with pytest.raises((ec.SchemaError, ec.LeakageError)):
        ec.run_ablation(rows, {}, "M0", plan, policy=ec.MetricPolicy(1, 1))


def test_r9_feature_sources_required_before_fit():
    rows = labelled(6, 6)
    plan = ec.make_splits(rows, "leave_study_out", n_folds=3)
    feats = {r["observation_id"]: {"m1_wc_pairs": 1.} for r in rows}
    calls = []
    def fit(train, test, cols):
        calls.append(1)
        return {r["observation_id"]: .5 for r in test}
    with pytest.raises(ec.SchemaError, match="source"):
        ec.run_ablation(rows, feats, "M1", plan, fit, policy=ec.MetricPolicy(1, 1))
    assert not calls


@pytest.mark.parametrize("path", ["%2570si25/a", "ｐｓｉ２５/a", r"feng\.\terminal/a",
                                  r"feng\safe\..\terminal/a", "%252570si25/a"])
def test_r10_encoded_unicode_and_separator_aliases_are_forbidden(path):
    with pytest.raises(ec.ForbiddenSourceError):
        ec.check_source_not_forbidden(path)


def test_r11_excluded_rows_do_not_clear_ci_cluster_minimum():
    rows = labelled(4, 4, n_unl=20)
    with pytest.raises(ec.UncertaintyRefused, match="8 independent cluster"):
        ec.bootstrap_metric_ci("auroc", rows, [1.] * 4 + [0.] * 24,
                               policy=ec.MetricPolicy(1, 1), n_boot=300, seed=1)


def test_r12_nonfinite_point_cannot_have_finite_ci():
    with pytest.raises(ec.UncertaintyRefused, match="point"):
        ec.cluster_bootstrap_ci(independent_rows(10),
            lambda idx: float("nan") if tuple(idx) == tuple(range(10)) else 1., n_boot=20)


@pytest.mark.parametrize("kwargs", [
    {"min_clusters": float("nan")}, {"min_clusters": 0}, {"min_clusters": 1.5},
    {"min_clusters": True}, {"n_boot": float("nan")}, {"n_boot": True},
    {"max_failed_fraction": float("nan")}, {"max_failed_fraction": -1},
    {"max_failed_fraction": 1.1}, {"alpha": True},
])
def test_r12_invalid_ci_controls_are_rejected(kwargs):
    with pytest.raises(ec.SchemaError):
        ec.cluster_bootstrap_ci(independent_rows(10), lambda idx: 1., **kwargs)


def test_r13_pu_ablation_records_k_outside_eligible_universe():
    rows = labelled(4, 3, n_unl=6, n_unassayed=5, n_censored=2)
    plan = ec.make_splits(rows, "leave_study_out")
    result = ec.run_ablation(rows, {}, "M0", plan, pu=True,
                             policy=ec.MetricPolicy(1, 1), top_k=(15,))
    assert result.pooled[0].pu.n_universe == 10
    assert "universe of 10" in result.pooled[0].blocked["pu_recall_at_top_k[15]"]


def test_r4_duplicate_assay_records_are_reported_and_never_independent():
    rows = [row("copy1"), row("copy2")]
    table = ec.validate_table(rows)
    assert table.duplicate_assay_records == (("copy1", "copy2"),)
    assert len(set(ec.effective_clusters(table, merge_axes=()))) == 1


def test_r4_sequence_bridge_connects_unit_counts_transitively():
    rows = labelled(3, 3)
    for i, seq in enumerate(("A" * 18, "A" * 16 + "CC", "A" * 14 + "CCCC")):
        rows[i]["aso_sequence"] = seq
    with pytest.raises(ec.MetricBlocked, match="1 independent measured_positive"):
        ec.compute_metric("auroc", rows, [1.] * 3 + [0.] * 3, policy=POLICY)


@pytest.mark.parametrize("kind", ["config", "groups", "folds"])
def test_r8_recomputes_plan_before_fit_even_when_ids_and_hash_are_unchanged(kind):
    rows = labelled(6, 6)
    plan = ec.make_splits(rows, "leave_study_out", n_folds=3)
    if kind == "config":
        plan = replace(plan, config=ec.GroupingConfig(window_padding=1))
    elif kind == "groups":
        plan = replace(plan, group_ids={**plan.group_ids, "sequence": tuple(map(str, range(12)))})
    else:
        plan = replace(plan, folds=plan.folds[:-1])
    feats = features_for(rows, random.Random(2))
    def fit(*args):
        pytest.fail("invalid plan reached fitter")
    with pytest.raises(ec.LeakageError):
        ec.run_ablation(rows, feats, "M1", plan, fit, **feature_lineage(feats))


@pytest.mark.parametrize("bad", ["missing_cell", "empty_paths", "unlisted", "forbidden", "no_manifest",
                                 "no_inventory"])
def test_r9_manifest_covers_every_feature_cell_and_checks_sources_before_fit(bad):
    rows = labelled(6, 6)
    plan = ec.make_splits(rows, "leave_study_out", n_folds=3)
    feats = features_for(rows, random.Random(2))
    lineage = feature_lineage(feats)
    oid = rows[0]["observation_id"]
    col = next(iter(feats[oid]))
    if bad == "missing_cell":
        del lineage["feature_sources"][oid][col]
    elif bad == "no_manifest":
        del lineage["feature_sources"]
    elif bad == "no_inventory":
        lineage["source_paths"] = ()
    else:
        lineage["feature_sources"][oid][col] = {
            "empty_paths": (), "unlisted": ("other.json",), "forbidden": ("%2570si25/data",)
        }[bad]
    def fit(*args):
        pytest.fail("invalid manifest reached fitter")
    with pytest.raises((ec.SchemaError, ec.ForbiddenSourceError)):
        ec.run_ablation(rows, feats, "M1", plan, fit, **lineage)


def test_r11_excluded_rows_do_not_change_resampling_or_cluster_counts():
    rows = labelled(6, 6, n_unl=20, n_unassayed=3, n_censored=2)
    scores = [1.] * 6 + [0.] * 31
    kwargs = dict(policy=ec.MetricPolicy(1, 1), n_boot=100, seed=1)
    full = ec.bootstrap_metric_ci("auroc", rows, scores, **kwargs)
    eligible = ec.bootstrap_metric_ci("auroc", rows[:12], scores[:12], **kwargs)
    assert full == eligible


def test_r12_ci_result_itself_rejects_nonfinite_point():
    with pytest.raises(ec.UncertaintyRefused):
        ec.CIResult(float("nan"), 0., 1., .05, 10, 100, 0, 1, (), "test", "derived")


def test_r12_all_failed_replicates_refused_even_if_tolerance_is_one():
    with pytest.raises(ec.UncertaintyRefused, match="degenerate"):
        ec.cluster_bootstrap_ci(independent_rows(10),
            lambda idx: 1. if tuple(idx) == tuple(range(10)) else float("nan"),
            n_boot=20, max_failed_fraction=1.)


@pytest.mark.parametrize("n_pos,n_neg", [(10, 0), (0, 10)])
def test_r2_both_classes_required_at_lowest_valid_policy(n_pos, n_neg):
    for metric in ec.SUPERVISED_METRICS:
        with pytest.raises(ec.MetricBlocked):
            ec.compute_metric(metric, labelled(n_pos, n_neg), [.5] * 10,
                              policy=ec.MetricPolicy(1, 1), threshold=.5)


@pytest.mark.parametrize("bad", [0, float("nan"), True])
def test_r2_pu_minimum_requires_a_positive_integer(bad):
    with pytest.raises(ec.SchemaError):
        ec.pu_report(labelled(3, 0, n_unl=3), [.5] * 6, min_positives=bad)


@pytest.mark.parametrize("change", [{"seed": 17}, {"min_independent_groups": 2}])
def test_r8_fingerprint_binds_split_controls_even_when_folds_would_be_identical(change):
    rows = labelled(6, 6)
    plan = ec.make_splits(rows, "leave_study_out")
    with pytest.raises(ec.LeakageError):
        ec.run_ablation(rows, {}, "M0", replace(plan, **change))


def test_r11_excluded_bridge_cannot_split_preexisting_clusters():
    rows = labelled(6, 6)
    rows.append(row("bridge", "unlabeled", seq="A" * 18, start=1000, end=1517))
    with pytest.raises(ec.UncertaintyRefused, match="11 independent cluster"):
        ec.bootstrap_metric_ci("auroc", rows, [1.] * 6 + [0.] * 7,
                               policy=ec.MetricPolicy(1, 1), min_clusters=12, n_boot=30)


# F3-FIX2: re-review reproductions and boundary variants (synthetic data only).
@pytest.mark.parametrize("axes", [(), ("sequence",), ("window",), ("study",)])
def test_fix2_ci_cannot_disable_correlation_axes(axes):
    rows = independent_rows(12)
    for r in rows:
        r.update(aso_sequence="A" * 18, start=1000, end=1017)
    with pytest.raises(ec.UncertaintyRefused):
        ec.cluster_bootstrap_ci(rows, mean_score(list(range(12))), merge_axes=axes, n_boot=20)


@pytest.mark.parametrize("field", ["sequence_threshold", "window_padding",
                                  "min_containment_length", "max_unique_sequences"])
@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -1, 0.5, True, "2", None])
def test_fix2_grouping_integer_controls(field, bad):
    with pytest.raises(ec.SchemaError):
        ec.GroupingConfig(**{field: bad})


@pytest.mark.parametrize("field", ["min_containment_length", "max_unique_sequences"])
def test_fix2_grouping_positive_controls(field):
    with pytest.raises(ec.SchemaError):
        ec.GroupingConfig(**{field: 0})


@pytest.mark.parametrize("field", ["containment_is_identical", "strand_aware"])
@pytest.mark.parametrize("bad", [0, 1, None, "false"])
def test_fix2_grouping_boolean_controls(field, bad):
    with pytest.raises(ec.SchemaError):
        ec.GroupingConfig(**{field: bad})


@pytest.mark.parametrize("bad", [None, "study", ([],), ("study", "study")])
def test_fix2_grouping_relax_is_validated(bad):
    with pytest.raises(ec.SchemaError):
        ec.GroupingConfig(relax=bad)


def test_fix2_grouping_relax_is_detached():
    axes = ["chemistry"]
    cfg = ec.GroupingConfig(relax=axes)
    axes.append("sequence")
    assert cfg.relax == ("chemistry",)


@pytest.mark.parametrize("left,right", [("chr5", "5"), ("chrX", "X"),
                                        ("chrY", "Y"), ("chrM", "MT"), ("M", "chrMT")])
def test_fix2_chromosome_aliases_never_straddle(left, right):
    rows = independent_rows(8)
    rows[0]["chrom"] = left
    rows[1].update(chrom=right, start=rows[0]["start"], end=rows[0]["end"])
    for mode in ec.MODES:
        plan = ec.make_splits(rows, mode)
        assert plan.group_ids["window"][0] == plan.group_ids["window"][1]
        for fold in plan.folds:
            assert not ({0, 1} & set(fold.train_idx) and {0, 1} & set(fold.test_idx))
    assert ec.effective_clusters(rows)[0] == ec.effective_clusters(rows)[1]


@pytest.mark.parametrize("chrom", ["chr0", "chr23", "05", "NC_000005.10", "unknown"])
def test_fix2_unrecognized_chromosomes_refused(chrom):
    with pytest.raises(ec.SchemaError, match="chrom"):
        ec.axis_groups([row("a", chrom=chrom)], axes=("window",))


@pytest.mark.parametrize("entry", ["outer", "nested_outer", "nested_inner"])
@pytest.mark.parametrize("bad", [float("nan"), float("inf"), 0, -1, 1.5, True, "3", None])
def test_fix2_split_minima_are_positive_integers(entry, bad):
    rows = independent_rows(8)
    with pytest.raises(ec.SchemaError):
        if entry == "outer":
            ec.make_splits(rows, "leave_study_out", min_independent_groups=bad)
        else:
            kw = {"min_outer_groups" if entry == "nested_outer" else "min_inner_groups": bad}
            ec.make_nested_splits(rows, outer_mode="leave_study_out", **kw)


@pytest.mark.parametrize("metric", ["auroc", "auprc", "accuracy", "precision", "recall",
                                   "fdr", "brier", "calibration_error", "pr_curve", "ci"])
def test_fix2_exact_assay_copies_cannot_reweight_metrics(metric):
    rows = labelled(10, 10)
    scores = [.5] * 10 + [1.] + [0.] * 9
    assert ec.compute_metric("auroc", rows, scores).value == pytest.approx(.9)
    rows += [{**rows[10], "observation_id": f"copy{i}", "independent_unit_id": f"u{i}"}
             for i in range(20)]
    scores += [1.] * 20
    with pytest.raises(ec.MetricBlocked, match="duplicate assay"):
        if metric == "pr_curve":
            ec.precision_recall_curve(rows, scores)
        elif metric == "ci":
            ec.bootstrap_metric_ci("auroc", rows, scores, n_boot=20)
        else:
            ec.compute_metric(metric, rows, scores, threshold=.5)


@pytest.mark.parametrize("state", ["measured_positive", "unlabeled", "unassayed"])
def test_fix2_pu_refuses_duplicate_candidates(state):
    rows = labelled(10, 0, n_unl=10, n_unassayed=2)
    source = next(r for r in rows if r["label_state"] == state)
    rows.append({**source, "observation_id": "copy", "independent_unit_id": "copy_unit"})
    with pytest.raises(ec.MetricBlocked, match="duplicate assay"):
        ec.pu_report(rows, list(range(len(rows))), top_k=(5,), include_unassayed=True)


def test_fix2_generic_ci_refuses_duplicate_contributors():
    rows = independent_rows(12)
    rows.append({**rows[0], "observation_id": "copy", "independent_unit_id": "copy_unit"})
    with pytest.raises(ec.UncertaintyRefused, match="duplicate assay"):
        ec.cluster_bootstrap_ci(rows, mean_score(list(range(13))), n_boot=20)


@pytest.mark.parametrize("metric", ["auroc", "auprc", "accuracy", "precision", "recall",
                                   "fdr", "brier", "calibration_error", "pr_curve", "ci"])
def test_fix3_annotations_cannot_disguise_duplicate_assays(metric):
    rows = labelled(10, 10)
    scores = [.5] * 10 + [1.] + [0.] * 9
    assert ec.compute_metric("auroc", rows, scores).value == pytest.approx(.9)
    rows += [{**rows[10], "observation_id": f"copy{i}", "independent_unit_id": f"u{i}",
              "note": f"irrelevant {i}",
              "provenance": {"source": "synthetic:test", "note": f"import {i}"},
              "negative_evidence": {**NEG_EVIDENCE, "note": f"annotation {i}"}}
             for i in range(20)]
    scores += [1.] * 20
    with pytest.raises(ec.MetricBlocked, match="duplicate assay"):
        if metric == "pr_curve":
            ec.precision_recall_curve(rows, scores)
        elif metric == "ci":
            ec.bootstrap_metric_ci("auroc", rows, scores, n_boot=20)
        else:
            ec.compute_metric(metric, rows, scores, threshold=.5)


@pytest.mark.parametrize("state", ["measured_positive", "unlabeled", "unassayed"])
def test_fix3_pu_refuses_annotated_copies(state):
    rows = labelled(10, 0, n_unl=10, n_unassayed=2)
    source = next(r for r in rows if r["label_state"] == state)
    rows.append({**source, "observation_id": "copy", "independent_unit_id": "copy_unit",
                 "note": "different annotation"})
    with pytest.raises(ec.MetricBlocked, match="duplicate assay"):
        ec.pu_report(rows, list(range(len(rows))), top_k=(5,), include_unassayed=True)


def test_fix3_generic_ci_refuses_annotated_copies_before_statistic():
    rows = independent_rows(12)
    rows.append({**rows[0], "observation_id": "copy", "independent_unit_id": "copy_unit",
                 "note": "different annotation"})

    def statistic(indices):
        pytest.fail("duplicate contributors must be refused before resampling")

    with pytest.raises(ec.UncertaintyRefused, match="duplicate assay"):
        ec.cluster_bootstrap_ci(rows, statistic, n_boot=20)


@pytest.mark.parametrize("field,value", [("assay_id", "assay2"), ("study_id", "S2"),
                                         ("aso_sequence", "T" * 18),
                                         ("chemistry_family", "chemZ"), ("chrom", "chr6"),
                                         ("start", 1001), ("end", 1018), ("strand", "-")])
def test_fix3_distinct_assay_defining_fields_preserved(field, value):
    table = ec.validate_table([row("a"), row("b", **{field: value})])
    assert table.duplicate_assay_records == ()


@pytest.mark.parametrize("metric", ["accuracy", "precision", "recall", "fdr"])
@pytest.mark.parametrize("threshold", [None, float("nan"), float("inf"), -float("inf"),
                                       True, False, "0.5", [], {}])
def test_fix3_classification_requires_explicit_finite_numeric_threshold(metric, threshold):
    rows = labelled(10, 10)
    with pytest.raises(ec.SchemaError, match="threshold"):
        ec.compute_metric(metric, rows, [.9] * 10 + [.1] * 10, threshold=threshold)


@pytest.mark.parametrize("threshold", [0, 1, .5, -1.])
def test_fix3_finite_numeric_thresholds_are_valid(threshold):
    rows = labelled(10, 10)
    result = ec.compute_metric("accuracy", rows, [.9] * 10 + [.1] * 10, threshold=threshold)
    assert result.value == (1. if threshold == .5 else .5)


@pytest.mark.parametrize("limit", [-1, -.001, float("inf"), -float("inf")])
def test_fix3_invalid_detection_limits_refused(limit):
    evidence = {**NEG_EVIDENCE, "detection_limit": limit}
    with pytest.raises(ec.NegativeEvidenceError, match="detection_limit"):
        ec.validate_table([row("n", "verified_negative", negative_evidence=evidence)])


def test_fix3_units_do_not_authorize_unsupported_signed_detection_limits():
    evidence = {**NEG_EVIDENCE, "detection_limit": -1, "unit": "log10",
                "convention": "signed"}
    with pytest.raises(ec.NegativeEvidenceError, match="detection_limit"):
        ec.validate_table([row("n", "verified_negative", negative_evidence=evidence)])


@pytest.mark.parametrize("limit", [0, .05])
def test_fix3_nonnegative_detection_limits_accepted(limit):
    evidence = {**NEG_EVIDENCE, "detection_limit": limit}
    assert len(ec.validate_table([row("n", "verified_negative", negative_evidence=evidence)])) == 1


_FIX4_ALIASES = [
    ("assay_id", lambda s, i: s + "\u200b" * (i + 1), 20),
    ("study_id", lambda s, i: s + "\u200b" * (i + 1), 20),
    ("assay_id", lambda s, i: "\u0430ssay1", 1),
    ("study_id", lambda s, i: "\u040510", 1),
    ("assay_id", lambda s, i: s + " " * (i + 1), 20),
    ("assay_id", lambda s, i: " ＡＳＳＡＹ１ ", 1),
    ("study_id", lambda s, i: " Ｓ１０ ", 1),
    ("aso_sequence", lambda s, i: "".join(
        c.lower() if (i + 1) & (1 << j) else c for j, c in enumerate(s)), 20),
    ("aso_sequence", lambda s, i: s.lower().replace("t", "u"), 1),
    ("chrom", lambda s, i: "5", 1),
    ("chemistry_family", lambda s, i: s.upper(), 1),
    ("chemistry_family", lambda s, i: "2'-O-methoxyethyl", 1),
]


@pytest.mark.parametrize("field,alias,n_copies", _FIX4_ALIASES)
@pytest.mark.parametrize("consumer", ["auroc", "pr_curve", "metric_ci", "generic_ci", "pu"])
def test_fix4_identity_alias_copies_cannot_reweight_statistics(field, alias, n_copies, consumer):
    rows = labelled(10, 10)
    rows[10]["chemistry_family"] = "moe"
    scores = [.5] * 10 + [1.] + [0.] * 9
    assert ec.compute_metric("auroc", rows, scores).value == pytest.approx(.9)
    if consumer == "pu":
        rows[10]["label_state"] = "unlabeled"
    rows += [{**rows[10], field: alias(rows[10][field], i),
              "observation_id": f"copy{i}", "independent_unit_id": f"copy_unit{i}"}
             for i in range(n_copies)]
    scores += [1.] * n_copies
    if consumer == "generic_ci":
        def statistic(indices):
            pytest.fail("aliased copies must be refused before resampling")
        with pytest.raises(ec.SchemaError if field.endswith("_id") else ec.UncertaintyRefused,
                           match=field if field.endswith("_id") else "duplicate assay"):
            ec.cluster_bootstrap_ci(rows, statistic, n_boot=20)
    else:
        with pytest.raises(ec.SchemaError if field.endswith("_id") else ec.MetricBlocked,
                           match=field if field.endswith("_id") else "duplicate assay"):
            if consumer == "auroc":
                ec.compute_metric("auroc", rows, scores)
            elif consumer == "pr_curve":
                ec.precision_recall_curve(rows, scores)
            elif consumer == "metric_ci":
                ec.bootstrap_metric_ci("auroc", rows, scores, n_boot=20)
            else:
                ec.pu_report(rows, scores, top_k=(5,))


@pytest.mark.parametrize("field,alias,n_copies", _FIX4_ALIASES)
def test_fix4_identity_aliases_agree_in_validation_grouping_and_clusters(field, alias, n_copies):
    source = row("original", chem="moe", study="S10")
    copy = {**source, field: alias(source[field], 0),
            "observation_id": "copy", "independent_unit_id": "copy_unit"}
    if field.endswith("_id"):
        with pytest.raises(ec.SchemaError, match=field):
            ec.validate_table([source, copy])
        return
    for rows in ([source, copy], ec.validate_table([source, copy])):
        table = ec.validate_table(rows)
        assert table.duplicate_assay_records == (("original", "copy"),)
        assert table.rows[0][field] == table.rows[1][field]
        groups = ec.axis_groups(rows, axes=("sequence", "chemistry", "study", "window"))
        assert all(g[0] == g[1] for g in groups.values())
        clusters = ec.effective_clusters(rows, merge_axes=())
        assert clusters[0] == clusters[1]


def test_fix4_genuine_distinct_assay_replicates_still_contribute():
    rows = labelled(10, 10)
    rows.append({**rows[10], "assay_id": "replicate2", "observation_id": "replicate",
                 "independent_unit_id": "replicate_unit"})
    assert ec.validate_table(rows).duplicate_assay_records == ()
    clusters = ec.effective_clusters(rows, merge_axes=())
    assert len(set(clusters)) == 21
    assert ec.compute_metric("auroc", rows, [.5] * 10 + [1.] + [0.] * 9 + [1.]).value \
        == pytest.approx(9 / 11)


@pytest.mark.parametrize("field,value", [("aso_sequence", "ACGX"), ("aso_sequence", ""),
    ("chemistry_family", "MOE/LNA"), ("chrom", "chr23"), ("start", "1000"),
    ("start", True), ("end", 999), ("strand", "unknown")])
def test_fix4_invalid_identity_rejected_before_hashing(field, value):
    with pytest.raises(ec.SchemaError):
        ec.validate_table([row("a", **{field: value})])


_ID_FIELDS = ("assay_id", "study_id", "independent_unit_id", "observation_id")
_BAD_IDS = ("", " ", "id ", " id", "i d", "id\t", "id\n", "id\r", "id\x00",
            "id\x1f", "id\x7f", "id\u200b", "id\u200d", "id\ufeff", "id\u202e",
            "\u0430ssay1", "\uff21ssay1", "id\u00a0", "caf\u00e9", "id/1", "id@1",
            "a" * 129, None, 123, True)


@pytest.mark.parametrize("field", _ID_FIELDS)
@pytest.mark.parametrize("value", _BAD_IDS)
def test_fix5_invalid_identifiers_rejected(field, value):
    with pytest.raises(ec.SchemaError, match=field):
        ec.validate_table([row("o", **{field: value})])


@pytest.mark.parametrize("field", _ID_FIELDS)
@pytest.mark.parametrize("value", ["AazZ09._:+-", "a" * 128])
def test_fix5_identifier_vocabulary_preserved_exactly(field, value):
    table = ec.validate_table([row("o", **{field: value})])
    assert table.rows[0][field] == value


@pytest.mark.parametrize("field", _ID_FIELDS)
@pytest.mark.parametrize("entry", ["validate", "groups", "split", "nested", "metric", "pr",
                                    "pu", "clusters", "ci", "metric_ci", "verify", "ablation"])
def test_fix5_public_consumers_revalidate_forged_table_identifiers(field, entry):
    rows = labelled(10, 10)
    plan = ec.make_splits(rows, "leave_study_out")
    rows[0][field] += "\u200b"
    table = ec.ObservationTable(tuple(rows))
    scores = [.5] * 20

    def statistic(indices):
        pytest.fail("invalid identifiers must be rejected before statistic")
    calls = {
        "validate": lambda: ec.validate_table(table),
        "groups": lambda: ec.axis_groups(table),
        "split": lambda: ec.make_splits(table, "leave_study_out"),
        "nested": lambda: ec.make_nested_splits(table, outer_mode="leave_study_out"),
        "metric": lambda: ec.compute_metric("auroc", table, scores),
        "pr": lambda: ec.precision_recall_curve(table, scores),
        "pu": lambda: ec.pu_report(table, scores),
        "clusters": lambda: ec.effective_clusters(table),
        "ci": lambda: ec.cluster_bootstrap_ci(table, statistic),
        "metric_ci": lambda: ec.bootstrap_metric_ci("auroc", table, scores),
        "verify": lambda: ec.verify_no_leakage(plan, table),
        "ablation": lambda: ec.run_ablation(table, {}, "M0", plan),
    }
    with pytest.raises(ec.SchemaError, match=field):
        calls[entry]()


@pytest.mark.parametrize("field", _ID_FIELDS)
@pytest.mark.parametrize("suffix", [" ", "\u200b", "\u0430"])
@pytest.mark.parametrize("extension", ["json", "jsonl"])
def test_fix5_read_table_rejects_invalid_identifiers(tmp_path, field, suffix, extension):
    path = tmp_path / f"observations.{extension}"
    record = row("o")
    record[field] += suffix
    path.write_text(json.dumps([record] if extension == "json" else record), encoding="utf-8")
    with pytest.raises(ec.SchemaError, match=field):
        ec.read_table(path)


@pytest.mark.parametrize("field", ["assay_id", "study_id"])
def test_fix5_case_distinct_assays_are_genuine_replicates(field):
    rows = labelled(10, 10)
    rows[10][field] = "AssayA"
    rows.append({**rows[10], field: "assaya", "observation_id": "Replicate",
                 "independent_unit_id": "ReplicateUnit"})
    table = ec.validate_table(rows)
    assert table.rows[10][field] == "AssayA"
    assert table.rows[-1][field] == "assaya"
    assert table.duplicate_assay_records == ()
    assert len(set(ec.effective_clusters(table, merge_axes=()))) == 21
    assert ec.compute_metric("auroc", table, [.5] * 10 + [1.] + [0.] * 9 + [1.]).value \
        == pytest.approx(9 / 11)
    if field == "study_id":
        groups = ec.axis_groups(table, axes=("study",))["study"]
        assert groups[10] == "study:AssayA"
        assert groups[-1] == "study:assaya"


def test_fix5_case_distinct_observation_and_unit_ids_remain_distinct():
    rows = independent_rows(2)
    for r, oid, unit in zip(rows, ("ObsA", "obsa"), ("UnitA", "unita")):
        r.update(observation_id=oid, independent_unit_id=unit)
    assert ec.validate_table(rows).ids == ("ObsA", "obsa")
    assert len(set(ec.effective_clusters(rows, merge_axes=()))) == 2


def test_fix5_explicit_responder_only_study_ids_are_opaque():
    record = row("n", "verified_negative", study="StudyA")
    with pytest.raises(ec.NegativeEvidenceError):
        ec.validate_table([record], responder_only_studies=["StudyA"])
    assert len(ec.validate_table([record], responder_only_studies=["studya"])) == 1
    with pytest.raises(ec.SchemaError):
        ec.validate_table([], responder_only_studies=["StudyA "])


def test_fix5_permuted_reference_scores_validate_and_preserve_ids():
    with pytest.raises(ec.SchemaError):
        ec.permuted_reference_scores(["Obs\u200b"])
    assert set(ec.permuted_reference_scores(["ObsA", "obsa"])) == {"ObsA", "obsa"}
