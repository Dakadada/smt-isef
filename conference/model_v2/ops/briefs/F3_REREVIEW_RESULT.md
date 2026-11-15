# F3 re-review — DO-NOT-ADVANCE
Targeted gate: 209 passed, 0 failed, 0 skipped (2.50 s); no jobs or sealed data accessed.
Original R1–R13 repros rerun; synthetic chemistry aliases registered as in test fixture where needed.
R1–R3: invalid evidence/policy and responder-only source refused; R3 absence-basis variant refused.
R4: duplicate independent-unit inflation blocked, but duplicate score weighting remains (finding 4).
R5–R7: no unit/nested bridge straddle; MOE aliases group and purge correctly.
R8–R10: changed plan and missing feature manifests refused; all three encoded/Unicode/separator paths refused.
R11–R12: excluded rows do not inflate CI N; nonfinite point and NaN min_clusters refused.
R13: k=15 against PU universe 10 is recorded blocked (original pytest.raises now fails as expected).

- HIGH `conference/model_v2/evaluation_core.py:1128`: `merge_axes=()` gives a 12-cluster CI for 12 identical-sequence, overlapping-window rows; default correctly refuses N=1. Make sequence/window merging mandatory for uncertainty, or refuse weakened axes.
- HIGH `conference/model_v2/evaluation_core.py:406`: `GroupingConfig(sequence_threshold=float('nan'))` accepts NaN; two one-mismatch sequences split across train/test. Require finite nonnegative integer threshold/padding and validate all grouping controls.
- HIGH `conference/model_v2/evaluation_core.py:497`: overlapping `chr5` and `5` windows form separate groups and enter opposite train/test folds. Canonicalize GRCh38 chromosome names or reject noncanonical names before clustering.
- MED `conference/model_v2/evaluation_core.py:915`: duplicating one misranked verified negative 20 times changes AUROC 0.9 to 0.3 while independent-negative N remains 10. Reject/deduplicate exact assay copies or use explicit independent-unit weights for metrics and PU statistics.
- MED `conference/model_v2/evaluation_core.py:665`: `min_independent_groups=float('nan')` permits a two-group split despite the minimum-N guard. Require a positive integer at outer and inner split entry points.

All findings reproduced with synthetic in-memory rows; none fixed in this read-only review. Feature manifests still depend on truthful caller declarations, as noted in F3_FIX_RESULT. Task-specific pytest gate overrides the unrelated feng full-suite instruction.
