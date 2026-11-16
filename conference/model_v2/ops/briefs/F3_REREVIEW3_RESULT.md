# F3-REREVIEW3 — DO-NOT-ADVANCE

- MED `conference/model_v2/evaluation_core.py:342`: identity hashes raw values, so aliases evade duplicate detection and reweight metrics. Synthetic 10 positives/10 negatives: AUROC 0.9; append 20 renamed copies of one negative with distinct trailing spaces in `assay_id`, or distinct sequence letter casing, and AUROC becomes 0.3 with zero duplicate groups. Canonicalize/validate identity fields consistently before hashing; include chromosome and chemistry aliases.
- Repro: a single copy with `chrom='5'` instead of `chr5`, or uppercase chemistry instead of its accepted lowercase spelling, also evades detection (AUROC 0.9 → 0.81818). Field-order reversal is correctly detected as duplicate.
- Prior MED annotation repro: 20 renamed copies with differing notes/provenance/evidence annotations now raise `MetricBlocked`. Prior MED threshold repro: NaN and -inf, plus boolean and string thresholds, raise `SchemaError`.
- Detection limits NaN, +inf, -1 and True raise `NegativeEvidenceError`; 0 accepted. Signed/log-scale limits are explicitly unsupported.
- Replicates with distinct assay IDs remain distinct; same study/assay ID and identical remaining identity fields are flagged even with distinct independent-unit IDs. This follows the documented contract but requires callers to assign unique assay IDs to genuine replicates; reused IDs can falsely block them.
- Targeted pytest: 373 passed, 0 failed, 0 skipped (4.30 s). No source edits, submissions, sealed data, or full feng suite; task-specific gate used.
