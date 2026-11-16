# F3-REREVIEW2 — DO-NOT-ADVANCE

- MED `conference/model_v2/evaluation_core.py:336`: duplicate identity hashes every metadata field; adding a different irrelevant `note` to each of 20 renamed copies of one verified negative bypasses the guard and changes AUROC from 0.9 to 0.3 while independent negative N stays 10. Define assay identity from canonical assay-defining fields, or require explicit deduplication/independent-unit weighting before all statistics, including PU and CIs.
- MED `conference/model_v2/evaluation_core.py:899`: `accuracy` with `threshold=float('nan')` returns 0.5 on 10 positive and 10 verified-negative rows because every comparison is false; `-inf` likewise yields FDR 0.5. Require an explicit finite, non-boolean numeric threshold before classification metrics.

Prior five repros rerun with synthetic rows: weakened CI axes refused; NaN sequence threshold refused; chr5/5 grouped; exact renamed copies refused (baseline AUROC 0.9); NaN split minimum refused. The first finding above extends the fourth fix.
NaN/inf detection limits refused; negative detection limit -1 is accepted (its physical validity depends on assay units and was not specified). Case, dot-segment, repeatedly encoded and Unicode forbidden paths refused in probes. Manifests still depend on truthful caller declarations; arbitrary in-memory source provenance cannot be authenticated.
Targeted gate: 304 passed, 0 failed, 0 skipped (4.02 s). No source edits, jobs, sealed data, or full feng suite; task-specific gate used. Findings remain unfixed in this read-only review.
