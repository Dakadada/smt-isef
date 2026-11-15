F3-FIX2 — ADVANCE (software gate only): all five re-review findings fixed; none deferred.
HIGH 1 — `evaluation_core.py:1172`: CIs require both sequence/window merge axes; weakened axes refuse before the statistic callback.
HIGH 2 — `evaluation_core.py:392`: finite integer threshold/padding and positive size controls, strict booleans, validated/detached relax axes; public clustering revalidates config.
HIGH 3 — `evaluation_core.py:505`: canonical primary-chromosome names merge chr5/5, X/Y, and M/MT aliases before window clustering in every split/CI path.
MED 4 — `evaluation_core.py:970`: exact renamed assay copies among contributors block supervised metrics/PR curves, PU quantities, and CIs; copies cannot reweight scores.
MED 5 — `evaluation_core.py:667`, `:743`: positive integer minimum-group validation at outer and both nested entry points, plus the shared split builder.
Regression tests: `tests/test_evaluation_core.py:1117` onward adds 95 cases; pre-fix targeted run: 91 failed, 5 passed, 208 deselected; all five findings reproduced.
Existing bug-encoding test corrected: `tests/test_evaluation_core.py:586` previously expected merge_axes=() to inflate overlapping-window CI N to 16; now expects refusal. Other 208 cases unchanged.
Final task gate: **304 passed, 0 failed, 0 skipped (3.45 s)**; all 209 existing cases retained, including the corrected assertion.
Command: `PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 /home/dashii/orcd/pool/smt-isef-env/bin/python -m pytest conference/model_v2/tests/test_evaluation_core.py -q`.
Caller restrictions: resolve exact assay copies before statistics; CI merge_axes must include sequence/window; non-primary/unknown chromosome names (including accession aliases) fail closed.
Limitations: duplicate detection covers exact metadata copies excluding observation/unit IDs, not arbitrary near-duplicate assays; source manifests still rely on truthful declarations.
Only the two owned files and this requested report edited; no new dependencies, jobs, commits, or sealed/outcome data access. Unfixed findings: none.
