# F3-FIX3 — fixed; targeted gate PASS
- `conference/model_v2/evaluation_core.py:202,342`: duplicate identity uses only study/assay IDs, sequence, chemistry and target window; notes, provenance annotations, labels and evidence metadata cannot disguise copies.
- `conference/model_v2/evaluation_core.py:905`: classification thresholds require explicit finite int/float values, excluding booleans; NaN, ±inf and nonnumbers raise SchemaError.
- `conference/model_v2/evaluation_core.py:316`: detection limits must be finite and nonnegative; reject all negative limits, including claimed signed/log units, because no validated signed convention exists here.
- `conference/model_v2/tests/test_evaluation_core.py:1227`: 69 added cases cover annotated 20-copy AUROC repro across metrics/PR/CI, PU, CI callback refusal, distinct assay fields, thresholds and detection limits.
- TDD: new cases before implementation: 49 failed, 20 passed, 304 deselected; failures reproduced the missing guards.
- Requested gate, smt-isef-env interpreter and specified thread/PYTHONPATH settings: 373 passed, 0 failed, 0 skipped (4.29 s); all 304 existing cases preserved; none encoded a bug.
- Unfixed findings: none. Limitation: callers must assign distinct assay_id values within a study for genuine replicate assays or differing experimental conditions; provenance remains caller-declared.
- git diff --check passed. Task-specific gate only; no full feng suite. No commits, submissions, jobs or sealed data access.
