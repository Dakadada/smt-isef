**F3-FIX — ADVANCE (software gate only).** All R1–R13 fixed; none deferred. Changes limited to the two owned files plus this explicitly requested report.
Final gate: **209 passed, 0 failed, 0 skipped** (2.51 s): 122 existing cases retained + 87 added cases.
Command: `PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 /home/dashii/orcd/pool/smt-isef-env/bin/python -m pytest conference/model_v2/tests/test_evaluation_core.py -q`.
Regression evidence: `tests/test_evaluation_core.py:844` onward covers all 13 repros; initial run reproduced 62 failures, with additional duplicate/config/CI-bridge tests also observed failing before their fixes.

- R1 HIGH — `evaluation_core.py:341`: every public table consumer revalidates; validated rows and nested metadata are detached, immutable mappings/tuples.
- R2 HIGH — `evaluation_core.py:759`: positive integer policy/PU minima, rejecting bool/NaN/inf/fractions/zero; supervised quantities always require both classes.
- R3 HIGH — `evaluation_core.py:247`: canonical source/evidence categories; built-in Holgersen S4 study/provenance/path recognition; responder-only and unknown negative bases refused.
- R4 HIGH — `evaluation_core.py:1105`: supervised/PU minima use unit+sequence+window connected components; exact renamed assay copies reported in `duplicate_assay_records` and always clustered.
- R5 HIGH — `evaluation_core.py:517`: independent-unit identity is a mandatory, non-relaxable axis for every split and nested level.
- R6 HIGH — `evaluation_core.py:708`: inner folds inherit outer cluster identities, preserving transitive groups after bridge purging; inner/outer leakage checks enforced.
- R7 HIGH — `evaluation_core.py:376`: explicit chemistry alias registry unifies MOE/2'-MOE/2-MOE/2'-O-methoxyethyl; ambiguous/unknown families rejected.
- R8 HIGH — `evaluation_core.py:613`: SHA-256 data/config fingerprints include grouping and split controls; groups and complete folds are recomputed before fitting.
- R9 HIGH — `evaluation_core.py:1352`: feature baselines require nonempty source inventories and per-cell source manifests; coverage and every declared path checked before fitting; features snapshotted.
- R10 MED — `evaluation_core.py:106`: bounded fixed-point URL decoding/NFKC/separator normalization, dot-segment and symlink checks; non-convergence fails closed.
- R11 MED — `evaluation_core.py:1189`: CI counts/resamples only metric-contributing observations while retaining original cluster identities through excluded bridges.
- R12 MED — `evaluation_core.py:1126`: finite point/bounds, positive integer cluster/replicate controls, finite alpha/failure tolerance; all-degenerate draws refused even at tolerance=1.
- R13 MED — `evaluation_core.py:1324`: ablation top-k bounded by positive+unlabeled universe; out-of-range quantities recorded as blocked, including folds with insufficient positives.

Caller change: `feature_sources[observation_id][column] = (source_path, ...)` must cover provided/consumed features, with each path in nonempty `source_paths`; feature-free M0 remains exempt.
Caller changes: rebuild old plans; supply recognized negative-evidence `basis`; use canonical chemistry aliases or explicitly extend `CHEMISTRY_ALIASES`; metric inputs now require sequence/window metadata for independence checks.
Legacy bug-encoding input changed: `tests/test_evaluation_core.py:502` used `MetricPolicy(0,0)`; now `(1,1)` while retaining the PU-block assertion. No existing test removed.
Fixture updates: synthetic chemistry families explicitly registered; feature-ablation tests now declare source manifests. Existing numerical/grouping assertions retained.
Limitations: manifests authenticate neither file contents nor caller truthfulness; exact-duplicate reporting excludes only row/unit IDs, retains rows, and does not deduplicate point-estimate weights. Other correlations collapse through configured components.
Scope: no other source edits, dependencies, submissions, commits, sealed/outcome data access. Task-specific interpreter/test command overrides the standing feng gate; unrelated suites were not run.
