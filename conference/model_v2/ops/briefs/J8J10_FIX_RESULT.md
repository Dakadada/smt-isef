# J8J10-FIX result
F1–F6 fixed; ready for adversarial rereview. Paths below are relative to conference/model_v2/.
F1 — j8_feature_join.py:249,336: indexed MAX(end-start) supplies a validated lower start bound; no chromosome-prefix scan; inclusive overlaps preserved across ingestion/restart.
F2 — j10_robustness.py:164: discovery-selected shift distributions retain observed/decoy scores but p_value=None with an explicit selection/exchangeability refusal.
F3 — j10_robustness.py:331: failed context/offset/channel yields blocked nulls carrying upstream reasons; candidates remain in rank coverage.
F4 — j8_feature_join.py:66,214,344,403; j10_robustness.py:33,318: guard resolved sidecars/catalogs/manifests before opening, including SQLite recovery files and output markers.
F5 — j10_robustness.py:407: verify J3 .done and J9 manifest rows/SHA before feature reads; contract records both marker paths and hashes.
F6 — j8_feature_join.py:361; j10_robustness.py:360: exact present J8 contract required; J10 schema + canonical payload SHA required; damaged same-contract results rebuild.
Regression coverage — tests/test_j8_feature_join.py:246; tests/test_j10_robustness.py:205: 33 new cases from repros, including forbidden-target symlinks and rehashed incomplete results.
Benchmark — tests/test_j8_feature_join.py:290: synthetic full-density slice, 1,000 configs × 216 cells = 216,000 rows; 100× prefix growth gives 58→41 miss VM steps; 216-hit query 1,556.
Measured — 329–340 dense annotations/s; 9,534–10,161 misses/s; 2.18 s build; 67,117,056-byte index, 2,449,408-byte width index; process peak RSS 171,520 KiB.
Sizing — jobs/j8_full.sbatch:13: 11,002,701 joins + prior ingestion/IO costs project ~4.4–13.4 h; request 24 h; RAM ~2 GiB estimate/8G allocation; disk ~37 GB + 0.05 GB width index.
Job templates — fresh j8_bounded_v2/j8_genome_v2/j10_bounded_v2 destinations; old catalogs/contracts intentionally require fresh outputs.
Tests — requested root J8/J10 pytest gate: 73 passed, 0 skipped, 0 failed (5.21 s); original 40 retained plus 33 regressions; shell syntax and scoped diff checks pass.
Old tests corrected — shift-null test encoded the F2 bug; CLI fixture omitted J3 SHA and J9 commit manifest, now supplies authenticated synthetic commits.
Unfixed — none of F1–F6. Limits: synthetic density slice, not full cardinality; overlapping loci/shared IO/output growth may exceed estimates; full-genome runtime/RSS and feng-wide gate not measured/run.
Scope — only owned source/tests/jobs and this requested report written; no real wet-lab/outcome reads, submissions, commits, or /orcd traversal; concurrent unrelated edits untouched.
