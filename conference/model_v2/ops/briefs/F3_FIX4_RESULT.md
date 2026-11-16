# F3-FIX4 — PASS (task-specific gate)
- `conference/model_v2/evaluation_core.py:253,308`: canonicalize identity before hashing; NFKC/strip/casefold IDs, existing sequence/alphabet, chromosome and chemistry normalizers; validate coordinates/strand. Duplicate, grouping and CI paths share validated rows.
- `conference/model_v2/tests/test_evaluation_core.py:1328`: all four rereview repros plus Unicode/study, U/T and chemistry vocabulary aliases block AUROC/PR/metric-CI/generic-CI/PU copy reweighting.
- `conference/model_v2/tests/test_evaluation_core.py:1357,1371,1385`: grouping/cluster agreement, genuine distinct-assay replicate contribution, invalid identity rejection; 57 added cases.
- `conference/model_v2/tests/test_evaluation_core.py:241,320`: two existing expectations encoded raw-case study-group behavior; updated to canonical lowercase; all 373 existing cases retained.
- Requested root command: 430 passed, 0 failed, 0 skipped (5.57 s); scoped diff whitespace check clean.
- Unfixed: none within scope. Compatibility: validated identity values/group IDs become canonical; supplied malformed optional identity fields now fail earlier; genuine replicates still require distinct assay IDs.
- Only two owned code/test files and this brief edited; no commits, submissions, sealed data, or full feng suite.
