# F3-REREVIEW5 — ADVANCE
- LOW `conference/model_v2/evaluation_core.py:212`: callers must supply honest assay IDs: 20 identical negative content copies with new `A.0`…`A.19` IDs yield zero duplicate groups and AUROC 0.9→0.3; document explicitly that copies must retain the same assay ID and validate provenance upstream. Genuine replicates may have distinct IDs.
- `conference/model_v2/evaluation_core.py:22`: caller responsibility is only partially explicit: it requires distinct IDs for genuine replicates, but does not explicitly require identical IDs for copies; this is the only remaining finding.
- Prior MEDs fixed: `conference/model_v2/evaluation_core.py:195,263` rejects zero-width/non-ASCII IDs and preserves case-sensitive IDs; permitted distinct IDs such as `A.1` and `A_1` remain distinct.
- No further duplicate-assay, metric, leakage, or forbidden-source bypass established in this targeted review; guards depend on truthful caller-supplied IDs and source lineage.
- Targeted pytest: 639 passed, 0 failed, 0 skipped (6.51 s); bytecode and pytest cache disabled. Full feng suite not run (task-specific review gate).
- No source edits, submissions, or sealed data accessed.
