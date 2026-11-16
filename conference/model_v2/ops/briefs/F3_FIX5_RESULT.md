# F3-FIX5 — fixed; targeted gate PASS
- `conference/model_v2/evaluation_core.py:22,263`: documented/enforced opaque case-sensitive ASCII `[A-Za-z0-9._:+-]{1,128}` for all four IDs; no ID normalization.
- `conference/model_v2/evaluation_core.py:295,1386`: table entry points/read_table revalidate IDs; explicit responder-only study IDs compare exactly; M0 ID inputs validate too. Existing forbidden/known-source guards retained.
- `conference/model_v2/tests/test_evaluation_core.py:1332,1409`: trailing-space, 20-copy zero-width, Cyrillic/fullwidth, controls, length/type, JSON/JSONL, forged public-table and case-distinct replicate regressions; earlier sequence/chrom/chemistry/annotation guards retained.
- Updated old normalization tests: `test_fix4_identity_alias_copies_cannot_reweight_statistics`, `test_fix4_identity_aliases_agree_in_validation_grouping_and_clusters`: invalid ID aliases now SchemaError.
- Updated old casefold expectations: `test_trap_duplicate_aso_under_different_ids_never_straddles`, `test_leave_study_out_holds_out_whole_studies`: study groups preserve uppercase S.
- Exact requested pytest command: 639 passed, 0 failed, 0 skipped (6.44 s); owned-file diff whitespace check clean. New regressions first reproduced 158 failures before implementation.
- Unfixed: none in scope. Compatibility: legacy invalid IDs now reject; callers must supply consistent case and distinct replicate assay IDs. Full feng suite not run (task-specific gate overrides standing default).
- Only owned source/test files plus this requested report edited; concurrent unrelated changes left untouched. No commits, submissions, or sealed data.
