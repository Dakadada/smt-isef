J10-OFFTARGET — implementation ready; production summaries NOT run; jobs NOT submitted.
j10_offtarget_profile.py:131 — independent sequence/duplex/accessibility quantities; F18/F14/F10, mismatch bins 0..3.
j10_offtarget_profile.py:28 — fixed energy <= {-5,-10,-15,-20} kcal/mol and accessibility >= {.01,.1,.5}; no outcome fitting.
j10_offtarget_profile.py:89 — exact SMN2/SMN1 anchors and other SMN gene rows separately retained in JSON.
j10_offtarget_profile.py:152 — identical gene-specific site IDs collapse; conflicting duplicates fail; different genes retained.
j10_offtarget_profile.py:245 — hash-owned streaming partials; source SHA/rows/provenance validated; atomic parquet/JSON plus SHA256 .done.
j10_offtarget_profile.py:276 — merge requires all 256 partials, checksums/common provenance, and 11,002,701 raw rows.
Fractions: accessibility_pass/accessibility_evaluated; None for empty denominator; per-channel status categories retain failures/missing.
Labels: derived / model_assumption (unmodified RNA, chemistry not modelled) / proxy; no occupancy summation or collapsed score.
J3 window sensitivity: only W=120/L=60/flank=80 configuration; context clipped at boundaries; no invented variants.
J6 gene-only region may be not_available; gene biotype unavailable in J6 schema; reported explicitly, no external joins.
Sizing estimate: 1–5 min and <1 GiB per ~42,979-row shard; request 20 min/2 GiB/one CPU, array 0-255%64, mit_normal, requeue.
Basis: recorded J3 10,837-row size test 642.72 s / 1,231,396 KiB RSS / 9,814,128 filesystem-input blocks; includes J6 scan + computation, NOT isolated parquet read timing.
Targeted pytest: 17 passed, 0 skipped, 0 failed; model_v2: 771 passed, 1 skipped, 0 failed.
Synthetic oracle: 5 input / 4 unique / 1 duplicate; duplex 3 evaluated / 2 passing; accessibility 3 evaluated / 3 passing; joint 1.
Submit (NOT executed): J10_ARRAY_JOB=$(sbatch --parsable conference/model_v2/jobs/j10_array.sbatch)
Submit (NOT executed): sbatch --dependency=afterok:${J10_ARRAY_JOB%%;*} conference/model_v2/jobs/j10_merge.sbatch
Output: ~/orcd/scratch/feng/model_v2/j10_offtarget_v1/j10_profile.{parquet,json,done}; SMN site rows are in JSON.
Ownership conflict: requested ops/briefs/J10_RESULT.md is outside owner-only paths; this report stays in allowed jobs/j10_* instead.
feng gate: 581 passed, 5 skipped, 1 failed (230.89 s): test_cleanroom.py:23 rejects existing src/ text in feng/src/feng/source_auth.py; outside ownership, unfixed; global gate DO-NOT-ADVANCE.
