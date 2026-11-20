J10-OFFTARGET-FIX: fixed; no jobs submitted or commits made.
conference/model_v2/j10_offtarget_profile.py:213 — load holds the writer's partial lock across marker/hash/JSON/Parquet reads; merge uses this locked load; partial reuses its existing lock via _load_locked (avoids nested-lock deadlock).
conference/model_v2/tests/test_j10_offtarget_profile.py:90 — deterministic generation replacement regression covers load and merge; both failed before the fix, both pass afterward; all 17 existing tests retained.
conference/model_v2/jobs/j10_array.sbatch:20; conference/model_v2/jobs/j10_merge.sbatch:16 — explicit new output ~/orcd/scratch/feng/model_v2/j10_offtarget_v2; array 0-255%64, mit_normal, 20 min, 2G retained for both jobs.
Pytest (requested interpreter/environment/command): 19 passed, 0 skipped, 0 failed (1.17s); bash -n and git diff --check passed.
From /home/dashii/smt-isef, submit array: array_job_id=$(sbatch --parsable conference/model_v2/jobs/j10_array.sbatch)
Then submit merge: sbatch --dependency="afterok:${array_job_id%%;*}" conference/model_v2/jobs/j10_merge.sbatch
Unfixed: none in requested scope. Limitation: v1 output from 29254991 has stale code identity; regenerate J10 partials and merge in v2 before downstream use; no live output validation performed.
