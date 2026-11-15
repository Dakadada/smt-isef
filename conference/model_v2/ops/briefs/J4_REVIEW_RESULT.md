DO-NOT-ADVANCE — 1 HIGH, 1 MED; source unchanged.
Scope: gpu_perturbation_atlas.py, tests/test_gpu_perturbation_atlas.py, jobs/j4_sentinel.sbatch against SCIENTIFIC_LOCK.md and ops/briefs/J4.md.
Paths below are relative to conference/model_v2/.

HIGH — gpu_perturbation_atlas.py:502,508; gpu_perturbation_atlas.py:257 — Reference assembly is asserted, never authenticated.
Scenario: GRCh38 config coordinates plus an accidentally selected GRCh37/custom --fasta pass the config-only fence; parquet unconditionally declares GRCh38 / GENCODE v50. Smoke bypasses even the config declaration.
Reproduced with a synthetic non-reference FASTA and stubbed CLI: exit 0, status complete, metadata GRCh38; this was infrastructure verification, not model evidence.
Fix: require a verified GRCh38 reference identity for CLI FASTA overrides and smoke, bind it to manifest/metadata, and reject mismatches before constructing Pangolin; test incorrect-reference rejection.

MED — jobs/j4_sentinel.sbatch:37,38; gpu_perturbation_atlas.py:31,507 — Startup SIGTERM can prevent requeue.
Scenario: SIGTERM while numpy/pyarrow imports run precedes Python's first handler; forwarded TERM kills the payload, wait fails, and the sole requeue owner exits 1 without scontrol requeue.
Reproduced through the actual module with a delayed numpy import and mocked module/scontrol commands: shell exit 1, zero requeues, "payload checkpoint flush failed".
Fix: establish a startup/readiness protocol that permits safe requeue before checkpoint handling is ready while preserving failure refusal after readiness; extend the shell test beyond its already-installed fake handler.

Checked: exon-7 readout 70076574, five 18-base smoke masks including F18, minus-strand mask/index reversal, masked-minus-reference sign, weight-group/channel order (1/4/7/10), cache identities, failure nulls/evidence labels, OOM record preservation, checkpoint snapshots, and shell-only requeue ownership after handler installation.
Pangolin channel mapping agrees with the installed model.py/pangolin.py and feng scorer; per-weight hashes and model/source hashes are recorded. No further confirmed defect in these checks.
Pytest: 26 passed, 1 skipped, 0 failed (1.21s); specified interpreter, target, and thread limits; bytecode/cache writes disabled for read-only review.
Skip: real CPU/GPU tolerance requires a Slurm GPU allocation. No real inference, GPU use, or Slurm submission; weight loading and device tolerance remain unverified.
Unfixed: both findings, as requested for this read-only review. Scope/test command follow the task's conference-specific override of AGENTS.md's feng defaults; no feng mutation.
