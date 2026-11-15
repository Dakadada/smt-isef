# AGENTS.md — Codex working agreement for smt-isef

Standing instructions for any Codex session in this repo. Read before acting. The
human partner drives; you implement, investigate, and review. When a task prompt
conflicts with this file, the task prompt wins — but flag the conflict.

## Repo layout

- `feng/` — the active **prospective** arm (Arms A–E + appendix). All current work
  lives here. This is where you operate unless told otherwise.
- `src/`, `pre/` — the older committed program. Separate effort. Do **not** touch it
  while working `feng/`. Its handoff is `docs/HANDOFF.md`.
- Entry point for state: **`feng/HANDOFF.md`** — read its top `LIVE CHECKPOINT` block
  first every session; it holds the current live job/monitor/thread IDs and the
  next action. `feng/PLAN.md` + `feng/TASKBOARD.md` are the frozen design.

## Environment (exact)

- Python interpreter: `/home/dashii/orcd/pool/feng-env/bin/python`
- Source path: run with `PYTHONPATH=feng/src`
- Test suite (the gate): `PYTHONPATH=feng/src /home/dashii/orcd/pool/feng-env/bin/python -m pytest feng/tests/ -q`
- GPU env: `~/orcd/pool/feng-gpu-env` (torch 2.6.0+cu124 / openmm 8.2 / pangolin 1.0.2).
- Binaries in feng-env: `IntaRNA`, `RNAplfold`, `samtools` (pysam absent → samtools
  fallback is the live path). Do not assume `sklearn` — it is **not** installed.
- Scratch: `$SCRATCH=/home/dashii/orcd/scratch`; pool: `/home/dashii/orcd/pool`.

## HPC guard rails (shared, quota-limited node)

- Never scan/walk/`find`/`du`/`ls -R` `/orcd/scratch` or `/orcd/pool`, or traverse
  `/orcd` broadly. Scope every filesystem op to a known subdir with bounded depth.
- Read quota from `~/orcd/.quota` (cached), not `du`. Read `/etc/motd.d/*` for
  maintenance/policy before substantive work.
- Heavy compute goes to Slurm on `mit_preemptable` with `--requeue` + checkpointing;
  do not run long jobs in the session shell.

## Scientific discipline — non-negotiable

1. **Pre-FREEZE.** No `feng-freeze*` git tag is cut. The prereg seal is UNOPENED and
   the sealed Phase-9 test is one-shot and terminal. **Never** cut a freeze tag, lock
   a real seal, set `FENG_PHASE9_UNLOCK`, or read sealed ΔPSI unless the human
   explicitly says so in that session. These are fail-closed by design — keep them so.
2. **Never manufacture a GO.** Do not move a pre-specified criterion, substitute a
   weaker endpoint, or tune until a result appears. A clean NO-GO is a real result;
   report it honestly. If an endpoint is under-powered or a curve incomplete, return
   indeterminate/None — never a verdict.
3. **Never fabricate a score.** Missing/failed computations return `None` or raise —
   never a fabricated `0.0`. Follow the honest 3-state convention already in the
   scorers. Recompute features leakage-free; never copy a label-derived value into a
   feature column.
4. **Coordinates.** Genome is GRCh38 / GENCODE v50. hg19/GRCh37 inputs must be lifted
   (`feng.rbp.liftover`) before any join; assembly fences are hard — do not weaken them.
5. **Freeze-safety of edits.** Additive, TDD, control-only where matching requires it.
   No new ΔPSI/rank/permutation readout outside the guarded Phase-9 path.

## Token / usage economy — you share one limited account

Codex here is driven by Claude Code and runs on a shared account with a **5-hour
rolling usage limit** that has stalled work mid-task. Be economical so more real
work fits per window:

- **Trust the prompt.** Claude hands you the findings, file:lines, and context you
  need. Don't re-explore the repo, re-read HANDOFF end-to-end, or re-derive what the
  prompt already states. Go straight to the named files.
- **Targeted reads.** Open the specific functions/lines you're changing, not whole
  large files or trees. No broad `grep`/`find` sweeps when a path is given.
- **Don't re-run the full suite gratuitously.** Run the pytest gate once when your
  edits are complete (and only the relevant subset while iterating). Reuse a green
  result the prompt already states rather than reproducing it.
- **Batch the work.** Fix all findings in a pass before the test run; don't
  round-trip per finding.
- **Return tight.** Structured `file:line` + one line each + counts. No prose walls,
  no re-quoting the prompt back, no restating unchanged context.

## Workflow

- **TDD.** Write/extend a test for each change; the full suite must pass before you
  call a task done. Report exact pass/skip/fail counts, never "should pass."
- **Verify on disk**, don't assert from intent — check counts, nulls, mtimes, actual
  file contents.
- **Do not commit, `sbatch`, or launch anything** unless the task explicitly says to.
  Leave unrelated untracked files (e.g. Arm C rank/readout artifacts) unstaged.
- **Reviews are adversarial.** Assume defects exist; return a clear `ADVANCE` /
  `DO-NOT-ADVANCE` verdict with `file:line` + concrete failure scenario + fix per
  finding, ranked HIGH/MED/LOW. Do not raise cosmetic nits at a launch gate.
- On any deviation from a brief (missing dep, wrong assumption, a limitation you
  couldn't fix), flag it explicitly in your return summary rather than hiding it.

## Return format

End substantive tasks with: (1) per-finding/change `file:line` + one line each;
(2) anything unfixed + why; (3) exact pytest counts; (4) new limitations a reviewer
or the human should know.
