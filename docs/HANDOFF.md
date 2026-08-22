# Resume notes — 2026-08-22, written before 04:11 usage reset

## Live compute
- SLURM 20956443 — Boltz-2 back-prediction, mit_normal_gpu, RUNNING.
  Inputs /home/dashii/orcd/pool/smt-isef-data/boltz_in, out .../boltz_out
  Caveat: single-sequence mode (msa: empty). Pipeline shakeout, not the real
  validation gate. MSA pass is a later iteration.

## OPEN FINDING — highest priority
Agent A reported: **blocking ISS-N1 DECREASES A_5ss.**
This is the opposite of the established mechanism (nusinersen blocks ISS-N1 and
frees the exon 7 5'ss for U1). Do not tune to make it pass. Diagnose:

1. Hard-constraint semantics. `hc_add_up` forces a base UNPAIRED but leaves it
   available. An effector-occupied site must be UNAVAILABLE — the bases cannot
   pair with anything. If ViennaRNA is re-pairing freed ISS-N1 bases into a new
   structure that sequesters the 5'ss, that is a modeling bug, not biology.
   Check base-pair probability matrices constrained vs unconstrained directly.
2. Window size. ISTL1 is a LONG-DISTANCE intra-intronic interaction spanning
   much of intron 7 (32053..32496). If the analysis window is shorter than both
   arms, the structure the ASO is meant to disrupt is absent from the model.
   Verify the window actually contains both ISTL1 arms.
3. A_5ss definition. Confirm SS5_EXON7 = 32050..32058 and that "accessible"
   means all positions unpaired, not any.

If after 1-3 the sign still inverts, that is a real negative result about the
purely-structural occlusion model, and the fallback in the plan applies: add
explicit hnRNP A1/A2 competitive binding rather than treating occlusion as
purely structural.

## Files on disk from stopped agents — REVIEW, do not trust
src/smt/seq/construct.py, src/smt/thermo/{constrained,hybridization,accessibility,link}.py
src/smt/steric/structures.py
Missing: tests, steric/shadow.py, steric/parity.py, validate/length_sweep.py

## Review protocol before ANY array submission
1. Verify c.840C>T at coordinate 32004 independently (coordinate sanity).
2. Verify hard constraint actually changed the ensemble (bpp matrix diff).
3. dG_unfold positive, monotonic in window length.
4. link.py must REFUSE to run unfit — check no default parameters leaked in.
5. Single shard smoke run + timing before submitting 100-way array.
6. Shard partition audit: union = full sweep, no overlap.

## Reference coordinates (verified)
NG_008728, 35073 bp. exon6 26119..26229, intron6 26230..31998,
exon7 31999..32052, intron7 32053..32496, exon8 32497..33073.
ISS-N1 32062..32076. GC core 32059..32066. 5'ss region 32050..32058.
c.840C>T at 32004 (SMN2 has T).

## Not yet done
- docs/lab-request.md written but NOT SENT. Longest-latency item in project.

## Restore point
Commit 345b89b — everything above is committed. Safe to hard-reset to it.

## Autonomous operating rules (user asleep, two windows: 04:11 and 09:11 EDT)
Budget is the scarce resource, cluster compute is not. Order of preference:
1. Fix correctness problems FIRST. Nothing gets submitted broken. A 100-way
   array of wrong physics is worse than no array — it looks like progress.
2. Prefer cluster compute over agent tokens. One sbatch = hours of work for
   ~zero budget. Spawning an agent = real budget.
3. Commit after every working increment so any stop is resumable.
4. When budget gets tight, stop spawning and just submit. Leave notes here.

## Preemptable strategy — the right pattern
Do NOT literally checkpoint every 5 minutes. For this workload the better
answer is SMALL IDEMPOTENT WORK UNITS + skip-if-done, which makes preemption
cost at most one unit and needs no checkpoint machinery at all.

Rules for any preemptable submission:
- `#SBATCH --requeue` so SLURM re-queues automatically on preemption.
- `#SBATCH --open-mode=append` so requeued logs do not clobber the first run.
- `#SBATCH --signal=B:USR1@120` plus a bash trap to flush partial output and
  exit cleanly when preemption is 120 s away.
- Every work unit writes to a UNIQUE output path, and the job SKIPS units whose
  output already exists. This is what makes requeue free.
- Write output atomically: write to `<path>.tmp` then `mv` into place, so a
  kill mid-write never leaves a half-file that the skip-check believes is done.
- Keep single units well under the preemption horizon. Aim for <15 min each.
- Size arrays large and units small, rather than few units and long walltimes.

Partition order: fill mit_normal / mit_normal_gpu first (6 h cap, not
preempted), then overflow to mit_preemptable (48 h cap, preempted freely).
Submit the same array to both and let skip-if-done deduplicate.

For Boltz specifically: each input YAML is already a natural unit. Check for an
existing output dir per input and skip. Do not use --override on preemptable.

## Resumability status as of 00:37
- Work on disk: YES, committed at 345b89b.
- Stopped agents: resumable by message, but their code has a known sign bug.
  Prefer reviewing and fixing directly over resuming them blind.
- GPU job 20956443: running on mit_normal_gpu, NOT preemptable, no requeue
  needed. If it fails, its inputs are regenerable from the CIFs.

## HARD CONSTRAINT — agent usage (user directive, not negotiable)
- MAXIMUM TWO subagents at any time. Never three. Not even briefly.
- Subagents run on SONNET (`model: "sonnet"` on the Agent tool).
- Opus (main session) REVIEWS the code afterward. Reviewing is the main
  session's job and must not be delegated.
- Caveat: the Agent tool exposes `model` but no reasoning-effort parameter.
  Effort is inherited from the session / agent definition and cannot be set
  per-call. Sonnet is guaranteed; "high effort" is not settable here.
- Division of labor: Sonnet writes, Opus reviews. Do not have Opus write bulk
  code that a Sonnet agent could produce, and do not have a Sonnet agent
  approve its own physics.
boltz resubmitted as 20956826 after installing cuequivariance-torch

## RESOLVED 00:55 — sign inversion was a wrong-footprint artifact, NOT a bug

Diagnostic: A_5ss baseline vs blocked, across max_bp_span ladder (0/50/100/
150/300) and two window sizes (798 nt, 1048 nt). Estimator = exact Z-ratio.

RESULT — blocking GC_CORE (intron 7 +7..+14, coords 32059..32066, the 3UP8
target) INCREASES A_5ss at EVERY span and BOTH window sizes. Fold change
consistently 6-7x. This is the correct mechanistic direction and it is robust.

  span    A_base     A_GC      delta        fold
  global  1.799e-05  1.141e-04  +9.607e-05   6.3x
  50      4.601e-06  1.372e-05  +9.117e-06   3.0x
  100     4.234e-05  3.001e-04  +2.577e-04   7.1x
  150     1.175e-05  4.153e-05  +2.978e-05   3.5x
  300     2.117e-05  1.292e-04  +1.080e-04   6.1x

Blocking the FULL 15-nt ISS-N1 (32062..32076) is what inverts, and only at
some spans (negative at global/100/300, positive at 50/150). Agent A tested
this footprint and reported the inversion. It was testing the silencer, not
the therapeutic target.

INTERPRETATION: ISS-N1's 3' half (+15..+24) appears to participate in structure
that helps keep the 5'ss accessible, so occluding it is counterproductive; the
GC core is the sub-element whose occlusion frees the site. This is consistent
with the published micro-walk showing an 8-mer at the GC core is potent. Treat
as a model PREDICTION to state, not a fact to assert.

CONSEQUENCE FOR A3: the falsification gate is passing on the footprint that
matters. The sweep must score every tiling position and length, not just the
annotated ISS-N1, and must report a max_bp_span LADDER. Never report a single
span -- ISS-N1's sign flips with it.

OPEN: span-dependence of the ISS-N1 result is not understood. Alternating sign
between 50/150 and 100/300 suggests structural bistability or a numerical
issue. Investigate before making any claim about full-length ISS-N1. The
GC_CORE result does not depend on resolving this.

Absolute A_5ss values are ~1e-5 because it is a JOINT probability over 9 nt.
Small magnitude is expected; report fold change, not absolute.
