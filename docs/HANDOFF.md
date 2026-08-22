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
