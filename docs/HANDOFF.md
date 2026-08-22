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

## 01:15 — Boltz back-prediction FAILED THE GATE (single-sequence mode)

Job 20956826 COMPLETED in 8:33 on mit_normal_gpu. All 4 structures predicted.
Pipeline works end to end. The SCIENCE does not pass.

  PDB   conf    ptm   protein_RMSD(A)
  6E9F  0.325  0.274   30.75
  8FTI  0.344  0.306   41.42
  9LZU  0.346  0.290   33.89
  9M8Q  0.320  0.273   33.59

30-41 A protein RMSD = wrong global fold placement, not a near miss. Boltz's own
confidence (0.32-0.35) and pTM (0.27-0.31) agree the predictions are bad.
RNA chains superpose at 0.7-4.8 A but over only 4-27 residues, which is too few
to be meaningful.

CAUSE: `msa: empty`. A ~900-residue Cas13 cannot be folded single-sequence.
This was flagged as a caveat at submission and is now confirmed quantitatively.

DO NOT report A4 structure validation as passing. It has not.

NEXT STEP for A4: real MSAs required.
  - Preferred: `boltz predict --use_msa_server` (ColabFold MMseqs2). Requires
    outbound internet FROM THE COMPUTE NODE. Login node has internet (curl to
    RCSB/NCBI works); compute nodes are UNTESTED. Test before relying on it.
  - If compute nodes are air-gapped: generate MSAs on the login node, write
    them alongside the YAML inputs, and run prediction offline on GPU.
  - Do not use --override on any rerun; keep completed outputs for comparison.
Re-run the RMSD comparison after MSAs and only then judge the gate.

## 01:35 — A4 GATE PASSES with real MSAs (job 20957171, 4:41 on mit_normal_gpu)

Compute nodes DO have outbound internet (nettest: github 200, colabfold API
reachable; its "/" returns 404 but endpoints are live). So --use_msa_server
works from compute nodes. No login-node MSA pipeline needed.

  PDB   family     single-seq  with-MSA   conf   ptm    nres
  6E9F  Cas13d       30.75 A     4.63 A  0.855  0.944   864
  9LZU  Cas13d       33.89 A     5.71 A  0.715  0.901   853
  9M8Q  Cas13d       33.59 A     5.33 A  0.678  0.861   904
  8FTI  Cas13bt3     41.42 A    23.82 A  0.710  0.693   737

VERDICT: gate passes for the Cas13d subfamily (6E9F, 9LZU, 9M8Q) at 4.6-5.7 A
over 850-900 residues. That is FOLD-LEVEL accuracy, not atomic. dCas13x/CasRx
is Cas13d, so this is the relevant family. 8FTI is Cas13bt3, a different
subfamily, and out of scope.

WHY 5 A IS ACCEPTABLE HERE, stated explicitly so it is not overclaimed:
the steric shadow uses a ~40 A probe (U1 snRNP sized). A 5 A coordinate error
is ~8x smaller than the probe radius, so excluded-volume estimates are not
sensitive to it. This argument does NOT license using these models for
anything requiring atomic detail (contacts, H-bonds, base-specific readout).

USAGE RULE: compute shadows from DEPOSITED structures wherever they exist.
Boltz is only needed for the novel dCas13x-gRNA-SMN2 complex, and the
back-prediction above is what justifies trusting it there.

Outputs: .../boltz_out (single-seq, keep for comparison), .../boltz_out_msa.
Rerun recipe: slurm/boltz_msa.sbatch <in_dir> <out_dir>, inputs must NOT
contain `msa: empty`.

## 04:30 — Opus code review COMPLETE, all checklist items pass

  [1] c.840C>T at 32004 = T                      PASS
      Independently confirmed: intron 7 opens GTAAGT (canonical 5'ss),
      ISS_N1 = CCAGCATTATGAAAG, GC_CORE = CTGCCAGC. Coordinates sound.
  [2] hard constraint reaches the ensemble       PASS
      GC_CORE total pairing prob 1.4967 -> 0.0 when blocked; global bpp
      matrix L1 change 20.73. The constraint is real, not silently dropped.
  [3] dG_unfold positive and monotonic           PASS  0.91 -> 4.29 kcal/mol, L=4..20
  [4] link.py refuses to run unfit               PASS  no defaults, fit provenance required
  [5] single-shard smoke run                     PASS  97 s / 3 footprints (~32 s each)
  [6] shard partition audit                      PASS  1065 total = union = unique, 10-11 per shard

## A3 SWEEP ARRAY SUBMITTED — job 20961463, 400 shards, mit_normal
1065 footprints (starts 32040..32110, lengths 6..20) x 5-span ladder x 2 PF
= 10650 partition functions, ~9.5 core-hours. --requeue + --open-mode=append;
each shard skips if its output exists and writes via .tmp rename, so preemption
costs at most one shard. Output: /home/dashii/orcd/pool/smt-isef-data/sweep

GRADING (do only after the array completes, and do it honestly):
the published answer -- an 8-mer over the GC-rich core (intron 7 +7..+14) --
is NOT used in scoring. Columns overlaps_gc_core / overlaps_iss_n1 are
annotations for grading only. Gate passes if top-ranked footprints by delta_A
concentrate on the GC core at short lengths, across the span ladder.

## KNOWN LIMITATION found during smoke run — do not ignore when grading
Footprint (32040, len 6) sits INSIDE exon 7 and scores fold_change 1218x.
Structurally sensible: it blocks one arm of TSL2, freeing the other. But a real
ASO there would occlude the exon 7 SRSF1 ESE and likely SUPPRESS inclusion.
The accessibility model cannot see ESE disruption -- it measures structure only.
Exonic footprints will therefore be systematically over-ranked. This is the
same positional-polarity problem A5 must solve (intronic = silencer, exonic =
enhancer). When grading A3, either restrict to intronic footprints or report
exonic hits separately with this caveat stated.

## GitHub
https://github.com/Dakadada/smt-isef (private). Remote `origin`, branch main.
