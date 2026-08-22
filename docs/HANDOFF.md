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

## A3 GATE PASSES — job 20962114, 200/200 shards, 0 failed

795 footprints x 5-span ladder = 3975 rows. Graded on the 780 strictly
intronic footprints (start > 32058); ranked by MEDIAN delta_A across the span
ladder, which is robust to the ISS-N1 sign flips noted earlier.

TOP INTRONIC FOOTPRINTS (median delta_A):
  start  len   end     medDelta   GC  ISS
  32059    9  32067   9.629e-05   1    1
  32059    8  32066   9.607e-05   1    1   <-- EXACTLY GC_CORE (intron7 +7..+14)
  32059    7  32065   9.320e-05   1    1
  32060    8  32067   8.715e-05   1    1
  32061    6  32066   8.447e-05   1    1
  32090    8  32097   3.873e-05   0    0   <-- secondary site, NOT annotated

GRADE: 18/20 of top-20 intronic overlap the GC core, against a base rate of
0.154 across all intronic footprints. Top-20 lengths are 6-13, concentrated at
6-9. The model recovered BOTH the position and the short-length regime of the
published micro-walk result, blind -- overlaps_gc_core is an annotation column
and never enters scoring.

This is the go/no-go from the plan. It passes. Downstream workstreams (A4
parity, A5 sieve, A7 cryptic exons, A8 kinetics) are now justified.

NEW PREDICTION worth following up: footprint (32090, len 8) ranks 11th with no
GC-core or ISS-N1 overlap. Either a real secondary accessibility site or an
artifact of the compact window. Check it before claiming it.

CAVEAT REPEATED: exonic footprints were excluded from grading because the
accessibility model cannot see ESE disruption and systematically over-ranks
them. Any published version must state this.

## NEXT
A5 sieve (transcriptome scan + positional polarity), then GEO validation
against the published high-concentration ISS-N1 ASO RNA-seq (PMC8395096).
Still unwritten: steric/shadow.py, steric/parity.py, sieve/*, kinetics/*.

## A4 STERIC SHADOW + PARITY LANE COMPLETE (commit 4ec5341)

structures.py and shadow.py were already on disk (agent-drafted, marked
UNREVIEWED in the 345b89b commit). Reviewed both line by line: classification
is residue-composition + base-pairing-geometry driven, not chain-ID
special-cased, and correctly handles both deposited architectures (separate
guide/target chains in 6E9F/9LZU/9M8Q; 8FTI's covalently fused
crRNA-linker-target chain, split at the spacer-duplex apex). All four
structures classify with guide+target confidently assigned, zero AMBIGUOUS
flags. One depositor-metadata oddity worth flagging: 8FTI's protein entity is
literally labelled "Integrase" in the mmCIF (`_entity.pdbx_description`),
carrying catalytic-dead mutations R84A/H89A/R739A/H744A -- a Cas13bt3
annotation quirk in the deposited file, not a bug in our classifier (which
never reads entity names to decide role).

Chain classification (guide / target / protein):
  6E9F  protein A (864 aa) | guide B (52 nt) | target C (27 nt, footprint 7..26, 20 bp duplex)
  8FTI  protein A (737 aa) | fused chain B (97 nt): guide=3' arm (66 nt) target=5' arm (-67..-37, 30 bp duplex, footprint -66..-37)
  9LZU  protein A (853 aa) | guide C (51 nt) | target B (14 nt, ENTIRELY hybridized, footprint == full target, 13 bp duplex)
  9M8Q  protein A (904 aa) | guide B (55 nt) | target C (25 nt of a 30-mer construct, footprint 3..26, 17 bp duplex)

U1 snRNP probe radius: R = (3M/(4*pi*rho*N_A))^(1/3), M=250 kDa,
rho=1.35 g/cm^3 (large-protein partial-specific-volume asymptote, Fischer/
Polikarpov/Craievich 2004) -> R = 41.87 A (>= the "roughly 40 A" sanity
target in the brief). Flagged as an upper bound: U1 snRNP is ~1/3 RNA
(denser, ~1.7 g/cm^3) and is not spherical, so this likely over-estimates
slightly.

RNP vs bare-duplex shadow width at the U1 probe (threshold s(i)>=0.5,
n_points=300):
  PDB    rnp_occluded/target   duplex_occluded/target   wider by
  6E9F        26 / 27                18 / 27              +8 nt
  8FTI        31 / 31                19 / 31             +12 nt
  9LZU        14 / 14                14 / 14               0 nt  (see below)
  9M8Q        25 / 25                19 / 25              +6 nt

Central claim (RNP shadow strictly wider than the bare guide:target duplex
at large probe) HOLDS on 6E9F, 8FTI, 9M8Q. It does NOT distinguish on 9LZU,
but this is a data-availability artifact, not a physics failure: 9LZU's
crystallized target is exactly 14 nt and the hybridized footprint IS the
entire resolved target (no flanking density at all survived deposition), so
neither state has any flanking nucleotide on which to show an asymmetric
extension -- both saturate to 14/14 because there is nothing left to
saturate past. tests/test_steric.py encodes this honestly: the weak claim
(RNP >= duplex) is asserted on all four; the strict widening claim is
skipped on 9LZU with the reason stated, not loosened.

ASO length reaching parity with the RNP (occluded-nt count, U1 probe, swept
8-20 nt, capped by each structure's resolved duplex length):
  6E9F: L=18 (best available; still short of parity, gap 18 nt)
  8FTI: L=20 (best available; still short of parity, gap 23 nt)
  9LZU: L=14 (max resolvable; gap 6 nt)
  9M8Q: L=20 (best available; still short of parity, gap 16 nt)
No ASO length in the 8-20 nt therapeutic window reaches full parity with the
RNP shadow on any structure -- the sweep always bottoms out at the longest
testable length, i.e. even a maximal-length ASO leaves a measurable gap.
That gap, not a crossover point, is the headline number.

Unexpected finding, NOT tuned away (see test_steric.py module docstring for
the full derivation): per-nucleotide shadow s(i) is an area-weighted average
over a residue's atoms, and while each individual ATOM's own accessible
fraction is provably monotonic non-increasing in probe radius (proved from
the exact hard-sphere legality condition), the residue-level weighted
average is NOT guaranteed monotonic -- weights (each atom's own free-state
exposure) shift with radius, a Simpson's-paradox-style reversal. Confirmed
at up to 6000 points/atom (not quadrature noise): individual nucleotides in
8FTI (both states) and the bare-duplex state generally (fewer atoms, less
averaging) show local dips up to ~30% swept nucleotide-by-nucleotide across
the probe ladder. Averaged over the WHOLE target strand this cancels and the
mean is cleanly monotonic everywhere checked -- that mean is what
test_steric.py asserts monotonicity on, with the per-nucleotide caveat
stated rather than hidden.

Also unwritten still: sieve/*, kinetics/*, wiring shadow.py into
constrained.py's occlusion footprints (not requested by A4's scope, noted
for whoever picks up A5/A6).

## A5 TRANSCRIPTOME-WIDE SIEVE COMPLETE (sieve/transcriptome.py, scan.py, polarity.py)

Data: GENCODE human release 44 (GRCh38.p14), "basic" GTF (29,570,410 bytes)
+ primary_assembly genome FASTA (844,691,642 bytes gz / 3,151,417,447 bytes
decompressed), both MD5-verified against GENCODE's own MD5SUMS on fetch.
Under /home/dashii/orcd/pool/smt-isef-data/transcriptome/. Indexed to
~20,035 canonical (Ensembl_canonical), protein_coding transcripts --
one/gene, on canonical chromosomes only (see transcriptome.py module
docstring, "Scope decisions", for why: full comprehensive+all-isoform+
full-intron coverage is not tractable here). Introns represented only within
200 nt of each splice site (donor+acceptor flanks), not in full -- a
disclosed tractability assumption, not a measurement.

SANITY CHECK the index gives for free: SMN2's canonical transcript
(ENST00000380743.9) has 9 exons, matching this project's own
SMN2_EXONS numbering (1,2a,2b,3,4,5,6,7,8) exactly, and its 8th
transcript-order intron (444 nt) is byte-for-byte the same intron 7 already
pinned in smt.seq.reference/construct (same length, ISS-N1 and GC_CORE
sequences fall at the same offsets). Independent confirmation the GENCODE
coordinates and this project's RefSeqGene-based coordinates agree.

Scan (nusinersen, TCACTTTCATAATGCTGG): stage-1 seed+extend (tiling 9-mers,
<=3 mismatches) over all 20,035 transcripts took 2.9 s, 37 candidates in 34
genes. Stage 2 (IntaRNA E_hybrid + RNAplfold dG_unfold) on all 37 took 15 s
total. RESULT: SMN1 and SMN2 tie for RANK 1 by ddG_net (-22.85 kcal/mol
each, 0 mismatches, exactly ISS-N1, intron 7 +10..+27) -- the sieve recovers
its own intended on-target site as the single best transcriptome-wide hit,
found blind (the scan has no special-casing of SMN1/SMN2). Polarity module
gives SMN1/SMN2 the single most POSITIVE signed_effect in the whole run
(+16.37, correctly predicting INCLUSION, matching nusinersen's real
mechanism) -- also the top-ranked score overall, not just top of its
category.

POLR2H does NOT appear among the 37 stage-1 candidates. Diagnostic
brute-force (no seed filter) found POLR2H's canonical mRNA has an 18-nt
window at 4/18 mismatches (78% identity, genomic chr3:184364973-184364990);
scored via the same stage-2 pipeline it gets ddG_net=-17.11, which would
rank ~9th of 38 if included. WHY THE SEED STAGE MISSED IT (real, disclosed
limitation, not a bug): the seed design's documented pigeonhole guarantee
only covers <=1 mismatch; this candidate's 4 mismatches are spread so that
NEITHER 9-mer tile is exact (nor would a 7-mer tile be -- the longest exact
run in this alignment is 5 nt, verified by brute force). Also found: this
diagnostic window sits in EXON 3 of GENCODE's canonical POLR2H transcript,
not "exon 2" -- an open discrepancy with the literature's exon numbering
(possibly a different reference transcript/isoform; not resolved here).

GEO validation (full writeup: docs/geo-validation.md): PMC8395096 (Ottesen
et al. 2021, Int J Mol Sci 22:8378, PMID 34445083) has NO GEO/SRA accession
-- checked directly against the article's own Data Availability Statement,
which states all data is in the article's supplementary files only (Tables
S1-S17, itemized in the doc). POLR2H is not mentioned anywhere in that
paper's text. The POLR2H exon-2/F18MOE claim traces instead to a newer,
related GEO series, GSE307431 (deposited 2026-01-06, same lab, no PMID
linked yet), whose abstract explicitly describes POLR2H exon 2 skipping via
direct exonic annealing of an 18-mer MOE ASO. Comparison protocol (rank
correlation + sign accuracy against PMC8395096's Tables S3-S6, POLR2H case
study against GSE307431 once downloaded) is written up in
docs/geo-validation.md but NOT executed -- no RNA-seq or supplementary PSI
data was downloaded or used in polarity.py's scoring, per the "no
experimental PSI data, no placeholder calibration" rule.

tests/test_sieve.py: 24 tests, all passing, self-contained (synthetic
Transcript objects, no dependency on the 3 GB downloaded index) --
classification correctness (UTR/exonic/intronic/proximal, both strands),
seed-stage recovery of a planted perfect match + the pigeonhole 1-mismatch
guarantee, ddG_net = dG_hybrid + dG_unfold via the existing thermo modules
(no reimplemented energies), the exonic/intronic sign flip at matched
distance and ddG_net, and an exact shard-partition audit. Full suite: 101
passed, 1 skipped (pre-existing, unrelated).

slurm/sieve_array.sbatch written (200-way array, mit_general/mit_normal,
--requeue + --open-mode=append, skip-if-exists + atomic .tmp rename per
shard) but NOT submitted -- the direct single-process run above already
covers the real candidate volume (37 total) faster than array overhead
would; the array script is provided for reproducibility / larger future
oligos or relaxed seed thresholds where candidate volume could grow.

Environment oddity worth flagging: this session observed docs/HANDOFF.md and
several src/smt files (steric/shadow.py, steric/parity.py, the A4 sections
above) grow on disk during this session without this agent's own edits,
apparently from a concurrent process. Not investigated further -- outside
this task's scope -- but noted here since "why did HANDOFF.md change
underneath me" is a fair question for whoever reads this next.

## NEXT
kinetics/* is the only unwritten lane left per the original plan. A5's own
loose ends: (1) resolve the POLR2H exon-2-vs-exon-3 numbering discrepancy
against the specific transcript GSE307431/the literature actually used; (2)
decide whether to add a smaller/approximate seed pass specifically to widen
sensitivity beyond the 1-mismatch guarantee (the POLR2H diagnostic shows the
default 7-10 nt seed range, used as specified, has a real blind spot); (3)
execute the GEO comparison protocol in docs/geo-validation.md once
PMC8395096's supplementary tables and/or GSE307431's processed data are
actually downloaded.

## CORRECTION 05:20 — the sieve validation dataset was misidentified

Earlier plan text and HANDOFF entries claimed PMC8395096 (Ottesen et al. 2021,
Int J Mol Sci 22:8378) provides transcriptome-wide RNA-seq for calibrating A5,
"available today". BOTH claims were wrong:
  - PMC8395096 has NO GEO/SRA deposit. Its Data Availability statement puts the
    data in Supplementary Tables S1-S17 only.
  - PMC8395096 never mentions POLR2H anywhere.

The POLR2H exon-2 off-target and the usable RNA-seq come from GSE307431
(same lab, deposited 2026-01-06, no linked PMID yet). VERIFIED via NCBI eutils:
20 samples, "Expression profiling by high throughput sequencing", titled
"A therapeutic antisense oligonucleotide encompassing 2'-O-methoxyethyl
modification triggers unique perturbation of the transcriptome". Its abstract
states skipping of POLR2H exon 2 is triggered by direct annealing of F18MOE.

Anywhere the plan or docs still cite PMC8395096 as the A5 calibration set,
that is stale. Use GSE307431.

## A5 BUILT — GENCODE v44, 20035 transcripts, but the seed stage MISSES POLR2H

Built: sieve/transcriptome.py, scan.py, polarity.py, tests/test_sieve.py,
slurm/sieve_array.sbatch, docs/geo-validation.md. Suite 101 passed, 1 skipped.

Transcriptome: GENCODE release 44 (GRCh38.p14) basic, MD5-verified against
GENCODE MD5SUMS. 3.0 G in .../transcriptome/. 20,035 canonical protein-coding
transcripts, one per gene via Ensembl_canonical tag -- flagged in code as a
tractability ASSUMPTION, not a measurement. Introns stored as 200 nt flanks
around splice sites only -- also flagged as an assumption.
Index sanity check that passed: SMN2 canonical has 9 exons, and its 8th intron
(444 nt) matches this project's independently pinned intron 7 length and
ISS-N1/GC_CORE positions exactly.

STRONG RESULT: scanning nusinersen's 18-mer (seed_k=9, <=3 mismatches) over all
20,035 transcripts yields 37 candidates across 34 genes in 2.9 s. SMN1/SMN2
ranks #1, ddG_net = -22.85 kcal/mol, zero mismatches, on ISS-N1, and polarity
assigns signed_effect = +16.37 -- the highest positive value in the entire
table, meaning "promotes inclusion", matching nusinersen's real mechanism.
Produced blind; there is no SMN special-case in the code.

PROBLEM -- DO NOT CALL A5 VALIDATED YET: POLR2H is NOT among the 37. Its 4
mismatches are distributed across both 9-mer tiles, so no 9-mer (nor even any
7-mer) matches exactly. The two-tile seed design only guarantees detection at
<=1 mismatch. Brute-force search without seed filtering finds the site at
ddG_net = -17.11, which would rank ~9/38. A safety sieve that misses the
field's canonical off-target does not yet function as a safety sieve.

FIX REQUIRED before the sieve array is worth submitting: increase seed
sensitivity. Options are spaced/gapped seeds, more tiles with smaller k, or a
pigeonhole scheme sized for the intended mismatch budget (for 18 nt and <=4
mismatches, 5 tiles of ~3-4 nt guarantees one exact tile). Whatever is chosen,
state the guaranteed mismatch bound in the docstring, and re-run to confirm
POLR2H appears.

UNRESOLVED: the POLR2H window falls in exon 3 of the GENCODE canonical
transcript, not exon 2 as the literature says. Could be transcript-isoform
numbering, or a real discrepancy. Resolve before publishing anything about it.

NOTE: a previous agent observed files changing under it -- that was the other
agent working in the same repo concurrently, not corruption.

## 04:20 — plan file corrected, A8 dispatched
Plan at ~/.claude/plans/jiggly-hugging-biscuit.md now cites GSE307431 for A5
validation. PMC8395096 appears there only inside the correction note.
Sieve array 20963760 submitted, PENDING in queue. A8 kinetics agent dispatched.

## A5 SIEVE COMPLETE — job 20964129, 200/200 shards, 0 failed

741 scored sites across 705 genes. Mismatch distribution: 0mm=2, 2mm=7,
3mm=90, 4mm=642.

RANKING by ddG_net (stronger = more negative):
  SMN2  rank 1/741   ddG -22.85  intron:8:donor  0 mm
  SMN1  rank 1/741   ddG -22.85  intron:8:donor  0 mm   (tied, paralogues)
  POLR2H rank 34/741 ddG -17.11  mrna            4 mm
SMN1/SMN2 are the ONLY zero-mismatch sites in 140.6 Mnt. Nearest competitor is
USP28 at -21.69 with 3 mismatches in mRNA.

POSITIONAL POLARITY — the sieve's central claim, validated on both known cases
blind (no gene is special-cased anywhere in the code):
  SMN2   category=intronic_proximal  signed_effect=+16.37 -> predicts INCLUSION
  POLR2H category=exonic_proximal    signed_effect=-10.38 -> predicts SKIPPING
Both match the literature: nusinersen restores SMN2 exon 7 inclusion, and
induces POLR2H exon skipping. Same oligo, opposite sign, decided by position.
SMN2 is the single most positive of all 741 sites.

BE PRECISE WHEN REPORTING THIS: POLR2H's SIGN is correct, but its magnitude is
mid-pack -- the most negative site is GGA1 at -14.39, not POLR2H at -10.38. The
validated claim is the DIRECTION of the polarity call, not a rank-1 skipping
prediction. Do not overstate it.

Category counts: intronic_core 242, three_utr 187, exonic_core 168,
exonic_proximal 92, intronic_proximal 42, five_utr 10.
Merged output: /home/dashii/orcd/pool/smt-isef-data/sieve/annotated.csv

## TWO THRESHOLD BUGS CAUGHT IN THIS LANE — both silent, neither raised an error
1. Seed heuristic guaranteed only <=1 mismatch; POLR2H (4 mm, split across both
   tiles) was invisible. Fixed by exhaustive scanning (~14 s for 140.6 Mnt).
2. DEFAULT_MAX_MISMATCHES was 3; POLR2H needs 4. First array returned a clean
   99-row table with POLR2H simply absent. Fixed, budget now 4.
Both would have produced publishable-looking output that quietly omitted the
one off-target the sieve exists to find. Watch for this failure mode.

## NEXT: GEO validation against GSE307431
20 samples RNA-seq. Compare predicted signed_effect against measured splicing
changes: rank correlation plus sign accuracy, POLR2H held out as case study.
Not started. Downloading and aligning 20 RNA-seq samples is the largest compute
job remaining -- size it before submitting.

## A8 COMPLETE — C2B34 kinetics, three-hypothesis discrimination

The writing agent was killed mid-task by a session limit, AFTER writing all
code. src/smt/kinetics/{alu,model,hypotheses}.py and tests/test_kinetics.py
all survived; 22 tests pass. Opus wrote docs/c2b34-discriminating-experiment.md
and committed. Nothing was lost.

INVERTED REPEATS, real search on SMN2 (not asserted -- measured, with a null
control). Top pairs flanking the C2B34 region:
  rank  upstream       downstream     E(kcal/mol)  z_null  fwd  revcomp  inverted
  0     9465-9744      23240-23516    -476.9       -52.6   18   209      True
  1     16205-16491    23225-23516    -431.5       -33.9   14   184      True
  2     13668-13944    23252-23523    -410.2       -35.7   15   170      True
Upstream hits sit in intron 1 (5245..18892), downstream in intron 4
(22925..24712) -- exactly the introns flanking exons 2a-4. The revcomp-vs-
forward asymmetry (209 vs 18) is the inverted architecture backsplicing needs.
z_null = -52.6 means this is far outside the null distribution, not just a
large raw energy.

SEPARABILITY POWER (6 doses x 8 timepoints, 200 reps):
  sigma_rel   H1_feedback  H2_direct  H3_stability
  0.02        1.000        1.000      0.995
  0.05        1.000        0.935      0.935
  0.10        1.000        0.800      0.780
  0.20        1.000        0.620      0.655
  0.40        0.970        0.620      0.580
H1 is trivially separable at ANY realistic precision because of its lag.
H2 vs H3 is the hard pair and needs sigma <= 5% for >90% power -- beyond gel
densitometry, within ddPCR/qPCR range. This is the deliverable: not a claimed
mechanism, but the experiment that decides between three, with its required
precision.

Rate constants are symbolic. illustrative_rates() raises if fields are omitted
and is labelled explicitly as not a measurement, per the no-placeholder rule.

## PHASE A STATUS — all four pillars have results
  A3 gate    PASS  blind recovery of published GC-core 8-mer at 32059..32066
  A4 parity  DONE  no ASO in 8-20 nt reaches steric parity with a dCas13 RNP
  A5 sieve   DONE  SMN2 +16.37 inclusion / POLR2H -10.38 skipping, both blind
  A8 kinetics DONE three-hypothesis discriminator + required precision
REMAINING: GEO validation of A5 against GSE307431. A6 minigene compression and
A7 STMN2/UNC13A retarget are unstarted (both compressible per the plan).

## GEO VALIDATION, STAGE 1 (counts proxy) — partial pass, one known case FAILS

GSE307431 provides ONLY gene-level counts (GSE307431_Counts.csv.gz, 5.8 MB).
No splicing tables. Design: GM03813 SMA patient cells, 20 samples, 3 reps,
comparing untreated / scrambled / ISS-N1-targeting ASOs in OMe, MOE and PMO
chemistries. F18MOE is nusinersen's chemistry; ScrMOE is its scrambled control.

Gene counts can only test skipping that is NMD-coupled. This is a PROXY, and a
weak one. Result, F18MOE vs ScrMOE, CPM-normalised, expressed genes only:

  predicted SKIPPING   median log2FC -0.0417  (n=152)
  predicted INCLUSION  median log2FC -0.0053  (n=191)
  background           median log2FC -0.0015  (n=13374)
  Mann-Whitney skip < background : p = 6.67e-05
  Mann-Whitney skip < inclusion  : p = 3.31e-03

AGGREGATE: the polarity sign carries real information. Genes predicted to be
skipped are significantly down relative to both background and predicted-
inclusion genes.

BUT THE KNOWN CASE FAILS: POLR2H log2FC = +0.150 (slightly UP), against a
predicted-skipping call of -10.38. SMN1 -0.171, SMN2 +0.060.

DO NOT REPORT A5 AS VALIDATED ON POLR2H. It is not. Possible reasons, none
resolved here:
  - POLR2H exon 2 skipping may not be efficiently NMD-coupled, so gene-level
    counts would not move even if the skipping is real
  - the skipped isoform may be stabilised
  - 3 replicates is low power for a single gene at these effect sizes
Effect sizes throughout are tiny (-0.04 log2FC is a ~3% expression change), so
the aggregate result is statistically significant but biologically small. State
it that way; do not inflate it.

## GEO VALIDATION, STAGE 2 (real splicing) — NOT STARTED, sized
Raw reads: SRA study PRJNA1321660, 20 paired runs, ~3.1 GB per sample,
~62 GB total. SCRATCH has ~884 GB free, so storage is not a constraint.
Needs STAR (or similar) + rMATS to get per-exon PSI; neither is installed.
This is the only way to test POLR2H exon 2 directly, and the only way to
validate the sieve on SPLICING rather than on an NMD proxy.
Run accessions are in the ENA filereport for PRJNA1321660.

## A7 COMPLETE — cryptic exon retarget, and a PROJECT-WIDE SCOPE BOUNDARY

Suite 138 passed, 1 skipped (was 123/1). Commit dfa8f2d. New files:
src/smt/seq/targets.py, src/smt/thermo/occlusion.py,
src/smt/validate/cryptic_sweep.py, slurm/cryptic_sweep.sbatch,
tests/test_cryptic.py. cryptic_sweep NOT submitted (17,820 footprints x 5 spans).

### THE SCOPE BOUNDARY — read this before claiming anything about ASO therapy
`hc_add_up` forces a base UNPAIRED-AND-AVAILABLE. That correctly represents
"this base can no longer participate in inhibitory structure" -- an INDIRECT,
STRUCTURAL effect where the effector sits at site X and the readout is a
DIFFERENT site Y. That is exactly the SMN2/ISS-N1 case, and it is why the A3
gate passed.

It does NOT represent "this base is physically occupied by an annealed ASO".
Forcing an occupied site unpaired makes it read as trivially ACCESSIBLE --
backwards for occlusion. So a footprint that overlaps its own readout window is
REFUSED (a_5ss already raises; cryptic_sweep flags direct_overlap=1), not faked.

CONSEQUENCE, state it openly: real cryptic-exon ASO therapy frequently works by
DIRECTLY occupying the cryptic splice site. This framework cannot model that
modality. It models indirect structural modulation only. This limitation
applies to the WHOLE project, not just A7, and it is the honest answer to the
concern raised early on that hc_add_up "frees bases instead of occupying them".
The A3 result stands because SMN2's mechanism genuinely is indirect.

### Coordinates pinned (no RefSeqGene exists for either gene -- verified, 0 hits)
STMN2:  NC_000008.11:79,611,116-79,666,157, len 55042.
        Cryptic exon 2a Interval(5707,5933) = chr8:79,616,822-79,617,048
        (Melamed 2019, PMC6348009). HIGH confidence: 227 nt matches paper, AG
        acceptor confirmed at 5705-5706, AUUAAA PAS at 5910.
        TDP-43 GU-repeat site Interval(5759,5786) is NOT literature-pinned --
        derived by GT-repeat scan, flagged UNVERIFIED-PRECISE in code.
UNC13A: NC_000019.10:17,601,335-17,688,353, MINUS strand, len 87019,
        NM_001080421.3 (44 exons, matches NCBI Gene).
        CE-128 Interval(45813,45940), CE-178 Interval(45763,45940), shared
        donor (Brown 2022, PMC8891019). HIGH: lengths match, GT donor and both
        AG acceptors confirmed, sits between canonical exons 20 and 21.
SNPs:   rs12608932 local 46474 T>G, rs12973192 local 45924 G>C. Positions HIGH
        (dbSNP SPDI + direct base check after strand-complementing). WHICH
        allele is "risk" is carried from paper notation only -- MEDIUM, not
        independently re-derived.

### Signed figure of merit
FOM = direction.sign * delta_A; sign=+1 FREE (SMN2), -1 OCCLUDE (cryptic).
Pure sign flip on the same physical quantity, same units. Zero new
thermodynamics -- reuses a_5ss's ss5_feature parametrization unchanged.
Regression-guarded against the SMN2 GC_CORE result.

### Results, including the negative ones
STMN2  best (5759-5778, len 20) over the TDP-43 GU-repeat region.
       median FOM +5.54e-3, desired direction 5/5 spans. Independently lands on
       the same region Baughn et al. 2023 mask therapeutically.
UNC13A best (45987-46006, len 20), ~40 nt downstream of the shared donor.
       median FOM +5.42e-2, desired 4/5 spans.
       NEGATIVE RESULT: footprints immediately flanking the donor/acceptors came
       back robustly WRONG-DIRECTION at 4-5/5 spans. Mirrors the SMN2 lesson
       that the full silencer inverts while a specific sub-element works.
rs12608932  fold ~1.000 at every span -- essentially no structural effect.
       Consistent with a motif-disruption mechanism this purely structural model
       cannot see. NOT evidence the SNP does nothing.
rs12973192  span- and site-dependent, no clean sign. Against the two acceptors
       at span=100 it gives 6.34x and 4.96x INCREASED accessibility, matching
       the clinical risk-increases-inclusion direction -- but NOT consistently
       across the ladder, and the donor readout decreases at 4/5 spans.
       UNRESOLVED. Same class as the ISS-N1 span flip. Do not cherry-pick the
       span that agrees with expectation.
