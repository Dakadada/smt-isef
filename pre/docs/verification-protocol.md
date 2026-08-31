# Verification protocol — what actually catches errors here

Twelve silent failures have been caught in this project. Every one returned a
**value of the correct type**: a float where a float belonged, a count where a
count belonged. Nothing crashed. There is no type error for "reference
coordinates used in construct space" or "two conformations of one molecule
summed", so ordinary code hygiene does not see them.

This file records what has actually worked, measured across all twelve, and
adds the one gate that was never armed.

## Why this domain is unusually hostile

**Coordinates have many plausible wrong readings.** Reference vs construct
space, 0- vs 1-based, strand, 3'->5' vs 5'->3'. The lab supplied ASOs written
3'->5' and a naive comparison against a 5'->3' target produces a confident
number, not an error.

**The biology is full of near-identical entities.** SMN1/SMN2 at 99.9%. Two
ligand copies in an asymmetric unit. Alternate conformations of one residue. An
inverted repeat with one arm antisense by construction. Each carries a
count-once/count-twice failure that looks fine either way.

**The tools return values instead of raising.** ViennaRNA segfaults on a bad
index rather than throwing. IntaRNA absent produces no output, which reads as
"no duplex forms". The aligner tie-breaks silently. Inherited design, absorbed
consequences.

## The only three detectors that have ever worked

Measured across all twelve catches. Nothing else has caught anything, and
"it looked reasonable" has a **0% detection rate**.

### 1. Contradiction with a previously recorded number
- W53 read 39 heavy atoms against its own 32-atom SMILES -> altlocs summed.
- `UGGCAG` reported REJECTED when this repo's own docstring records it binding
  at -9.09 kcal/mol -> IntaRNA missing from PATH.
- POLR2H counts +0.150 (proxy) vs -0.920 (rMATS) -> the proxy was wrong.

**Rule:** any computed quantity that has a previously recorded value gets a
regression assertion, not an eyeball.

### 2. Contradiction with the literature
- G1's "NO GOLDILOCKS POINT" stood committed for a day while the field's own
  papers on this exact drug report the opposite: shorter ASOs carry FEWER
  hybridization off-targets, because longer ones tolerate mismatches. G1 had
  ranked on relative margin; its own table shows absolute off-target affinity
  moving the other way at p = 0.0165.

**This gate was never armed.** It caught a committed conclusion the first time
it was used.

### 3. A deliberately broken input
- A ligand displaced 10 A must FAIL the pose grader.
- A mis-signed workbook must raise.
- A 90%-downstream rule must be checkable against a majority-class base rate.

**Rule, already in force:** a grader that cannot fail is not a gate.

## THE NEW GATE — literature check before commit, not after

> **Every headline conclusion is checked against published work on the same
> system BEFORE it is committed.** Not the method, not the code — the CLAIM.
> If the field has measured the same quantity, find it and compare. If the
> field disagrees, that is either a finding or a bug, and it must be resolved
> before the conclusion enters the record.

Cost: minutes. Yield on first use: one corrected committed conclusion. On
present evidence this has the highest yield per minute of any check run here.

Applies especially to a result phrased as an absence ("no optimum", "no
correlation", "no effect"), because an absence is exactly what a wrong axis or
an underpowered test produces.

## A thirteenth category, new as of 2026-08-25

**A correct number reported on an axis that does not answer the question
asked.** G1's margins are arithmetically right. The conclusion drawn from them
is not, because relative specificity is not off-target burden. Harder to catch
than a bug: nothing is wrong except the framing, so every regression test
passes.

The only defense found so far is gate 2 — someone else measured the same thing
and phrased it differently.

## At authoring time

1. State what the number would look like if the code were wrong. If the answer
   is "the same", stop and add a discriminating check first.
2. Never let absence and failure share a return value. Raise on the
   environment; return `None` only for the measured thing.
3. Any new "no result" path gets a case in `tests/test_no_silent_failures.py`.
4. Report the detectable effect size beside every p-value. "No effect" and "no
   power" are different claims and have been written the same way here.
