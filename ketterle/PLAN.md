# Ketterle — design plan (DRAFT, back burner)

**Status: DRAFT / BACK BURNER (2026-09-07).** Shape agreed with user; to be fleshed out
later. NOT approved for implementation. No code until this is reviewed + sealed.
Owes a Codex/opus-5 adversarial structural review before it becomes plan-of-record
(deferred — Codex out until 2026-09-10).

Standalone in-silico ASO **specificity** arm (Wolfgang Ketterle). Nusinersen / SMN2
exon-7 anchored. Runs on the existing toolset: Pangolin (splice-site usage under steric
N-masking), IntaRNA (duplex ΔG + target accessibility ED), RNAfold/ViennaRNA.

---

## 1. Why this arm exists / what changed

- The original bivalent coincidence / AND-gate idea is **DEAD on SMN2** (see [HANDOFF](HANDOFF.md)):
  structural vise — ISTL1 arms 279 nt apart (coincidence common); TSL2 8-bp stem too weak
  (a single 9-mer opens it, no split-weak pair). Confirmed by feasibility gate.
- Pivot (user, fork 1+2): abandon bivalent-on-SMN2; build a **different in-silico specificity
  mechanism** that yields a POSITIVE, preregisterable result (a documented null is NOT the goal).

## 2. What Ketterle delivers (three nested layers, one mechanism)

- **A — SCORE:** a biophysically-grounded off-target-discrimination score that ranks ASO
  off-target risk better than a naive edit-distance / seed-count baseline.
- **C — LAW:** how off-target burden scales with ASO length/chemistry (the specificity side of
  "shorter-is-better").
- **B — PRINCIPLE:** design rule = minimize off-target burden at fixed on-target potency.

These nest: the score reveals the law, which implies the principle. One mechanism, not three.

## 3. Scope decisions (locked with user)

1. **Method-first, then length.** Build the specificity method as a standalone benchmark
   (provable against a baseline), THEN apply to the length question.
2. **Ketterle's genuinely-new contribution = OFF-TARGET specificity.** The on-target /
   shorter-is-better work already substantially exists in `pre/` (`smt.validate.length_sweep`
   blind ISS-N1 length sweep reproduces "8-mer over the GC core is potent"; `smt.thermo.
   constrained.a_5ss` 5′ss-accessibility scorer). Do NOT rebuild it clean-room.
3. **Two layers, bridge mandatory (not optional):**
   - **Core** — feng clean-room, public data only, predictions sealed first.
   - **Bridge** — separate, walled, REQUIRED; validates against wet-lab / `pre/` reality.

## 4. Core layer (feng clean-room)

**Mechanism = causal gate** (revised after the Holgersen spike — the naive "thermo beats
edit-distance" is pre-empted by Holgersen's own near-null result, so the lever changes):

1. **Hybridization filter** — enumerate near-complement sites (edit-distance/seed baseline,
   as the existing Arm E scan does), keep those with a favorable **IntaRNA duplex ΔG including
   target accessibility (ED)**. Necessary but not sufficient (Holgersen: off-target splicing IS
   hybridization-driven, but binding alone doesn't predict effect).
2. **Splice-consequence discriminator (PRIMARY signal)** — score each surviving site by
   **Pangolin ΔΨ** at that site. Holgersen's decisive gain came from a splicing-effect model,
   NOT binding energy → this is the load-bearing lever.

This is a **defined causal chain** (accessible → stable → splice-consequential), each stage an
ablation — NOT a blind Pangolin+IntaRNA ensemble (satisfies the honesty rule against ensemble
stacking without a defined causal role + locked ablation).

**Benchmark:** Holgersen 2021 (`feng/holgersen-…pdf`) — 81 PS-MOE steric-block oligos 16–20 nt
(= nusinersen class), per-(SBO, exon) labels, LOASO-ready 70:15:15 disjoint splits.
- **SPLICING labels ONLY** (hybridization-driven). **Drop expression labels** (NOT
  hybridization-driven → a complementarity score can't beat baseline on them).
- Label: per-(SBO, exon) pair, exon body ± ≤200 bp intron, edit distance ≤5; positive =
  bootstrap-significant |dPSI| ≥ 0.5 (strict) or ≥ 0.2 (model); negative = ≥15-read coverage,
  no change. Prefer the non-targeted/promiscuous SBO subset (avoids selection inflation).

**Novelty (since the naive version is pre-empted):** (a) a **training-free** biophysics+Pangolin
pipeline that matches/beats Holgersen's *trained* CNN (reproducible, no fitting); (b) the
**length / off-target-vs-length** application + the wet-lab bridge — neither of which Holgersen did.
"Better binding energy" is explicitly NOT the claim.

**Falsifiable core:** the causal-gate score beats the edit-distance/seed-count baseline on
Holgersen splicing labels by a preregistered margin, LOASO. Prereg the score, margin, splits,
and label thresholds before reading.

## 5. Bridge layer (separate, walled, mandatory — may reference `pre/`)

- One-directional: bridge consumes the core's **sealed** predictions + reads on-target PSI.
  **Core never sees PSI.** Blind-grade only (reuse `pre/`'s discipline: ground truth touched
  only at grading — `smt.validate.length_sweep`, `smt.sieve.ottesen`).
- **On-target / shorter-is-better ground truth available NOW:** Singh 2009 (RNA Biology, 3UP8
  8-mer length series) on the classic **pCI-SMN2 minigene** (444-nt intron 7 → ISS-N1 + ISTL1
  intact); Singh 2006 ISS-N1 walk. Reconstructable now from the papers' fragment description
  against the genomic reference (`pre/` already builds this style via `smt.minigene.spec`).
- **Super-minigene = eventual, lab-pending extension.** Ottesen 2024 `SMN2Sup3445` (full
  intron 7 → ISTL1 intact) is the target variant; truncated variants (757T/371T/2BM4) are NOT
  ISTL1-intact. Exact sequence blocked on the lab (Ask 1 / Suppl S2 accession).
- Reuse from `pre/` (deliberate — prior arm by same author): `smt.thermo.constrained.a_5ss`,
  `smt.seq.construct` coords, `smt.minigene.spec`. Must stay walled from: feeding any PSI into
  the scorer / fitting loop.

## 6. Honesty architecture

- Additive `candidate` score defs; all params + the benchmark protocol frozen before reading.
- Causal-chain gate with locked ablations (not an accuracy-driven ensemble).
- Core predictions SEALED before the bridge reads any PSI (predict-then-confront).
- Blind grading in the bridge; no post-hoc threshold picking; freeze discipline intact.

## 7. Honest gaps / risks (stated up front)

1. **Off-target-vs-length LAW (deliverable C) has NO wet-lab validation** — Holgersen is 16–20 nt
   only; super-minigene has no ASO data. Off-target-vs-length stays **in-silico**; only
   *on-target*-vs-length is wet-lab-validated (Singh 2009 / `pre/`). Do not overclaim C.
2. **Headline-hypothesis risk:** Holgersen tested ΔG vs edit-distance and found no clear thermo
   win (ceiling ~99% FDR — a hit must also overlap a regulatory element). The pivot to
   splice-consequence + training-free framing is the response; the arm must beat the baseline on
   the *splice-consequence* signal, not binding energy.
3. **Reproducibility ceiling:** only ~19% of Holgersen splicing events reproduce between repeats
   — use the large-effect (|dPSI| ≥ 0.5) labels.

## 8. Hard dependencies (blockers before build)

- **Holgersen Supplementary Tables S1 (sequences) + S4–S6 (per-pair labels)** — NOT in repo, raw
  reads proprietary. Core benchmark is blocked without them. **Obtain / confirm access first.**
- **Super-minigene:** exact sequence + `SMN2Sup3445` confirmation — lab-pending (Ask 1). Bridge
  can start on Singh 2009 pCI without it.

## 9. Open items to flesh out later (the "we'll flesh it out more")

- Exact IntaRNA accessibility (ED) settings vs Holgersen's RNAcofold/RNAplex ΔG — the one genuine
  methodological opening; specify precisely.
- Exact Pangolin readout at off-target sites (which tissue tracks, which usage delta).
- The preregistered baseline definition + margin + LOASO protocol details.
- How the off-target-vs-length prediction combines with `pre/`'s on-target length result into a
  single therapeutic-window statement.
- Clean-room boundary spec: what exactly the bridge may import from `pre/` + a manifest test that
  fails if the CORE touches `pre/` or any PSI.
- Adversarial structural review of this plan (deferred; Codex out until 09-10).

## Spike record (2026-09-07, opus 5)

- **Holgersen benchmark:** PARTIAL-GO. Right chemistry/labels (splicing only, LOASO-ready) BUT the
  naive thermo-beats-editdistance claim is pre-empted → pivot to splice-consequence primary +
  training-free framing. Suppl tables S1/S4–S6 required (not in repo).
- **Super-minigene / bridge:** PARTIAL. On-target shorter-is-better ground truth exists NOW (Singh
  2009 pCI-SMN2); super-minigene has no ASO data + sequence lab-blocked; ISTL1 intact only in
  SMN2Sup3445. `pre/` already has a reusable on-target stack + blind-grading discipline.
