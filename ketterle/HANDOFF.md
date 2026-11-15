# HANDOFF — ketterle/ (standalone in-silico specificity arm)

**Created 2026-09-06. Spun out of the feng program into its own top-level folder.**
Named for Wolfgang Ketterle. Sibling to `feng/`, not under it. Clean-room + honesty
architecture + freeze discipline inherited from the feng program (see `feng/HANDOFF.md`
and `feng/CLAUDE.md`).

---

## ★ LIVE STATE — 2026-09-06

**Goal of the arm:** an *in-silico specificity* result for splice-switching ASOs (nusinersen /
SMN2 exon-7 anchored), computable on the existing feng toolset (Pangolin steric-mask splice
scoring; IntaRNA duplex ΔG; RNAfold/ViennaRNA structure). Must produce a POSITIVE, honest,
preregisterable result — a documented null is explicitly NOT the goal (user, 2026-09-06).

**Current status: PLAN DRAFTED — BACK BURNER (2026-09-07).** Bivalent AND-gate DEAD; replacement =
an off-target specificity method (causal gate: IntaRNA-accessibility filter → Pangolin
splice-consequence), design drafted in **[PLAN.md](PLAN.md)**. NOT implementing; to flesh out later.
Two feasibility spikes done (Holgersen benchmark, super-minigene/bridge) — see PLAN.md §Spike record.

### DEAD END #1 — bivalent coincidence / AND-gate on SMN2 exon-7 (NO-GO)
The idea: a pair of short oligos that must BOTH bind a compact on-target to act; if the two
required sites are close (small window W), off-target co-occurrence is rare → specificity gain.
Killed by a structural **vise** — the two real SMN2 exon-7 targets fail at OPPOSITE ends:

- **ISTL1** (intron-7 long-distance stem): arms **279 nt apart** (Singh 2013 NAR 41:8144).
  W≈279 → two near-complement sites co-occur within that window all over the transcriptome →
  no specificity gain. (Whole ISTL1/2/3 family shares the +275–297 3′ strand → all large-W.)
- **TSL2** (5′ss terminal hairpin, exon7/intron7 junction): tiny W (~3–19 nt) but the **8-bp stem
  is too weak**. Feasibility gate (RNAfold 2.7.2 + IntaRNA 3.4.1, feng-env):
  - TSL2 = `AUUCCUUAAAUUAAGGAGU` (19 nt; last 17 nt exon 7 + first 2 nt intron 7; chr5:70076558–
    70076576 +; canonical GT donor confirmed). Fold `((((((((...))))))))` MFE **−5.60**, ens −6.02.
  - A single 9-mer vs the 3′ arm binds **−4.39 kcal/mol** → opens the whole (cooperative) stem
    alone. NO split-weak pair exists (3′ arm never sub-threshold; robust δ∈[2,4]). Real ASO
    chemistry (MOE/PS/LNA) binds tighter → reinforces NO-GO.
  - Gate script: `scratchpad/gate.py` (session scratchpad, not committed).

**The vise (the real finding):** compact enough for rare coincidence ⟺ stem weak enough that ONE
oligo opens it (bivalency pointless); strong enough to need two oligos ⟺ arms far apart (large W,
common coincidence). Bivalent AND-gate specificity on SMN2 exon-7 is not viable. Other candidates
(bipartite ISS-N1 = inside nusinersen's footprint; 3′ss tandem ESS = one oligo spans; exon+6 /
intron+100 = W~140; URC1/2 = enhancers/wrong direction; TSL1+TSL2 = not coupled) all fail too.

### Archived design (for reference / negative result)
The full seal-ready v3 design for the (now-dead) bivalent arm — `H-off-pair` predicate with
inter-site complementarity, within-reference positional null, distinct-site burden, TSL2 primary +
ISTL1 large-W negative control, monotonic-in-W prediction, W-titration — lives in this session's
transcript. It passed two adversarial-review rounds (DO-NOT-ADVANCE → revised → pivoted) before the
feasibility gate killed the on-target. Preserve as a documented negative if ever needed; do NOT
resurrect without a target that beats the vise.

### NEXT (all deferred — BACK BURNER)
- Flesh out PLAN.md open items (§9): IntaRNA-ED settings, Pangolin off-target readout, prereg
  baseline+margin+LOASO protocol, off-target-vs-length ↔ pre/ on-target combination, clean-room
  boundary + manifest test.
- **Obtain Holgersen Supplementary Tables S1 (sequences) + S4–S6 (labels)** — hard blocker for the
  core benchmark (not in repo, raw reads proprietary).
- Adversarial structural review of PLAN.md (Codex out until 09-10 → opus 5 fallback) → seal
  candidate score def → build. None of this until user takes it off the back burner.

## Model tiering (inherited; Codex OUT until 2026-09-10)
Dev/design/review all on the **opus 5** Claude fallback lane while Codex is exhausted (banked reset
NOT spent — nothing so far was Codex-only). Monitoring = haiku. Main chat = opus 4.8. Do NOT spend
the banked reset unless write-heavy Codex-only work or a hard pre-10th gate appears.

## Tools confirmed (feng-env)
`/home/dashii/orcd/pool/feng-env/bin/` has `RNAfold` (ViennaRNA 2.7.0/2.7.2) + full `IntaRNA` 3.4.1
suite. Interim pytest runner: `/home/dashii/orcd/pool/smt-isef-env/bin/python`.
