# TODO — state as of 2026-08-24

Derived from plan rev 5. Every item below is either open, blocked, or optional.
Anything not listed here is done and recorded in `docs/HANDOFF.md`.

Cluster state at time of writing: **queue empty, nothing running.**

---

## OPEN — can start now, nothing blocking

### T4b / G1 — Goldilocks off-target panel   [CPU, small]
The only remaining half of the Goldilocks study. The on-target axis is done
(`docs/g-ontarget-axis.txt`, `src/smt/validate/goldilocks.py`) and its own
printed conclusion flags this gap:

> "Off-target scoring (G1, targeted panel) attaches to these winners; until it
> runs, this is the on-target axis ALONE and implies nothing about a trade-off."

Work: for each length 6-20, take that length's winning span-robust on-target
footprint, build the ASO as the reverse complement of its window, and score it
against a targeted panel with `smt.sieve.scan.score_candidate` /
`seed_scan_transcript`.

Panel: POLR2H, MAP1S, ANKRD28, CCSER2 (the four validated hits) plus a matched
control set drawn from `sieve/annotated.csv`. Mark nusinersen's 18-mer
(`scan.py:90`) and the published 8-mer on the curve.

Pre-stated expectation, to be honest about before running: shorter ASOs may win
on off-target simply because shorter seeds match less often. **If the curve is
monotonic in length with no optimum, there is no Goldilocks point** — report
that, do not re-bin until one appears.

Also folds in the reframed "polymerase hot spots" idea as a covariate: does the
ASO hybridize more promiscuously in repeat-rich transcripts?

### T-suite — full pytest re-run   [CPU, minutes]
Last recorded baseline is **197 passed / 1 skipped** (job 21065504), taken
*before* `structure/`, `goldilocks.py`, and `kinetics/alu_identity.py` landed.
The current number is unknown. Run `slurm/pytest.sbatch` and record it.

---

## OPTIONAL — GPU is free; these were planned but never gated as required

### T6 / G2 — full-transcriptome sieve for the winning length
Explicitly conditional: **only if G1 shows a real trade-off.** Do not run before
G1 answers.

### S2 — is the accessibility ceiling visible in 3D?
Unbound Boltz predictions of the four cryptic targets. Pre-specified test: their
3D accessibility must rank in the same order as their `A_5ss` baselines (STMN2
most exposed, UNC13A_DONOR most buried). A disagreement is a real falsification
signal and goes in the record either way. ~40 predictions.

### S3 — hnRNP A1 ternary complex
RNA alone / RNA+UP1 / RNA+UP1+ASO. This is the **spec bullet that is currently
unmet** — "synthetic chemistries and biological guide RNAs as identical steric
blockades against hnRNP A1/A2." A4 measured only geometric shadows from
deposited dCas13 structures and never modeled the protein being competed with.
~75 predictions. Of the optional GPU work this is the one with the clearest
claim attached.

### S4 — C2B34 Alu-arm duplex in 3D
The inverted-repeat duplex that loops the backsplice into proximity. ~25
predictions. Note P came back negative, so this is now a *structure* figure, not
a primate-specificity mechanism figure — describe it accordingly.

### T7 — write-up
Lead with the **accessibility ceiling** (A7). Never lead with the sieve.
Blocked in practice on G1 finishing so the Goldilocks section is whole.

---

## BLOCKED — waiting on the wet lab

### T9 / Phase B1-B4 — the arm that decides the project
Blocked on partner densitometry. Request sent 2026-08-23; not yet returned.

- **B1** ingest, written against the *first real file* — spec already delivered
  in `docs/densitometry-spec.md` so nothing has to be built when it lands
- **B2** fit the `A_5ss -> PSI` link; functional form and fit/held-out split
  declared *before* fitting
- **B3** held-out prediction — report held-out error, not fit error
- **B4** risdiplam as a separate U1 on-rate parameter, never folded into the
  occlusion term

**Standing rule, unchanged:** nothing consumes PSI until real PSI exists. No
synthetic calibration, no placeholder link fit, no illustrative numbers.

Phase B's arrival must **not** quietly re-open S1's or P1's pre-registered
decision rules.

---

## CLOSED THIS ROUND — for the record, not for re-doing

| | outcome |
|---|---|
| T1 pooled rMATS + enrichment | A5 pooled **corrected p = 0.0253** — the only result surviving full correction |
| T2 S-gate | **PASSED, F1 0.971** — and my stated expectation that RNA accuracy would break it was wrong |
| T3 / S1 | **NULL, p = 0.4545, 16/20 usable.** Deep run made the pre-registered criterion worse, not better; kept as primary anyway |
| T4a G on-target | done, plus a scan-edge artifact caught |
| T5 P Alu identity | **NOT ESTABLISHED** — one arm of four, and the wrong one. Verdict is the pre-declared fallback: "Alu-scale inverted repeats, subfamily unresolved" |
| T8 Phase B prep | `docs/densitometry-spec.md` written before any data exists |

Two of my own hypotheses failed this round. Both stay in the record.
