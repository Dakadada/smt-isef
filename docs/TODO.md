# TODO — state as of 2026-08-24 (evening)

Derived from plan rev 5. Everything listed as OPEN in the morning revision is
now closed. See `docs/HANDOFF.md` for the durable record.

Cluster state: idle except job 21133941 (full pytest, new baseline pending).

---

## BLOCKED — waiting on the wet lab

### T9 / Phase B1-B4 — the arm that decides the project
Blocked on partner densitometry. Request sent 2026-08-23; not yet returned.

- **B1** ingest, written against the *first real file* — spec already delivered
  in `docs/densitometry-spec.md`, so nothing has to be built when it lands
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

## CONDITIONAL — gated on a result that did not arrive

### T6 / G2 — full-transcriptome sieve for the winning length
Explicitly conditional on G1 showing a real trade-off. **G1 found no interior
optimum**, so the trigger condition is not met and G2 does not run. Recorded as
a decision, not an omission.

---

## DELIBERATELY NOT RUN

### S4 — C2B34 Alu-arm duplex in 3D
Assessed and declined, with reasons, rather than left ambiguous:

1. **Ungateable at the required size.** The arms are 280 and 292 nt with a
   250-bp duplex core. Boltz's largest validated input in this project was
   356 nt of RNA, and there is no deposited RNA duplex at that scale to gate
   against. The S-gate precedent is that ungated predictions carry no claim.
2. **The claim it was meant to illustrate did not survive.** S4 was planned as
   the figure for the primate-specificity mechanism. P is corrected but still
   NOT ESTABLISHED, so there is no mechanism for the figure to depict.
3. What remains is decoration, and the project does not need another
   uninterpretable picture.

A defensible reduced version exists — predict a 40-60 nt sub-duplex of the
paired core and gate it on a deposited RNA duplex — if a figure is later wanted
for the poster. It would be illustration, explicitly not evidence.

---

## CLOSED THIS ROUND

| | outcome |
|---|---|
| T1 pooled rMATS + enrichment | A5 pooled **corrected p = 0.0253** — the only result surviving full correction |
| T2 S1-gate | **PASSED, F1 0.971** (my prediction that RNA accuracy would break it was wrong) |
| T3 / S1 | **NULL, p = 0.4545**, 16/20 usable |
| T4a G on-target | done; scan-edge artifact caught |
| **T4b / G1 off-target** | **DONE.** No Goldilocks point in 6-18; margin rises monotonically with length. Feasibility pre-gate added — A3 length 6/7 winners form no duplex at all. Free validation: the nusinersen 18-mer has exactly two exact transcriptome sites, SMN2 and SMN1. |
| T5 / P Alu identity | **CORRECTED.** Strand bug found: upstream arm z -0.54 -> +2.57 with a perfect Pol III promoter. Verdict still NOT ESTABLISHED (both arms just under the pre-declared 2.58), but the earlier asymmetry was an artifact. |
| **S2 ceiling in 3D** | **INCONCLUSIVE**, rho = -0.200, p = 0.4583 at n=4. No new GPU — the unbound predictions already existed. |
| **S3 hnRNP A1** | **UNRESOLVED.** Gate passed (F1 0.857) but both gate RNAs are UAG-rich and ISS-N1 has no UAG, so the target is outside the gate regime. Not reported as a negative. |
| T8 Phase B prep | `docs/densitometry-spec.md`, written before any data exists |

---

## REMAINING HOUSEKEEPING

- **Full pytest baseline** — job 21133941 running. Previous 197 passed / 1
  skipped was re-confirmed today but predates the 7 new G1 tests; expect ~204.
- **T7 write-up** — lead with the **accessibility ceiling** (A7) and the
  **feasibility gate** (G1). Never lead with the sieve.
