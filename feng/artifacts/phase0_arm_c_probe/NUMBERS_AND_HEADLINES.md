# Arm C — Numbers, Results, Headlines

Working reference. Every actionable number from the Arm C bind→act gate (job `21881182`,
`GATE_RESULT.md`). **Honesty tags are load-bearing** — do not strip them when pulling a line
into a slide or the paper.

- `[PREREG]` = the locked, blind-set criterion. Untouchable.
- `[RAW]` = true statement about the data, prereg-neutral. Safe to say out loud.
- `[REFRAME]` = same data, metric appropriate for a 0.001-prevalence task. Legit, but label
  it as "the right metric for rare events," NOT as clearing the prereg bar.
- `[PENDING]` = not yet computed; fill from the DeLong/bootstrap fork.

> ⚠️ **CODEX SOL ADVERSARIAL REVIEW (structural) — 3 HIGH blockers. This is NOT publishable as
> the prespecified Arm C verdict yet.** The metric *implementations* are sound (orientation,
> no leakage, sklearn-convention AUC/AP, real IntaRNA), but: (1) the +0.05-on-both bar is a
> SUBSTITUTED, weaker endpoint than PLAN.md:134's locked one (nested LOAO, combined baseline,
> AUCPR gain ≥0.10 or ≥2×, bootstrap p<0.01, FDR, external replication); (2) negatives are
> RESPONDER-ONLY (S4) — "silent binding" rows include exons never assayed/expressed, so the
> 1-in-1000 promiscuity claim is confounded; (3) 16 low-complexity SBOs were dropped pre-scan,
> so only **65/81** SBOs were analyzed (headline said 81). Until resolved, this is an
> **exploratory Phase-0 probe**, not the locked prereg verdict. See §8.

---

## 1. The verdict (exploratory probe — NOT the locked endpoint)

| | value |
|---|---|
| **Prereg GO bar** `[PREREG]` | beat best baseline by **+0.05 on BOTH** AUC & AP |
| **Verdict** `[PREREG]` | **NO-GO** — cleared AUC by +0.041, AP by +0.002; both under margin |
| Freeze status | candidate / paper_supported; no sealed ΔPSI read → **freeze-safe** |

The NO-GO stands until a prespecified significance test says otherwise. Do not lower the bar.

## 2. Headline metrics

| metric | Model | edit_distance | duplex_dg | note |
|---|---|---|---|---|
| **AUC** | **0.838** | 0.674 | 0.797 | model tops both |
| **AP** | **0.042** | 0.040 | 0.004 | model tops both |

- `[RAW]` **Model tops both baselines on both POOLED metrics.** (Per-fold, individual baselines
  still beat the model on some SBOs — e.g. duplex_dg > model for DG1063, DG132. "Clean sweep"
  is pooled only, not per-SBO.) Margin under the +0.05 bar.
- `[RAW]` AUC **0.838** absolute — a genuinely good classifier on a transcriptome-wide task.
- `[RAW]` Per-fold dG AUC frequently **0.8–0.99** — on individual SBOs, discrimination often
  near-perfect.

## 3. The multipliers (rare-event framing)

- `[REFRAME]` **AP 0.042 vs prevalence 0.00114 = ~37× enrichment over random.** Model surfaces
  candidate responders ~37-fold better than chance. (Was mis-stated as 42× from rounded inputs.)
- `[REFRAME]` **Model AP 0.042 vs duplex_dg AP 0.004 = ~10.5×.** On precision, binding beats
  energy by an order of magnitude. NOTE (Codex-corrected): the +0.05 unit bar is *demanding* at
  0.00114 prevalence but **not mathematically unavailable** — AP can reach 1.0 at any prevalence,
  so do NOT claim the unit-delta is "structurally near-zero." Ratio framing is the appropriate
  lens for rare events, but that's a judgment, not a math necessity.
- `[RAW]` AUC margin over sequence baseline: 0.838 − 0.674 = **+0.164** (landslide in units).
- `[RAW]` AUC margin over energy baseline: 0.838 − 0.797 = **+0.041**.
- `[RAW]` AP margin over sequence baseline: 0.042 − 0.040 = **+0.002**.

## 4. The finding nobody's leading with

- `[RAW]` **Reported-responder rate = 217 / 190743 = 0.00114** among predicted binding
  candidates (sequence/energy similarity from the scanner — NOT assayed binding events).
- ⚠️ **CONFOUND (Codex HIGH-2):** the 190,526 "non-responder" rows are candidates *absent from
  S4's responder-only list* — that includes exons never assayed or not expressed, not just true
  non-responders. So "binding ≠ acting" is **not yet established**; the honest phrasing is
  *"predicted binding candidates not reported as responders."* Needs coverage-aware,
  power-analyzed negatives before the promiscuity headline can stand.
- **Headline (conditional):** *IF properly-powered negatives confirm it,* the story is binding is
  promiscuous, action is rare, specificity is the game → feeds Arm E. **Not bankable until the
  negative set is fixed.**

## 5. Scale / credibility

| | value |
|---|---|
| Bound (SBO,exon) rows | **190,743** |
| Positives (actionable) | **217** |
| Prevalence | **0.001137** (217/190743) |
| SBO groups | 65 (**40 evaluable**; 25 degenerate 0-pos/0-neg held-out — pool but yield no per-fold metric) |
| SBOs scanned | 81 (GENCODE-v50 genome scan) |
| Model | LOSO logistic, binding features (numpy, deterministic) |

## 6. Mechanistic read (bio)

- `[RAW]` **duplex_dg** (binding energy) = dominant *global* ranker (AUC 0.797), near-useless on
  *precision* (AP 0.004). Tighter duplex → more likely to act. Steric block needs stable
  occupancy.
- `[RAW]` **edit_distance** (on-target sequence identity) = weak ranker, strong *precision*
  signal (AP 0.040). Hitting the intended site is what concentrates real hits.
- `[RAW]` **Model** = combines both. No hidden "binding grammar" beyond *bind tight + bind the
  right place*. Not a failed experiment — a demonstration that steric-block action follows
  thermodynamics + target identity.

## 7. Significance test (DONE — DeLong + cluster bootstrap)

Predictions **recomputed** (not persisted) by replaying the deterministic 65-fold LOSO loop over
`armc_features.parquet` (190743×15) — numpy, no RNG, no cluster/IntaRNA re-run. Pooled AUC: model
0.8377 / dg 0.7972 / edit 0.6743. Pooled AP: model 0.0420 / dg 0.0036 / edit 0.0396.

- `[RAW]` **DeLong (iid), model AUC vs duplex_dg:** diff **+0.0405**, z=3.801, **p=0.000144**,
  95% CI [+0.0196, +0.0613]. ⚠️ **iid ignores LOSO clustering (pooled scores span 40 fold-models)
  → OVERSTATES significance. Do not quote this p as the result.**
- `[RAW]` **Cluster bootstrap (B=2000, resampled by SBO group) — the honest number:**
  - **AUC diff (model − duplex_dg):** median **+0.0404**, 95% CI **[+0.0003, +0.0782]**, boot
    **p = 0.049**. Significant, but **barely** — CI lower bound grazes zero (+0.0003). Fragile.
  - **AP diff (model − edit_distance):** median **+0.0042**, 95% CI **[−0.0055, +0.0175]**, boot
    **p = 0.399**. **Not significant** — CI spans zero.
- **Verdict — PARTLY:** binding significantly out-*ranks* duplex energy (AUC p≈0.049) but the edge
  is real-yet-fragile; on *precision* it's flat (AP p=0.40) at 0.001 prevalence. **Consistent with
  the prespecified NO-GO** — the +0.041 AUC edge is ~significant but below the +0.05 effect-size
  bar, and AP is noise.
- **Honest headline you CAN say:** *"Binding adds a statistically detectable but small ranking edge
  over duplex energy (cluster-bootstrap p≈0.05, AUC +0.04); the edge does not reach the
  prespecified effect-size bar and there is no detectable gain in precision."* Earned, not spin.
- **Do NOT say** "p<0.05, binding beats energy" as a clean win — the honest p is 0.049 with a CI
  kissing zero, not the iid p=0.0001.

---

## 8. Codex sol adversarial review — findings (2026-09-03)

Read-only structural audit. **Overall: NO-GO metrics are sound, but NOT publishable as the
prespecified verdict.** Blockers must be resolved or the result relabeled as an exploratory probe.

**HIGH (blockers):**
1. **Substituted endpoint** (`armc_gate.py:51-54,266-272`) — the +0.05-on-both bar is NOT
   PLAN.md:134's locked endpoint (nested LOAO, combined multi-feature baseline, AUCPR gain
   ≥0.10 or ≥2×, bootstrap p<0.01, sensitivity/FDR, external replication). Fix: relabel as
   exploratory Phase-0 probe, OR evaluate the actual locked endpoint.
2. **Responder-only negatives** (`README.md:22`, `armc_build.py:770-781`, admitted `:46-56`) —
   non-S4 candidates labeled "bound-not-acted" without coverage-aware power analysis;
   unexpressed exons become false negatives. Fix: "not reported as responders," or build
   assayed/expressed/powered negatives + label-stability analysis.
3. **Undisclosed 16-SBO exclusion** (`armc_build.py:598,788-794`) — headline says 81 scanned;
   only 65 analyzed. Fix: disclose 65/81, include the 16 with a scalable method or prespecify +
   sensitivity-test the exclusion.

**MED:** micro-pooling weights by candidate count & compares different fold-models' probabilities
(`armc_gate.py:204`) — report group-aware/macro + cluster-bootstrap (done, §7); hand-rolled
optimizer has no convergence/parity test (`armc_gate.py:108`); merge readiness counts parquet
files not ID-matched sentinels (`armc_gate.sbatch:49`); 42×→~37× + false "structurally near-zero"
(fixed §3); "no baseline beats it anywhere" false per-fold (fixed §2).

**CLEAN (attacked, held):** scan orientation (revcomp once), minus-strand exon orientation,
1-based interval→exon_key conversion, positive-locus injection (no boundary under-labeling), no
S4-energy leakage, all 190743 rows real IntaRNA, sklearn-convention AUC/AP, pdfplumber parse,
samtools fallback, DG111/DG2386 duplicate.

---

_Source: `GATE_RESULT.md`, job 21881182. DeLong/bootstrap fork DONE (§7). Codex sol structural
review DONE (§8). Standing gate: findings logged; HIGH blockers open pending user decision._
