# Phase 0 — Arm C (Bind→Act gate) data/power probe

**Status: AMBER → GO-leaning (updated 2026-09-02 after parsing Holgersen S4/S5).**
The clustering fear did NOT materialize. See "S4/S5 extraction" below.
Author: Claude (general-purpose subagent) for the initial probe; extraction + revised
verdict by Claude main thread. Read-only; no heavy compute. Date 2026-09-02.

## S4/S5 extraction (2026-09-02, Claude main thread — actual parse)

Both supplementary PDFs (`feng/sj-pdf-4-...pdf` = S4, `feng/sj-pdf-5-...pdf` = S5) parsed
with pdfplumber. Findings (all `paper_supported`):

- **Independence unit = ~75 SBOs** (S4 splicing) / 64 (S5 expression), spread across
  DIVERSE genomic loci (chr1/3/5/8/16/22/…), NOT the feared ~16 clustered genes. This
  clears the proposed ≥25 LOAO-fold bar with headroom. The Arm-A clustering/poolability
  failure mode does **not** apply.
- **Both axes are present in one harmonized source, no cross-study pooling.**
  - Binding: RNAplex ΔG + RNAcofold ΔG + edit-distance / longest-exact-match categories.
  - Action: S4 = per-exon Case vs Control PSI (splicing, 724 scored sites,
    median 6 sites/SBO, |ΔPSI| median 0.58); S5 = per-gene Log fold change
    (expression, 30,493 rows, 22,515 sequence-matched, |LFC| median 1.5).
- **Both tables are HIT/RESPONDER lists only** — S4 lists fired exons; S5 is a
  pre-thresholded DE-gene list (every matched row |LFC|≥1). **Neither tabulates the
  bind-but-no-act negative class** the two-gate model needs.

**Revised decisive quantity:** no longer "is N≥25?" (yes — 75 SBOs). It is now
**"does a constructed bind-but-no-act negative set yield a contrast a held-out model can
beat an edit-distance/energy baseline on?"** The negatives are constructible with feng's
existing seed-extend scanner (enumerate per-SBO sequence-matched transcriptome sites,
label those absent from the S4/S5 hit lists as non-responders) — no new data, no heavy
compute, before any 5k–20k core-hr commit.

**Verdict: AMBER → GO-leaning, conditional on the negative-set construction above.**

---

### Original probe body (pre-extraction, AMBER) follows.

## The question

Arm C's primary endpoint needs a held-out-ASO actionability model (duplex/accessibility
features → active/inactive label) beating an edit-distance/energy baseline under **nested
leave-one-ASO-out**, with an ASO-cluster bootstrap on the AUCPR difference. That requires
enough *independent* ASOs, each anchored to a *comparable* binary/continuous action label,
to (a) fit something, (b) leave folds out without leaking regional sequence context, and
(c) get a stable AUCPR-gain estimate. Arm A died on exactly this shape of problem
(poolability + power, not access) after apparent overlap looked fine on paper. This probe
asks the same question of Arm C before its 5,000–20,000 CPU-h build starts.

**Preregistered Arm C threshold check:** `feng/PLAN.md` gives Arm C an AUCPR-gain success
bar (≥0.10 absolute, bootstrap p<0.01), an external-replication bar (GSE307431 + Scharner,
directional only, "neither used for tuning"), and a kill bar (>30% of negative labels flip
under reasonable coverage/power re-thresholding, or external AUCPR gain <0.05) — but,
unlike Arm A's explicit "≥25 independent recoverable ASO windows," **Arm C has no
preregistered minimum effective-N.** Arm C also has no `**Data:**` line at all in the plan
(Arms A/B/D each do) — the training corpus is implied only by the "Why it is open"
paragraph's citations to Scharner and Holgersen. Both gaps are addressed below (finding 3
proposes a threshold; the inventory below fills the missing Data line).

## Data inventory — binding axis × action axis, same ASOs

| Source | ASOs (nominal) | Binding axis | Action axis | Access | Role for Arm C |
|---|---|---|---|---|---|
| **Holgersen et al. 2021**, *Nucleic Acid Ther.* 31(6):392-403, doi:10.1089/nat.2020.0921 (bioRxiv 10.1101/2020.09.03.281667; PMC8713556) | **81 SBOs** (40 designed + 31 "promiscuous" + 5 RBP-motif + 3 non-targeting + 2 sense) across **16 genes**, 105 RNA-seq runs (HepG2 n=48, HEK293T n=53, PXB hepatocytes n=4), all PS-MOE, 16–20 nt | **Computed**: edit distance, longest exact match, ViennaRNA ΔG (RNA-cofold/RNA-plex) | **Measured**: RNA-seq PSI, custom spliced-read counting, ΔPSI; 724 large splicing changes (\|ΔPSI\|>0.5) reported genome-wide | Raw FASTQ **proprietary** ("only to not-for-profit groups not competing with Deep Genomics" — PLAN.md line 283 already flags this). Processed per-site values in **Supp. Tables S4 (715 KB) / S5 (18.4 MB)**, both PDFs — schema/machine-readability **not yet inspected**, open risk. | Only plausible **primary training corpus** for the nested LOAO fit. |
| **Scharner et al. 2020**, *Nucleic Acids Res.* 48(2):802-816, PMID 31802121 | **10 ASO labels**, but these collapse to **~2 base target sequences** ("ASO1", "ASO2") plus one nusinersen-adjacent SMN control, each run in multiple **chemistry/conjugate variants** (MOE, cEt/DNA, MOE/DNA, GalNAc3) | **Computed** (RNAhybrid) **+ measured Tm** for a subset of duplexes (Supp. Table S3) — the one genuine *measured*-binding data point in this inventory | **Measured**: radioactive RT-PCR / autoradiography, % exon inclusion, ΔPSIδ; 140 candidates → 108 valid amplicons → 17 validated off-targets (>10% splicing change) | No GEO/SRA. Values mostly in figures/S1-S3; not fully tabulated. | **External replication only** (already so scoped in PLAN.md: "GSE307431 and Scharner, neither used for tuning"). The apparent N=10 is a **chemistry-variant inflation trap** on ~2-3 true target sequences — structurally the same trap the plan explicitly forbids for F18MOE/F18OMe/F20PMO. Correctly excluded from LOAO training by the plan's own design. |
| **GSE307431 / PRJNA1321660** — Ottesen et al. 2026, *NAR Mol. Med.* 3(1):ugag002, PMC12805893 | **3 chemistry variants at ONE target site** (F18MOE, F18OMe, F20PMO = nusinersen-identical/near-identical sequence at ISS-N1) | **Computed** baseline (edit distance ≤3 → 14% sensitivity at 99% FDR, per PLAN.md's own citation) | **Measured**: RNA-seq, 8 orthogonally-validated off-target exons (incl. POLR2H exon 2) out of a genome-wide candidate universe | **Public**, raw reads available (GEO/SRA). | **External replication only** — again already correctly scoped that way in PLAN.md. Also the plan's explicit "no chemistry model on F18MOE/F18OMe/F20PMO" do-not applies directly: 1 effective target-site cluster, not 3 independent ASOs. |
| **Hua 2007 / Hua 2008 / Singh 2010 SMN2 walks** (already harmonized in `feng/artifacts/phase0_arm_a_overlap/harmonized_aso_positions.csv`) | 125 records / 115 distinct target intervals, all at **one locus** (SMN2 exon 7 ± flanking intron) | **Computable** (sequence on file; duplex/accessibility energy not yet computed) | **Measured/mixed**: 65/108 rows have a stated numeric or directional effect; 43/108 are `figure_reported` only (no safe numeric value — an Arm A finding, applies identically here) | Public (open-access supplements), already extracted. | Structurally usable as *on-target* positional-tiling analogue of the same binding→action question, but single-locus: Arm A already showed the **defensible independent-window count for this exact ASO set is 11** (pairwise non-overlapping), not 115. Same clustering discount applies if used here. |
| **"Determining off-target effects..." (2025)**, *Hum. Mol. Genet.* 34(22):1912-, PMC12581821 / bioRxiv 10.1101/2025.01.09.631923 — APP/ATXN3/HTT AONs in neuronal iPSC | 4 ASOs (APP, ATXN3 ×2, HTT), MOE/PS, 3 iPSC lines | **Computed** (GGGenome, ≤3 mismatches) | **Measured**: RNA-seq, rMATS + Whippet ΔPSI | Raw FASTQ **destined for EGA, no accession yet published** — currently inaccessible. Authors report only **~2% of detected splicing changes correlate with a predicted binding site** — i.e. most action events have no candidate binding site to pair with. | Not usable now (embargoed). If it lands, it would add genes outside SMN2 for cross-locus generalization, but the 2% overlap number is itself evidence *for* Arm C's core claim (binding ≠ action) and *against* naive site-matching yielding many usable pairs. |

**Distinguishing measured vs. computed binding, explicitly:** across every source above, the
*action* axis is always measured (RNA-seq PSI or RT-PCR quantification). The *binding* axis
is almost always **computed** (edit distance / mismatch count / ViennaRNA or RNAhybrid
energy) — Scharner's Tm subset is the sole example of *measured* duplex binding tied to
splicing outcome, and it covers only a handful of duplexes, not a usable training set on its
own.

## Count + poolability

**Nominal count:** 81 (Holgersen) + 10 (Scharner) + 3 (GSE307431) + 125 (Hua/Singh) + 4
(2025 iPSC, inaccessible) = well over 200 "ASO" labels on paper.

**Defensible count is far smaller, and the discount mechanism is the same one Arm A hit:**

1. **Chemistry-variant inflation.** Scharner's 10 and GSE307431's 3 are chemistry/conjugate
   variants of ~2-3 and 1 target sequence respectively — not independent binding→action
   observations. PLAN.md's own do-not list ("Do not train a chemistry model on
   F18MOE/F18OMe/F20PMO") already forbids treating these as distinct points, and its Arm
   C/B text already scopes both sources as *external replication only*, never LOAO training.
   This is correct as written — flagging it here only to make explicit that the "13
   external-source ASOs" figure is really ~3-4 independent target sequences.
2. **Locus/gene clustering inside the one usable training corpus.** Holgersen's 81 SBOs
   span only **16 genes** — roughly 5 ASOs per gene. ASOs targeting the same gene share
   local sequence composition, splice-site architecture, RBP-motif density, and often the
   same nearby candidate off-target exons. A LOAO fold that holds out one ASO but trains on
   four siblings targeting the same gene leaks regional context into "held-out" evaluation
   — structurally identical to Arm A's finding that ASO-walk positions were confounded with
   study/region, not exchangeable. **The honest independence unit is plausibly gene/locus,
   not ASO**, which would put Holgersen's effective-N near **16**, not 81. This has not been
   verified against the actual gene list (not yet extracted) — it is a risk bound, not a
   confirmed number.
3. **The 31 "promiscuous" SBOs are selected, not random.** Holgersen's paper says these were
   *chosen* to have many near-matches — a positive-outcome-biased design for roughly 38% of
   the corpus. This inflates apparent hit rate/variance structure relative to a randomly
   drawn ASO population and needs explicit stratification (design-source as a covariate or
   held fixed across folds), or it silently changes what "active" means across the corpus.
4. **Label-definition mismatch across studies, mirroring Arm A's non-poolable-outcomes
   finding exactly.** Holgersen: |ΔPSI|>0.5 via RNA-seq spliced-read counting. Scharner:
   >10% splicing change via radioactive gel RT-PCR. GSE307431: rMATS/other, own threshold.
   None of these are the same measurement on the same scale; a shared "active/inactive"
   binary needs a **defined, preregistered per-study recalibration** (e.g., an
   assay-specific noise floor, not a single global ΔPSI cut), or the external-replication
   check is comparing incommensurable quantities the way Arm A's cross-study Spearman would
   have.
5. **Confounds, named:** gene/locus (co-varies with ASO within Holgersen); assay method
   (RNA-seq vs. radioactive RT-PCR, co-varies with study); chemistry/conjugate (co-varies
   with "ASO identity" in Scharner/GSE307431 unless explicitly collapsed); design intent
   (promiscuous vs. designed, co-varies with expected positive rate). All four are the same
   *shape* of confound Arm A hit (study confounded with region/assay/chemistry), just
   relocated to a different arm.

## Power

Treating **gene/locus** as the honest LOAO-cluster unit (per finding 2), Holgersen's
corpus gives an effective-N somewhere in **[16, 81]**, with 16 the defensible floor absent
further stratification and 81 only defensible if within-gene ASO independence can be shown
(e.g., low correlation of residual features after removing gene-level effects — untested).
AUCPR estimation with an imbalanced positive rate (724 large changes reported against a much
larger candidate-site denominator; exact ratio not yet extracted) needs enough *positive*
folds specifically — a LOAO fold contributes no information to the pooled AUCPR estimate if
its held-out ASO has zero active sites, which the "median 5 (range 0–140)" ΔPSI-changed-exon
count in the source implies is plausible for some ASOs (0 at the low end of that range).

**No preregistered Arm C effective-N threshold exists** (see "The question" above). Proposed
threshold, mirroring Arm A's ≥25 windows but adapted to Arm C's LOAO-cluster unit: **≥25
independent gene/locus clusters, each contributing at least one active and one
candidate-inactive site**, before committing the 5,000–20,000 CPU-h build. Under this
proposed bar, the gene-level floor (~16) fails and the ASO-level ceiling (~81) passes —
the actual number is unknown until the supplement is parsed, which is why this probe cannot
resolve GO/NO-GO on its own.

## Access vs. quantity/comparability

**Not primarily an access problem**, mirroring Arm A's own conclusion. Every action-axis
and binding-axis value needed to build Arm C's features and labels is published in
open-access papers or supplements; only Holgersen's *raw* FASTQ is restricted, and Arm C
does not need raw reads — it needs the processed per-site table, which is public alongside
the article. The one live access gap is the 2025 iPSC dataset (EGA accession not yet
issued) — not currently usable, immaterial to the current build decision since it was never
part of the plan's cited corpus.

**It is a quantity/comparability problem**, structurally identical to Arm A: the nominal
ASO count is inflated by chemistry-variant duplication (finding 1) and the true
LOAO-independence unit is very likely coarser (gene/locus, finding 2) than the naive ASO
count, on top of a genuine label-definition mismatch across the studies meant to cross-check
each other (finding 4).

## Verdict: AMBER

**The single decisive gating quantity:** the number of distinct **genes/loci** (not ASO
labels) represented among Holgersen's 81 SBOs, and the per-gene ASO count distribution —
neither has been extracted from Supp. Tables S4/S5 yet. If that count is ≥25 with
enough within-cluster positives/negatives, Arm C clears a defensible power bar comparable to
Arm A's; if it is closer to the ~16-gene floor implied by "16 genes total," Arm C inherits
Arm A's exact failure mode at a smaller scale.

**Why AMBER and not GO:** unlike Arm A, there is no confirmed sub-threshold number yet —
just a plausible risk band that straddles a reasonable analogue of Arm A's bar. Calling GO
now would repeat Arm A's original mistake (asserting overlap/power before extraction).

**Why AMBER and not NO-GO:** unlike Arm A (defensible count 11, hard below 25), nothing
here is yet confirmed to fail. 81 ASOs and 16 genes are both extractable facts sitting on
either side of a reasonable bar, and the extraction cost is two PDF supplementary tables —
cheap relative to the 5,000–20,000 CPU-h this gates.

**What would revive/resolve it:** a bounded, cheap (non-sbatch) extraction pass on
Holgersen Supp. Tables S4/S5 to get (a) the actual gene list and per-gene ASO counts, (b)
whether per-site genomic coordinates/sequences are present (required for feature
computation — not yet confirmed machine-readable), and (c) the true positive/negative site
ratio per ASO. That single extraction converts this AMBER into a defensible GO or NO-GO
without touching sealed data, chemistry modeling, or any forbidden ensemble, and costs
orders of magnitude less than the Arm C build it gates.

## Honesty note

All effect and count values above are quoted from the source articles/abstracts as reported
(`paper_supported`); no sealed `terminal/` or Phase-9 ΔPSI data was read. No numeric AUCPR,
sensitivity, or independence estimate above is a measured result of this project — all are
either literature-reported (`paper_supported`) or explicit risk bounds pending extraction
(`candidate`), per `feng/src/feng/honesty/states.py`'s allowed states. No chemistry model
was invoked or implied; the two chemistry-variant sources (Scharner, GSE307431) are treated
strictly as external, direction-only replication, consistent with PLAN.md's existing
do-nots — this memo did not need to relax them.

## Sources

- Holgersen et al. 2021, *Nucleic Acid Therapeutics* 31(6):392-403. doi:10.1089/nat.2020.0921 ; https://dx.doi.org/10.1089/nat.2020.0921 ; bioRxiv 10.1101/2020.09.03.281667 ; PMC8713556
- Scharner et al. 2020, *Nucleic Acids Research* 48(2):802-816. PMID 31802121 ; https://academic.oup.com/nar/article/48/2/802/5658446
- Ottesen et al. 2026, *NAR Molecular Medicine* 3(1):ugag002. PMC12805893 ; https://academic.oup.com/narmolmed/article/3/1/ugag002/8415575 ; GEO GSE307431 ; SRA PRJNA1321660
- "Determining off-target effects of splice-switching antisense oligonucleotides using short
  read RNAseq in neuronally differentiated human induced pluripotent stem cells", 2025,
  *Human Molecular Genetics* 34(22):1912-. PMC12581821 ; bioRxiv 10.1101/2025.01.09.631923
- Hua 2007 (doi:10.1371/journal.pbio.0050073), Hua 2008 (PMC2427210), Singh 2010
  (rnajournal.cshlp.org/content/16/6/1167) — as harmonized in
  `feng/artifacts/phase0_arm_a_overlap/harmonized_aso_positions.csv`
- `feng/PLAN.md` — Arm C section (lines 114-136), Arm C build gate (line 183), do-nots
  (lines 275-291)
