# Phase 0 — Arm A (StericScope) overlap feasibility probe

**Status: FINAL — RED under the defensible independence definition (11 independent
recoverable windows < 25 kill threshold). The earlier exact-interval GREEN call was
retracted by adversarial review 2026-09-01 as claim inflation; see the "Adversarial review"
log at the end of this file. Per user decision (2026-09-01) the RED was **CONTESTED** — and
the contest LOST honestly: pairwise-non-overlap is indeed the wrong independence unit (the
tiling walk has effective dimension ~12, not 11), but a proper spatial effective-N is only
**~11 (< 25 frozen bar)** and the held-out transfer Spearman has ~19% power at the frozen
ρ=0.5 target (needs an implausible ρ≈0.85 for 80%), with **0 analysis-ready numeric ASO
outcomes on disk**. See `CONTEST.md`. **RED STANDS → Arm A is DESCRIPTIVE** (a cross-effector
comparison), not a held-out transfer test, on the available overlap. The only path back:
add genuinely new, independent numeric ASO data.**
Author: Claude (Codex hit its usage limit mid-probe; this is the Claude-side completion of
the source-scoping half). Date 2026-09-01.

> NOTE (2026-09-01, post-review): the exact-interval GREEN interpretation is retracted and
> superseded by the RED verdict. The logged `+9` coordinate-frame error and incomplete
> Hua-2008 intron-6 microwalk have been corrected below; see resolved findings 4 and 5.

## The question

Arm A's headline claim needs the SpliceRUSH dCas13d guide-tiling field and the historical
SMN2 ASO walks to co-locate on ONE coordinate frame densely enough to yield **≥25
independent recoverable ASO windows under guide-covered field**. GREEN → Arm A buildable;
thin overlap → Arm A degrades to a descriptive cross-effector comparison up front.

## Data sources (verified, with availability)

| Source | What it gives | Accession / locator | Availability |
|---|---|---|---|
| **SpliceRUSH** (Recinos et al. 2024, *Nat Commun*) | dCas13d gRNA tiling of SMN2 incl. distal introns; per-guide effect on exon-7 inclusion; **also ~3 ASOs at the distal SRE ISE-D1** | Paper PMC11076525 / doi:10.1038/s41467-024-47140-8; raw reads **SRA PRJNA1052769**; Supplementary Data 1–4 (.zip) | Public. **BUT per-guide effect scores are NOT tabulated** — only raw read counts (Suppl Data 4) + figure plots. gRNA + ASO sequences are in **Suppl Data 1**. |
| **Hua 2007** (*PLoS Biology*) | ASO **exon walk** across SMN2 exon 7 (MOE), inclusion/skipping/complex effects | doi:10.1371/journal.pbio.0050073 | Open access (PLoS). Positions/effects in figures + supplement → digitization needed. |
| **Hua 2008** (*AJHG*) | ASO **intron tiling** of intron 6/7 flanks; found the ISS (drug site) | PMC2427210 | Open access. Supplement tables → digitization needed. |
| **Singh 2010** (*RNA*) | ASO **microwalk** across ISS-N1 (single-nt-stepped, incl. 8–15-mers), MOE | rnajournal.cshlp.org/content/16/6/1167 | Open access. Positions/effects in figures/supplement → digitization needed. |

## Key feasibility findings (why this is non-trivial)

1. **SpliceRUSH per-guide effects are not a table.** They must be RE-DERIVED from the raw
   SRA reads (PRJNA1052769) through a guide-count → enrichment/Z-score pipeline, or
   digitized from figures. This is exactly the Phase-3 board item "SpliceRUSH guide-level
   effects separated from reporter enrichment" — it is real work, not a download.
2. **Coordinate frames differ three ways.** SpliceRUSH reports guide/ASO positions
   **relative to the splice reporter construct** ("N nt from the 5′ ss") and used **hg19**
   for RNA-seq; the historical walks use their own SMN2 exon-7 / intron-7 numbering (e.g.
   ISS-N1 = intron-7 positions 10–24). All must be harmonized onto ONE frame (recommend
   SMN2 genomic, GRCh38 chr5 or RefSeqGene NG_008691), handling SMN1/SMN2 near-identity
   carefully so guide/ASO positions are not mis-assigned to the paralog.
3. **The ≥25 windows come from the historical walks, not SpliceRUSH's ASOs.** SpliceRUSH
   itself reports only ~3 ASOs (ISE-D1, inclusion 81%→36–48.5%). The ASO count must come
   from Hua 2007 (exon) + Hua 2008 (intron) + Singh 2010 (microwalk), all of which need
   supplement digitization.

## UPDATE (Claude, 2026-09-01) — dCas13d field EXTRACTED; verdict upgraded to GREEN-leaning

The "per-guide effects aren't tabulated" concern was **wrong** — Dataset S4 directly
provides, per guide: `OligoID`, guide seq, `Minigene` (SMN1/SMN2), `start pos`/`end pos`/
`length` (minigene coords), `Region` (E6/I6/E7/I7/boundaries), and a **precomputed combined
`Z-score`** (+ p-value, N, FDR) plus per-replicate Dual-IN/Dual-EX Top/Bot-vs-unsorted
enrichments. No enrichment computation or raw-SRA reprocessing required.

**Extracted** to `$SCRATCH/feng/phase0/splicerush_dcas13d_field.csv` (1,937 guides;
**1,625 SMN2**). Region coverage:

| Region | guides | Region | guides |
|---|---|---|---|
| I7 (intron 7) | **842** | E6I6 | 45 |
| I6 | 353 | E7I7 | 45 |
| E6 | 195 | I7E8 | 35 |
| E7 | 62 | (None) | 312 |
| I6E7 | 48 | | |

**1,032 guides tile the E7/I7 neighborhood** (minigene pos 296–835) — exactly the region
the historical ASO walks cover (Hua 2007 exon 7; Hua 2008 intron 6/7; Singh 2010 ISS-N1 in
intron 7). Guide coverage over the ASO-walked field is therefore DENSE. Given 842 intron-7
guides alone, ≥25 ASO windows falling under guide-covered field is very likely to clear once
the ASO positions are mapped in.

**Remaining before a final GREEN:** digitize the Hua/Singh ASO walk positions+effects and
harmonize them to the S4 minigene coordinate frame (map exon-7/intron-7 numbering ↔ minigene
`start pos`), then count independent ASO windows under guide coverage. This is the only open
step; the dCas13d half is done.

## Original verdict (pre-extraction): AMBER — superseded by the GREEN-leaning update above

- The two modalities DO share the SMN2 exon-7 ± flanking-intron region, and SpliceRUSH
  tiles it (plus distal introns), so guide coverage over the ASO-walked region is very
  likely.
- The combined historical ASO inventory (exon walk + intron tiling + ISS-N1 microwalk) is
  *plausibly* well above 25 positions — but this is a **literature-based expectation, not a
  measured count**. The number of *independent, non-redundant* ASO windows that fall under
  dCas13d-guide-covered field CANNOT be asserted without doing the extraction + mapping.
- **No number is claimed here.** Honest status: the overlap is promising enough to proceed
  to extraction; GREEN (≥25 confirmed) is unproven until the steps below run.

## Data staged (Claude, 2026-09-01) — under `/home/dashii/orcd/scratch/feng/phase0/`

Downloaded from Springer static-content (the CDN URLs a fetch suggested 404'd; canonical
host is `static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-47140-8/MediaObjects/`):
- `moesm4.zip` → `supp_1-4/Dataset S1–S4.xlsx` (Suppl Data 1–4).
  - **Dataset S1** → sheet `Oligos` (44 rows): SMN2 plasmid parts + gRNA + ASO sequences
    (5′-3′, Note, Length); sheet `Illumina_Primers`.
  - Dataset S2 = gRNA library oligos; S3 = cell/coverage counts; **S4 = read-count
    summaries (1.2 MB)** — the substrate for per-guide enrichment.
- `moesm6.xlsx` (36 MB) = **Source Data** (values behind every figure incl. Fig 6 dCas13d
  + ASO validation) — likely holds the plotted per-guide effects; large, parse offline.
- **`openpyxl` was missing from feng-env** — added to `env/environment.yml` and installed
  into the live env; required to read these xlsx files.

## Execution TODOs (hand to Codex on usage reset ~9:25 PM, or resume Claude-side)

1. [x] Download SpliceRUSH Suppl Data 1–4 + Source Data → staged in `$SCRATCH/feng/phase0/`
   (see above). Structure confirmed.
2. Extract per-guide effect scores. **Refined (Claude inspected the xlsx):** the Source
   Data (`moesm6.xlsx`, 9 sheets: Fig 2/3/6 + supp) are **validation gel quantifications**
   (Fig 2b = nusinersen ASO dose-response 0–80 nM; Fig 6b/c = dCas13d+ASO validation at
   ISE-D1), NOT the 1,937-guide screen scores. So the dCas13d **field must be re-derived
   from `Dataset S4.xlsx`** (single sheet, **1941×70** = ~1,937 guides × count columns for
   the FACS-sorted bins) via a sorted-bin enrichment calc — a standard computation, **no
   raw-SRA reprocessing required.** ASO seqs/positions come from Dataset S1 `Oligos`
   (44 rows: plasmid parts + gRNAs + the ~3 ISE-D1 ASOs). Keep provenance; clean-room
   (no `pre/`/`src/` logic).
   - Cross-effector validation data DOES exist in Source Data (dCas13d vs ASO at matched
     ISE-D1 sites, Fig 6) but is small-n — useful as a sanity check, not the ≥25 windows.
3. Digitize the Hua 2007 exon walk, Hua 2008 intron tiling, and Singh 2010 microwalk ASO
   tables (position + measured effect + chemistry) from the open-access supplements.
4. Harmonize all positions onto one SMN2 genomic frame (state the reference); map guides and
   ASOs; resolve SMN1/SMN2.
5. Count independent ASO windows under guide-covered field; report count + overlap span +
   density vs the ≥25 threshold → finalize GREEN / AMBER / RED.

## Honesty note

Positions, counts, and the verdict above are source-scoping and literature-based
expectations. No overlap count has been measured. Any final Arm-A go/no-go must cite the
computed count from the TODO steps, not this memo's expectation.

## Sources

- SpliceRUSH: https://pmc.ncbi.nlm.nih.gov/articles/PMC11076525/ ; SRA PRJNA1052769
- Hua 2007: https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0050073
- Hua 2008: https://pmc.ncbi.nlm.nih.gov/articles/PMC2427210/
- Singh 2010: https://rnajournal.cshlp.org/content/16/6/1167

## ASO-walk harmonization — 2026-09-01

### Result

The harmonized inventory contains **125 published ASO records**, comprising **115 distinct
exact target intervals** after removing ten exact cross-paper duplicates. All **115/115**
distinct intervals are fully contained in at least one SMN2 guide interval in
`splicerush_dcas13d_field.csv`. These exact intervals are an inventory, not an independence
count. The Phase-0 overlap gate remains **RED** because the maximum pairwise-non-overlapping
count is **11**, below the required 25 independent recoverable windows.

The row-level artifact is
[`harmonized_aso_positions.csv`](./harmonized_aso_positions.csv). It was intended for
`$SCRATCH/feng/phase0/harmonized_aso_positions.csv`, but this session's sandbox permits
reading, not writing, `/home/dashii/orcd/scratch`; the repository artifact is the complete
machine-readable deliverable and remains uncommitted.

### Coordinate frame and transforms

The S4 field coordinate frame includes a 9-nt 5′ prefix (`ATGGCTAGC`) before the Dataset S1
exon-6 sequence, followed by 111 nt exon 6, 199 nt shortened intron 6, a 55 nt exon 7
containing an extra **G at construct exon position 37**, then 444 nt intron 7. Consequently
exon 7 starts at S4 field position 320 and intron 7 starts at position 375. The +9 offset was
confirmed by direct target-sequence alignment for guides 1001, 1340, and 1384. The
historical papers use local splice-feature coordinates, not hg19 coordinates.

| Published frame | Explicit transform to inclusive S4 minigene coordinates | Check |
|---|---|---|
| Hua 2008 intron 6, with `-1` immediately upstream of exon 7 | `m = 320 + q`; e.g. `-60..-46 -> 260..274`, `-15..-1 -> 305..319` | S1 intron-6/exon-7 junction plus the 9-nt S4 prefix |
| Hua 2007 natural 54-nt exon 7, `+1..+54` | for `p <= 36`, `m = 319 + p`; for `p >= 37`, `m = 320 + p` | The S1-only inserted G is minigene 356; a historical ASO crossing natural 36/37 maps as two segments with a one-base hull |
| Hua 2008 or Singh 2010 intron 7, `+1` immediately after exon 7 | `m = 374 + q`; e.g. `+1..+15 -> 375..389`, `+10..+29 -> 384..403` | S1 exon-7/intron-7 junction and reverse-complement sequence alignment, plus the 9-nt S4 prefix |
| hg19 | **No transform used.** The SpliceRUSH hg19 reference applies to its RNA-seq analysis, whereas S4 `start pos` is reporter-construct-relative. | Avoids SMN1/SMN2 paralog ambiguity and an erroneous genomic detour |

For targets spanning the artificial exon-7 G, the CSV records both the inclusive coordinate
hull and the two native target segments. Guide containment is evaluated on the more
conservative hull, so the inserted base cannot inflate the overlap call.

### Published ASO positions and effects

Effects below are transcribed only where the article text or caption states a direction or
magnitude. `figure_reported` in the CSV means that the published gel/diagram contains the
per-ASO result but the indexed figure does not expose a safe numeric value; no value was
invented. Effects are descriptive context and do not enter the overlap gate.

#### Hua 2007, PLoS Biology — exon-7 walk

Source: [Hua et al. 2007](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0050073),
Table 1 for every target interval; Figures 1B/2 for the initial walk and Figures 3B-C,
5A-D, and 9 for the microwalk effects.

| Walk | Original exon-7 target intervals | S4 span | Published effect summary |
|---|---|---|---|
| Initial, 9 ASOs | 01-15, 06-20, 11-25, 16-30, 21-35, 26-40, 31-45, 36-50, 40-54 | 320-374 | 06-20 and 36-50 increase inclusion; 01-15/16-30/21-35/26-40 strongly decrease; 11-25/31-45 decrease; 40-54 neutral. Figure 3 gives 45% inclusion for 06-20 at 100 nM. |
| Region A microwalk, 17 ASOs | 08-22, 07-21, 05-19, 04-18, 03-17, 02-16, 05-20, 12-23, 11-22, 10-21, 09-20, 08-19, 07-18, 06-17, 05-16, 04-15, 03-14 | 321-342 | 07-21 reaches 62% inclusion at 100 nM (EC50 16 nM). Figures 5/9 record the remaining directions; documented exceptions and assay discordance are retained row-wise in the CSV. |
| Region B microwalk, 22 ASOs | 36-53, 35-52, 34-51, 33-50, 32-49, 39-53, 38-52, 37-51, 35-49, 34-48, 33-47, 32-46, 42-53, 41-52, 40-51, 39-50, 38-49, 37-48, 36-47, 35-46, 34-45, 33-44 | 355-373 (native segments; hulls through 374 where applicable) | 34-48 robustly increases inclusion (EC50 13 nM); all five 18-mers inhibit to varying degrees; 39-53/38-52 are neutral; several targets are assay-discordant as recorded in the CSV. |

#### Hua 2008, AJHG — flanking-intron walks

Source: [Hua et al. 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2427210/),
Table 1 and Figure 1 for every interval, Figure 2A-C for the initial effects, and Figure 7
for both intron microwalks.

| Walk | Original target intervals | S4 span | Published effect summary |
|---|---|---|---|
| Intron 6 initial, 10 ASOs | -60..-46, -55..-41, -50..-36, -45..-31, -40..-26, -35..-21, -30..-16, -25..-10, -20..-6, -15..-1 | 260-319 | -55..-41 nearly doubles inclusion; the four junction/polypyrimidine/branch-point targets -30..-16 through -15..-1 strongly inhibit. Other lane values remain figure-reported. |
| Intron 6 optimization, 17 ASOs | -51..-37, -52..-38, -53..-39, -54..-40, -56..-42, -57..-43, -58..-44, -59..-45; -48..-37, -49..-38, -50..-39, -51..-40, -52..-41, -53..-42, -54..-43, -55..-44, -56..-45 | 261-283 | Figure 7 reports the direction and relative strength of every microwalk ASO; where indexed text does not safely expose a per-ASO value, the CSV uses `figure_reported`, matching the intron-7 evidence standard. |
| Intron 7 initial, 10 ASOs | +1..+15, +6..+20, +11..+25, +16..+30, +21..+35, +26..+40, +31..+45, +36..+50, +41..+55, +46..+60 | 375-434 | +11..+25 gives the strongest initial increase; +16..+30 also increases; +1..+15 strongly inhibits. |
| Intron 7 optimization, 21 ASOs | 08-25, 09-26, 10-27, 11-28; 15-29, 14-28, 13-27, 12-26, 10-24, 09-23, 08-22, 07-21; 18-29, 17-28, 16-27, 15-26, 14-25, 13-24, 12-23, 11-22, 10-21 | 381-403 | Figure 7 identifies all four 18-mers as strongest, with 10-27 slightly best; 09-23 is the best 15-mer and 10-21 the best but considerably weaker 12-mer. |

#### Singh 2010, RNA — ISS-N1 ultra-refined microwalk

Source: [Singh et al. 2010](https://rnajournal.cshlp.org/content/16/6/1167),
Table 1 for ASO sequences and Figure 1A-C for coordinates and effects. Target coordinates
were recovered by reverse-complement alignment to the numbered intron-7 sequence in
Figure 1A, not inferred from ASO names alone.

| ASOs | Original intron-7 target intervals | S4 span | Published effect summary |
|---|---|---|---|
| Anti-N1; F15; 1UP15, 2UP15, 3UP15; 1DN15, 2DN15, 3DN15 | 10-29; 10-24; 9-23, 8-22, 7-21; 11-25, 12-26, 13-27 | 381-403 | Figure 1B defines green as exon-7 skipping <=45%, gray as no effect, and pink as skipping >=55%. ASOs sequestering intron-7 10C are stimulatory controls; downstream-only ASOs are non-stimulatory. |
| F14, L14; F13, M13, L13; F12, M12F, M12L; F11, F10, F9 | 10-23, 11-24; 10-22, 11-23, 12-24; 10-21, 11-22, 12-23; 10-20, 10-19, 10-18 | 384-398 | F14 strongly increases inclusion; the one-base-shifted L14 produces the strongest skipping phenotype. L13/M13/M12F/M12L have reduced or no negative effect relative to L14. |

### Overlap counting method

1. Keep only published, sequence-addressable ASOs from the three walks; exclude unrelated,
   mismatched, mutant-specific, and chemistry-only replicate controls.
2. Harmonize to inclusive S4 coordinates with the transforms above.
3. Deduplicate exact `(region, natural start, natural end)` target intervals across papers.
   Overlapping but nonidentical stepped targets remain distinct perturbation windows, which
   is the intended meaning of a walk and of the gate's “distinct target sites.” This yields
   48 exon-7 + 27 intron-6 + 40 intron-7 = **115** windows.
4. Call covered only when the complete ASO hull is contained in at least one SMN2 S4 guide
   interval (`guide.start_pos <= ASO.start` and `guide.end_pos >= ASO.end`). Result:
   **115/115 covered**. Base-level union coverage has no gaps over minigene 260-434.

A stricter sensitivity analysis that requires every counted ASO interval to be pairwise
non-overlapping gives only **11** windows (intron 6: 4; exon 7: 3; intron 7: 4) and would be
RED. That rule discards the single-nucleotide and five-nucleotide steps that define an ASO
walk. Programmatic recomputation after both corrections reconfirms **11** (intron 6: 4;
exon 7: 3; intron 7: 4), so the verdict remains **RED**.

### Assumptions and direction of risk

| Assumption | Could it flip the exact-interval inventory count? | Direction if wrong |
|---|---|---|
| “Independent/distinct” means unique exact target intervals, not pairwise-disjoint intervals. | **Yes.** Pairwise-disjoint sensitivity is 11. | Decreases the 115-window inventory to 11 independent windows; verdict RED. |
| The S1 exon-7 extra G is after natural exon position 36, as stated and sequence-aligned. | Not with observed dense coverage; all conservative hulls remain contained. | A wrong insertion placement could only decrease defensible exon-window overlap near the insertion. |
| Hua intron-6 negative positions count backward from the exon-7 acceptor; Singh/Hua intron-7 `+1` is the first intronic base. | A large frame/orientation error could; one-base errors cannot because coverage is gapless. | Misorientation would decrease or invalidate mapped windows, never increase the defensible count. |
| hg19 is not an intermediate frame for S4. | Forcing hg19 could create an SMN1/SMN2 paralog assignment error. | Decreases the defensible count through ambiguous or rejected mappings. |

**Final Arm-A verdict: RED.**
The defensible count of *independent* recoverable windows is **11** (pairwise non-overlapping;
intron 6: 4, exon 7: 3, intron 7: 4), below the ≥25 kill threshold. 115 counts 1-nt/5-nt tiling
shifts of the same 15-mers as separate observations — that is inventory, not independence.
See the adversarial review log below.

---

## Adversarial review — 2026-09-01 (Codex gpt-5.6-sol, structural)

Verdict: **RED**. Standing count to defend: **11 independent recoverable coordinate windows**
(itself an upper bound — pairwise non-overlap prevents nucleotide reuse but does not prove
statistical independence among measurements from the same short locus/study/assay).

Findings (severity-tagged; logged, not all yet fixed):

1. **[blocker] Exact intervals are an inventory count, not an independence count.** The gate requires
   "≥25 *independent* recoverable ASO windows"; finalization silently redefined independence
   as unique exact endpoints. 1-nt-shifted 15-mers share 14/15 bases and probe the same local
   field. Recomputed max pairwise-disjoint packing = 11 (I6:4, E7:3, I7:4). 11 < 25 → kill
   criterion met. *(Resolved in verdict above.)*
2. **[blocker] "Preregistered distinct-site interpretation" is unsupported / post-hoc.**
   Neither HANDOFF nor the frozen plan defines independence as exact-interval uniqueness; the
   plan says only "independent recoverable ASO windows" and schedules formal preregistration
   for month 8. The exact-interval rule first appears here, where it flips RED→GREEN. That is
   definition-shopping — an honesty-layer violation. *(Corrected: claim retracted above.)*
3. **[blocker] "Recoverable" effects are not recovered enough for the planned test.** Arm A's
   endpoint needs within-study standardized effects + a Spearman correlation. CSV holds
   qualitative text: 43/108 rows are `effect_direction=figure_reported`; effects declared
   "descriptive context" that don't enter the gate. So even 11 establishes only *coordinate*
   recoverability, not 11 usable standardized outcomes. **OPEN.**
4. **[major, RESOLVED] S4 coordinate transform was off by +9 nt.** Field CSV has a 9-nt
   prefix (`ATGGCTAGC…`) omitted from the listed S1 exon-6 sequence. Direct sequence
   alignment confirmed guide 1001 contains the prefix before S1 exon 6, guide 1340 aligns at
   S1-concatenated 331–355 versus field 340–364, and guide 1384 aligns at 375–399 versus
   field 384–408. All ASO coordinates, segments, and keys were shifted into the S4 field
   frame; recomputed exact-interval coverage is 115/115.
5. **[major, RESOLVED] Hua 2008 intron-6 microwalk completed.** Added all 17 published
   intron-6 microwalk ASOs from Table 1 / Figure 7 using the same figure-level evidence
   standard as the intron-7 microwalk. Inventory is now 125 records / 115 exact intervals;
   programmatic interval scheduling reconfirms the pairwise-disjoint maximum is 11 and the
   verdict remains RED.
6. **[major] Exact-window dedup hides incompatible study-specific outcomes.** 10 coords appear
   in both Hua 2008 and Singh 2010; at `I7:+11..+25` Hua labels strongest stimulation while
   Singh's 1DN15 is non-stimulatory/inhibitory (different chemistries/contexts). An exact
   interval is not one exchangeable observation — Arm A needs a preregistered study-specific
   outcome/aggregation model, not coordinate-only dedup. **OPEN.**
7. **[minor] Spot-checked source positions and headline directions otherwise supported** (Hua
   2007 initial + microwalks; Hua 2008 included intervals; Singh 2010 reverse-complement
   alignments; the artificial-G placement). The 43 `figure_reported` directions remain
   unverified rather than confirmed.

Implication for the plan: Arm A's overlap kill-gate is **RED**. Per the plan's pre-flight
logic ("thin overlap → Arm A degrades to a descriptive cross-effector comparison up front"),
StericScope cannot be built as a held-out transfer-law fit on this overlap. Arm D (locked
nusinersen benchmark) remains the guaranteed floor.
