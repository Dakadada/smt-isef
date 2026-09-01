# Phase 0 — Arm A (StericScope) overlap feasibility probe

**Status: PARTIAL (data-source landscape + execution plan complete; the actual ≥25-window
overlap count is NOT yet computed — it requires the extraction work described below).**
Author: Claude (Codex hit its usage limit mid-probe; this is the Claude-side completion of
the source-scoping half). Date 2026-09-01.

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

## Preliminary verdict: **AMBER (leaning green on raw feasibility, GATED on extraction)**

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
