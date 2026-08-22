# A5 validation target — GEO accession lookup, and what did NOT pan out

Written while building the A5 transcriptome-wide sieve
(`src/smt/sieve/{transcriptome,scan,polarity}.py`). This documents exactly
what was found, including a discrepancy in the task's own framing, rather
than papering over it.

## The requested paper: PMC8395096

- **Title**: "High Concentration of an ISS-N1-Targeting Antisense
  Oligonucleotide Causes Massive Perturbation of the Transcriptome"
- **Authors**: Ottesen EW, Luo D, Singh NN, Singh RN
- **Journal**: *International Journal of Molecular Sciences* 22(16):8378, 2021
- **DOI**: 10.3390/ijms22168378
- **PMID**: 34445083 / **PMCID**: PMC8395096
- Design: SMA patient fibroblasts treated with the ISS-N1-targeting ASO
  ("anti-N1") at two doses: **HiA** = 100 nM (the "high concentration" of the
  title), **LoA** = 5 nM. RNA-seq, differential expression + differential
  splicing analysis (rMATS-style event calling: exon inclusion, exon
  skipping, intron retention, increased intron removal, altered 5'/3'
  splice sites, mutually-exclusive exons).

### GEO accession: THERE ISN'T ONE

Checked directly (`Bash: curl` on the raw PMC HTML, not just the WebFetch
summarizer, because the first pass truncated the article and missed this):
the paper's own **Data Availability Statement** reads, verbatim:

> "All data generated or analyzed during this study are included in this
> published article and its supplementary information files."

No GEO series, no SRA/BioProject accession, nothing on NCBI's GDS database
links to this PMID (checked via `elink dbfrom=pubmed db=gds id=34445083`:
empty linkset). The processed splicing calls live ONLY in the article's own
supplementary materials (a single 2.3 MB zip attached to the PMC record),
specifically:

| Table | Contents |
|---|---|
| S1 / S2 | Up/downregulated genes, HiA treatment |
| S3 / S4 | Exon inclusion (EIN) events, HiA / LoA |
| S5 / S6 | Exon skipping (ESK) events, HiA / LoA |
| S7 / S8 | Intron retention (IRT) events, HiA / LoA |
| S9 / S10 | Increased intron removal (IRM) events, HiA / LoA |
| S11 / S12 | Altered 5' splice site (A5S) events, HiA / LoA |
| S13 / S14 | Altered 3' splice site (A3S) events, HiA / LoA |
| S15 / S16 | Mutually exclusive exon (MXE) events, HiA / LoA |
| S17 | Primers used for RT-PCR validation of alternative splicing |

**Tables S5/S6 (exon skipping) and S3/S4 (exon inclusion) at the HiA (100 nM)
dose are the ones that would grade this sieve** (see comparison plan below).

### POLR2H is not in this paper

Also checked directly: the string "POLR2H" does not appear anywhere in
PMC8395096's full text. This paper reports thousands of genome-wide
splicing events in aggregate (per the table list above) but does not name
POLR2H as a case study. **The task description's specific claim -- "it
induces exon 2 skipping in POLR2H by binding an exonic splicing enhancer" --
does not come from this paper.**

## Where the POLR2H claim actually comes from: a newer, related GEO series

Searching GEO directly (`esearch db=gds`, not just PubMed-linked records)
for the same lab's more recent work turned up:

- **GSE307431** — "A therapeutic antisense oligonucleotide encompassing
  2'-O-methoxyethyl modification triggers unique perturbation of the
  transcriptome." Homo sapiens, 20 samples, deposited 2026-01-06. No linked
  PMID in GEO yet as of this writing (2026-08-22) -- likely ahead of print
  or not yet cross-indexed.
- Its abstract states, verbatim: *"Using minigenes, we showed that skipping
  of POLR2H exon 2 and seven other exons is triggered by direct annealing of
  F18MOE to target exonic sequences, and the chemistry-specific effects of
  F18MOE on exon skipping are primarily due to enhanced mismatch tolerance of
  MOE."* F18MOE = an 18-mer, 2'-MOE/phosphorothioate ASO targeting ISS-N1 --
  i.e. nusinersen's own chemistry and length. This is almost certainly the
  correct source for the POLR2H claim in the task description, from the same
  Singh lab, comparing three ISS-N1-targeting ASO chemistries (F18OMe,
  F18MOE, F20PMO) rather than one HiA/LoA dose series.

**This is disclosed as a finding, not resolved further**: this task did not
download or align either dataset's RNA-seq (out of scope, per the task), so
GSE307431's exact POLR2H fold-change / PSI numbers and their supplementary
table locations are not yet catalogued here. That is the natural next step
if POLR2H specifically (rather than the aggregate HiA/LoA event lists) is
the grading target.

## The comparison that would actually grade this sieve

Two comparisons, at two levels of rigor, against the two datasets above:

### 1. Aggregate rank correlation (PMC8395096, Tables S5/S6 + S3/S4, HiA)

1. Download Supplementary Tables S3-S6 (HiA columns) from the PMC8395096
   supplementary zip. Each row is a gene/exon with an rMATS-style
   inclusion-level difference (deltaPSI) and significance.
2. Run this sieve's full pipeline (`smt.sieve.scan` then
   `smt.sieve.polarity`) for the SAME ASO sequence and dose regime this
   paper used (the ISS-N1-targeting "anti-N1" ASO -- confirm its exact
   chemistry/sequence against the paper's Methods before treating it as
   equivalent to nusinersen; sequence may match but chemistry differs from
   the approved 2'-MOE drug).
3. Restrict to genes/exons this sieve actually found candidate sites for
   (this sieve found 37 transcriptome-wide at default seed threshold for
   nusinersen -- see the run log/report from this task -- so N will be
   small; that is a REAL constraint on statistical power, not a design
   flaw to paper over).
4. Compute Spearman rank correlation between this sieve's `signed_effect`
   (from `smt.sieve.polarity`) and the paper's measured deltaPSI, across
   the matched gene/exon set.
5. **Sign accuracy**: of the matched set, what fraction have
   `sign(signed_effect) == sign(deltaPSI)`? Report this SEPARATELY from
   rank correlation -- a sieve that gets every sign right but ranks
   magnitudes poorly is a different (and for this project's purpose, more
   important) result than one with good rank correlation but frequent sign
   errors, given the whole point of A5 is the sign, not the magnitude.

### 2. POLR2H case study (GSE307431, once its processed tables are located)

1. Confirm F18MOE's exact sequence against nusinersen's (this project's
   `smt.sieve.scan.NUSINERSEN`) -- expected to match exactly, since "F18MOE"
   names an 18-mer, 2'-O-methoxyethyl, ISS-N1-targeting ASO, which is
   nusinersen's own description.
2. Locate GSE307431's processed splicing-event table (once downloaded) for
   the POLR2H exon 2 skipping event specifically, and its effect size.
3. Report this sieve's `signed_effect` AND rank (among all transcriptome-wide
   candidates for nusinersen) for POLR2H specifically, held out as a single
   case study distinct from the aggregate comparison above. **State the
   result honestly whether it succeeds or fails** -- per this task's own
   diagnostic run (see the task report), the DEFAULT seed-and-extend stage-1
   filter did NOT recover a POLR2H candidate for nusinersen (a supplementary
   brute-force check found a 4-mismatch exonic window, which the tiling-seed
   heuristic's documented pigeonhole guarantee only covers up to 1
   mismatch), so this case study starts from a known, disclosed limitation,
   not a clean slate.

## Explicitly NOT done in this task (by design)

- No RNA-seq from either GEO series or the PMC8395096 supplementary zip was
  downloaded, aligned, or otherwise processed.
- No PSI/deltaPSI values were used anywhere in `smt.sieve.polarity`'s
  scoring. This project has no experimental PSI data and placeholder
  calibration against any is explicitly forbidden by the task brief; the
  comparison plan above is written for whoever does that download next, not
  performed here.
