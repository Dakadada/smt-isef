# Data & Construct Request — Computational Splicing Model

Sending this before any bench work starts so the computational model is built
against the exact molecules you assay, and so the data lands in a form that
needs no rework later.

Split into four asks. **Ask 1 and Ask 2 unblock work immediately.** Ask 3 is the
format contract for data that arrives later. Ask 4 is optional and only matters
if we pursue the circRNA arm.

---

## Ask 1 — Which construct, and its exact sequence

The model computes RNA secondary structure and effector binding energetics
directly from sequence. A construct with truncated introns is a different
molecule from the genomic locus, so predictions have to be run on the actual
plasmid sequence or they will not match your gels.

Please confirm:

1. **Which construct will be used for the ASO experiments?**
   - `SMN2Sup3445` (or another `SMN2Sup` variant — `SMN2Sup757T`, `SMN2Sup-CMV`, etc.)
   - the classic pCI-based exon 6–8 minigene
   - something else

2. **The full sequence of that construct**, in whichever form is easiest:
   - the NCBI GenBank accession — Ottesen et al. 2024 (NAR 52:3547) states super
     minigene sequences were deposited and lists accessions in Supplementary
     Table S2; I have not been able to retrieve S2, so the accession itself is
     enough
   - a SnapGene/GenBank/FASTA file of the plasmid
   - or the fragment junction coordinates plus restriction sites, which is
     enough to reconstruct it against the genomic reference

3. **Confirmation on intron 7.** The paper states "the entire introns 3 and 7
   were used." I am assuming intron 7 is full-length and therefore that the
   ISTL1 long-distance interaction is intact. Please confirm — if intron 7 is
   truncated in the variant you actually run, the structural model changes
   substantially.

4. **Reference SMN2 allele/haplotype** used to build the construct, if it differs
   from the genome reference at any position in exons 6–8 or introns 6–7.

---

## Ask 2 — Any existing PSI data you can share now

Model calibration does not need new experiments. Any prior or published
inclusion data on this construct — ASO dose-response, length series, mutant
panels — lets calibration start immediately and keeps the bench off the
critical path.

Specifically useful:

- ASO length series at ISS-N1 (including the 8-mer / 3UP8 series if available)
- Nusinersen dose-response on the same construct
- any mutant minigene panel with measured exon 7 inclusion
- untreated / mock baseline inclusion for the construct

Published is fine — I mainly need to know which figures correspond to which
construct, since that mapping is often ambiguous from the paper alone.

---

## Ask 3 — Format for PSI data as it is generated

Rather than converting gel images or summary percentages after the fact, one
row per lane, in a CSV or spreadsheet. Requested columns:

| Column | Meaning |
|---|---|
| `experiment_id` | Free-text ID grouping one gel / one experiment |
| `gel_id` | Identifier of the physical gel |
| `lane` | Lane number on that gel |
| `construct_id` | Construct name, matching Ask 1 |
| `cell_line` | e.g. HeLa, SH-SY5Y, patient fibroblast line |
| `effector_id` | ASO/compound name, or `none` for untreated |
| `effector_sequence` | 5'→3' sequence, or blank for small molecules |
| `effector_chemistry` | e.g. `MOE-PS`, `OMe-PS`, `mixed`, per-position if mixed |
| `effector_length` | nt |
| `concentration_nM` | Final concentration; `0` for untreated |
| `treatment_hours` | Time from treatment to harvest |
| `replicate` | Biological replicate number |
| `technical_replicate` | Technical replicate number, if any |
| `band_included_intensity` | Raw densitometry, exon-7-included band |
| `band_skipped_intensity` | Raw densitometry, exon-7-skipped band |
| `background_intensity` | Local background for that lane |
| `saturated` | `TRUE`/`FALSE` — any band pixel-saturated |
| `notes` | Anything unusual about the lane |

**Raw band intensities matter more than the computed percentage.** With raw
values plus background, uncertainty can be propagated properly and the
normalization can be re-derived if conventions change. A pre-computed PSI
column is welcome in addition, but not as a replacement.

Also needed once, not per row:

- **Background subtraction method** — local rolling-ball, fixed rectangle, or none
- **Gel percentage** and whether polyacrylamide or agarose
- **Imager and software** used for densitometry
- **Whether intensities are linearized** or raw 8/16-bit pixel values
- **Primer sequences** for the RT-PCR readout, and the expected amplicon sizes
  for included and skipped products
- **RT priming strategy** — random hexamer, oligo-dT, or gene-specific
- **Cycle number**, and whether the reaction was checked to be in the linear range
- **Transfection reagent and method** (affects effective intracellular concentration)

A single filled-in example row for one lane would resolve most ambiguity faster
than any amount of specification.

---

## Ask 4 — circRNA arm (only if we pursue it)

Relevant only if the C2B34 / backsplicing work goes forward:

- **Divergent primer sequences** used to detect the backspliced junction
- **RNase R treatment** — used or not; if used, the digestion conditions
- Confirmation that RT priming avoids oligo-dT for these measurements, since
  circles are not polyadenylated
- Whether C2B34 quantification is from the same RNA prep as the linear
  measurements, or a separate one
- Any existing linear-to-circular ratio measurements, treated and untreated

---

## What is being built, for context

A biophysical model that predicts, from sequence alone: how strongly an effector
binds a target site, how much that binding opens or blocks the exon 7 5′ splice
site, and how that maps to inclusion. It is calibrated against measured
inclusion, then used to predict effects at sites that have not been tested and
to screen for off-target splicing changes elsewhere in the transcriptome.

The model can be built and tested against published data before any new bench
data exists. Nothing here is blocking on new experiments — Ask 1 is the only
item that gates real progress, since running the model on the wrong sequence
would invalidate everything downstream.
