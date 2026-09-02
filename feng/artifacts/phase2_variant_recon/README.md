# Phase-2 SMN2 variant recon: two-day KILL-GATE

**Status (2026-09-02): SUPERSEDED — resolved to `NO-GO` by the staged HapSMA
extraction.** See `EXTRACTION_HAPSMA.md` (this directory) for the authoritative
verdict: 0 qualifying donor-independent SMN2 target-site ALTs on the staged public
HapSMA workbook (GO requires ≥2 unrelated donors). That is an evidence-gate result,
explicitly NOT a biological zero-variant, conservation-fraction, or patient-response
claim. The `NEEDS-EXTRACTION` narrative below is the prior recon state, retained for
the record.

**(prior, 2026-09-01): `NEEDS-EXTRACTION`.** Public, non-controlled,
copy-resolved material exists, but this recon did not measure the number of
qualifying non-reference SMN2 target-site haplotypes. No target-allele count is
claimed below.

## Claim boundary: what the Appendix could actually say

The Appendix arm asks whether a natural, copy-resolved SMN2 allele inside the
18-nt nusinersen-binding footprint changes the *predicted* action of the ASO over
and above the allele's predicted cis effect. For a fixed prediction function
`f(sequence, ASO condition)`, it would compute the four cells of a 2 x 2 design:

```text
cis effect  = f(ALT, -ASO) - f(REF, -ASO)

interaction = [f(ALT, +ASO) - f(ALT, -ASO)]
            - [f(REF, +ASO) - f(REF, -ASO)]
```

A surviving arm could therefore claim only that a naturally observed,
copy-resolved target allele is a **candidate** for allele-dependent ASO action in
the specified model. It could decompose model output into baseline cis and
variant-by-ASO terms. It could not claim a biological interaction, differential
patient response, clinical effect, or resistance to nusinersen unless the same
alleles were measured experimentally in a variant x ASO design. No public paired
variant +/- nusinersen measurements were verified in this recon.

This is kill-gated because without a natural target-site ALT on at least two
independent public SMN2 haplotypes, the four-cell exercise has no population
substrate. It would be a synthetic perturbation study, not the promised natural
variant Appendix. Even if the haplotype gate passes, the result initially remains
`candidate`; a literature-reported observation can be `paper_supported`. Nothing
in this recon qualifies as `clean_reanalysis_reproduced`.

## Qualifying evidence specification

Before looking for an ALT, the extraction must lock the footprint by all of the
following: the published 18-nt nusinersen sequence and orientation; its exact
ISS-N1 target sequence; and a strand-aware coordinate interval in every reference
frame used (at minimum the Paraphase GRCh38/SMN1-realignment frame and the HapSMA
T2T-CHM13 frame). A hit qualifies only if:

1. the ALT base or indel overlaps one of those 18 target nucleotides;
2. the complete target interval is callable on that haplotype, rather than absent
   from a variants-only table;
3. the allele is assigned to an SMN2 copy by read-backed long-read phasing across
   paralog-discriminating positions; and
4. independence is at the donor level: two SMN2 copies from one person count as
   one donor, not two independent observations.

Haplogroup labels, copy number, a locus-level VCF, a short-read call, an unphased
SMN1/SMN2 mixture, or an assembly segment with uncertain paralog identity does not
qualify.

## Public source inventory

Availability was checked at the locators below. “Usable” means usable for the
target-haplotype audit after coordinate/QC extraction, not proof that a target ALT
exists.

| Source and locator | What is public now | Copy resolution | Gate use / limitation | Honesty state |
|---|---|---|---|---|
| Chen et al., 2023, *AJHG*, DOI [10.1016/j.ajhg.2023.01.001](https://doi.org/10.1016/j.ajhg.2023.01.001), PMCID [PMC9943720](https://pmc.ncbi.nlm.nih.gov/articles/PMC9943720/); [SMN_phased_data](https://github.com/xiao-chen-xc/SMN_phased_data) | Public per-sample IGV sessions and `igv_sessions/data/<sample>.bam`, `.bam.bai`, `.bed.gz`, and `.bed.gz.tbi` haplotype-visualization files. The repository is directly downloadable without an application. | Paraphase gathers reads from both paralogs, realigns them to one SMN-family reference, phases them, and emits haplotagged bamlets; the paper describes SMN1/SMN2 classification and phased variant calling. | **Usable for extraction and non-patient conservation.** The BAM/BED frame is a Paraphase realignment frame, so the 18-nt interval must be mapped and every putative ALT checked in the haplotagged reads. There is a cohort-label discrepancy: the repository README calls the sessions “107 Coriell samples,” while the paper's data-availability statement identifies downloadable bamlets for GIAB and HPRC samples and the filenames are `HG...`. Do not treat this deposit as an SMA-patient cohort until sample provenance is reconciled. | `paper_supported` availability and method; target-site status unmeasured |
| Chen et al., 2023, public article supplements, especially Data S1 (`mmc2.xlsx`) linked from [PMC9943720](https://pmc.ncbi.nlm.nih.gov/articles/PMC9943720/#appsec1) | Public haplogroup/sample summary tables. | Derived from Paraphase calls, but the workbook is a summary rather than a guaranteed base-by-base callable target matrix. | **Supportive, not sufficient alone.** Use it to reconcile sample/haplotype labels; use read-backed BAM/BED evidence for the footprint decision. | `paper_supported` |
| Zwartkruis et al., 2025, *Genome Medicine*, DOI [10.1186/s13073-025-01448-2](https://doi.org/10.1186/s13073-025-01448-2), PMCID [PMC11927269](https://pmc.ncbi.nlm.nih.gov/articles/PMC11927269/), Additional file 5 `13073_2025_1448_MOESM5_ESM.xlsx` | The publisher supplement is public and is described as Table S3, “an overview of all identified variants per haplotype” in T2T-CHM13 coordinates. Additional file 3 gives paths for the public HPRC long reads; Additional files 4, 6, and 7 provide SMN1/2-specific positions and coordinate aids. Code is public at [HapSMA](https://github.com/UMCUGenetics/HapSMA) and [ManuscriptSMNGeneConversion](https://github.com/UMCUGenetics/ManuscriptSMNGeneConversion). | HapSMA performs polyploid long-read phasing; the paper classifies phased copies as SMN1 or SMN2 using the exon-7 paralogous sequence variant (`c.840C>T`) and reports per-haplotype Clair3 calls with depth/allele-fraction filtering. | **Best public tabular extraction candidate.** It may permit a footprint query without controlled raw patient reads, but this recon has not parsed the workbook or verified complete callability across the 18 bases. T2T-to-target coordinate reconciliation is mandatory. HPRC rows are controls, not patients. | `paper_supported` availability and method; target-site status unmeasured |
| Public HPRC ONT/HiFi reads listed in Zwartkruis Additional file 3 and the [Human Pangenome open-data registry](https://github.com/human-pangenomics/HPP_Year1_Data_Freeze_v1.0) (HPRC Genome Sequencing BioProject `PRJNA701308`) | Public long reads for the paper's selected control donors; file paths are supplied in the supplement. | Can be reprocessed by HapSMA/Paraphase for read-backed, copy-resolved haplotypes. | **Usable only as non-patient conservation/technical evidence.** HPRC assemblies are not patients and must never support a patient-response statement. Prefer the published phased table/bamlets within a two-day gate; whole-genome reprocessing is a fallback only if it fits the gate. | `paper_supported` source; any new call would require clean reanalysis |
| Zwartkruis et al. SMA-patient ONT reads; DataverseNL DOI [10.34894/G7YG0V](https://doi.org/10.34894/G7YG0V) | Metadata and an access procedure are public; raw patient long reads are legally restricted and require a request/data-sharing agreement. Some per-haplotype variant summaries are published in the article supplements. | The study's HapSMA output is copy-resolved, but the raw patient reads are not open. | **Raw reads are not usable for this gate.** Do not apply or wait. Public supplement rows may be used if they expose the necessary allele, SMN2-copy assignment, donor independence, and target callability; otherwise they cannot be rescued by the controlled data. | public summaries `paper_supported`; raw data unavailable |
| Blasco-Perez et al., 2021, *Human Mutation*, PMCID [PMC8252042](https://pmc.ncbi.nlm.nih.gov/articles/PMC8252042/) (“Beyond copy number”) | A public method paper reports SMN2-targeted long-range PCR/NGS and states that no described modifier in its study occupied ISS-N1. A later thesis reports a larger internal series, but no verified public per-copy target-haplotype deposit or accession was found in this recon. | The assay is SMN2-targeted, but this recon did not verify phasing of separate SMN2 copies through the 18-nt footprint. | **Not gate-qualifying.** A negative cohort statement is contextual only under the hard rule that every variant observation must be copy-resolved. It cannot supply either a qualifying ALT count or a definitive copy-resolved conservation count. | contextual `paper_supported`, not qualifying evidence |

### Explicitly excluded

- gnomAD, ordinary 1000 Genomes short-read calls, dbSNP, ClinVar, and similar
  locus-level catalogs: no target-region call is accepted unless it is assigned to
  a specific SMN2 copy by qualifying evidence. Donor overlap with HPRC does not
  upgrade an ordinary short-read call.
- HPRC assemblies as patient observations: disallowed. Public HPRC long reads or
  published phased outputs may support a control conservation audit only.
- The controlled Zwartkruis patient raw reads, EGAS50000001222, dbGaP/EGA data,
  proprietary reads, and author-mediated access: unavailable within this gate;
  “application submitted” is not data.
- Any reported target-site absence from an unphased mixture of multiple SMN2
  copies: not copy-resolved and therefore not decisive under this recon's rules.

## Two-day decision procedure and hard KILL CRITERIA

### GO — all conditions are required

1. The 18-nt footprint is sequence-verified and mapped unambiguously, with strand
   and reference-frame crosswalks recorded.
2. At least **two unrelated donors** in material that is public and downloadable
   now each have a non-reference allele within that footprint on a callable,
   read-backed, copy-resolved SMN2 haplotype.
3. For each qualifying donor, the audit records public locator, sample ID,
   haplotype/copy ID, REF, ALT, coordinate frame, depth/read support or equivalent
   published QC, and the paralog-discriminating evidence used to call it SMN2.
4. Manual read review (or an equivalently auditable published haplotype record)
   does not show the hit to be a mapping, strand, homopolymer, or paralog-assignment
   artifact. Same-donor duplicate copies and the same allele copied between files
   are deduplicated.
5. The Appendix is presented as a computational `candidate` decomposition unless
   actual paired variant +/- ASO measurements are found. No HPRC donor is called a
   patient and no modeled interaction is called a biological interaction.

If all five conditions hold by the end of the two-day extraction, the Appendix is
**GO** for the narrow model-decomposition arm.

### NO-GO — any one condition kills the arm

- After complete extraction of the qualifying public sources above, fewer than two
  unrelated donors have a qualifying copy-resolved SMN2 ALT inside the footprint.
- The apparent hits cannot be assigned to SMN2, cannot be separated across SMN2
  copies, are not callable over the full footprint, or disappear on read/QC review.
- The threshold can be met only by counting two copies from one donor, HPRC
  assemblies as patients, generic short-read calls, controlled/author-mediated
  data, or a source that cannot be downloaded within the gate.
- The target sequence, strand, or coordinate crosswalk remains ambiguous enough
  that overlap cannot be adjudicated reproducibly within two days.

On NO-GO, publish only a target-site conservation audit with its callable
denominator and limitations, then stop. Do not replace natural alleles with in
silico mutagenesis and retain the original claim.

### NEEDS-EXTRACTION — the present provisional verdict

This is the only defensible current verdict. The public Paraphase BAM/BED deposit
and the public HapSMA per-haplotype workbook appear capable of answering the gate,
but their 18-nt target rows have not been extracted and reconciled. Literature
statements suggesting ISS-N1 conservation do not substitute for the required
copy-resolved audit, and no qualifying-ALT count has been measured here.

The **single biggest uncertainty** is whether the public, per-haplotype HapSMA table
and/or Paraphase haplotagged bamlets contain even two donor-independent,
read-supported SMN2 ALTs inside the correctly crosswalked 18-nt footprint. That is
the first and decisive extraction; everything else is secondary.

## Freeze-safe execution boundary

This memo is recon only. It was prepared without reading `feng/terminal/` and
without reading, importing, moving, or path-referencing code from the excluded
legacy/source packages. A later extraction must remain self-contained under this
artifact directory or use newly re-derived logic; it must not inspect terminal
results. No implementation or result-producing code is introduced here.

## Extraction result — 2026-09-02

**Resolved verdict: `NO-GO — insufficient public evidence retrievable and
auditable within the gate`.** This supersedes the provisional
`NEEDS-EXTRACTION` status above. It is an evidence-availability verdict, not a
claim that natural target-site variation is absent.

### Footprint and coordinate reconciliation

The locked drug sequence is ASO 10-27/nusinersen
`5'-TCACTTTCATAATGCTGG-3'`. Its plus-strand SMN pre-mRNA target is the reverse
complement `5'-CCAGCATTATGAAAGTGA-3'`, intron 7 positions +10 through +27.

The GENCODE v50 GRCh38 primary-assembly FASTA already staged under the feng
scratch area was queried directly. The canonical SMN2 exon-7 interval ends at
GRCh38 `chr5:70076574`; therefore +10..+27 is the 1-based inclusive interval
`chr5:70076584-70076601`. The directly extracted flank was:

```text
GRCh38 SMN2 chr5:70076575-70076610
GTAAGTCTGCCAGCATTATGAAAGTGAATCTTACTT
         |------------------|
         chr5:70076584-70076601
```

Paraphase realigns both paralogs to the GRCh38 SMN1 frame. The paralog anchor
PSV13 (`c.840C>T`) is at SMN2 `chr5:70076526` and SMN1
`chr5:70951946`, a local offset of **+875,420 bp**. Applying that offset gives
the Paraphase query interval `chr5:70952004-70952021`. Direct extraction of
`chr5:70951995-70952030` from the same FASTA produced the identical 36-nt flank
and the identical target sequence. This is a sequence-checked local paralog
transform, not an assumption that the full loci have a single constant offset.

HapSMA reports variants after directing SMN1- and SMN2-derived reads to the
T2T-CHM13 SMN1 phasing frame. Its public PSV liftover table records PSV13 as
0-based BED `chr5:71408733-71408734` (1-based `chr5:71408734`). The target begins
58 bases and ends 75 bases after the PSV13 base in the GRCh38 SMN1 frame.
Applying those local offsets yields the 1-based HapSMA query interval
`chr5:71408792-71408809` (0-based half-open
`chr5:71408791-71408809`). This interval lies upstream of the same table's PSV14
at 1-based `chr5:71408882`. The transform is anchored and reproducible, but a
T2T reference slice or the workbook itself could not be retrieved in this
execution, so no donor call was accepted on coordinate inference alone.

### Retrieval and QC outcome

No whole genome was downloaded. No controlled or author-mediated source was
used. The exact public HapSMA workbook URL resolved to
`13073_2025_1448_MOESM5_ESM.xlsx`, but binary retrieval failed before a file
could be staged. The public Paraphase repository exposed the per-sample
`.bam`, `.bam.bai`, `.bed.gz`, and `.tbi` paths, but those binary artifacts also
could not be retrieved. Direct shell networking was disabled, the mandated
`$SCRATCH/feng/phase2/` path was read-only to the execution sandbox, the Slurm
controller was unreachable for a tightly scoped staging job, and the available
web fetcher returned cache-miss/no-content results for both binary sources.

The public article and repositories were sufficient to verify the methods and
coordinate anchors, but not to inspect any donor's target bases or read support.
The controlled SMA-patient reads remained excluded. HPRC samples were treated
only as non-patient controls and were not relabeled as patients.

### Per-donor evidence table

The table intentionally has no donor rows: no donor-level target record was
retrieved, so populating a donor ID, ALT, depth, or QC result would invent
evidence.

| Donor ID | Source | Copy-resolution method | Target ALT and read support | QC |
|---|---|---|---|---|
| *No qualifying donor-level record available for inspection* | HapSMA workbook and Paraphase BAM/BED sources attempted | Published methods are copy-resolved; no individual target record was locally auditable | Not measured | **FAIL / not evaluable** |

Source-level disposition:

| Public source | Tight target | Outcome | Gate disposition |
|---|---|---|---|
| HapSMA Additional file 5, Table S3 | T2T-CHM13 `chr5:71408792-71408809` | Workbook locator verified; binary table not retrievable/stageable in this execution; no sample/haplotype row or DP/AF inspected | Does not count |
| Paraphase `SMN_phased_data` haplotagged BAM/BED files | GRCh38 SMN1-realignment `chr5:70952004-70952021` | Public paths verified; binary records not retrievable/stageable in this execution; no HP-tagged reads inspected | Does not count |
| Zwartkruis SMA-patient raw ONT reads | Same T2T target after HapSMA remapping | Controlled/data-sharing access; not requested and not used | Excluded |

### Measured count and decision

- **Qualifying donor-independent target-site ALTs evidenced by this extraction:
  0.** This is a count of qualifying, auditable evidence records found, not a
  measured allele count and not evidence of target-site conservation.
- **Callable donor denominator: not measured.** No donor may be called REF merely
  because it is absent from an inaccessible variants-only file.
- **Threshold:** GO requires at least two unrelated donors with public, callable,
  read-backed, copy-resolved SMN2 target ALTs.
- **Decision:** `NO-GO — insufficient public evidence`. The threshold was not met,
  and the memo explicitly kills the arm when the source needed to meet the
  threshold cannot be downloaded and audited within the gate. The variant x ASO
  Appendix stops here. A target-site conservation claim is also not made because
  the callable denominator was unavailable.

Honesty state: the footprint definition, public-source availability, and the
published copy-resolution methods are `paper_supported`; the local coordinate
sequence checks are supporting audit facts. There is no `candidate` ALT and no
`clean_reanalysis_reproduced` variant result.
