# Arm D control-only detectable-exon universe (GENCODE v50)

> Claim boundary: descriptive background construction only. Allowed
> result states are `candidate` and `paper_supported`. No treatment
> effect is read or reported here.

## What this is

The background set of annotated transcript-internal exons that are
adequately measured in the GSE307431 baseline libraries, built from the
three scrambled MOE (ScrMOE) control alignments only. It supplies the matched
background pool for the locked nusinersen off-target ranking; it does not
contain or depend on any ASO-treated measurement.

## Control samples used [paper_supported]

| run_accession | condition | chemistry | control_type |
| --- | --- | --- | --- |
| SRR35298704 | control | MOE | scrambled |
| SRR35298699 | control | MOE | scrambled |
| SRR35298715 | control | MOE | scrambled |

Total control libraries pooled: 3 (all ScrMOE).

## Method [candidate]

- Candidate identities and transcript flanks are seeded from every
  transcript-internal exon in the GENCODE v50 annotation. Detectability
  is therefore not gated by whether rMATS enumerates a skipped-exon event.
- rMATS-turbo SE events are retained as source-event provenance. They were
  enumerated from GENCODE v50 with `--task both`
  and `--statoff`; both rMATS groups are controls. `--statoff` disables
  statistical testing, but rMATS still emits a dummy-control
  `IncLevelDifference` in some source rows. The builder ignores and
  explicitly rejects all contrast fields before parquet publication.
  Neither rMATS group supplies the published exon metrics; those metrics
  are recomputed from the three ScrMOE STAR junction tables and never
  combine the three ASO chemistries.
- Alternate annotated transcript flanks are collapsed by the genomic exon key
  `chrom:exonStart_0base-exonEnd:strand`; each parquet row is one unique
  exon interval. Source annotation IDs and multiplicity remain in audit columns.
  If one interval is annotated to multiple genes, all gene IDs/symbols and
  per-gene control CPM values are retained in aggregate provenance columns;
  `multi_gene` makes that ambiguity explicit.
- Inclusion level (PSI) uses the rMATS JC estimator applied directly to
  annotated exon flanks and the three ScrMOE STAR junction tables; an
  exon's junction depth is the mean of inclusion-junction plus
  skipping-junction uniquely counted reads across control samples.
- Internal exons only: candidates come from allowed GENCODE transcript
  classes with at least three annotated exons. Sorting each transcript
  by genomic position and retaining only positions between its first and
  last exon excludes the 5'-most and 3'-most exons by construction.
- Gene expression is mean control CPM from STAR reverse-strand gene
  counts (fr-firststrand). GC fraction and exon length come from the
  reference exon coordinates.
- Frozen detectability policy `armD-detectability-scrmoe-v1`: {"migration_from":"pooled-11-control minimum 6 of 11","min_mean_total_junction_reads":20.0,"min_samples_with_psi":2,"policy_id":"armD-detectability-scrmoe-v1","rationale":"With exactly three prespecified ScrMOE replicates, require PSI in at least two independent replicates (a strict majority); retain the prespecified mean-junction-depth threshold of 20."}

## Detectability thresholds [candidate]

- mean total junction reads >= 20
- inclusion level estimable in >= 2 control samples

## Counts [candidate]

- Unique internal exon intervals after alternate-flank collapse: 287731
- Detectable background exons: 147002
- Distinct genes among detectable exons: 9632
- Detectable baseline PSI quartiles (Q1/median/Q3): 0.996 / 1.000 / 1.000

### Detectable exons by chromosome (density) [candidate]

| chrom | detectable_exons |
| --- | --- |
| chr1 | 14973 |
| chr2 | 10900 |
| chr17 | 9129 |
| chr3 | 8818 |
| chr11 | 8759 |
| chr19 | 8674 |
| chr12 | 8582 |
| chr6 | 7027 |
| chr16 | 6869 |
| chr5 | 6786 |
| chr7 | 6532 |
| chr9 | 6110 |
| chr10 | 5714 |
| chr15 | 5113 |
| chr4 | 5100 |
| chr14 | 4988 |
| chr8 | 4939 |
| chrX | 4060 |
| chr20 | 3859 |
| chr22 | 3728 |
| chr13 | 2519 |
| chr18 | 2126 |
| chr21 | 1558 |
| chrY | 130 |
| KI270711.1 | 8 |
| GL000194.1 | 1 |

## Provenance [candidate]

- SE source table: `/orcd/scratch/orcd/011/dashii/feng/bench/detectable_universe/SE.MATS.JC.control-only.txt`
- SE source SHA-256: `0537c969b5d8f966e0f908649d98bfad30cc96a884f57dd2035e157b9b09820c`
- Reference FASTA: `/orcd/scratch/orcd/011/dashii/feng/gencode/v50/GRCh38.primary_assembly.genome.fa`
- Reference FASTA SHA-256: `e49b92b3e4f321bf254c042f25b726d9931c4d74c7523e8b6bb530e63b0cfd4b`
- Junction read length: `100`
- GENCODE GTF: `/orcd/scratch/orcd/011/dashii/feng/gencode/v50/gencode.v50.primary_assembly.annotation.gtf`
- GENCODE GTF SHA-256: `04af7379ad6c214f49758db270938a38c41bebbf7b9dc8e24065d8f3a9f08d0d`
- ScrMOE STAR junction SRR35298704: `/orcd/scratch/orcd/011/dashii/feng/phase3/star/SRR35298704/SRR35298704.SJ.out.tab`
- ScrMOE STAR junction SRR35298704 SHA-256: `060e1a2d10441173ed332661717b4d9643739255652ffa168454bcaaa4f07332`
- ScrMOE STAR junction SRR35298699: `/orcd/scratch/orcd/011/dashii/feng/phase3/star/SRR35298699/SRR35298699.SJ.out.tab`
- ScrMOE STAR junction SRR35298699 SHA-256: `f1562a82db5b05c5d07110b5170699e7dd0abf224aad83ffec0ab9ab1a290c42`
- ScrMOE STAR junction SRR35298715: `/orcd/scratch/orcd/011/dashii/feng/phase3/star/SRR35298715/SRR35298715.SJ.out.tab`
- ScrMOE STAR junction SRR35298715 SHA-256: `c0527f908454df24e90d1cb1ab3c98042942f0b4162f112ca17c94ecff79a3ad`
- ScrMOE STAR gene count SRR35298704: `/orcd/scratch/orcd/011/dashii/feng/phase3/star/SRR35298704/SRR35298704.ReadsPerGene.out.tab`
- ScrMOE STAR gene count SRR35298704 SHA-256: `54aadba6384d9e565d5c5d93e32dec5108475938aab278d6ff6974fce146e3b6`
- ScrMOE STAR gene count SRR35298699: `/orcd/scratch/orcd/011/dashii/feng/phase3/star/SRR35298699/SRR35298699.ReadsPerGene.out.tab`
- ScrMOE STAR gene count SRR35298699 SHA-256: `f9bda0fa8e67ce698817470f6a0e13367ac5735b680bb9d778957c4f5a6987cc`
- ScrMOE STAR gene count SRR35298715: `/orcd/scratch/orcd/011/dashii/feng/phase3/star/SRR35298715/SRR35298715.ReadsPerGene.out.tab`
- ScrMOE STAR gene count SRR35298715 SHA-256: `bc3306f670c3390e9183ef83edf6e5d1b594568585467282908601a3c4cfcd22`
- Conditions table: `/home/dashii/smt-isef/feng/artifacts/phase3/sample_conditions.tsv`
- Conditions table SHA-256: `e14f80ba3341ff52d088e19cd016d4aeff885237804d3ba986d0d76533b9e564`
- rMATS b1 control list: `/home/dashii/smt-isef/feng/artifacts/bench/provenance/b1.txt`
- rMATS b1 control list SHA-256: `c94242ea007b43c55a2ef85bd273336bd8a86f3710454e0fd588bb5950225eda`
- rMATS b2 control list: `/home/dashii/smt-isef/feng/artifacts/bench/provenance/b2.txt`
- rMATS b2 control list SHA-256: `04e09a69cfeb585482497b5d1229ce0f5f32ea94b0ac3ec84bed6e700e9b3790`

## Downstream gate note

The locked benchmark requires at least 50 eligible background exons per
positive off-target exon. That gate is evaluated once the positives
table lands; this file reports the detectable count and per-chromosome
density so the gate can be checked against the matched strata.
