# Arm D control-only detectable-exon universe (GENCODE v50)

> Claim boundary: descriptive background construction only. Allowed
> result states are `candidate` and `paper_supported`. No treatment
> effect is read or reported here.

## What this is

The background set of internal cassette exons (rMATS SE events) that are
adequately measured in the GSE307431 baseline libraries, built from the
scrambled and untreated control alignments only. It supplies the matched
background pool for the locked nusinersen off-target ranking; it does not
contain or depend on any ASO-treated measurement.

## Control samples used [paper_supported]

| run_accession | condition | chemistry | control_type |
| --- | --- | --- | --- |
| SRR35298717 | control | OMe | scrambled |
| SRR35298701 | control | PMO | scrambled |
| SRR35298704 | control | MOE | scrambled |
| SRR35298706 | control | PMO | scrambled |
| SRR35298709 | control | OMe | scrambled |
| SRR35298711 | control | OMe | scrambled |
| SRR35298718 | control | untreated | untreated |
| SRR35298699 | control | MOE | scrambled |
| SRR35298708 | control | untreated | untreated |
| SRR35298713 | control | PMO | scrambled |
| SRR35298715 | control | MOE | scrambled |

Total control libraries pooled: 11 (9 scrambled + 2 untreated).

## Method [candidate]

- rMATS-turbo SE events enumerated from GENCODE v50 with `--task both`
  and `--statoff`; both rMATS groups are controls. `--statoff` disables
  statistical testing, but rMATS still emits a dummy-control
  `IncLevelDifference` in some source rows. The builder ignores and
  explicitly rejects all contrast fields before parquet publication.
  Junction counts and inclusion levels from both
  groups are pooled into one baseline set per exon. This baseline pool
  never combines the three ASO chemistries.
- Alternate-flank rMATS events are collapsed by the genomic exon key
  `chrom:exonStart_0base-exonEnd:strand`; each parquet row is one unique
  exon interval. Source event IDs and multiplicity remain in audit columns.
  If one interval is annotated to multiple genes, all gene IDs/symbols and
  per-gene control CPM values are retained in aggregate provenance columns;
  `multi_gene` makes that ambiguity explicit.
- Inclusion level (PSI) is the rMATS junction-count inclusion level; an
  exon's junction depth is the mean of inclusion-junction plus
  skipping-junction uniquely counted reads across control samples.
- Internal exons only: an SE event always has a defined flanking exon on
  each side, so the 5'-most and 3'-most exons of a transcript are
  excluded by construction.
- Gene expression is mean control CPM from STAR reverse-strand gene
  counts (fr-firststrand). GC fraction and exon length come from the
  reference exon coordinates.

## Detectability thresholds [candidate]

- mean total junction reads >= 20
- inclusion level estimable in >= 6 control samples

## Counts [candidate]

- Unique internal exon intervals after alternate-flank collapse: 170584
- Detectable background exons: 81621
- Distinct genes among detectable exons: 9407
- Detectable baseline PSI quartiles (Q1/median/Q3): 0.106 / 0.990 / 0.999

### Detectable exons by chromosome (density) [candidate]

| chrom | detectable_exons |
| --- | --- |
| chr1 | 8415 |
| chr2 | 6006 |
| chr3 | 5056 |
| chr17 | 5002 |
| chr19 | 4822 |
| chr11 | 4815 |
| chr12 | 4627 |
| chr7 | 3962 |
| chr16 | 3862 |
| chr5 | 3765 |
| chr6 | 3696 |
| chr9 | 3354 |
| chr10 | 3269 |
| chr14 | 2889 |
| chr8 | 2790 |
| chr15 | 2741 |
| chr4 | 2679 |
| chr20 | 2268 |
| chr22 | 2117 |
| chrX | 2030 |
| chr13 | 1305 |
| chr18 | 1240 |
| chr21 | 857 |
| chrY | 54 |

## Provenance [candidate]

- SE source table: `/home/dashii/orcd/scratch/feng/bench/detectable_universe/SE.MATS.JC.control-only.txt`
- SE source SHA-256: `0226cacb2efd9521651ca540771edf85bb6bd25e2923bb147bbd276bfc7de2bd`
- Reference FASTA: `/home/dashii/orcd/scratch/feng/gencode/v50/GRCh38.primary_assembly.genome.fa`
- Conditions table: `feng/artifacts/phase3/sample_conditions.tsv`
- Conditions table SHA-256: `e14f80ba3341ff52d088e19cd016d4aeff885237804d3ba986d0d76533b9e564`
- rMATS b1 control list: `feng/artifacts/bench/provenance/b1.txt`
- rMATS b1 control list SHA-256: `0c278cbef3f77e610b58940864eeca351ebb8ca91c6a41f5f36c3b070521ea08`
- rMATS b2 control list: `feng/artifacts/bench/provenance/b2.txt`
- rMATS b2 control list SHA-256: `14efb75911938b83495fd48ff7bf9e472ea7075e421e12d8e62b18cf683b6fcd`

## Downstream gate note

The locked benchmark requires at least 50 eligible background exons per
positive off-target exon. That gate is evaluated once the positives
table lands; this file reports the detectable count and per-chromosome
density so the gate can be checked against the matched strata.
