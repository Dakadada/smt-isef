# Item 3 — GRCh38 off-target distinct-gene counts `[candidate]`

These are sequence-derived aligner predictions, not measured off-targets or
clinical/safety claims.

## Inputs and methods

- Genome: staged GENCODE-v50 GRCh38 primary assembly.
- Annotation: matching `gencode.v50.primary_assembly.annotation.gtf.gz`.
- Bowtie 1: end-to-end FASTA, exhaustive `-a -v 3`; no best-stratum filter.
- BLAST+: `blastn-short`, E-value threshold fixed at 1000, `-dust no`, both
  strands. This is a different search model; it is not forced to mimic `-v 3`.
- Primary count: distinct GTF `gene` spans overlapped on the transcript-compatible
  strand. Any-strand distinct count and raw alignment count are audit columns.
- Alignments stream into a mapper that emits only `query × overlap-class ×
  gene-id` keys. GNU `sort -u` performs disk-backed deduplication with a fixed
  2 GiB RAM budget, and a streaming reducer counts the sorted keys using
  O(number-of-queries) memory. SAM/BLAST hit dumps are not retained. Indexes are built once under
  `/orcd/scratch/orcd/011/dashii/conference_smn2_aso/`.
- F18 and F14 run as separate Bowtie partitions with four threads. Each unique
  F10 sequence runs serially in its own one-thread Bowtie process, so only one
  combinatorially large 10-mer enumeration can occupy memory at a time. The
  Slurm allocation is 220 GiB as a safety ceiling, not as a substitute for the
  bounded pipeline.

## Output tables

- `item3_query_identity_map.tsv`: persisted one-to-many mapping from each
  `Q_<hash>` to ASO, real/scramble state, scramble index, and sequence.
- `item3_real_aso_gene_counts.tsv`: clean one-row-per-real-ASO table for F18,
  F14, and F10, with Bowtie and BLAST counts side by side.
- `item3_headline_f18_vs_f10.tsv`: one-row F18-versus-F10 headline comparison.
- `item6_alignment_scramble_summary.tsv`: separate per-ASO/per-aligner scramble
  n, mean, population SD, q05, median, and q95.
- The longer `item3_gene_counts.tsv`, `item6_alignment_null_raw.tsv`, and
  `item6_alignment_null_summary.tsv` remain audit tables.

## ENERGY-FILTERED extension `[candidate]`

**Predeclared prediction:** at the primary 10-fold equilibrium-weakening
threshold, F10 would collapse toward near-perfect variants and have fewer
energy-filtered transcript-compatible genes than F18.

**Result: not confirmed.** F10 had fewer viable variant sequences than F18 at
the primary threshold (22 versus 29), but those short sequences occurred in far
more genes: **F18 = 2 versus F10 = 12,532** transcript-compatible genes. For
comparison, the unfiltered Bowtie `-v 3` counts were F18 = 1,518 and F10 =
71,148. Energy filtering reduced F10's count, but it did not establish F10
specificity over F18. The exact F10 target sequence by itself occurred in 7,778
genomic alignments overlapping 2,785 transcript-compatible gene spans, so a
relative thermodynamic filter cannot make this 10-mer locus-unique.

| ASO | Maximum weakening vs on-target | ΔΔG margin (kcal/mol) | Viable variants | Energy-filtered genes | Raw `-v 3` genes |
|---|---:|---:|---:|---:|---:|
| F18 | 1× | 0.000 | 4 | 2 | 1,518 |
| F18 | 3× | 0.677 | 15 | 2 | 1,518 |
| F18 | **10× primary** | **1.419** | **29** | **2** | **1,518** |
| F18 | 30× | 2.096 | 60 | 5 | 1,518 |
| F18 | 100× | 2.838 | 149 | 11 | 1,518 |
| F14 | 1× | 0.000 | 1 | 19 | 37,713 |
| F14 | 3× | 0.677 | 26 | 222 | 37,713 |
| F14 | **10× primary** | **1.419** | **64** | **421** | **37,713** |
| F14 | 30× | 2.096 | 107 | 844 | 37,713 |
| F14 | 100× | 2.838 | 223 | 1,671 | 37,713 |
| F10 | 1× | 0.000 | 1 | 2,785 | 71,148 |
| F10 | 3× | 0.677 | 7 | 9,325 | 71,148 |
| F10 | **10× primary** | **1.419** | **22** | **12,532** | **71,148** |
| F10 | 30× | 2.096 | 47 | 27,746 | 71,148 |
| F10 | 100× | 2.838 | 174 | 40,488 | 71,148 |

The threshold was fixed before scoring as
`ΔG_variant <= ΔG_on-target + RT ln(fold)` at 310.15 K, with a sweep of 1×,
3×, 10×, 30×, and 100× and 10× as primary. The exact target and every variant
were scored identically: ASO as IntaRNA query, isolated genomic target variant
as IntaRNA target, default RNA:RNA model at 37 °C, and the best interaction's
`E` column. This is IntaRNA's total interaction free energy, matching the
existing `intarna_on_target` wrapper and item-1 score; it is not the separate
raw `E_hybrid` component. Exact on-target values reproduced item 1: F18
-27.00, F14 -18.89, and F10 -12.80 kcal/mol.

Every unique target variant within at most three substitutions was scored
(F18 23,464; F14 10,690; F10 3,676). IntaRNA no-interaction results (`0.0` plus
the explicit no-duplex flag) were always non-viable. Only the union passing the
loosest 100× cutoff (546 variants) was reverse-complemented and aligned to
GRCh38 with Bowtie `-f -a -v 0`; each hit contributed to every stricter cutoff
that variant passed. Gene IDs were unioned across variants within each
ASO/cutoff before counting. 3UP8 remained excluded as predeclared because its
8-mer search is saturated.

Viable-variant counts include the exact intended target (therefore off-target
variant sequences equal the displayed count minus one). Gene counts likewise
retain every transcript-compatible matched gene, including the intended SMN2
gene if overlapped, matching the existing raw item-3 counter; they are not
post-hoc target-subtracted.

New outputs:

- `item3_energy_variant_scores.tsv`: all 37,830 scored variants, ΔΔG,
  equilibrium-weakening ratio, mismatch count, orientation, and sweep flags.
- `item3_energy_viable_100x.fa`: the candidate-tagged 546-query exact-alignment
  input.
- `item3_energy_exact_v0_counts.tsv`: raw exact alignments and unioned gene
  counts per ASO/cutoff.
- `item3_energy_filtered_sweep.tsv`: the full sweep with viable-variant counts,
  energy-filtered genes, and the separate raw `-v 3` comparator.
- `item3_energy_filtered_headline_f18_vs_f10.tsv`: the predeclared 10× headline
  and explicit `not_confirmed` verdict.
- `item3_energy_provenance.txt`: IntaRNA 3.4.1, Bowtie 1.3.1, reference paths,
  commands, threshold, orientation, and exclusion assumptions.

The equilibrium-ratio interpretation assumes ideal thermodynamics and does not
model concentration, kinetics, transcript abundance, full-transcript target
accessibility, cellular uptake, chemical modifications, or measured splice
effects. Gene spans include introns because ISS-N1 is intronic. All values are
sequence-derived candidate predictions, not wet-lab, efficacy, or safety
measurements.

## BLAST threshold flag — human decision required

`blastn-short -evalue 1000` is a **maximally permissive** short-query setting.
It produces thousands of raw/gene hits for some ASOs and is materially looser
than the approximately 445-gene biological framing benchmark. It should be
treated as a permissive sensitivity analysis, not the manuscript headline.
The defensible primary headline is Bowtie end-to-end `-v 3`, because its error
model is fixed and directly interpretable across the prespecified F18-versus-F10
comparison. A tighter BLAST E-value may be scientifically reasonable, but it
must be selected prospectively with human sign-off; this implementation does
not retune the threshold after seeing results.

The random-genome preflight estimate is
`2*(G-L+1)*sum[C(L,k)*3^k/4^L, k=0..3]`, with `G=3,088,269,832`; the factor two
represents both genomic orientations. The fixed saturation threshold is ten
million expected hits. Saturated nonmandatory queries have no numeric count and
are labeled `saturated / non-specific by construction - count uninformative`.
F18 and F10 families are mandatory exhaustive searches; 3UP8 is saturated and
is not aligned. F14 runs below threshold.

## Assumptions and limitations

- The primary fixed-mismatch headline is F18 versus F10. Applying three allowed
  mismatches to every length is intentionally retained and is less stringent
  proportionally for short ASOs.
- Gene-span overlap includes introns because ISS-N1 itself is intronic and the
  biological substrate is pre-mRNA.
- A Bowtie/BLAST query aligned to the reverse genomic strand is compatible with
  a `+` GTF gene; a forward query alignment is compatible with a `-` gene.
- Repetitive/multimapping loci are exhaustive unless preflight saturation is
  declared. No cap, silent truncation, or numeric substitute is used.
- Completed family/query count files carry `.done` sentinels. SIGTERM stops the
  active process group and requeues; completed partitions are skipped after
  restart, while incomplete external-sort temporary files are replaced.
- **Human confirmation required:** 3UP8's inferred target is context positions
  7–14 (`CUGCCAGC`), three bases upstream of F18.
- Software versions and exact live commands are written to
  `item3_software_versions.txt` by the completed cluster job.
