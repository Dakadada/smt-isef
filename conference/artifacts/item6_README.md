# Item 6 — dinucleotide-preserving null `[candidate]`

These are sequence-derived null-model predictions, not measurements or
wet-lab outcomes.

## Outputs and methods

- `item6_thermo_null_raw.tsv`: 100 draws per ASO through IntaRNA-default and
  Turner-2004 NN/TBE models.
- `item6_thermo_null_summary.tsv`: mean, population SD, linear quantiles,
  midrank percentile, and predeclared one-sided empirical enrichment tail.
- `item6_alignment_queries.tsv`: real and null alignment identities.
- `item6_alignment_null_raw.tsv` and `item6_alignment_null_summary.tsv`: created
  by the Item-3 cluster workflow.
- `item6_alignment_scramble_summary.tsv`: clean alignment-null distribution
  (n, mean, population SD, q05, median, q95) for F18/F14/F10, separately for
  Bowtie and BLAST.

RNG seed is **20260926**. Each draw uses a randomized Euler trail through the
sequence's directed dinucleotide multigraph, preserving the exact dinucleotide
multiset and length. Draws are with replacement. Duplicate sequences retain
their sampling multiplicity but are scored/aligned once computationally.
Scrambled ASOs are scored against the original biological target window, so the
test measures loss of sequence-specific affinity while preserving nearest-
neighbor composition.

For ΔG, more-negative is the predeclared enrichment direction. IntaRNA's `0.0`
no-favorable-duplex convention is stored as `intarna_no_duplex` and excluded
from numeric ΔG summaries; it is never interpreted as measured zero energy. A
separate all-draw empirical tail treats no-duplex outcomes only as censored
weaker draws (denominator contribution, no assigned ΔG). For gene counts,
larger is the predeclared nonspecific-hit enrichment direction.

## Assumptions and limitations

- Exact dinucleotide constraints create small null supports for short ASOs:
  F10 has 8 unique sequences among these draws and 3UP8 has 6. Identity draws
  are retained rather than post-hoc rejected.
- A null with many IntaRNA no-duplex outcomes has a smaller numeric summary;
  both requested and numeric sample sizes are reported.
- 3UP8 alignment nulls are all marked saturated/noninformative rather than
  assigned fabricated gene counts. F10-family nulls retain the mandatory
  exhaustive override so the requested 10-mer baseline remains estimable.
- The null establishes pipeline behavior only. It cannot prove or disprove
  biological specificity.
- **Human confirmation required:** 3UP8's inferred original target is
  `CUGCCAGC` at supplied-context positions 7–14.
