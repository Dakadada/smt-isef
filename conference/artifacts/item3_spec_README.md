# Item 3 physical-specificity metrics

All outputs are **[candidate]** computational estimates, not wet-lab results.
They exclude the intended GENCODE-v50 SMN2 locus
`ENSG00000205571.17` from off-target load while retaining SMN1 and every other
genomic locus. Genomic multiplicity is preserved.

## Declared parameters

- Temperature: 310.15 K (37 C).
- Gas constant: `1.987204e-3 kcal mol^-1 K^-1`.
- Concentrations: 1 nM, 10 nM, 100 nM, and 1 uM.
- Site occupancy: `theta = [A] / ([A] + exp(dG/(R*T)))`, with concentration
  and Kd in molar units.
- Absolute dG cutoffs: -8, -10, -12, -14, and -16 kcal/mol; a site passes
  when `dG <= cutoff`.
- Accessibility cutoff: mean per-base RNAplfold unpaired probability `>= 0.5`.
  The reused `feng.titration.evaluate.rnaplfold_accessibility` defaults are
  70-nt flank padding, RNAplfold `-W 80`, and `-L 40`.
- Random seed: none. Enumeration, exact alignment, and calculations are
  deterministic.

## Inputs and state handling

Variant energies come from `item3_energy_variant_scores.tsv` (IntaRNA 3.4.1).
The Slurm replay uses Bowtie1 exact (`-a -v 0`) alignment. Occupancy uses the
published 100x-viable variants. The absolute sweep adds already-scored variants
passing the weakest absolute cutoff; it does not re-enumerate or re-score.

The cited public-RNA-seq `artifacts/judge-demo/.../evidence.json` contains
event-level SMN2 PSI observations, not a transcriptome-wide expressed-gene or
expressed-transcript set. Therefore `expression_state=unavailable` and no
expression values are inferred. The gated table is accessibility-only and says
so explicitly. RNAplfold failures remain failed/missing, distinct from loci
below the cutoff and loci without transcript-compatible gene overlap.

3UP8 remains `saturated`; no numeric specificity value is emitted for it.
There is no equivalence margin, parity classification, or pass/fail verdict.

## Tables

- `item3_spec_occupancy.tsv`: raw locus count and derived `theta_on`,
  `sum_theta_off`, `theta_on/sum_theta_off`, and
  `theta_on/(theta_on+sum_theta_off)` at each concentration. Convenience
  F10/F18 and F14/F18 burden and specificity ratios are derived columns.
- `item3_spec_gated.tsv`: the same occupancy quantities after the declared
  accessibility gate, alongside ungated load, evaluated/failed counts, and
  gated distinct transcript-compatible genes. Expression is unavailable.
- `item3_spec_absolute.tsv`: raw off-target locus count and derived distinct
  transcript-compatible gene count under each shared absolute dG cutoff.

Until the Slurm replay completes, tables contain `pending_slurm` rather than
invented numeric values. The job atomically replaces them with computed rows.

