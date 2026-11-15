# Item 1 — hybridization free energy `[candidate]`

These are sequence-derived model predictions, not measurements, wet-lab
outcomes, efficacy estimates, or clinical/safety claims.

## Outputs

- `item1_dg_raw.tsv`: one row per ASO. IntaRNA-default, raw Turner-reference,
  salt adjustment, and salt-adjusted NN ΔG are separate fields.
- `item1_salt_sensitivity.tsv`: predeclared effective-monovalent grid.
- `item1_dg_ranks.tsv`: derived ranks only; it does not replace raw energies.
- `item1_provenance.json`: versions and fixed parameters.

## Methods

All sequences are modeled as RNA:RNA at 37 C. Model A reuses
`feng.titration.evaluate.intarna_on_target` and IntaRNA 3.4.1 defaults
(Turner 2004; no explicit salt). Model B uses RNAduplex 2.7.2 and its default
Turner-2004 RNA nearest-neighbor parameters. RNAduplex raw MFE is referenced to
its reported 1.021 M default. ViennaRNA 2.7.2 accepts `--salt`, but direct
checks showed that RNAduplex MFE output does not change across the salt grid;
therefore the output applies and separately records the standard monovalent
entropy adjustment
`ΔS = 0.368*(Npaired-1)*ln([M+]/1.021 M)` and `ΔGsalt = -TΔS/1000`.

The primary operational buffer (user-specified 2026-09-30) is **50 mM monovalent
+ 1.5 mM Mg2+**. The divalent is folded into a single effective monovalent via
the von Ahsen (2001) equivalence `[Na+]_eq(mM) = [mono](mM) + 120*sqrt([Mg2+](mM))`,
giving **0.197 M effective monovalent**, which then feeds the SantaLucia
correction above. This is a DNA-derived empirical rule applied to an RNA:RNA
duplex, so the Mg2+ term is a transparent heuristic, not a validated RNA
measurement. Per-ASO buffer ΔG is recorded in `nn_buffer_dg_kcal_mol` with its
`nn_buffer_{monovalent,divalent,effective_monovalent}_molar` basis columns.
The predeclared monovalent sensitivity grid is retained: 0.025, 0.050, 0.089,
0.100, and 0.150 M. (The earlier 0.089 M 1x-TBE Tris-proxy point is superseded
as the primary but remains representable on the grid.)

References: Mathews et al., *PNAS* 2004, 101:7287–7292,
doi:10.1073/pnas.0401799101 (Turner-2004 RNA parameters); SantaLucia,
*PNAS* 1998, 95:1460–1465, doi:10.1073/pnas.95.4.1460 (monovalent entropy
form); von Ahsen et al., *Clin Chem* 2001, 47:1956–1961 (Mg2+→Na+ equivalence);
Green & Sambrook, *Cold Spring Harb Protoc* 2019,
doi:10.1101/pdb.rec101246 (TBE recipe).

## Assumptions and limitations

- ASOs are interpreted 5'→3'; intended RNA targets are exact reverse
  complements found uniquely in the supplied context.
- F18/F14/F10 targets are context positions 10–27/10–23/10–19 (1-based).
- **Human confirmation required:** 3UP8 is inferred at positions 7–14,
  `CUGCCAGC`, three nucleotides upstream of F18's start.
- `Npaired` is the number of formed pairs in RNAduplex's selected structure;
  this makes the correction explicit for partial scramble duplexes.
- The SantaLucia monovalent form and the von Ahsen Mg2+→Na+ equivalence are
  DNA-derived approximations layered onto an RNA NN model. The grid exposes the
  monovalent sensitivity; neither the buffer ΔG nor the Mg2+ term is a calibrated
  RNA thermodynamic measurement.
- The primary chemistry is RNA:RNA. Code accepts a parameter file for reruns,
  but a DNA:RNA rerun requires a validated hybrid parameter set; Vienna's `DNA`
  preset must not be mislabeled as DNA:RNA.

