# J3 model card

Inputs: `/home/dashii/orcd/scratch/feng/model_v2/j3_genome_v1/gate_report.json` (SHA256 `f6f1f45e8d2c70fc83f9032ad7a9d17af21e29b5c9fcc61847b4998fb16921f9`); `/home/dashii/smt-isef/conference/model_v2/j3_scoring.py` (SHA256 `93614579945832ddde19443b7b15e8dbe0390cb6ce2002f7529c1fb2a1d5eb6f`).

Evidence class: derived / proxy / model_assumption.

Output: Sequence complementarity, RNA duplex energy and accessibility with explicit channel failures.

Labels: unlabeled; no outcome labels or verified negatives are read, joined, fitted or selected.

Assumptions: GRCh38 / GENCODE v50; sequence-derived SMN coordinates remain unconfirmed against GTF. Chemistry unknown; any duplex term assumes unmodified RNA. Canonical context and model tissue channels do not establish cell-specific ASO response.

Limitations: No potency or measured occupancy; failed rows remain failed.

Future independent wet-lab comparison must be made against this snapshot. Already-seen July gels are NOT prospective validation.
