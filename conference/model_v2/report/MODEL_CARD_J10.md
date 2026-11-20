# J10 model card

Inputs: `/home/dashii/orcd/scratch/feng/model_v2/j10_offtarget_v2/j10_profile.json` (SHA256 `3ab6e8787cea3ddda4bd4911a03e1120ae04d7465036516e0551d985a1af76ab`); `/home/dashii/smt-isef/conference/model_v2/j10_offtarget_profile.py` (SHA256 `e6163a53788c4288736e8390bbcf382558b38188bb9b4591858943d7cc1a86e8`).

Evidence class: derived / proxy / model_assumption.

Output: Production v2 site counts by mismatch; sequence counts include SMN anchors.

Labels: unlabeled; no outcome labels or verified negatives are read, joined, fitted or selected.

Assumptions: GRCh38 / GENCODE v50; sequence-derived SMN coordinates remain unconfirmed against GTF. Chemistry unknown; any duplex term assumes unmodified RNA. Canonical context and model tissue channels do not establish cell-specific ASO response.

Limitations: v1 is stale and excluded. v2 review is bounded; no risk or safety claim.

Future independent wet-lab comparison must be made against this snapshot. Already-seen July gels are NOT prospective validation.
