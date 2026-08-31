# Generative splice-switching design engine

This package is the prospective design surface for ranked steric effectors. Its
stable API is `design_for(target) -> DesignReport`; `target` may be a built-in
event ID or a `TargetDefinition` containing a strand-oriented transcript context,
physics event definition, bounded design interval, and frozen off-target
transcriptome. Reports contain ranked ASO, dCas13 gRNA, and U7 candidates, the
multi-objective Pareto front, uncertainty, and immutable-run provenance.

Every report begins with **PROSPECTIVE PREDICTION — NOT EXPERIMENTALLY
VALIDATED**. Without adapter-accepted typed occupancy evidence, efficacy is only a
robust direction and confidence category. No ΔPSI or effect magnitude is invented.
Adapter abstentions remain abstentions and cannot enter the Pareto front. The
off-target burden is recall-oriented and unsigned: it is not a prediction that a
hit changes splicing.

## Offline use

```python
from design_engine import design_for, render_json, render_text

report = design_for("SMN2_exon7_inclusion")
print(render_text(report))
payload = render_json(report)
```

The built-in SMN2 and DMD contexts and deterministic folder are CI/demo fixtures,
not biological models. They exercise all wiring without network access or weights.
SMN2 is explicitly labeled for exon inclusion, while DMD is explicitly labeled for
therapeutic exon skipping even though the measured event is named exon53 inclusion.

## Networked real-prediction path

In a networked, reproducible session, construct a `TargetDefinition` using the
version-exact full pre-mRNA sequence and the event coordinates derived from the
frozen GRCh38 annotation. Supply a production `Folder` implementation to
`DesignEngine` that performs the real constrained folding requests and returns the
adapter's required constraint receipts. Freeze and supply the complete expressed
pre-mRNA transcriptome, expression/localization metadata, and signed eCLIP or
perturbation evidence to the target for off-target enumeration.

Run SpliceAI and Pangolin with verified official weights and golden-output checks as
separate predictive evidence; do not silently substitute them for the physics
adapter or reinterpret their scores as validation. Record model/weight/reference
hashes in the deposited run manifest. Real folding plus verified SpliceAI/Pangolin
plus a real transcriptome enable real *prospective predictions*—they still do not
make candidates experimentally validated. Independent cellular and safety assays
remain required.

Configuration is deterministic and bounded through `GenerationConfig`. A frozen
clock can be injected for byte-identical test/deposition regeneration; ordinary
runs record a UTC timestamp. JSON is emitted with sorted keys for immutable storage.
