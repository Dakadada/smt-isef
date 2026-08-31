# Decisive cis-vs-trans benchmark

`DecisiveBenchmark` scores the frozen equilibrium physics adapter and every
nonstructural baseline on the same connected-component, identity-disjoint outer
folds from `baselines.evaluation`. Direction is primary. Magnitude remains absent
unless both the harmonized label and the adapter carry matched typed theta
evidence. Foldwise ridge EDF must be within one degree of freedom of physics.

## Real networked run

Do network work only in a networked/HPC session. Freeze and hash the reference,
event definitions, official SpliceAI/Pangolin weights and golden-output receipts;
run `rawread_reanalysis resolve`, `fetch`, `quantify`, `compute`, and `harmonize`
as described in the repository README. Build `FeatureBuilder` inputs with verified
cis backends, build typed adapter cases, then load the harmonized TSV with
`baselines.evaluation.load_harmonized`, run `DecisiveBenchmark`, and persist both
artifacts with `write_report`. Do not substitute archive metadata for verified
reads, and do not enable `directional_power_validated` without a separately frozen,
successful directional power simulation. At the current four-family data floor the
only possible positive verdict is `feasibility-only`.
