# Leak-proof nonstructural baseline arena

The arena predicts a schema-proven treated-minus-one-matched-control endpoint. Every
label must carry nonempty study, locus, intervention, chemistry, and biological-sample
cluster IDs. Outer folds are connected components of the graph induced by sharing
*any* identity. The harness asserts disjointness on all five dimensions and refuses
missing identities. Family is reserved for weighting and uncertainty, not leakage control.

Magnitude labels require `harmonized-treated-minus-matched-control/v1`, explicit
contrast type/control cluster/count/sign convention/units, and typed theta evidence.
The protocol minimum for decisive magnitude or CI output is eight families and cannot
be lowered by configuration.

## Strong frozen baselines

Affinity uses Sugimoto RNA/DNA nearest-neighbor ΔG37 terms, a logarithmic monovalent-salt
correction, and prospectively registered DNA/2'-OMe/MOE/LNA/PMO corrections. Missing,
unaligned, noncanonical, or non-complementary oligos fail. Motif features scan every
overlapping position using the complete four-matrix
[ESEfinder 2.0](https://esefinder.ahc.umn.edu/tools/ESE2/ESEmatrix.html)
functional-SELEX model and the experimentally supported SFmap hnRNP/PTBP silencer
[library](https://pmc.ncbi.nlm.nih.gov/articles/PMC2896136/). ESEfinder reports
678 SF2/ASF, 669 SRp40, 2,599 SC35, and 133 SRp55
above-threshold sequences, so this is not a short handpicked motif list.
Splice-distance features are signed and donor/acceptor plus inclusion/skipping typed.
The k-mer search is pre-registered at k={3,4,5}; a full regularized representation is
fit, with no variance shortlist.

`default_config.json` freezes both complete searches. Baselines search k and ridge
lambda only inside identity-component-disjoint training folds. Physics searches the
frozen theta grid and uses the frozen monotone-link knots and group ridge precision;
it receives the identical EDF budget. Search amendments require a schema/version
change, not a runtime argument.

Ridge capacity is the real weighted hat-matrix EDF,
`tr((X'WX + lambda P)^-1 X'WX)`, with an unpenalized intercept. Ridge and k searches
are training-only. The physics link and every baseline declare their measured EDF and
must fit the same frozen EDF budget; fold reports contain both values.

## Event-aware occlusion

Each splice site declares donor/acceptor and whether more propensity increases or
decreases the endpoint. Substitution changes remain site-specific, are signed into the
common endpoint, and are only then summed. The same `occlusion_importance` transform
must be used for every compared sequence model. Physics native treatment predictions
must be reported separately and never relabeled as derived occlusion sensitivity.

## Strict SpliceAI/Pangolin adapter contract

`VerifiedCisBackend` requires model name, package version, exact weight SHA-256,
required context, padding, coordinate convention, strand handling, donor/acceptor
channel set, reference/GENCODE convention, and nonempty golden outputs copied from an
official implementation. Construction fails unless loaded weights/reference/golden
observations match the frozen manifest. Scoring then requires both bounded channels at
the exact unpadded input length. A scalar callable, a merely length-matched vector, a
stub, or an unverified weight file cannot be called SpliceAI/Pangolin.

In a networked validation session, install the pinned official packages, record the
licensed weight hashes, generate fixed forward/minus-strand and boundary-padding cases
with the official CLI/API, store those values in `golden_outputs`, and run the contract
test. Offline tests deliberately verify the fail-closed placeholder until real hashes
and official golden values are frozen. RNAplfold remains a separately version-pinned
accessibility backend.

Uncertainty samples families with replacement and carries all their interventions and
replicate-cluster identities. Reports include unique family support, Kish effective
clusters, and leave-one-family-out sensitivity, avoiding fake precision near eight
families.
