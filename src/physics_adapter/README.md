# Effector-physics adapter

This package implements the Task-5 estimand from the locked architecture:
`P(productive | c) = (1-theta) P_U + theta P_B`. The primary result is the
equilibrium full-transcript direction. Co-transcriptional folding is secondary;
eCLIP/probing is an ablation/feature variant, never a third vote. Concordance is
reported only as a diagnostic.

Runtime is pinned to CPython 3.12 (tested with 3.12.13); the adapter itself has
no third-party dependencies. A production `Folder` may add a separately pinned
folding engine. Offline tests use only the standard library and `TableFolder`.
`default_config.json` records the frozen theta grid, empirical elongation-rate
distribution, and co-transcriptional span; load decoded configuration through
`ModelConfig.from_mapping`, which rejects unknown fields and grid changes.

Productive states are explicit joint events: all U1-contact positions and any
declared branchpoint/PPT positions are accessible, while every declared
silencer is free of repressive RBP. In the bound ensemble the target is an
explicit `IntermolecularComplexSpec`, not a constrained transcript region. The
spec fixes the coordinate convention, ASO sequence and chemistry, strand
identities, antiparallel register, and every intermolecular base pair. The target
is also unavailable to cellular factors; a bulky RNP's separately declared
shadow is unavailable as well.

`Folder` is the intentionally small integration boundary for a real partition
function implementation. It returns a `FoldingResult` containing the probability
and a constraint receipt: applied constraint type, coordinate convention,
sequence/strand/pairing/chemistry hashes, applied orientation and register,
backend version, and an explicit duplex-enforced assertion. The adapter compares
the receipt with the request and refuses the probability on any mismatch.
`TableFolder` is only an auditable deterministic stub for small offline tests.

The theta grid is frozen at `{0,.25,.5,.75,1}`. Magnitude authorization requires
an immutable `MagnitudeEvidence` whose input and protocol hashes match the
prediction and whose occupancy source and matching criterion are frozen. Without
it, callers receive theta-conditional probabilities/deltas and a conservative
interval-based sign only; the result expressly prohibits magnitude and
calibration claims. Theta is never inferred from endpoint PSI.

`CapacityConstrainedLink` freezes the monotone-link budget (global intercept and
slope, at most three knots, regularized study/chemistry intercepts), rejects
identity-like group columns and forbidden locus/oligo/interaction/label fields.
It also enforces observations per fitted degree of freedom, connected crossed
study-by-chemistry support, within-group label variation, and study-disjoint
outer folds with chemistry represented across folds. It supplies the identical
budget for each declared nonphysical baseline.

Assumptions/limits: coordinates are zero-based half-open; event annotations are
supplied by callers; no RNA folding engine is bundled; uncertainty is supplied
by the backend and mixed conservatively endpoint-by-endpoint; the empirical
elongation-rate distribution and 200-nt span are frozen defaults that must be
versioned if prospectively amended. The feature variant passes probing/eCLIP
data through the context but interpretation belongs to its backend.
