# Pre-registered steric-intervention benchmark protocol

**Lock date:** 2026-08-28. Primary scope is the SMN2 case study and public-data feasibility/negative-result benchmark. Multi-locus evaluation activates only after byte verification of >=3 independent study families across >=2 loci.

## Estimand and hierarchy

The primary feasibility endpoint is **direction only**: the sign of cross-fitted predicted treated-minus-control delta-PSI for the preregistered event. Evaluate it with directional balanced accuracy and directional concordance. The exact mapping, frozen primary equilibrium model, constrained link, theta proxy/grid downgrade, and eCLIP/probing ablation status are locked in `docs/architecture.md`. Delta-PSI magnitude is a conditional secondary endpoint, not a current primary estimand.

The hierarchy is **study family -> biological sample -> oligo/intervention -> event**. Shared controls, sources, laboratories, constructs, or near-duplicate oligos are grouped into one family. Samples are replicates, oligos/guides/cassettes are intervention units, and events are repeated outcomes that **never count as intervention replication**. One targeted event per intervention is primary.

Outer splits hold out whole families; inner selection uses training families only. With fewer than 3 eligible families there is no generalization/model-selection claim or terminal reservation.

## Eligibility and theta set

A strict family requires all endpoint runs downloaded and checksum-matched, exact intervention sequence/chemistry, resolved treated/control map, and a primary event with >=20 informative junction reads per condition. Replicate uncertainty is retained. Archive statements are counted separately.

No current candidate intervention has the matched external occupancy proxy in `docs/architecture.md`, so the theta-requiring magnitude estimand currently applies to zero interventions. Magnitude analysis activates conditionally and only for interventions with documented matched theta evidence. All others receive the frozen theta-grid curve and must have a direction unchanged over that grid; they enter the primary sign metrics but never magnitude, MAE, Spearman, calibration, occupancy, or mechanism analyses. Missing predictions are reported; outcomes cannot impute them.

## Baselines

SpliceAI/Pangolin are frozen unchanged-sequence feature extractors. Their untreated propensity and frozen in-silico occlusion importance, plus affinity, accessibility, splice-element distance/orientation, chemistry/dose/length, and k-mers form the measurable baselines in `docs/architecture.md`. They are not direct trans-effect competitors. Baselines get identical splits, metadata, coefficient/knot budget (within one degree of freedom), tuning budget, and coverage rule. Physics must beat the best baseline selected wholly inside each outer fold.

## Power simulation before thresholds

Before outcomes open, publish a seeded hierarchical simulation using the verified hierarchy and missingness: study-effect SD `{3,7,12}` PSI points; interventions/family `{1,2,4,8}`; biological replicates `{2,3,4}`; measurement SD `{3,7,12}`; true Spearman `{0,.2,.4,.6}`; MAE gains `{0,10%,20%,30%}`; manifest-matched locus imbalance; theta availability `{25%,50%,100%}`. Run 10,000 datasets/cell through the exact cross-fitting/resampling code.

For the conditional magnitude endpoints, fix a numeric success threshold only when the verified theta-supported design has **>=80% power** at rho=.4/20% MAE gain and **<=5% false-positive probability** under rho=0/no gain; otherwise the metric is descriptive. The primary direction endpoint is a feasibility report unless a separately frozen directional power simulation supports a numeric threshold. Freeze the report, minimum intervention count, and any stricter threshold before outcome access. It may tighten, never relax, the floors below. Do not compute a study-cluster CI with fewer than **8 families**; the confirmatory target is **10**. Below 8, show study estimates/ranges without cluster CI or inferential p-value.

## Staged feasibility tiers

`N` is the byte-verified independent family count.

| Tier | Diversity floor | Exactly claimable | Not claimable |
|---|---|---|---|
| **N=2** | >=2 interventions | paired pipeline reproducibility, per-case theta-grid/robust direction, observed effects, failures | cross-study performance, correlation/calibration, physics superiority, terminal test, locus/effector generalization |
| **N=3** | >=2 loci, >=3 interventions | mechanical three-fold LOSO **feasibility**, fold errors/coverage/direction, descriptive matched-baseline differences | validated benchmark/predictor, cluster CI, general multi-locus/mechanism claim |
| **N=4** | >=2 loci, >=5 interventions | descriptive four-fold LOSO stress test, feasibility/negative-result paper, per-study errors | confirmatory generalization, cluster CI, terminal claim, rank/calibration unless separate floors pass |

Above N=4, N=4 restrictions persist until metric-specific and simulation floors pass. A never-touched terminal study requires **>=10 families total**, leaving at least 9 for development; otherwise evaluation is cross-fitted and non-terminal.

## Metric floors and gates

- **Spearman:** >=12 independent interventions, >=2 loci, >=8 families, and >=4 interventions/locus. Confirmatory use requires 10 families and simulated power. If powered: rho >=.40 and study-cluster 95% lower bound >0.
- **Calibration:** >=20 independent interventions, >=2 loci, >=10 families, >=5 interventions/locus, with matched theta for every case. If powered: slope 0.75-1.25 with 95% CI containing 1; nominal 90% coverage must be 85-95%.
- **MAE value:** theta-supported >=12 interventions, >=2 loci, >=8 families; confirmatory use needs 10 families and power. If powered: >=20% reduction versus best matched-capacity baseline and cluster 95% lower bound >10%.
- **Direction (primary feasibility endpoint):** directional balanced accuracy and concordance use the frozen no-change `|delta-PSI|<5` (2/10 sensitivities). They are descriptive below 8 families; no study-cluster CI or inferential p-value is computed below that floor.
- **Magnitude (conditional secondary endpoint):** MAE and calibration slope activate only for interventions with documented matched theta evidence and remain subject to their intervention/locus/family and simulation floors above. With the current candidates, this endpoint has zero eligible interventions.

When confirmatory floors pass, all powered applicable components must pass; subgroup/composite results cannot rescue failure. Below floors, report feasibility/negative results rather than a hard `>=3 studies` failure.

## Weighting, uncertainty, leakage

Give each family equal total weight and each intervention within family equal weight. Propagate beta-binomial junction uncertainty. Resampling follows families, then biological samples, then interventions; events can describe secondary outcomes but never supply replication. Use 10,000 draws/frozen seed only above the 8-family CI floor.

Fingerprint sequences/windows; group shared controls/sources; exclude intervention samples from cis training; hide held-out FASTQ outcomes, tables, and figures until a signed prediction manifest exists. Never infer theta from endpoint PSI. Hash reference, event, strand, contrast, filters, code, container, model, link, grid, and baselines before evaluation.

## Reporting

Report strict/archive counts, loci, unique interventions, theta availability, every fold/case, exclusions, coverage, detectable-effect limits, baselines, and ablations. At the audited floor (0 byte-verified; about 4 candidates/2 loci/5 interventions), the deliverable is an SMN2 case study plus feasibility/negative-result benchmark and immutable designs. Do not claim a validated steric predictor or biological/therapeutic mechanism. Novelty/IP language cites `docs/prior-art-defense.md`: the benchmark, not biological discovery, is durable.
