# Arm A overlap verdict contest

**Date:** 2026-09-01  
**Scope:** contest the use of pairwise non-overlap as the independence unit, without
presuming that the contest passes. This memo does not change the overlap artifact or the
standing verdict.

## Decision in one paragraph

Pairwise non-overlap is a poor model of the information in a tiling experiment: two
14/15-overlapping ASOs are not independent, but their one-base contrast can still localize
an edge in the regulatory field. Independence is therefore continuous, not a binary
property of intervals. That objection is valid but immaterial to the gate. The published
figures support eventual recovery of a numeric fixed-dose outcome for **86/115 distinct
covered windows** (87 study-window records because one coordinate occurs in two studies),
and a numeric-only maximum disjoint packing still spans **11/11 blocks** (I6: 4, E7: 3,
I7: 4). However, no standardized numeric response vector is yet present on disk, the
outcomes cannot honestly be pooled without study stratification, and the effective sample
size of the 86-window numeric subset is only **10.9** under the pre-specified spatial
estimator below. An optimistic block-permutation power calculation gives a minimum
detectable Spearman correlation of approximately **0.85** at 80% power. The frozen target
effect is 0.50. **Contest verdict: RED stands; Arm A is descriptive on the available
overlap.**

## Pass bar, fixed before interpreting the calculations

This contest passes only if both conditions hold:

1. A pre-specified numeric ASO endpoint is available for at least **25 effective
   observations**, after spatial dependence is accounted for. This preserves the plan's
   frozen kill criterion rather than redefining “independent” after seeing the walk.
2. The held-out Spearman test has at least **80% power** for a true rho of **0.50** at the
   plan's one-sided permutation alpha of **0.01**.

These are necessary overlap/power conditions, not substitutes for the other frozen Arm A
success requirements (improvement over all three baselines and the nusinersen rank). Under
the usual Fisher-z approximation, condition 2 alone requires about **37 effective
observations**: `3 + [(z_0.99 + z_0.80) / atanh(0.50)]^2 = 36.26`, rounded up.

## 1. Gating dependency: the numeric ASO endpoint

Finding #3 is still a blocker. The CSV contains coordinates and prose fields
(`effect_direction`, `effect_magnitude`), not a numeric outcome column standardized within
study. A few percentages and EC50 values occur inside prose, but they do not constitute an
analysis-ready response vector. Thus the current on-disk count of analysis-ready,
standardized outcomes is **0/115**. The following is an inventory of values that can be
digitized from the cited primary figures without inventing measurements.

| Study | Numeric outcome recoverable | Distinct covered windows | What remains qualitative or unsuitable |
|---|---|---:|---|
| Hua 2007 | Percent exon-7 inclusion is printed beneath each lane for the nine-ASO initial walk and both microwalks. A consistent in-cell endpoint can be taken from Figure 2 for the initial walk and Figure 5 for the microwalks. The paper also reports EC50 values for 06-20, 36-50, 07-21, and 34-48, plus selected dose-point percentages. | **48** | EC50 is potency, not the same response as fixed-dose inclusion; it must not be mixed into the primary response. In-vitro and in-cell values sometimes disagree, so the assay must be chosen before digitization. |
| Hua 2008 | Figure 2 prints percent inclusion beneath each lane for the 20 initial intron-walk ASOs in three assay contexts. One context can be selected consistently for all 20. | **20** | Figure 7 encodes the 17 I6 and 21 I7 optimization ASOs only by direction, color, and bar thickness. Those **38** windows are not safely numeric. Mouse time-course values for selected leads are a different organism, tissue, dose, and time endpoint and are excluded. |
| Singh 2010 | Figure 1C reports percent exon-7 skipping at 20 nM for all 19 ASOs (three experiments; the figure reports the quantitative result). Inclusion can be defined as `100 - percent skipping` if fixed in advance. | **19** | Figure 1B's green/gray/pink bins alone are categorical. The Singh assay is not numerically exchangeable with either Hua assay merely because all outcomes are percentages. |

Sources: [Hua 2007 article and figures](https://journals.plos.org/plosbiology/article/figures?id=10.1371%2Fjournal.pbio.0050073),
[Hua 2008 Table 1 and Figures 2/7](https://pmc.ncbi.nlm.nih.gov/articles/PMC2427210/), and
[Singh 2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2874169/). The row-to-source mapping
and corrected S4 coordinates are in
[`harmonized_aso_positions.csv`](./harmonized_aso_positions.csv).

The raw study counts sum to 87 study-window outcomes. Hua 2008 and Singh share one numeric
coordinate, corrected S4 `I7:385..399` (published `+11..+25`), so the union is **86 distinct
windows**. A greedy interval-scheduling check restricted to those 86 windows returns 11
pairwise-disjoint intervals: 4 I6, 3 E7, and 4 I7. Thus all 11 coarse spatial units contain
at least one digitizable numeric outcome, but 29/115 exact windows do not.

### Poolability decision

A single naive pooled response is not defensible. The papers differ in study, assay,
cellular context, dose, and chemistry, and the artifact already records a direction conflict
at the shared `I7:385..399` coordinate. Worse, study is spatially confounded: Hua 2007 is
the exon-7 walk, Hua 2008 supplies I6 and broad I7, and Singh is concentrated in the narrow
ISS-N1 portion of I7. A study offset can therefore masquerade as a spatial transfer signal.

The least-bad confirmatory construction would be:

- choose one fixed-dose, in-cell percentage endpoint within each paper before extracting
  values;
- convert skipping to inclusion where necessary, then rank or standardize **within study**;
- retain duplicate coordinates as separate study-specific observations rather than average
  incompatible assays;
- compute a stratified score and shift/permutate outcomes only within contiguous
  study-by-region tracks.

That is a within-study analysis with a stratified summary, not one exchangeable 86-point
response. There are only three studies and four contiguous study-region tracks (Hua 2007
E7: 48; Hua 2008 I6: 10; Hua 2008 I7: 10; Singh I7: 19), so a random-effects synthesis is
not credible. Until the 87 figure values are actually digitized under this rule, there is
no correlation test at all.

## 2. Effective sample size: a continuous alternative to “11 or 115”

The defensible objection to pairwise non-overlap is that overlapping perturbations provide
partial information. If adjacent ASOs have correlation near one, their contribution is
almost redundant; if a one-base shift crosses a sharp field boundary, their contrast is
informative. Counting exact sequences assigns both cases a full degree of freedom, while
pairwise packing assigns the contrast none.

### Pre-specified estimator

I used the complete SMN2 dCas13d field as the only available numeric spatial track:

1. Restrict to guide starts 260--434, the corrected ASO span.
2. At each integer start, take the median `z_combined` across guide lengths. This gives a
   gapless 175-position track and prevents the 22-nt and 25-nt guides at one start from
   being counted as independent replicates.
3. Compute the empirical autocorrelation at lags 1--12 nt, stopping before the first
   negative value. The observed lag-1 autocorrelation is **0.790**; the positive-sequence
   integrated correlation time is **9.83 nt**.
4. Fit `r(d) = exp(-d / ell)` by unweighted least squares to those 12 autocorrelations.
   The fitted e-folding length is **ell = 5.61 nt**.
5. For ASO center coordinates `c_i`, set `R_ij = exp(-|c_i-c_j| / ell)`. For a
   correlation between two spatial tracks, use the product-ACF effective size
   `N_eff = N^2 / sum_ij(R_ij^2)`. This is the irregular-spacing analogue of the usual
   variance-inflation formula `N / [1 + 2 sum_k r_x(k) r_y(k)]`, under the stated
   assumption that the ASO response has the same correlation length as the dCas-derived
   predictor.

This yields:

| Input windows | Nominal N | Effective N |
|---|---:|---:|
| Complete coordinate inventory | 115 | **12.1** |
| Numerically recoverable distinct-window subset | 86 | **10.9** |

The first row vindicates the narrow contest claim: the tiling walk contains slightly more
effective dimension than the 11-interval packing. The second row is the number relevant to
the actual transfer test, and it does not pass. The drop occurs because the 29 nonnumeric
windows include spatial coverage that cannot enter a correlation.

This estimator is assumption-dependent. A sensitivity range of `ell = 3--10 nt` gives
effective N **16.9--7.7** for the numeric subset; even the favorable 3-nt assumption remains
below 25. Treating the ASO outcomes as spatially independent would increase N, but there is
no evidence in the artifact supporting that assumption, and it contradicts the reason for
using a local-autocorrelation-preserving null. Conversely, study effects or a longer ASO
correlation length would reduce N. The reported **10.9** is therefore a model-based estimate,
not a newly discovered count of independent experiments.

## 3. Power under a local-dependence-preserving null

Power was estimated with an optimistic block-permutation proxy after collapsing the spatial
design to 11 effective units (rounding 10.9 upward). For each candidate true Spearman rho,
3,000 bivariate Gaussian-copula datasets of size 11 were generated. Response ranks were
permuted among blocks 9,999 times per dataset; rejection required a one-sided permutation
`p < 0.01`. The seed was `20260901`. This preserves the premise that dependence has already
been retained inside each block and only effective blocks are exchanged. It omits outcome
digitization error, study heterogeneity, and kernel-estimation uncertainty, so it is
favorable to the contest.

| True Spearman rho | Estimated power |
|---:|---:|
| 0.30 | 0.06 |
| 0.40 | 0.11 |
| **0.50 (frozen target)** | **0.19** |
| 0.60 | 0.31 |
| 0.70 | 0.49 |
| 0.80 | 0.70 |
| 0.85 | 0.82 |
| 0.90 | 0.93 |

The simulated minimum detectable correlation at approximately 80% power is therefore
**rho about 0.85**. The Fisher-z calculation at `N_eff = 10.9` gives a slightly more
optimistic MDE of **0.81**; using **0.85** is the defensible number because it reflects the
discrete rank/permutation test. At the frozen effect of 0.50, power is only about 19%, far
below the pre-set 80% bar.

The final preregistration would still need to spell out circular-shift geometry separately
for each contiguous study-region track. Shifting across the I6/E7/I7 junctions or across
studies would not preserve the observed dependence structure. This memo does not claim an
exact power calculation for that not-yet-specified null; the block calculation is an
optimistic feasibility bound. Failure of the optimistic bound is sufficient for the gate.

## 4. Contest verdict and ranked failure modes

**Verdict: RED stands.** Pairwise non-overlap is the wrong conceptual independence unit,
but replacing it with a spatial effective-N calculation does not clear the frozen bar.
Arm A cannot honestly be represented as an adequately powered held-out transfer test on
this overlap and should remain descriptive unless genuinely new, independent numeric ASO
data are added.

1. **[blocker] No analysis-ready numeric response exists.** There are 86 digitizable
   distinct windows, but zero standardized numeric outcomes on disk. The 38 Hua 2008
   optimization windows remain categorical. Finding #3 alone prevents the correlation.
2. **[blocker] Effective N and power fail the pre-set bar.** The relevant effective N is
   10.9, not 86 or 115; estimated power at rho 0.50 is about 19%, and 80% power requires an
   implausibly large rho of about 0.85.
3. **[major] Study is confounded with position.** Within-study standardization is mandatory,
   but it cannot remove the fact that each paper samples a different spatial region and
   assay context. A pooled rank test can mistake study differences for a field effect.
4. **[major] The exact circular-shift null is not yet specified for four short,
   discontinuous tracks.** Permuting across studies or splice regions is invalid; keeping
   shifts within them reduces the effective randomization support and cannot improve the
   information bound used here.
5. **[minor, pro-contest but insufficient] Pairwise packing discards real contrasts.** The
   full 115-window design has effective dimension 12.1 rather than exactly 11. This defeats
   the claim that disjoint intervals are the uniquely correct unit, but the numerical gain
   is too small to alter the verdict.

## Reproducibility notes

Counts were re-derived from the corrected CSV by deduplicating `distinct_window_key` and
classifying as numeric-recoverable all Hua 2007 rows, all Singh 2010 rows, and only Hua
2008 rows whose locator names Figures 1 and 2A-C. The numeric-only disjoint count was
recomputed by greedy interval scheduling in increasing end-coordinate order separately for
I6, E7, and I7. The dCas autocorrelation used `z_combined` from the read-only field file
`/home/dashii/orcd/scratch/feng/phase0/splicerush_dcas13d_field.csv`; no ASO labels entered
the correlation-length estimate.
