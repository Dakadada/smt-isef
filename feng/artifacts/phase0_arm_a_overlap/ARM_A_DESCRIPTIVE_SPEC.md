# Arm A descriptive cross-effector field map specification

**Status:** specification only  
**Evidence boundary:** descriptive and hypothesis-generating  
**Allowed result states:** `candidate`, `paper_supported`  
**Coordinate frame:** corrected SpliceRUSH S4 SMN2 minigene coordinates  
**Primary span:** inclusive positions 260-434

## 1. Purpose and claim boundary

Arm A will produce a descriptive map showing where the SpliceRUSH dCas13d guide field and
published SMN2 ASO effects occupy the same exon-7 neighborhood. It will align both effector
types on the corrected S4 coordinate frame and expose local same-direction,
opposite-direction, mixed, and unresolved patterns without treating the historical ASO
walks as an inferential validation set.

This deliverable answers only:

> At each part of the shared SMN2 intron-6/exon-7/intron-7 field, what dCas13d signal was
> recorded, which published ASOs overlap it, and what direction did each paper report in
> its own assay context?

It does not estimate performance, establish a cross-effector law, or combine the three ASO
studies into one exchangeable outcome.

## 2. Verified inputs and immutable scope

| Input | Use | Constraints |
|---|---|---|
| `/home/dashii/orcd/scratch/feng/phase0/splicerush_dcas13d_field.csv` | Per-guide `z_combined`, guide interval, guide length, minigene, and annotated region | Read only. Retain only `minigene == SMN2` records whose inclusive guide interval intersects 260-434. Preserve the source Z-score sign and scale. |
| `feng/artifacts/phase0_arm_a_overlap/harmonized_aso_positions.csv` | Published ASO intervals, study identity, qualitative effect fields, and source locators | Read only. Contains 125 study records and 115 exact windows. Keep study records separate, including the ten cross-paper coordinate duplicates. Do not average duplicate coordinates across papers. |
| `feng/artifacts/phase0_arm_a_overlap/README.md` | Coordinate transform, source inventory, and adversarial-review constraints | Read only. The corrected frame places exon 7 at 320-374 and intron 7 at 375 onward. |
| `feng/artifacts/phase0_arm_a_overlap/CONTEST.md` | Effective-sample and power failure, numeric-outcome inventory, and non-poolability decision | Read only. Its RED decision is fixed input to this specification. |

No new ASO outcome digitization is part of this deliverable. The 86 windows described as
digitizable in principle remain a future data-recovery possibility, not numeric input to
this map. There are zero standardized numeric ASO outcomes on disk.

## 3. Output package

When implemented, the descriptive package must be written under a new directory
`feng/artifacts/arm_a_descriptive/` and contain exactly these scientific outputs:

| Output | Format | Purpose |
|---|---|---|
| `arm_a_shared_field_map.pdf` and `arm_a_shared_field_map.svg` | Vector figure | Main aligned coordinate map with the dCas13d field and one ASO lane per study. |
| `arm_a_coarse_block_map.pdf` and `arm_a_coarse_block_map.svg` | Vector figure | Eleven-block, study-separated direction summary. |
| `arm_a_window_alignment.tsv` | Machine-readable table | One row per ASO study record with local dCas13d descriptive summaries and direction eligibility. |
| `arm_a_block_summary.tsv` | Machine-readable table | One row per block-by-study combination, including explicit numerators, denominators, and unresolved counts. |
| `arm_a_descriptive_report.md` | Narrative report | Captions, methods, result-state labels, limitations, and a fixed claim-boundary box. |

PNG previews may be generated for convenience but are not authoritative. The vector figures
and TSV tables are authoritative.

## 4. Fixed coordinate and direction conventions

### 4.1 Coordinate convention

- Use inclusive S4 coordinates everywhere.
- Display the complete 260-434 span even if a panel contains no mark at a position.
- Mark feature boundaries explicitly: I6 260-319, E7 320-374, and I7 375-434.
- Use `transformed_segments` for ASOs crossing the artificial exon-7 G. Draw the true
  segments rather than filling the artificial base. Use the inclusive hull only for guide
  containment summaries, matching the overlap artifact.
- Use guide midpoint `(start_pos + end_pos) / 2` only for the per-guide scatter layer.
  Interval coverage and local summaries must use the full inclusive guide interval.

### 4.2 ASO direction vocabulary

Map the CSV values only as follows. Do not infer a direction from free text.

| CSV `effect_direction` | Display class | Eligible for same/opposite direction counts? |
|---|---|---|
| `increase` | inclusion increase | Yes |
| `decrease` | inclusion decrease | Yes |
| `neutral` | neutral | No; count separately |
| `discordant` | assay-discordant | No; count separately |
| `neutral_or_decrease` | mixed/uncertain | No; count separately |
| `figure_reported` | unresolved, not digitized | No; count separately |
| Any other or missing value | unresolved, schema exception | No; fail the implementation check and list the value |

The map must never convert `effect_magnitude` prose into a numeric value. It may reproduce
short source-grounded annotations for named examples, with the source locator visible.

### 4.3 dCas13d sign convention

The primary displays retain raw `z_combined`. Before using the words “inclusion increase”
or “inclusion decrease” for dCas13d, the implementation must verify the Z-score sign
semantics from the SpliceRUSH source metadata and record that source in the report. If the
sign semantics cannot be verified, label the two classes only `positive Z` and `negative Z`
and omit all same/opposite direction counts. Do not infer biological direction from rank,
FDR, or the sign of selected example guides.

## 5. Figure 1: aligned shared-field map

`arm_a_shared_field_map` must use one shared horizontal axis from 260 through 434 and the
following vertically aligned panels.

### Panel A: construct landmarks

- Draw I6, E7, and I7 as labeled segments with boundaries at 319/320 and 374/375.
- Mark the artificial exon-7 G at coordinate 356.
- Print S4 coordinates at 5-nt major ticks and 1-nt minor ticks.

### Panel B: all dCas13d guides

- Plot one point per retained guide at its midpoint with `z_combined` on the y-axis.
- Encode 22-nt and 25-nt guide lengths by shape, not color.
- Use a fixed diverging color scale centered at zero for the raw Z-score.
- Draw a horizontal zero line.
- Do not filter guides by p-value, rank, FDR, or effect magnitude.
- State the retained guide count in the panel subtitle.

### Panel C: base-position dCas13d summary

For each integer base `x` from 260 to 434, define `G_x` as every retained guide with
`start_pos <= x <= end_pos`. Display:

- median `z_combined(G_x)` as the central line;
- first and third quartiles as a ribbon;
- `|G_x|` as a separate coverage strip directly below the ribbon.

No smoothing, interpolation, rolling average, or fitted curve is allowed. Quartiles must
use one named quantile convention, recorded in the report and held fixed across positions.

### Panels D-F: ASO records by study

Use separate lanes for Hua 2007, Hua 2008, and Singh 2010. For every one of the 125 records:

- draw its `transformed_segments` as an interval mark;
- color by the fixed display class in Section 4.2;
- label or tool-tip `aso_id`, original coordinate, source locator, and effect-magnitude
  prose;
- use hatching for `assay-discordant`, `mixed/uncertain`, and `unresolved` so those classes
  remain visible in grayscale;
- vertically stack overlapping intervals within a study without changing their horizontal
  coordinates.

Cross-paper duplicates must appear in both study lanes. No deduplicated “consensus ASO
effect” lane is permitted.

### Required caption

The caption must state that dCas13d values are numeric source measurements while ASO colors
are source-reported qualitative classes, that the studies differ in assay and chemistry,
and that visual proximity or same-direction color is hypothesis-generating only.

## 6. Figure 2: eleven-block descriptive summary

The eleven blocks are display bins, not independent experiments. Fix them before viewing
any block summaries:

| Region | Block ID | Inclusive S4 interval |
|---|---|---:|
| I6 | I6-1 | 260-274 |
| I6 | I6-2 | 275-289 |
| I6 | I6-3 | 290-304 |
| I6 | I6-4 | 305-319 |
| E7 | E7-1 | 320-334 |
| E7 | E7-2 | 335-350 |
| E7 | E7-3 | 351-374 |
| I7 | I7-1 | 375-389 |
| I7 | I7-2 | 390-404 |
| I7 | I7-3 | 405-419 |
| I7 | I7-4 | 420-434 |

Assign each ASO study record to exactly one block by the midpoint of its inclusive hull.
Assign a half-integer midpoint to the block containing that numeric midpoint. The bins
cover the complete span, so no in-range ASO may be unassigned. Label these as “coarse
display blocks” everywhere; never call them replicates or independent observations.

`arm_a_coarse_block_map` must contain:

1. A top heat strip showing the block median of the 1-nt positional medians from Figure 1
   Panel C. Print the numeric median in each cell.
2. Three separate ASO rows per block, one for each study. Each cell shows the counts of
   `increase`, `decrease`, `neutral`, `assay-discordant`, `mixed/uncertain`, and
   `unresolved` records assigned to that block.
3. If and only if dCas13d sign semantics are verified, a final study-separated strip with
   the descriptive same-direction fraction `n_same / (n_increase + n_decrease)` and the
   opposite-direction count. Never pool these fractions across studies. If the denominator
   is zero, print `NA (0 eligible)`.

No uncertainty interval, hypothesis test, p-value, permutation, regression, or pooled
correlation belongs in this figure.

## 7. Tables and exact descriptive statistics

### 7.1 Window alignment table

`arm_a_window_alignment.tsv` must have one row per input ASO study record and these columns:

| Column | Definition |
|---|---|
| `source`, `source_locator`, `aso_id` | Copied without modification from the ASO CSV. |
| `region`, `start`, `end`, `segments`, `distinct_window_key` | Corrected S4 coordinate fields. |
| `block_id` | Fixed midpoint assignment from Section 6. |
| `display_direction` | Exact mapping from Section 4.2. |
| `direction_eligible` | `true` only for increase/decrease. |
| `n_guides_containing_hull` | Number of retained SMN2 guides with `guide_start <= start` and `guide_end >= end`. |
| `local_z_median`, `local_z_q1`, `local_z_q3`, `local_z_min`, `local_z_max` | Descriptive distribution of `z_combined` among those containing guides. |
| `local_z_positive_n`, `local_z_negative_n`, `local_z_zero_n` | Sign counts among containing guides. |
| `direction_relation` | `same`, `opposite`, or `not_evaluable`; use `same/opposite` only after sign semantics are verified. |
| `result_state` | `paper_supported` for copied publication effects; `candidate` for the local cross-effector juxtaposition. |

The table must preserve all 125 records. A separate manifest line in the report must state
`125 records / 115 exact windows / 115 guide-covered exact windows` and must not attach an
independence interpretation to 115.

### 7.2 Block summary table

`arm_a_block_summary.tsv` must have 33 rows: 11 blocks multiplied by 3 studies, including
zero-count combinations. Required fields:

- block ID, region, inclusive start/end;
- study;
- total ASO record count;
- distinct exact-window count within that study and block;
- counts for all six display-direction classes;
- eligible direction count;
- same count, opposite count, and same-direction fraction, or explicit `NA` fields when
  dCas sign semantics are unavailable;
- block dCas positional-median summary: median, Q1, Q3, minimum, maximum, and number of
  integer positions;
- median and range of per-base guide coverage counts;
- result state, fixed to `candidate` for derived block summaries.

### 7.3 Inventory statistics in the narrative report

The report must include these counts, always with raw numerator and denominator:

- ASO records and exact windows by study and region;
- display-direction counts by study and region;
- unresolved or non-evaluable direction count by study and region;
- retained dCas13d guide count by annotated region and guide length;
- guide-coverage count per base summarized by region using median, minimum, and maximum;
- exact windows with at least one fully containing guide;
- block-by-study same/opposite counts only when sign semantics are verified.

Percentages are allowed only for composition within one named study and a printed
denominator. No percentage may combine Hua 2007, Hua 2008, and Singh 2010 into one effect
estimate.

## 8. Honest handling of carried limitations

### Qualitative ASO outcomes

- The analysis uses only the categorical `effect_direction` values already on disk.
- `figure_reported` means unresolved, not neutral and not missing at random.
- `discordant` and `neutral_or_decrease` remain distinct uncertainty classes.
- EC50 and isolated percentage values embedded in prose are annotations only.
- No rank, standardized effect, dose-response summary, or numeric ASO field may be
  reconstructed in this deliverable.

### Study-region confounding

- Show each study in a separate lane and separate block row.
- Report study-by-region cells even when they are empty.
- Do not subtract a study mean, impute unsampled regions, or compare a Hua-only region to a
  Singh-only region as though study were held constant.
- Any narrative pattern must name the study and coordinate interval that generated it.

### Assay and chemistry non-poolability

- Do not average effect directions across papers at duplicate coordinates.
- Do not majority-vote conflicting studies.
- Do not compute an all-study concordance fraction.
- Label every ASO record with its source; carry assay and chemistry metadata when available.
- Treat a cross-paper conflict as an observed context dependence, not as noise to remove.

## 9. Honesty-layer compliance

Every reported statement must carry one of two states:

- `paper_supported`: a coordinate, direction, or magnitude statement explicitly reported
  by a cited paper and transcribed without reinterpretation.
- `candidate`: a co-location, local sign relation, block pattern, or hypothesis derived by
  placing the two source types on the shared frame.

No output from this deliverable may use `clean_reanalysis_reproduced`. No `candidate` result
may be promoted because a pattern looks visually strong.

The implementation must scan figure text, captions, table headings, filenames, and report
prose for the strings `predict`, `prediction`, `transfer`, `held-out`, `precision`, `FDR`,
and `p-value`. Those strings are forbidden in scientific results. They may appear only in
the fixed claim-boundary section titled “What Arm A no longer claims.” A failed scan blocks
publication of the descriptive package.

### What Arm A no longer claims

- Arm A no longer claims that dCas13d guide activity predicts ASO splice modulation.
- Arm A no longer claims a cross-effector transfer law or a kernel validated on held-out ASOs.
- Arm A no longer claims that a correlation, p-value, precision, or FDR from the historical walks establishes predictive performance.
- Arm A no longer claims that 115 exact ASO windows are 115 independent observations or satisfy the minimum effective-sample gate.
- Arm A no longer claims that Hua 2007, Hua 2008, and Singh 2010 form one poolable numeric response.
- Arm A no longer claims that descriptive same-direction co-localization establishes mechanism, causality, portability beyond this reporter field, or clinical utility.

## 10. Revival condition for confirmatory modeling

The current numerically recoverable subset has effective N 10.9. A future data package may
reopen confirmatory Arm A work only after all of the following are true:

1. **Minimum reopening threshold:** the locked spatial effective-N calculation reaches at
   least 25 using numeric ASO outcomes. This requires at least **14.1 net-new effective
   observations** beyond the current estimate, rounded operationally to **15 net-new
   effective spatial units**. Raw one-base tiling windows do not count one-for-one.
2. **Powered target:** because the contest estimated that the frozen rho 0.50 target needs
   about 37 effective observations for 80% power at alpha 0.01, full reinstatement requires
   **effective N at least 37 or a prospectively specified simulation demonstrating at least
   80% power**. The acquisition target should therefore be at least **27 net-new effective
   units**, with a planning buffer of **30**.
3. **Independent datasets:** the new units must come from at least **two genuinely new
   experimental datasets**, not redigitization or technical replication of Hua/Singh. Each
   dataset must cover at least two of I6, E7, and I7 so study identity is not identical to
   region.
4. **Numeric endpoint:** provide per-ASO delta-PSI or percent exon-7 inclusion, replicate
   values, sample size, dose, time point, and an untreated or matched-control baseline.
5. **Assay harmonization:** use the same reporter or cellular system, endpoint definition,
   dose schedule, and measurement protocol across new datasets. If systems differ, include
   a shared reference ASO panel in every system and pre-specify the bridging model before
   viewing dCas13d alignment.
6. **Chemistry control:** use one ASO chemistry/backbone for the primary dataset, or use a
   crossed design in which chemistry is not confounded with target region and every
   chemistry has the same reference panel.
7. **Spatial independence:** distribute new targets across the field and, preferably,
   additional independently assayed SMN2 constructs. Estimate effective N with the locked
   autocorrelation method; do not substitute exact-sequence count or pairwise packing after
   seeing outcomes.
8. **Prospective power check:** freeze exclusions, endpoint, spatial null, and success bar
   before aligning the new ASO outcomes with the dCas13d field.

Reaching 25 effective observations permits reevaluation; it does not by itself prove
adequate power. Until the powered target and all harmonization conditions are met, the
descriptive claim boundary remains in force.

## 11. Clean-room and freeze-safe constraints

- Do not read from `feng/terminal/`, directly or indirectly. The descriptive task does not
  require terminal contrasts, labels, ranks, or sealed outputs.
- Do not import, copy, execute, or path-reference code from repository-root `src/` or
  `pre/src/`.
- Do not import `design_engine`, `offtarget`, `physics_adapter`, or `design_panel`.
- Derive any implementation only from the public source files named in Section 2 and new
  clean-room code under `feng/`.
- The data CSVs are read only. Generated files must go only to the new descriptive output
  directory.
- The implementation must remain reproducible without unsealing `feng/terminal/`.

## 12. Acceptance criteria for a future implementation

1. Both figures share the inclusive 260-434 S4 coordinate axis and mark boundaries at
   319/320 and 374/375 plus the artificial G at 356.
2. Figure 1 contains every retained SMN2 guide and all 125 ASO study records in three
   separate study lanes.
3. The window table has exactly 125 rows and preserves all source records and duplicate
   coordinates.
4. The block table has exactly 33 rows and uses the fixed block boundaries in Section 6.
5. All six ASO direction classes are reported; no unresolved or discordant record is
   silently recoded or omitted.
6. All local dCas13d summaries use fully containing guides, while the positional track uses
   all guides covering each base. The two denominators are labeled separately.
7. No pooled cross-study effect statistic, correlation, hypothesis test, uncertainty
   interval, or inferential significance language appears.
8. Every result statement is labeled `candidate` or `paper_supported`.
9. The fixed “What Arm A no longer claims” section appears verbatim in the narrative
   report.
10. The forbidden-language scan passes outside that fixed section.
11. Input hashes, row counts, coordinate span, quantile convention, and software versions
    are recorded in the report.
12. A file-access audit confirms no read from `feng/terminal/` and no import or path
    reference from repository-root `src/` or `pre/src/`.

## 13. Explicitly out of scope

- Digitizing the 86 potentially numeric ASO windows.
- Fitting any field, shadow, kernel, regression, classifier, or dose-response model.
- Computing a pooled or study-specific inferential correlation.
- Permutation, bootstrap, p-value, FDR, confidence interval, or power analysis.
- Selecting ASO candidates or ranking ASO target windows.
- Updating the standing RED decision, README, CONTEST, HANDOFF, or either input CSV.
- Reading sealed terminal data or reusing code from the earlier arm.
