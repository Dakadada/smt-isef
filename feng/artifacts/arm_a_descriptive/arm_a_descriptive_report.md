# Arm A descriptive SMN2 field atlas

> Claim boundary: descriptive and hypothesis-generating only. Allowed result states in this package are `candidate` and `paper_supported`.

## Manifest and coordinate reconciliation

- [paper_supported] Input ASO inventory: 125 records / 115 exact windows / 115 guide-covered exact windows. The exact-window count has no independence interpretation.
- [paper_supported] The source audit found that the earlier ASO artifact was nine nucleotides below the SpliceRUSH S4 field frame because it omitted the S4 5-prime prefix `ATGGCTAGC`. This package uses the corrected input coordinates: I6 maps as `320 + q` (within I6 the offset q is negative, placing I6 positions below the E7 start at 320), natural E7 starts at 320 with the artificial G at 356, and I7 maps as `374 + q`.
- [paper_supported] Sequence-alignment evidence recorded in the source memo: guide 1001 contains the prefix before the Dataset S1 exon-6 sequence; guide 1340 aligns at S1-concatenated 331-355 and S4 340-364; guide 1384 aligns at S1-concatenated 375-399 and S4 384-408. The generator independently recomputes all 125 corrected hulls and segments from paper-local coordinates and stops on any mismatch.
- [candidate] Display span: inclusive corrected S4 positions 260-434; construct regions I6 260-319, E7 320-374, and I7 375-434.
- [paper_supported] dCas13d sign semantics are verified: positive combined Z denotes exon-7 inclusion activation and negative combined Z denotes inhibition. Source: Recinos et al. 2024, Nature Communications, Figure 5a and Methods equations 3-5; https://doi.org/10.1038/s41467-024-47140-8.
- [candidate] Quantiles use NumPy method `linear` at every base, window, and coarse display block. There is no smoothing or interpolation.
- [candidate] Input SHA-256: dCas13d CSV `8166c7609d7817a757cc8000978e001b0f0362d9c98368c2088748c725558db4`; ASO CSV `b69db791747cfa060d19fcf4163f86d56b26ffc00d2d7f3213aea5f0a109b549`.
- [candidate] Input rows: 1937 total dCas13d rows, 395 retained SMN2 guides intersecting 260-434, and 125 ASO study records.
- [candidate] Software: Python 3.11.16, NumPy 1.26.4, pandas 3.0.5, matplotlib 3.10.9.
- [candidate] File-access audit: the generator reads only the two named CSV inputs. It does not read sealed terminal material or reuse earlier-arm code.
- [candidate] Block-edge assumption: the specification gives integer-inclusive blocks but also requires unique assignment of half-integer hull midpoints. Blocks are therefore interpreted with base-centered, left-closed/right-open edges; a midpoint exactly at 334.5 enters E7-2 rather than E7-1, and likewise at every adjacent boundary. The final right edge is closed.

## Methods

The raw guide scatter uses interval midpoints only for point placement. At each integer coordinate, the central track and quartile ribbon use every guide whose full inclusive interval covers that base. ASO window summaries instead use only guides that fully contain the inclusive ASO hull. Targets crossing the artificial exon-7 G are drawn as their two native segments, while containment uses the conservative hull. ASO directions come only from the fixed categorical input field; magnitude prose is carried as annotation and is not converted to numbers.

Each ASO study record is assigned once to one of eleven fixed coarse display blocks by its hull midpoint. Study lanes and summary rows remain separate. Duplicate coordinates in different papers remain separate records. No effect direction is averaged or voted across papers.

## Inventory

### ASO records and exact windows by study and region [paper_supported]

| source | region | records | exact_windows |
| --- | --- | --- | --- |
| Hua 2007 PLoS Biology | I6 | 0 | 0 |
| Hua 2007 PLoS Biology | E7 | 48 | 48 |
| Hua 2007 PLoS Biology | I7 | 0 | 0 |
| Hua 2008 AJHG | I6 | 27 | 27 |
| Hua 2008 AJHG | E7 | 0 | 0 |
| Hua 2008 AJHG | I7 | 31 | 31 |
| Singh 2010 RNA | I6 | 0 | 0 |
| Singh 2010 RNA | E7 | 0 | 0 |
| Singh 2010 RNA | I7 | 19 | 19 |

### Direction classes and non-evaluable records by study and region [paper_supported]

| source | region | inclusion increase | inclusion decrease | neutral | assay-discordant | mixed/uncertain | unresolved, not digitized | non_evaluable | study_region_records_denominator |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Hua 2007 PLoS Biology | I6 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Hua 2007 PLoS Biology | E7 | 7 | 15 | 3 | 7 | 0 | 16 | 26 | 48 |
| Hua 2007 PLoS Biology | I7 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Hua 2008 AJHG | I6 | 1 | 3 | 0 | 0 | 0 | 23 | 23 | 27 |
| Hua 2008 AJHG | E7 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Hua 2008 AJHG | I7 | 9 | 1 | 0 | 0 | 0 | 21 | 21 | 31 |
| Singh 2010 RNA | I6 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Singh 2010 RNA | E7 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Singh 2010 RNA | I7 | 11 | 1 | 0 | 0 | 7 | 0 | 7 | 19 |

### Retained guides by source annotation and length [paper_supported]

| region | length | guides | retained_guides_denominator |
| --- | --- | --- | --- |
| E7 | 22 | 31 | 395 |
| E7 | 25 | 31 | 395 |
| E7I7 | 22 | 21 | 395 |
| E7I7 | 25 | 24 | 395 |
| I6 | 22 | 60 | 395 |
| I6 | 25 | 60 | 395 |
| I6E7 | 22 | 24 | 395 |
| I6E7 | 25 | 24 | 395 |
| I7 | 22 | 60 | 395 |
| I7 | 25 | 60 | 395 |

### Per-base guide coverage by construct region [candidate]

| region | positions | coverage_median | coverage_min | coverage_max |
| --- | --- | --- | --- | --- |
| I6 | 60 | 47.0 | 47 | 47 |
| E7 | 55 | 47.0 | 47 | 47 |
| I7 | 60 | 47.0 | 47 | 47 |

### Exact-window coverage [candidate]

- 115 / 115 distinct exact windows have at least one fully containing retained guide.

### Block-by-study descriptive sign counts [candidate]

The denominator is the number of `increase` or `decrease` records in that named study and block. `same_direction_fraction` is `same_n / eligible_direction_n`; zero-denominator cells are explicit `NA`.
Two distinct reference medians appear in this package and are not interchangeable: the *window-local containing-guide median* (`local_z_median` in the window-alignment table — the median combined Z over the guides that fully contain an individual ASO hull) governs each window's own `direction_relation`, whereas the *coarse-block positional median* (`block_z_positional_median` — the median of the per-base positional Z medians across a display block) governs the same/opposite counts in the block-by-study table below.
When the block median is exactly zero, eligible ASOs are in the denominator but are neither same nor opposite; the printed fraction is therefore zero and both sign counts are zero.

| block_id | study | eligible_direction_n | same_n | opposite_n | same_direction_fraction |
| --- | --- | --- | --- | --- | --- |
| I6-1 | Hua 2007 PLoS Biology | 0 | 0 | 0 | NA |
| I6-1 | Hua 2008 AJHG | 1 | 0 | 0 | 0.0 |
| I6-1 | Singh 2010 RNA | 0 | 0 | 0 | NA |
| I6-2 | Hua 2007 PLoS Biology | 0 | 0 | 0 | NA |
| I6-2 | Hua 2008 AJHG | 0 | 0 | 0 | NA |
| I6-2 | Singh 2010 RNA | 0 | 0 | 0 | NA |
| I6-3 | Hua 2007 PLoS Biology | 0 | 0 | 0 | NA |
| I6-3 | Hua 2008 AJHG | 1 | 0 | 0 | 0.0 |
| I6-3 | Singh 2010 RNA | 0 | 0 | 0 | NA |
| I6-4 | Hua 2007 PLoS Biology | 0 | 0 | 0 | NA |
| I6-4 | Hua 2008 AJHG | 2 | 0 | 0 | 0.0 |
| I6-4 | Singh 2010 RNA | 0 | 0 | 0 | NA |
| E7-1 | Hua 2007 PLoS Biology | 5 | 3 | 2 | 0.6 |
| E7-1 | Hua 2008 AJHG | 0 | 0 | 0 | NA |
| E7-1 | Singh 2010 RNA | 0 | 0 | 0 | NA |
| E7-2 | Hua 2007 PLoS Biology | 5 | 5 | 0 | 1.0 |
| E7-2 | Hua 2008 AJHG | 0 | 0 | 0 | NA |
| E7-2 | Singh 2010 RNA | 0 | 0 | 0 | NA |
| E7-3 | Hua 2007 PLoS Biology | 12 | 8 | 4 | 0.6666666666666666 |
| E7-3 | Hua 2008 AJHG | 0 | 0 | 0 | NA |
| E7-3 | Singh 2010 RNA | 0 | 0 | 0 | NA |
| I7-1 | Hua 2007 PLoS Biology | 0 | 0 | 0 | NA |
| I7-1 | Hua 2008 AJHG | 1 | 1 | 0 | 1.0 |
| I7-1 | Singh 2010 RNA | 5 | 0 | 5 | 0.0 |
| I7-2 | Hua 2007 PLoS Biology | 0 | 0 | 0 | NA |
| I7-2 | Hua 2008 AJHG | 9 | 9 | 0 | 1.0 |
| I7-2 | Singh 2010 RNA | 7 | 6 | 1 | 0.8571428571428571 |
| I7-3 | Hua 2007 PLoS Biology | 0 | 0 | 0 | NA |
| I7-3 | Hua 2008 AJHG | 0 | 0 | 0 | NA |
| I7-3 | Singh 2010 RNA | 0 | 0 | 0 | NA |
| I7-4 | Hua 2007 PLoS Biology | 0 | 0 | 0 | NA |
| I7-4 | Hua 2008 AJHG | 0 | 0 | 0 | NA |
| I7-4 | Singh 2010 RNA | 0 | 0 | 0 | NA |

## Figure captions

- [candidate] `arm_a_shared_field_map`: dCas13d values are numeric source measurements, while ASO colors are source-reported qualitative classes. The studies differ in assay and chemistry. Visual proximity or same-direction color is hypothesis-generating only.
- [candidate] `arm_a_coarse_block_map`: eleven fixed coarse display blocks summarize positional Z medians and six ASO direction classes in separate study rows. Blocks are display bins, not independent experiments.

## Limitations

- [paper_supported] ASO outcomes are categorical transcriptions. `figure_reported` is unresolved, not neutral and not assumed randomly missing; assay-discordant and mixed/uncertain records remain separate.
- [candidate] Study identity is confounded with sampled region and assay context. Empty study-region cells remain visible, and no cross-study effect estimate is made.
- [candidate] Local and block sign relations are coordinate juxtapositions, not mechanistic or causal evidence.

## What Arm A no longer claims

- Arm A no longer claims that dCas13d guide activity predicts ASO splice modulation.
- Arm A no longer claims a cross-effector transfer law or a kernel validated on held-out ASOs.
- Arm A no longer claims that a correlation, p-value, precision, or FDR from the historical walks establishes predictive performance.
- Arm A no longer claims that 115 exact ASO windows are 115 independent observations or satisfy the minimum effective-sample gate.
- Arm A no longer claims that Hua 2007, Hua 2008, and Singh 2010 form one poolable numeric response.
- Arm A no longer claims that descriptive same-direction co-localization establishes mechanism, causality, portability beyond this reporter field, or clinical utility.
