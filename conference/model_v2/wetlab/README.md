# Exon-3 exploratory workbook audit

Source: user-supplied `SMN2_Exon3_Complete_Data.xlsx`, received date recorded in
[`../ops/WETLAB_MANIFEST.json`](../ops/WETLAB_MANIFEST.json). The source is opened
read-only and never saved. Exon-7 data, computational features, model training,
tuning, p-values, AUC, and prospective claims are excluded.

The workbook reports human SMN2Sup super-minigene RT-PCR/PAGE densitometry,
18-mer ASOs at 5 nM, and one measurement per ASO. These are **observed workbook
statements**, not independently established experimental facts. Raw exports are
the workbook's transcriptions of PDFs; original PDFs/images were not supplied
for this audit. Workbook creation time is not the experiment or receipt date.

Reproduce using:

```bash
MPLCONFIGDIR=/tmp/wetlab-exon3-mpl /home/dashii/orcd/pool/feng-env/bin/python \
  conference/model_v2/wetlab/exon3_parse.py \
  /home/dashii/smt-isef/conference/model_v2/data/wetlab/SMN2_Exon3_Complete_Data.xlsx \
  --output /home/dashii/orcd/scratch/feng/model_v2/wetlab_exon3/
```

**Execution deviation:** the requested scratch directory could not be created
(`Read-only file system`). This session wrote the artifacts to the ignored
[`outputs/`](outputs/) directory instead. The same command with
`--output conference/model_v2/wetlab/outputs` reproduces the actual location.
No commits, Slurm submissions, or changes to running jobs were made.

Artifacts are `audit.json` (every nonempty source cell, metadata, raw counts,
source SHA256, evidence schema, audit and confirmation list), `positional.tsv`,
`psi_audit.tsv`, and `positional.png`/`.svg`. TSV blanks mean null/undefined,
never zero; consult status columns and the JSON evidence schema. Figures plot
the same evidence-classed table values; horizontal segments show labeled target
intervals, not uncertainty intervals. Controls with undefined fractions are
omitted. The schema assigns every table column an evidence class; all cells and
metadata retain their workbook source. MW values are proxies, not verified bp.

For each gel/lane, Band 1 and Band 2 volumes are independently grouped from raw
rows; ladder rows are inventoried but excluded from calculation. Both bands must
have finite positive volumes to calculate `100*I/(I+E)`. Missing, zero, negative,
or nonnumeric bands yield null, never an imputed endpoint. Sample band duplicates,
unknown band numbers and unmatched lane sets fail closed. A restricted arithmetic
reader evaluates workbook PSI references separately and compares formulas and
caches with the raw result (tolerance 0.000001 percentage points); exported
3-decimal ratios use 0.0005 pp rounding tolerance. Unsupported formulas are
reported, never executed. The source formulas/caches are preserved verbatim.

The two control lanes have no Band 2 export row. The workbook states that export
completeness was confirmed and uses zero, producing 100% PSI. Completeness does
not establish zero skipped transcript. This audit preserves the workbook's
100% values as workbook-only values and leaves both independent fractions null.
For the 11 ASOs with both bands reported, all 22 PSI cells (Quantification and
Summary) match independent recomputation, and their cached values agree. All
26 PSI caches, including the two controls in both sheets, match their formulas.
The 11 ASO exported Band 1 ratios agree within rounding tolerance.

Exploratory values below are **derived conditional signal fractions**, assuming
the supplied Band 1=inclusion/Band 2=skipping assignments; they are not validated
transcript-molecule inclusion estimates. ASO positions and gel labels are
**observed, unverified workbook annotations**. No control-normalized deltas are
estimated; gel 7-9 has no on-gel controls.

| ASO | Gel | Labeled exon-3 nt interval | Conditional signal fraction (%) |
| --- | --- | --- | ---: |
| A | 7-7 | 2–19 | 60.647 |
| B | 7-7 | 16–33 | 59.986 |
| C | 7-7 | 30–47 | 58.565 |
| D | 7-7 | 43–60 | 43.952 |
| E | 7-7 | 57–74 | 37.294 |
| F | 7-9 | 71–88 | 52.463 |
| G | 7-9 | 87–104 | 8.024 |
| H | 7-9 | 103–120 | 10.469 |
| I | 7-9 | 119–136 | 96.436 |
| J | 7-9 | 135–152 | 41.622 |
| K | 7-9 | 151–168 | 41.776 |

All 11 labeled intervals have an internally consistent inclusive length of 18;
that arithmetic does not verify their sequence, strand, or coordinate origin.
Band 1 migrates more slowly and has a larger reported MW than Band 2 in all
11 pairs, consistent with the supplied size order but insufficient to establish
isoform identity. Exported MW proxies range from 547.923–569.977 (Band 1) and
285.788–294.360 (Band 2) on 7-7, versus 587.008–609.070 and 322.032–338.495 on
7-9. Pair separations range from 260.200–277.344 across both gels. These are
**derived ranges of unverified sizing proxies**: both different product lengths
and cross-gel calibration shifts require investigation. Primer-derived expected
sizes are null because neither primer nor construct sequences are present.
No fragment-length correction is applied; integrated signal fractions can differ
from molecular fractions when products have different lengths.

`needs_human_confirmation` in the JSON and manifest lists the unresolved items:
ASO sequences/orientation/chemistry, construct/version/exon origin, primer binding
sites/expected products, band identities, size shifts, original images/PDFs,
detection limits, experimental-unit lineage and biological/technical replication,
gel 7-9 controls, and coordinate compatibility. Biological and technical replicate
IDs remain separate null fields. Gel/lane identifiers identify measurements only;
they are not assigned either kind of replicate. No mean, SD or error bars are
justified by the workbook's reported n=1.

Computational comparison: **not comparable**. The assay is reported as a
super-minigene, and construct-to-genome coordinates, assembly, sequence and
strand are unverified. No feature join was attempted.

Follow-up validation: obtain the complete construct and ordered ASO/primer
sequences; align and verify all binding sites and expected amplicons; confirm
products with junction sequencing and size anchors; inspect original gel images
for background, saturation and detection thresholds; resolve size shifts; repeat
ASOs with independent biological samples and separately identified technical
repeats, randomized gel placement and matched controls on every gel. Establish
assay linearity and a defensible length-correction protocol before quantitative
transcript PSI interpretation. None of these recommendations retroactively
validate the present measurements.

Validation: the synthetic parser suite passed **20 tests, 0 skipped, 0 failed**.
The required feng gate returned **581 passed, 5 skipped, 1 failed** in 242.36 s.
`feng/tests/test_cleanroom.py::test_implementation_has_no_legacy_reference[source_file42]`
rejects `src/` substrings in unchanged `feng/src/feng/source_auth.py` (lines 31,
314 and 400). That file is outside this task's ownership and was left untouched;
the repository-wide gate remains red.
