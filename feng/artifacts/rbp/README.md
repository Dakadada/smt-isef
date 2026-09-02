# Arm B — ENCODE RBP layer (provenance memo)

> **Claim boundary.** This directory indexes *processed, published* ENCODE data for
> the RBPShadow arm. It contains no interaction result and no concordance statistic;
> the sign-concordance test is a separate, post-freeze step. Provenance state for
> every source here is `paper_supported`. Nothing is `clean_reanalysis_reproduced`.

## What this layer is

The RBPShadow concordance test needs two aligned pieces per RNA-binding protein:

1. **where the RBP binds** — ENCODE eCLIP CLIPper peaks, and
2. **which way its loss moves splicing** — the delta-PSI direction from an RBP
   knockdown RNA-seq experiment (batch-corrected rMATS).

This memo records the resolved ENCODE files for both, a hash-verified curated
subset staged to scratch, and the loaders that parse them. The ASO delta-PSI side
of the concordance test comes from the sealed GSE307431 rMATS output and is not
touched here.

## Sources (all ENCODE `PublicationData`, Van Nostrand et al. RBP map)

| collection | role | files | on-disk index |
|---|---|---|---|
| `ENCSR456FVU` | eCLIP peaks (BED narrowPeak) | 669 | `encode_manifest.tsv` role `eclip_peaks` |
| `ENCSR870OLK` | RBP-KD **batch-corrected** rMATS differential splicing (tar) | 472 | `encode_manifest.tsv` role `kd_splicing_rmats` |

`ENCSR870OLK` (batch-corrected) is used in preference to `ENCSR413YAF` (the
non-batch-corrected rMATS/MISO/DESeq/CUFFDIFF collection) for the knockdown
splicing direction. The eCLIP alignment/read files and the KD gene-quantification
TSVs are intentionally excluded — this layer is peaks + splicing deltas only.

The full manifest lists **1,141** processed files (0.72 GB eCLIP + 1.73 GB rMATS).
Each row carries: `collection, accession, role, target_rbp, cell_line, file_format,
output_type, assembly, needs_liftover, md5sum, sha256, file_size, experiment, url,
curated_subset`. Target RBP and cell line were resolved from each file's ENCODE
experiment (`target.label`, `biosample_ontology.term_name`).

## Coverage

- **203** RBP × cell-line combinations carry *both* an eCLIP peak set and a
  KD-splicing table — the joinable universe for the concordance test.
- Cell lines: HepG2 and K562 (HepG2 is where 48 of Holgersen's ASO experiments ran).

## Curated staged subset (hash-verified on disk)

To prove the fetch + parse pipeline without pulling the full 2.45 GB, a curated set
of known splicing-regulatory RBPs was staged and verified. These are the RBP × cell
combos that (a) carry both modalities and (b) match a known nuclear splicing-RBP
list — it includes the SMN2 exon-7 / ISS-N1 regulators HNRNPA1, TIA1, and TIAL1.

- **152 files, 317.1 MB**, staged to `$SCRATCH/feng/encode/{eclip_peaks,kd_splicing}/<RBP>/<cell>/`.
- **152 / 152 md5-verified** against the ENCODE-reported `md5sum` (0 failures);
  `sha256` recorded as secondary metadata after md5 verification.
- 38 RBP × cell combos (114 eCLIP beds + 38 KD rMATS tars). RBPs: FUS, HNRNPA1,
  HNRNPC, HNRNPL, HNRNPM, HNRNPU, KHSRP, PCBP1, PCBP2, PTBP1, RBFOX2, RBM22, SF3B4,
  SRSF1, SRSF7, SRSF9, TARDBP, TIA1, TIAL1, U2AF1, U2AF2.

`curated_staged.tsv` lists exactly these staged files. The remaining manifest rows
are indexed but not downloaded; their URL/md5 can be used for staging, but a SHA-256
must be recorded in the `EncodeFile` before either coordinate reader will accept the
downloaded bytes.

## Assembly / liftover (open item, deliberately deferred)

**Every ENCODE file in this layer is `hg19`.** The rest of the arm works in the
GENCODE v50 (GRCh38) frame. Coordinate-level use of eCLIP peaks against v50 exon
models therefore requires a lift to GRCh38. Liftover is **not** applied here; the
`needs_liftover` flag is `yes` on every row. Loaded frames retain their true `hg19`
assembly and use an assembly-aware dataframe whose `merge`/`join` methods reject a
GRCh38 frame. They are never stamped as GRCh38 before liftover. rMATS event
coordinates are likewise hg19.

## Loader (`feng/src/feng/rbp/encode.py`)

- `load_manifest(path)` → typed, validated `EncodeFile` rows (unique accessions;
  role vocabulary enforced).
- `read_eclip_peaks(path, source=...)` → SHA-256-verified tidy narrowPeak frame;
  eCLIP semantics on the numeric
  columns (`signal_value` = log2 fold enrichment over size-matched input,
  `neg_log10_p` = enrichment significance).
- `read_kd_splicing_rmats(path, source=...)` → SHA-256-verified
  `{event_type: DataFrame}` for SE/A3SS/A5SS/MXE/RI;
  `IncLevelDifference` is the knockdown-minus-control delta-PSI whose **sign** is the
  direction the concordance test will later use.

The loader parses and exposes only. It computes no concordance and asserts no
interaction. Tests: `feng/tests/test_rbp_encode.py` (synthetic fixtures + manifest
invariants), green under `feng-env`.
