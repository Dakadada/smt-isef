# HapSMA per-haplotype extraction: nusinersen target-site variants

**Audit date:** 2026-09-02  
**Operational verdict:** **NO-GO** for the natural-variant Appendix on the staged HapSMA evidence.  
**Measured qualifying donor-independent target-site ALT count:** **0** (GO requires at least 2 unrelated donors).

This is an evidence-gate result, not a claim that every assayed SMN2 copy is reference across the footprint and not a biological claim that target-site variation is absent. The variants-only workbook reports no target-overlapping ALT, but it does not expose the per-base callability or read metrics needed to turn blank cells into reference calls.

## Inputs and scope

The audit read the staged files directly and did not use the network:

| File | Role | SHA-256 |
|---|---|---|
| `MOESM3.xlsx` | Public HPRC read-path/sample aid | `f6e151e74de9b2632bdde9c86ab63bcab314d65fc446c667339cc307e917fb40` |
| `MOESM4.xlsx` | Per-haplotype copy type and PSV states | `e6f15deb71991785ee3ae97a1b75ccc7c9179ed17d2cefb52eb1ad7ecead82eb` |
| `MOESM5.xlsx` | Table S3, variant-by-haplotype matrix; decisive table | `f2eac5310dd42d00ff91e3fe2269eea010c729f930ba842e9189d5e6e08bb20b` |
| `MOESM6.xlsx` | PSV/copy-specific coordinate and aggregate-count aid | `437f38bbe6930c0139bce21abd97f1285ce1c3d85c766bc2e99ec6b3f0e78c9a` |
| `MOESM7.xlsx` | GRCh38-to-T2T coordinate pairs | `6802b831aae9d3700e1e88ae56243af0b806fef557a1a744dfe631812aab6bea` |

The repository's stated `feng/CLAUDE.md` path was absent; the root `CLAUDE.md`, `feng/HANDOFF.md`, and the authoritative phase-2 README were read and followed. No terminal results or excluded legacy/source packages were inspected.

## Coordinate and workbook-semantic reconciliation

The locked drug ASO is `5'-TCACTTTCATAATGCTGG-3'`; its plus-strand target is `5'-CCAGCATTATGAAAGTGA-3'`, intron 7 +10..+27.

| Frame / evidence | PSV13 (`c.840C>T`) | 18-nt target, 1-based inclusive | Reconciliation |
|---|---:|---:|---|
| GRCh38 SMN2 | `chr5:70076526` | `chr5:70076584-70076601` | Locked sequence-verified SMN2 interval; target offsets from PSV13 are +58..+75. |
| Paraphase GRCh38 SMN1 frame | `chr5:70951946` | `chr5:70952004-70952021` | The same local +58..+75 offsets in the SMN1 realignment frame. |
| HapSMA T2T-CHM13 SMN1 frame | `chr5:71408734` | **`chr5:71408792-71408809`** | MOESM6 identifies `71408734` as PSV13; adding +58..+75 yields the locked target. Zero-based half-open form: `chr5:71408791-71408809`. |

The crosswalk is independently supported by the workbooks rather than inferred from one offset:

- MOESM7 pairs GRCh38 SMN1 `70951845` with T2T `71408633`, a local offset of +456,788. Applying that local offset to `70952004-70952021` gives exactly `71408792-71408809`. This is used only as a local transform; MOESM7 shows an offset of +456,789 at its downstream `70952493 -> 71409282` pair, consistent with an intervening indel and warning against a locus-wide constant offset.
- MOESM6 places PSV13 at T2T `71408734` and PSV14 at `71408882`. Thus the target is 58 bases after PSV13 and ends 73 bases before PSV14.
- MOESM5 uses `pos_ref_alt` identifiers. It contains the exact PSV rows `71408734_C_T` and `71408882_A_G`. Its `71408734_C_T` carrier columns match, with no extras or omissions, the MOESM4 haplotypes whose PSV13 state is `2` and whose copy type is SMN2: 93 SMA haplotypes and 49 control haplotypes. MOESM6 independently reports 49/49 control SMN2 haplotypes ALT and 0 NA at PSV13. This exact match confirms the MOESM5 coordinate/allele frame and connects its columns to the SMN2 classification.
- Every MOESM5 haplotype column has an exact same-named MOESM4 sample/haplotype row (104/104 in `SMA`; 116/116 in `1000G`). All nonblank genotype-matrix cells in MOESM5 are `1`, i.e. ALT-presence flags; the workbook has no explicit `0/0`, reference-genotype, or no-call encoding.

## Target-overlap extraction

Every MOESM5 row in both sheets was parsed as `POS_REF_ALT`. For SNVs and deletions, overlap was tested against the inclusive reference span `[POS, POS + length(REF) - 1]`; insertion anchors were tested at `POS`. This matters because the workbook contains long REF alleles (maximum observed REF length 7,119 nt in `SMA`), so querying only variants whose POS lies inside the target would be unsafe.

**Result: no MOESM5 variant row overlaps `chr5:71408792-71408809` in either sheet.**

Boundary checks:

| Sheet | Nearest upstream reported variant | Distance to target start | Nearest downstream reported variant | Distance from target end |
|---|---|---:|---|---:|
| `SMA` | `71408744_CA_C`, REF span `71408744-71408745`, carried by `SMA_09_hap1` | 47 nt | `71408882_A_G` (PSV14) | 73 nt |
| `1000G` | `71408734_C_T` (PSV13) | 58 nt | `71408862_T_G`, carried by `HG02615_hap1` and `HG02809_hap1` | 53 nt |

Neither the nearby patient deletion nor either control variant touches a target nucleotide.

## Per-donor evidence

There are no qualifying target-ALT records, so the qualifying-evidence table is intentionally empty rather than populated with inferred reference alleles:

| Donor ID | Copy | REF | ALT | Coordinate frame | DP / AF / QC | PSV evidence |
|---|---|---|---|---|---|---|
| *No qualifying donor* | — | Not asserted | No target-overlapping ALT row | HapSMA T2T-CHM13 SMN1 frame, target `chr5:71408792-71408809` | MOESM5 supplies no DP, AF, FILTER, phase-quality, or per-base callability fields | No target hit to validate; SMN2 mapping was nevertheless cross-checked through MOESM4 PSV13 as described above |

For audit completeness, the following table enumerates all SMA donors and their MOESM4-classified SMN2 copies. “None recorded” means only that no MOESM5 ALT row overlaps the footprint; it must not be read as a per-copy REF call.

| Donor ID | SMN2 copy/copies (MOESM4) | Target REF | Target ALT in MOESM5 | DP / AF / target callability | PSV13 evidence |
|---|---|---|---|---|---|
| SMA_01 | hap1, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on both copies |
| SMA_03 | hap1, hap2, hap4 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_04 | hap1, hap2, hap3, hap4 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_05 | hap2, hap3, hap4 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_06 | hap1, hap2, hap3, hap4 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_09 | hap1, hap2, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_109 | hap1, hap4, hap5 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_11 | hap1, hap2, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_12 | hap1, hap2, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_13 | hap1, hap2, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_14 | hap1, hap2, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_15 | hap1, hap2, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_153 | hap1, hap2, hap3, hap4 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_17 | hap1, hap2, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_18 | hap1, hap2, hap3, hap4 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_19 | hap1, hap2, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_20 | hap1, hap2, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_23 | hap1, hap2, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_27 | hap1, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on both copies |
| SMA_28 | hap1, hap2 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on both copies |
| SMA_283 | hap1, hap2, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_29 | hap1, hap2, hap3, hap4 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_30 | hap1, hap2 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on both copies |
| SMA_323 | hap1, hap3, hap4 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_381 | hap1, hap2, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_44 | hap1, hap2, hap4 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_61 | hap2, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on both copies |
| SMA_63 | hap1, hap2, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_71 | hap1, hap2 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on both copies |
| SMA_89 | hap1, hap2, hap3, hap4 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |
| SMA_97 | hap1, hap2, hap3 | Not asserted | None recorded | Not provided | `2` / `71408734_C>T` on all listed copies |

Donor independence does not alter the count: no donor has even one reported target-site ALT, so duplicated copies within a donor cannot contribute.

## Denominators, controls, and limitations

| Workbook cohort | Donors represented | MOESM4-classified SMN2 haplotypes | Other haplotypes | Callable target denominator |
|---|---:|---:|---:|---:|
| `SMA` | 31 | 93 | 2 SMN1; 9 unclassified | **Not determinable** |
| `1000G` / HPRC controls | 29 | 49 | 57 SMN1; 10 unclassified | **Not determinable** |

The 31 SMA donors and 93 classified SMN2 copies are a **representation denominator**, not a callable denominator. MOESM4's named haplotypes and `region`/`bed` phasing-method fields support the paper's copy-resolved long-read assignments, and its PSV13 observations support SMN2 assignment 58 bases upstream. However, neither MOESM4 nor MOESM5 establishes that every one of the 18 target bases passed a callable threshold on each copy. MOESM5 has no DP/AF/QC fields and uses blank cells both where an ALT is not reported and where no explicit reference/no-call distinction is available. MOESM6 supplies aggregate REF/ALT/NA counts at selected PSVs for the control cohort, not per-base target callability or patient DP/AF.

The `1000G` sheet sample IDs are `HG...`, and MOESM3 points to Human Pangenome/HPRC long reads. These 29 donors are controls, **not SMA patients**. They were audited for target overlap but are not counted as patients and cannot support a patient-response claim.

To obtain a callable conservation denominator or to test whether variants-only-table blanks conceal missed calls would require per-copy gVCF/callability output or read-backed bamlets spanning the full footprint, with DP, allele fraction/read support, mapping/strand review, and retained phasing to PSV13 or equivalent paralog-discriminating positions. The separately staged Paraphase bamlets could strengthen the HPRC control audit, but they remain non-patient evidence; the HapSMA patient raw reads are controlled. These missing data block a firm **reference/conservation** claim, but they do not create any qualifying positive ALT in the published workbook.

## Gate decision and honesty state

- Qualifying donor-independent SMN2 target-site ALTs: **0**.
- GO threshold: **at least 2 unrelated donors**.
- Decision: **NO-GO** for the narrow natural-variant model-decomposition arm on the staged HapSMA evidence. The complete variant-row extraction has no target-overlapping candidate to assign, validate, or deduplicate, so the positive-evidence threshold is not met.
- This is not `still-insufficient` for the operational GO/NO-GO gate: lack of a callable denominator cannot convert an absent ALT record into qualifying positive evidence. It **is** insufficient for a biological zero-variant statement, a target-site conservation fraction, or a patient-response claim.
- Honesty state: workbook availability, haplotype/copy assignments, coordinate anchors, and the published variant records are `paper_supported`. No target ALT candidate was found. No result is labeled `clean_reanalysis_reproduced`, and no modeled or biological interaction is claimed.
