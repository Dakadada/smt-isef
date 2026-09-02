# Arm D positives — the 8 locked nusinersen off-target exons (GENCODE v50)

> Claim boundary: this file locks the *coordinates* of the 8 experimentally supported
> off-target skipped exons for the Arm D benchmark. Allowed honesty states here:
> `paper_supported` (the event list) and `candidate` (my exon-to-v50 coordinate derivation).
> No `clean_reanalysis_reproduced` claim is made; no ranking or effect is computed here.

## Source [paper_supported]

Ottesen, Murzyn, Kaas, Bertrand, Payne, Singh. *"A therapeutic antisense oligonucleotide
encompassing 2′-O-methoxyethyl modification triggers unique perturbation of the transcriptome."*
NAR Molecular Medicine, 2026 (pub. 2026-01-06). PMID 41551726 / PMCID PMC12805893. GEO GSE307431.

The paper validated splicing of **eight randomly selected exons whose skipping was enhanced by
F18MOE** (F18MOE = sequence- and chemistry-identical to nusinersen). Source: **Figure 3B** and
accompanying text ("we validated the alternative splicing of eight randomly selected exons whose
skipping was enhanced by F18MOE"). The 8 events, verbatim from Fig 3B:

POLR2H exon 2 · PITHD1 exon 2 · SERPINB7 exon 4 · RTTN exon 4 · REV3L exon 9 · PRKRA exon 2 ·
GOLGA4 exon 4 · PAK1 exon 2.

POLR2H exon 2 is the previously known nusinersen off-target — its presence is expected; the Arm D
evaluation is **locked, not blind** (this positive is known to us), with a mandatory
POLR2H-excluded sensitivity analysis downstream (plan line 151).

## Coordinate derivation [candidate]

**The paper gives no genomic coordinates** for the individual exons (Fig 3B is a validation
panel; reference build is GRCh38). Coordinates were therefore derived, not transcribed:

- One bounded single pass over the locked annotation
  `$SCRATCH/feng/gencode/v50/gencode.v50.primary_assembly.annotation.gtf` extracting the 8
  `gene_name` records.
- For each gene, Figure 3B's cassette-exon length and immediate flanking-exon/intron layout were
  matched against rMATS v50 events, then reconciled to a GENCODE v50 transcript. Seven map to MANE
  Select transcripts: POLR2H paper exon 2 is MANE exon 3, REV3L paper exon 9 is MANE exon 8, and the
  other five retain the paper ordinal. GOLGA4 instead maps to non-MANE transcript ENST00001121839.1:
  its 66-nt exon 4 uniquely reproduces the displayed `90-15357-99-6595-66-8356-315` chain. Output
  records the v50 ordinal and preserves the paper label in `event_id` and `notes`.

## Coordinate table

| event_id | gene | v50 transcript | exon | chrom | start–end (1-based) | strand |
|---|---|---|---|---|---|---|
| OT01 | POLR2H | ENST00000456318.6 | 3/6 | chr3 | 184364966–184365049 | + |
| OT02 | PITHD1 | ENST00000246151.9 | 2/6 | chr1 | 23779438–23779481 | + |
| OT03 | SERPINB7 | ENST00000398019.7 | 4/8 | chr18 | 63793161–63793277 | + |
| OT04 | RTTN | ENST00000640769.2 | 4/49 | chr18 | 70201894–70201983 | − |
| OT05 | REV3L | ENST00000368802.8 | 8/32 | chr6 | 111388001–111388085 | − |
| OT06 | PRKRA | ENST00000325748.9 | 2/8 | chr2 | 178450242–178450411 | − |
| OT07 | GOLGA4 | ENST00001121839.1 (non-MANE) | 4/16 | chr3 | 37273536–37273601 | + |
| OT08 | PAK1 | ENST00000356341.8 | 2/15 | chr11 | 77392331–77392541 | − |

## Map-to-v50 gate (plan line 153)

**PASS — 8/8.** Every positive maps unambiguously to a GENCODE v50 exon: seven on MANE Select
transcripts and GOLGA4 on the non-MANE primary-annotation transcript ENST00001121839.1. No positive
is unmappable, so the KILL condition ("any positive unmappable to the locked GENCODE universe") is
**not** triggered.

## Identity verification against Figure 3B

Figure 3B reports the cassette and neighboring exon/intron sizes. The v50 mappings reproduce its
eight cassette lengths exactly: POLR2H 84 nt, PITHD1 44 nt, SERPINB7 117 nt, RTTN 90 nt, REV3L 85
nt, PRKRA 170 nt, GOLGA4 66 nt, and PAK1 211 nt. The displayed immediate intron layouts also match.
This check resolves the POLR2H and REV3L transcript-numbering differences and corrects GOLGA4 from
an unrelated 48-nt MANE exon to the Figure 3B-anchored non-MANE 66-nt exon; the other five rows are
unchanged.

The paper does not print genomic coordinates, so the numeric intervals remain `candidate` mappings.
Coordinates are 1-based inclusive GTF coordinates on the primary assembly (chr* naming), matching
the STAR index and v50 universe used elsewhere in feng.

## Downstream

These 8 positives feed Arm D matched-stratum construction (`feng/src/feng/bench/strata.py`), which
must build ≥50 eligible background exons per positive from the **control-only** detectable universe
(no treatment data) — that is the separate ≥50/positive gate (plan line 153), checked once the
control-only universe lands.

**Coordinate-frame contract (read before any positive↔universe join).** The `start`/`end` columns
here are **1-based inclusive** (interval length = `end − start + 1`). The detectable universe keys
each exon in the rMATS **0-based half-open** frame (`exon_start_0base`), so its `exon_key` is
`chrom:(start−1)-end:strand`. Any coordinate join between this table and the universe **must** pass
through `feng.bench.strata.positive_interval_to_exon_key`, which drops exactly one base at the start.
Joining on the raw 1-based `start` shifts every positive by one base and silently yields zero
matches against the ≥50-background gate.
