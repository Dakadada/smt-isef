# Arm C gate — bind→act negative-set construction: BLOCKED at Task 1

**Status: NO-GO (blocked on missing input, not on a modeling result).**
Author: Claude (general-purpose subagent). Read-only. Date 2026-09-02.
No sealed `terminal/` or Phase-9 ΔPSI data read. No chemistry model, no
SpliceAI/Pangolin/SPLAM/IntaRNA ensemble invoked.

## Task 1 result: the 75 SBO sequences are NOT present in any file under `feng/`

```
$ ls feng/*.pdf
feng/sj-pdf-4-nat-10.1089_nat.2020.0921.pdf   (S4, 17 pages, 715 KB)
feng/sj-pdf-5-nat-10.1089_nat.2020.0921.pdf   (S5, 744 pages, 18.4 MB)
```

These are the **only** Holgersen et al. 2021 (doi:10.1089/nat.2020.0921) supplement
files present anywhere under `feng/`. No S1, S2, or S3, and no main-text PDF, exist
in this repo. (Confirmed against the prior Arm-C probe's own data inventory
(`feng/artifacts/phase0_arm_c_probe/README.md`), which already notes Scharner's
binding data lives in "S1-S3, not fully tabulated" for a *different* paper — no
equivalent Holgersen S1-S3 file is staged here either.)

**Both files were parsed directly (pdfplumber, `extract_text` + `extract_tables`,
full document, not a sample) to check for sequence content:**

- **S4** (17 pages): table schema is `ID | Exon interval | Case PSI | Control PSI |
  Edit distance exon body | Mismatches exon body | Edit distance with 200bp padding
  | Edit distance with flanking intron | Longest exact match | RNA-plex delta G |
  RNA-cofold delta G`. The `ID` column holds **75 unique numeric codes** (range
  3–228, non-sequential — these are the paper's opaque per-SBO IDs, confirming the
  "~75 SBO" count from the prior probe), not sequences. Full-text scan of all 17
  pages for the string `sequence` (case-insensitive): **0 hits**.
- **S5** (744 pages): schema is `ID | Gene | Log fold change | Edit distance
  frameshift exon | Edit distance 5' UTR | Edit distance 3' UTR | Edit distance
  full gene`. Same `ID`-code scheme, no sequence column. Full-text scan of all 744
  pages for `sequence` (case-insensitive) or the token `SBO`: **0 hits**.

**Both tables give edit-distance *categories* and ΔG energies per site — never the
oligo sequence itself, never a genomic coordinate for the SBO's binding site (only
for the *off-target exon*), and never a legend/footer page mapping numeric `ID` →
nucleotide sequence.** There is nothing in either file — table body, header, or
any of the 761 combined pages — that would let a scanner reconstruct what
sequence to search the genome for per SBO.

## Why this is the decisive finding, not a soft gap

Task 2 (build the bind-but-no-act negative class) requires, per SBO: enumerate
sequence-matched candidate sites in GENCODE v50 using feng's seed-extend scanner.
The scanner needs a **query sequence** to seed-extend against the transcriptome.
S4/S5 supply the *edit-distance category between the SBO and each already-found
off-target site*, not the SBO sequence that produced it — i.e. they report the
output of a match operation, not either of its inputs. Without the 75 SBO
sequences (or, as a weaker substitute, the target-gene binding-site coordinates
with strand, from which the antisense sequence could be read off GENCODE), there
is no query to seed-extend with, so no candidate universe can be enumerated, so no
bind-but-no-act negative label can be constructed, so no bind→act feature table
exists, so no LOAO test (task 4) is possible. Tasks 2–5 are not attempted — they
are gated by this single missing input, exactly as the task brief anticipated
("If the sequences are not in the files present under feng/, say so explicitly and
stop").

## Where the sequences would have to come from (not fetched — out of scope here)

Per the paper's own structure (per the prior probe's literature read, not
re-verified against a PDF that isn't present): Holgersen et al. 2021's SBO design
table (oligo ID → sequence → chemistry → target gene) is reported in the main text
(likely Table 1 / Table 2) and/or a Supp. Table S1–S3 that is **not staged in this
repo**. Acquiring it (open-access article + its own supplement bundle, distinct
from the S4/S5 files already downloaded here) is a prerequisite for re-attempting
task 2, and is outside this probe's scope (no fetch/download performed; this pass
is read-only over files already present under `feng/`, per instructions).

## VERDICT: NO-GO (for the negative-set construction as specified against files
currently in `feng/`)

**The one decisive number: 0 — hits for "sequence" across all 761 combined pages
of S4+S5, and 0 SBO sequences of any kind recoverable from either file.**

This is not a claim that Arm C's underlying scientific question is unanswerable —
it is a claim that it cannot be answered **from the two files currently in
`feng/`**, and that constructing the negative class as scoped (seed-extend against
GENCODE v50 per-SBO) requires an input (the SBO sequences, or equivalently strand
+ coordinate for each binding site) this repo does not yet hold. Re-scope options,
not pursued here:
1. Acquire the Holgersen main-text design table / actual S1-S3 supplement
   (external fetch, not performed by this read-only pass).
2. Fall back to a coarser negative construction that does not need the SBO
   sequence directly — e.g., if S4/S5 rows carried the *off-target exon's own*
   genomic interval plus a strand-consistent implied binding window, a negative
   set could in principle be built from "candidate off-target exons genome-wide
   that were scored but not listed as hit," using the off-target coordinates
   already in S4/S5 rather than the SBO sequence. **Not evaluated here** — it
   would require re-deriving, per SBO-ID, the full candidate-scoring universe the
   paper searched (not just the hits it reports), which S4/S5 alone do not
   describe either (no "sites considered but rejected" table, no genome-wide
   candidate count per SBO). This is a second, independent blocker on top of the
   sequence gap and is flagged for the record, not solved.

## Honesty note

All counts above (75 unique S4 IDs, page counts, 0 sequence-text hits) are
`paper_supported`/directly measured from the two files in `feng/` by this pass
(pdfplumber `extract_text`/`extract_tables` over the full document, not a
sample). No sealed `terminal/` or Phase-9 data was read. No chemistry model or
disallowed ensemble tool was invoked — this pass never reached feature/model
construction because task 1 did not clear.
