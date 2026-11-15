# ASO computational deliverables design

Status: approved by the human partner on 2026-09-26. Outputs are sequence-derived
model predictions and must carry the literal state `[candidate]`.

## Scope

Build three independent but interoperable pipelines for the F18, F14, F10, and
3UP8 ASOs: (1) two-model on-target hybridization free energy with a predeclared
1x-TBE salt sensitivity grid, (3) GRCh38/GENCODE-v50 distinct-gene off-target
counts from Bowtie 1 and BLASTn, and (6) 100 fixed-seed dinucleotide-preserving
scrambles through both pipelines. Item 2 is only a typed, refusing interface for
future replicate-level PSI data. No Arm E sealed/frozen path is read or changed.

## Sequence and thermodynamic model

ASOs are normalized from DNA letters to RNA. Each intended target is derived by
reverse-complement lookup in the supplied ISS-N1 context; uniqueness is required.
This places F18, F14, and F10 at 1-based context positions 10-27, 10-23, and
10-19. It places 3UP8 at positions 7-14 (`CUGCCAGC`), exactly three bases
upstream of F18's start. The 3UP8 placement is an explicit candidate assumption
requiring human confirmation.

Model A calls the existing `feng.titration.evaluate.intarna_on_target` wrapper,
using IntaRNA 3.4.1 defaults at 37 C, no explicit salt, and Turner 2004 RNA
parameters. Model B calls ViennaRNA 2.7.2 `RNAduplex` at 37 C with its Turner
2004 nearest-neighbor parameters. Direct verification showed that this
RNAduplex MFE does not change when only `--salt` changes, so the standard
monovalent entropy adjustment is applied explicitly and stored separately from
the raw 1.021 M reference energy. The operational 1x-TBE effective monovalent
concentration is 0.089 M, set equal to the documented
89 mM Tris/borate formulation; because this is an approximation rather than a
chemical Na+ measurement, the predeclared grid is 0.025, 0.050, 0.089, 0.100,
and 0.150 M. Raw model energies and derived ranks remain in separate columns.
RNA:RNA is the primary chemistry; the command layer accepts an alternate
Vienna parameter file (including `DNA`) for later DNA:RNA-oriented reruns, but
the current outputs make no DNA:RNA claim.

Scrambled ASOs are scored against the original biological target window, not a
new perfectly complementary target. This tests loss of sequence-specific
affinity while retaining dinucleotide stacking composition. IntaRNA and
RNAduplex may select partial/gapped interactions; no failed score becomes zero.

## Off-target model

The staged GENCODE-v50 GRCh38 primary-assembly FASTA and matching GTF are the
only reference inputs. Bowtie 1 uses exhaustive end-to-end `-f -a -v 3`
without best-stratum filtering and suppresses unused output fields. BLAST uses `-task blastn-short -evalue 1000
-dust no -strand both`, the documented short-query preset threshold. Alignments
are streamed through a coordinate-overlap mapper into disk-backed GNU
`sort -u` and a constant-memory reducer; SAM/BLAST hit dumps are never
persisted. The primary metric is distinct GTF gene-span overlaps whose
gene strand is compatible with an ASO binding the transcript. An any-strand
audit count is retained separately.

Before alignment, each query receives an analytically estimated random-genome
hit count under an equal-base iid model. A predeclared 10,000,000 expected-hit
threshold marks it `saturated / non-specific by construction - count
uninformative`; such queries are never sent to an aligner and receive no numeric
gene count. F18 and F10 are mandatory exhaustive headline queries despite the
preflight estimate, per human instruction. F14 is run when below threshold;
3UP8 is expected to be saturated. Unique sequences are aligned once and mapped
back to all 100 null draws, preserving the null multiplicities.

F18 and F14 are separate four-thread Bowtie partitions. Every unique F10
sequence is a separate serial one-thread Bowtie partition, preventing multiple
explosive 10-mer searches from co-residing. External sort has a 2 GiB memory
ceiling and the Slurm allocation has a 220 GiB safety ceiling. Indexes and database live under
`/orcd/scratch/orcd/011/dashii/conference_smn2_aso/`. A CPU
`mit_preemptable` Slurm job builds each once, uses stage checkpoint sentinels,
forwards SIGTERM to its payload, and requeues. Final compact TSV/JSON artifacts
are copied atomically into `conference/artifacts`; large indexes stay in
scratch.

## Artifact and failure contracts

Every row records `[candidate]`, model/aligner parameters, sequence identity,
and whether it is real or a null draw. Missing tools, failed commands,
unparseable records, saturation, and unavailable results remain distinct states;
none is replaced with `0` or an invented value. README files cite assumptions,
software versions, model parameters, and references. Item 2 validates only the
future schema (`aso_id`, `length`, `replicate`, `PSI_percent`) and raises a
clear unavailable-data exception rather than fitting anything.

## Tests

Tests are per production file. They cover reverse-complement target placement,
literal thermodynamic parsing, exact dinucleotide preservation and determinism,
saturation math and mandatory-query overrides, strand-compatible GTF overlap,
stream-parser behavior, and PSI-stub refusal. External-tool tests use the real
installed IntaRNA/RNAduplex binaries only where fast; no full Feng suite is run.
