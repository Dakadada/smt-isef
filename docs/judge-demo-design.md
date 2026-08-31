# #17 judge-facing demo: real SMN2 evidence to prospective design

**Status:** started; rendering is blocked by the #16 SMN2 real-context acceptance
gate below. Locking/publication is additionally blocked by the final adversarial review.
DMD and all other loci may remain visibly pending. This document intentionally records both the disputed #16 exit gate and the
#17 consumer contract so ownership cannot drift between them.

## Product decision

#16 is locked only after one end-to-end production path for SMN2 uses a real, versioned
context whose exact bytes and genomic-to-transcript mapping are hash-committed in the
project registry. That gate is now implemented; DMD remains pending.

The SMN2 sequence is a #16 responsibility, not a #17 demo fixture. `design_panel`
owns the coordinate registry, sequence digest, transcript identity, mapping digest,
therapeutic polarity, and call into `design_engine`. #17 must consume the deposited
panel artifact. It must not construct or inject a private context that the checked-in
panel cannot reproduce.

#16 may lock when all of the following are true:

1. The exact SMN2 context is bound to NCBI RefSeq NM_017411.4. NM_000344.4 is SMN1 and
   is rejected as the wrong gene identity. The intronic ISS-N1 interval is reconstructed
   from the published Anti-N1 sequence and committed source/mapping bytes.
2. The context boundaries cover the declared design interval and every regulatory
   interval used by the engine, with an explicit GRCh38 genomic-to-transcript map.
3. `config/event_junctions.GRCh38.tsv` commits the exact context SHA-256; provenance
   commits source, retrieval date, transcript version, strand handling, and mapping
   SHA-256.
4. The production loader builds `RealTranscriptContext` without calling
   `offline_target` or using the ACGT fixture.
5. An end-to-end test runs the checked-in panel and proves: at least one SMN2 design
   report is emitted, its manifest binds the committed sequence and mapping hashes,
   its machine-readable prediction status remains prospective, and a one-base context
   mutation fails closed. #17 owns validation and display of the exact human-readable
   banner at every rendering boundary.
6. A frozen SMN2 panel artifact is deposited successfully. DMD remains pending until
   its equivalent context is curated; that pending state is shown, not hidden.

## Judge experience

The demo is one local command:

```bash
make judge-demo RUN_DIR=/path/to/completed-stageA-run PANEL=/path/to/frozen-panel.json \
  BUILD_TIMESTAMP=2026-09-01T12:00:00Z
```

For the current MIT deployment, `RUN_DIR` is
`/orcd/scratch/orcd/011/dashii/smt-stageA-real`. `PANEL` is the exact no-overwrite
artifact deposited by #16; it is never discovered by “latest file” heuristics.

It creates a timestamped, no-overwrite directory under `artifacts/judge-demo/` with:

- `index.html`: the primary judge-facing artifact, self-contained and viewable without
  a server;
- `evidence.json`: the selected real-result rows plus hashes of their source tables;
- `prospective-design-panel.json`: the deposited #16 panel report, copied byte-for-byte;
- `provenance.json`: input hashes, code identities, command, supplied build timestamp,
  and claim-boundary metadata;
- `README.txt`: the exact command and a text fallback.

The command refuses to render if `RUN_DIR/STAGE_A_COMPLETE` is absent, the run config
schema is not `1.0`, a required input or selected row is missing, a selected row has
`coverage_pass != true` or fewer than 20 informative reads, an expected control is
absent, the panel schema is not `prospective-design-panel-report/2`, the SMN2 panel
contains zero designs, or any displayed design lacks the exact banner:

> PROSPECTIVE PREDICTION — NOT EXPERIMENTALLY VALIDATED

Every observation must also satisfy its declared `expected_direction`: `increase` is
positive dPSI, `decrease` is negative dPSI, and `no_change` requires `|dPSI| <= 2`.
Any mismatch fails rendering.

## Executable input and evidence contract

`config/judge_demo_evidence.json` is the declarative selector. It exhaustively lists
every eligible treated arm and matched control in the three accepted SMN2 families;
the demo may not select a favorable subset. It also lists excluded arms with a
machine-readable reason. Each observation names:

- `source`: `work/harmonized/tidy_splicing.tsv`;
- the composite key `family`, `sample`, `intervention`, `event`, `condition`, and
  `pair_set_id`;
- the expected direction (`increase`, `decrease`, or `no_change`), required comparator
  kind, unique `comparator_observation_id`, and minimum informative reads;
- whether the frozen numeric `magnitude` must match exactly, using its canonical TSV
  decimal string rather than binary floating-point equality;
- study/project and run accession, biological/technical replicate classification,
  group size, aggregation rule (`none`: points remain run/sample observations), and
  the primary publication supporting the expected direction.

The demo presents a frozen completed run, not a tolerant recomputation. Once Stage A is
complete and reviewed, the selected numeric values and source-file SHA-256 are committed
to that evidence manifest. A rerun with different bytes is a new evidence version and
must be reviewed and frozen separately. “Near zero” means absolute dPSI no greater than
2 percentage points; “coverage” means `coverage_pass=true` and `informative_reads >= 20`.
“Byte-verified” means the corresponding fetch checkpoint records the expected archive
MD5 and hashes the final FASTQ output. “Qualified endpoint evidence” means an eligible
RNA-seq row passing those integrity and coverage gates.

The builder reads every input once as bytes, hashes those bytes, and parses the same
buffer. It never hashes a path and later rereads changed contents. Required inputs are:

- `RUN_DIR/STAGE_A_COMPLETE`;
- `RUN_DIR/stageA.targeted.json` (`schema_version=1.0`) and
  `RUN_DIR/stageA.targeted.json.sha256`;
- `RUN_DIR/all_runs.frozen.tsv` and `RUN_DIR/all_runs.frozen.tsv.sha256`;
- `RUN_DIR/work/harmonized/tidy_splicing.tsv`;
- `RUN_DIR/work/counts/all_counts.tsv` and exactly
  `RUN_DIR/work/checkpoints/fetch/<run_accession>.json` for every accession enumerated
  by the evidence manifest; no checkpoint-directory scan is allowed;
- the explicit `PANEL` artifact;
- `config/judge_demo_evidence.json`.

The normative schemas are checked in as `config/schemas/judge-demo-selector-1.json`
for the selector, `config/schemas/judge-demo-evidence-1.json` for generated evidence,
`config/schemas/judge-demo-provenance-1.json` for provenance, and
`config/schemas/judge-demo-artifact-manifest-1.json` for publication. `evidence.json` uses schema
`judge-demo-evidence/1` and contains the selected rows,
stable observation IDs, units, selector keys, comparator IDs, source-byte digest, frozen
run-manifest digest, and archive-integrity references. `provenance.json` uses schema
`judge-demo-provenance/1` and contains the digest of every input, the digest of the
unchanged panel bytes, the evidence digest, renderer source identity (Git commit plus
an explicit sorted map of renderer/config relative paths to SHA-256; dirty bytes are
therefore represented directly), build command,
supplied UTC build timestamp, claim text, and a machine-readable list of prohibited
claims. Each plotted point links by observation ID to its evidence record and accession
integrity references.

The panel validator requires the exact report schema and prospective status, plus typed
candidate sequence, effector class, chemistry, transcript coordinates, uncertainty,
off-target fields, sequence/mapping/report SHA-256 values, and internal equality between
context hashes and report-manifest hashes. Distinct provenance fields are named
`run_manifest_sha256` and `panel_sha256`. `artifact_payload_sha256` appears only in
`artifact-manifest.json` and the frozen-artifact registry, avoiding a self-digest cycle.

Prohibited-claim validation uses `config/judge_demo_claim_policy.json`, whose nonempty
patterns are schema-validated literal phrases matched as normalized substrings. It scans only the headline, observation captions, candidate
cards, and demonstrated/prospective summary fields after Unicode normalization and
case-folding. It does not scan provenance or the explicit “not demonstrated” boundary
panel, where forbidden phrases must remain visible as disclaimers.

`index.html` embeds the validated JSON and uses no network resources. Every plot has a
text table equivalent, visible axes and percentage-point units, keyboard-focusable
points, and accessible expansion controls. A `judge-demo-preflight` target performs all
read-only dependency checks and reports every missing item before rendering.

Rendering is deterministic for identical input bytes plus `BUILD_TIMESTAMP` (strict UTC
RFC 3339). Publication writes and fsyncs all payload files in a sibling directory under
`artifacts/judge-demo/`, builds canonical JSON as UTF-8 with sorted keys and no
insignificant whitespace containing payload entries sorted by UTF-8 relative-path bytes
(`path`, byte `size`, `sha256`), and hashes those canonical JSON bytes. File modes are
excluded. It writes/fsyncs `artifact-manifest.json`, then
atomically renames the directory to `<BUILD_TIMESTAMP>-<payload-digest-prefix>`. A
parent-directory fsync follows the rename. A collision is accepted only after validating
the existing manifest and every payload hash: an identical build returns the existing
path and removes the temporary sibling; any mismatch fails.
A frozen artifact is one whose full digest is added by reviewed commit to
`artifacts/judge-demo/registry.tsv` (`timestamp`, `artifact_payload_sha256`, `path`,
`review_status`). Signing is optional and not treated as immutability.

Preflight diagnostics are staged: path/integrity checks run first, schema checks second,
evidence selection third, and cross-file provenance checks last. Each check reports
`PASS`, `FAIL`, or `BLOCKED_BY` so one broken prerequisite does not hide unrelated
missing inputs.

## What the judge sees

### 1. Real observations

A compact direction plot shows only byte-verified SMN2 endpoint results:

- three chemistry-specific targeting ASO arms: +25.5, +12.1, and +20.0 inclusion dPSI;
- F18MOE replication: +14.3 inclusion dPSI;
- two dCas13d D[333,357] replicates: -32.1 and -21.8 inclusion dPSI;
- chemistry-matched scramble and non-targeting-guide controls centered near zero.

Each point expands to accession, sample, intervention, comparator, informative-read
count, PSI, dPSI, approximation label, and source-table hash. These are observations,
not training performance and not evidence that a new candidate works.

### 2. One coherent locus-level evidence chain

The visual connects observations from two independent SMN2 BioProjects and their
constituent biological replicate libraries to intervention
placement and direction without claiming a learned
general rule: masking ISS-N1 with steric ASOs increases exon-7 inclusion; blocking the
downstream ISE-D1 region with dCas13d decreases inclusion. Controls remain near zero.
This establishes directional assay feasibility at one locus across distinct steric
effector implementations.

### 3. Prospective handoff

The demo then displays the real-context-bound SMN2 Pareto set produced by #16. Every
candidate card begins with the prospective banner and shows sequence, effector class,
chemistry, transcript coordinates, objective direction, uncertainty/abstention state,
off-target burden, and the sequence/mapping/report hashes needed for a future blind
assay. No candidate card displays an experimentally measured magnitude.

### 4. The boundary is part of the result

A final panel says plainly:

- demonstrated: directional recovery of known SMN2 splice responses from public raw
  reads, with clean matched controls;
- prospective: newly generated SMN2 candidate rankings;
- not demonstrated: a powered general benchmark, cross-locus generalization, clinical
  benefit, a validated molecule, or a new biological mechanism;
- pending: DMD and every other locus lacking a committed real context or qualified
  endpoint evidence.

## Headline claim

> **A compute-native SMN2 pipeline recovered the expected splice direction across three distinct steric implementations—ASO1, F18MOE, and dCas13d D[333,357]; genuine scramble-MOE and non-targeting-guide null arms stayed near zero; and the sequence-committed SMN2 target yielded frozen prospective designs for follow-up testing.**

This is the strongest defensible sentence because every verb maps to a different
artifact:

- “public raw reads” maps to checksum-verified FASTQs and the frozen manifests;
- “recovered the expected SMN2 splice direction” claims directional concordance at one
  named locus, not predictive accuracy or mechanism discovery;
- the three implementations are enumerated rather than inflated into plural ASO species;
- the null clause applies only to scramble-MOE and non-targeting-guide arms, never to
  untreated ASO1 replicate scatter;
- “sequence-bound prospective designs” is permitted only after the #16 context hashes
  are committed and the real panel emits designs;
- “froze ... before follow-up testing” describes the timestamped commitment. It does
  not imply that follow-up or experimental blinding has occurred.

Before `STAGE_A_COMPLETE`, the rendered sentence is prefixed `INCOMPLETE STAGE-A
SNAPSHOT:` and followed by `Unresolved arms, including F20PMO, are excluded.` The final
unprefixed headline cannot ship until Stage A completes.

The hostile-judge answer is direct: this is not presented as a general benchmark win.
It is a complete feasibility chain at SMN2: raw intervention reads, directionally
coherent positive and negative perturbations, matched null controls, and an immutable
prospective output that creates a falsifiable next experiment.

## Implementation slices

1. Finish the #16 SMN2 real-context gate and deposit one frozen panel artifact.
2. Add an offline, reproducibly installable demo builder that selects evidence through
   `config/judge_demo_evidence.json` rather than hard-coded row numbers.
3. Add fail-closed claim validation and per-design banner validation.
4. Render the static HTML/JSON/text artifact and add `make judge-demo`.
5. Test missing completion, changed hashes, missing controls, missing SMN2 designs,
   forbidden claim language, banner stripping, deterministic rendering, and
   no-overwrite deposits.
6. Run the demo against the completed real run, freeze the artifact, and conduct the
   required adversarial review before locking #17.
