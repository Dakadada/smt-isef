# J7/J9/J10 production v2 adversarial re-review — 2026-11-16

**J7: DO-NOT-ADVANCE. J9: ADVANCE with the bounded-verification limits below. J10: ADVANCE with the bounded-verification limits below.** These are review verdicts on the named artifacts and implementation, not validation of the scientific outcomes.

## Findings

- **HIGH J7 — `conference/model_v2/j7_context_features.py:140`: exon-map builder omits its implementation SHA.** `canonical_exons.json` records input hashes and version only; the merged marker's current implementation SHA identifies the shard/merge code, not the code that built the exon map. A different builder can produce a self-consistent map, whose bytes the current shards then faithfully commit. This matters especially for the unreviewed submission. Record the builder implementation SHA and verify it at shard/merge time, then rebuild or independently reproduce and authenticate the map before advancing. This is a provenance gap, not proof that the map is incorrect.
- **MED J9 — `conference/model_v2/j9_smn2_prediction.py:1275`: publication is not generation-atomic.** A rerun overwrites features and ranks before replacing the manifest at line 1311; interruption or concurrent reading can expose a mixture of generations. Write a staged generation and atomically publish a pointer or completion marker after all hashes are verified; consumers must verify hashes. This is a structural risk, not evidence of corruption in the observed build.

## Original eight findings

1. J7 producer partition/reference/universe: addressed at `j7_context_features.py:37,334`: marker/schema agreement, expected partition, inputs/reference/parameters, and expected J3 shard identities are checked at merge.
2. J9 wetlab allowlist: addressed at `j9_smn2_prediction.py:201,227`: exact repo files and restricted roots, lexical and resolved-path exclusions.
3. J9 false M5 evidence: addressed at `j9_smn2_prediction.py:694,705,812`: marker/schema provenance and genuine evidence are required; incomplete evidence is unavailable.
4. J9 endpoint mixing: addressed at `j9_smn2_prediction.py:776,798`: coordinates, strand, endpoint and tissue are checked, duplicates rejected, four channels required, endpoints kept separate; legacy values use donor only.
5. J9 incomplete stability: addressed at `j9_smn2_prediction.py:1134`: all alternatives and full identical ASO membership required, otherwise null/unavailable.
6. J9 reference identity: addressed at `j9_smn2_prediction.py:887,1013`: reference path, size, mtime and FAI hash captured and rechecked; this is bounded identity, not a FASTA content hash.
7. J7 stale shard cache: addressed at `j7_context_features.py:289,357`: implementation SHA participates in shard cache and merge checks. This does not resolve the separate exon-builder finding above.
8. J10 zero-shard gate: addressed at `j10_offtarget_profile.py:286,346`: invalid counts rejected and production CLI requires 256 shards.

## Production evidence and limits

- J7 `j7_context_v2/shards/j7_context.done`: 11,002,701 rows and 256 shard hashes; common implementation SHA `cb78f86c7dad168fc7eaa861839e347d35506e1eafbfac1b9e323e2d9727e9b8` matches current source. `canonical_exons.json` SHA `edfe503c7ac3ae8f7add8fe0a0e28b66a4cfc38b2cfb477d1bade8054c77ddcc` matches marker; records 20,033 transcripts and 209,760 exons; both recorded J2 manifest hashes match current repo files. Marker records RBP BLOCKED. Merge code checks shard integrity, uniqueness, and row total; this review did not rerun the heavy gate.
- J9 `j9_smn2_v2/manifest.json` and README: current J9 SHA `ebfa1d6b46ca33688a4651a668fba3d56fdaa296d52e9f486d33a7176b21074d`; curated `j4_atlas_v1/SMN2_exon7_issn1_fine`, committed atlas with 1,064 rows, run ID and shard SHA, recorded reference identity, 233 feature and 132 rank rows. The job log reports the same output counts/hashes. Current code reads and verifies atlas before feature assembly (`j9_smn2_prediction.py:947`); M5 is not read from the folding cache. Four seconds alone does not establish staleness; reuse of correctly keyed J3 folding cache is possible. Manifest is self-reported provenance; manifest-only restriction prevents independent verification of atlas marker/data or output bytes, and cannot prove historical execution against tampering or source changes during execution.
- J10 `j10_offtarget_v2/j10_profile.done`: current implementation SHA `e6163a53788c4288736e8390bbcf382558b38188bb9b4591858943d7cc1a86e8`; 256 distinct source indices, identical recorded configs, source rows sum to expected 11,002,701; output marker records 2,292 aggregate rows and JSON/parquet hashes. Current merge checks each partial's hashes, row counts and configuration (`j10_offtarget_profile.py:222,303`). Marker-only restriction precludes independently rehashing payloads or rerunning the merge gate.
- No production data tables, J3 inputs, atlas files, FASTA or GTF were read. J7 merged and J10 profile hashes and completeness are reported by their markers, not independently recomputed. J7 source-universe provenance is internally consistent as reported by the merge, not independently verified from J3. Canonical map and local J2 manifests were hashed as permitted.

## Tests and remaining work

Focused per-file tests with bytecode and pytest cache disabled: J7 29 passed, J9 58 passed, J10 19 passed, regression 43 passed, J10 robustness 37 passed: **186 passed, 0 skipped, 0 failed**. No full suite was run (explicit task constraint overrides the repository gate). No code fixed; J7 provenance blocker and J9 atomic-publication risk remain. No job launched, commit or push.
