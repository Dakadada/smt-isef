# Final launch audit — current production refresh

Audit time: 2026-11-20T20:17:04.830172+00:00

Campaign T0: **2026-11-15T06:52:49Z**; +48h deadline: **2026-11-17T06:52:49Z**. The campaign window had already elapsed before this refresh. This audit does not retroactively place completion within that window.

Overall: descriptive artifacts available; no scientific efficacy GO or campaign-wide certification. Review approval applies only to the named version and its stated limits. Scheduler states below are recorded evidence, not independently queried live states.

| Family | Status | Slurm IDs (recorded history, may include superseded runs) | Rows / units | Review verdict and file | Evidence class |
|---|---|---|---|---|---|
| J2 | Historical stage completion ledger | 29159943, 29162809, 29162822, 29163480, 29166102, 29166103 | phase1 59,450,148; sentinel 3,762; phase2 total not_available | UNREVIEWED: no standalone J2 review in named inputs; `ops/JOBS.jsonl` | derived |
| J3 | j3_genome_v1 gate | 29166104, 29181966, 29186280, 29190031 | 11002701 | DO-NOT-ADVANCE historical implementation result; current data gate passed, no new certification; `ops/briefs/J3_RESULT.md` | derived / proxy / model_assumption |
| J4 | j4_atlas_v1 deep gate | 29165682, 29167204, 29169335, 29186247, 29186278, 29186279, 29204876, 29207884 | 20395 complete / 20395 configs; 4 selected donor reference rows | DO-NOT-ADVANCE historical review; fixes and deep gate recorded, no fresh independent verdict; `ops/briefs/J4_REVIEW_RESULT.md; ops/briefs/J4_FIX_RESULT.md; ops/briefs/j4_gate.json` | derived / proxy |
| J5 | j5_j6_v4 gate ok:true; embeddings BLOCKED / not_available | 29266905, 29274598, 29274599, 29274600, 29280193, 29280194, 29280195, 29280507, 29280508, 29280509, 29292355, 29292356, 29292357, 29292358, 29292359 | 11002701 rows / 52424 valid windows | DO-NOT-ADVANCE pending profile/gate in historical review; v4 gate now ok, no rewritten review verdict; `ops/briefs/J5_REREVIEW3_RESULT.md` | derived / not_available |
| J6 | j6_genome_v1 gate | 29164482, 29166984, 29168354, 29168376, 29168377, 29168378 | 11002701 | ADVANCE (bounded code review); `ops/briefs/J6_REREVIEW_RESULT.md` | derived |
| J7 | j7_context_v3 merged; marker evidence | 29271719, 29271720, 29271721, 29303724, 29303726, 29303727, 29383227, 29383229, 29383235, 29383243, 29794068, 29794069, 29799681 | 11002701 / 256 shards | REVIEWED / ADVANCE with concurrency limitation; `ops/briefs/J7J9_V3_REVIEW_RESULT.md` | derived / unlabeled / not_available |
| J8 | j8_genome_v2 COMPLETED (user campaign record; not fresh sacct) | 29242485, 29253489, 29261102 | 11002701 / 52424 windows | ADVANCE (implementation review); `ops/briefs/J8J10_REREVIEW_RESULT.md` | derived / proxy |
| J9 | j9_smn2_v3 read_dataset verified both payloads | 29198432, 29223600, 29228637, 29230958, 29303731, 29383229 | 233 features / 132 rank/stability rows | REVIEWED / ADVANCE with bounded FASTA identity; `ops/briefs/J7J9_V3_REVIEW_RESULT.md` | derived / proxy / model_assumption |
| J10 | j10_offtarget_v2 merged; v1 excluded | 29242486, 29254990, 29254991, 29266568, 29266569, 29303729 | 11002701 input rows / 2292 aggregate rows / 256 shards | ADVANCE for v2, bounded verification; `ops/briefs/J7J9J10_REREVIEW_V2_RESULT.md` | derived / proxy / model_assumption |
| J11 | Bounded inline refresh; 6 GiB virtual-memory ceiling; no Slurm job submitted | none (inline) | 4 figures / 4 model cards / report + lock + manifests | UNREVIEWED refresh; focused tests 26 passed, 0 failed, 0 skipped; `report/J11_RESULT.md` | derived / proxy / model_assumption; report synthesis |

## Exact identity and verification scope

Each metadata file below is byte-hashed; its nested claims are producer-recorded unless explicitly stated. Large J5/J7/J8 outputs, J3/J6 shards, genome FASTA, atlas payload collection and external model weights were not reopened or rehashed. No directory-wide payload audit was performed.

- J2: `/home/dashii/smt-isef/conference/model_v2/ops/JOBS.jsonl` — SHA256 `a1883cd5fa5f75b66d0981dad3ecc3549186c7376f437ac2b3e07a3a638fcc3c`.

- J3: `/home/dashii/orcd/scratch/feng/model_v2/j3_genome_v1/gate_report.json` — SHA256 `f6f1f45e8d2c70fc83f9032ad7a9d17af21e29b5c9fcc61847b4998fb16921f9`.

- J4: `/home/dashii/smt-isef/conference/model_v2/ops/briefs/j4_gate.json` — SHA256 `786a35b281151b2e2373765ef1b70404581aeb819c6faa888e3f3da974b952a7`.

- J5: `/home/dashii/orcd/scratch/feng/model_v2/j5_gate_v4.json` — SHA256 `624e6140cee8e457642c2ec69f5b4ba2db411a62bf7be9e6a7822560c59be035`.

- J6: `/home/dashii/orcd/scratch/feng/model_v2/j6_genome_v1/gate_report.json` — SHA256 `31bc599df1c8cb20ff6b657464669deffa24990fb096a701a419584439de59b4`.

- J7: `/home/dashii/orcd/scratch/feng/model_v2/j7_context_v3/shards/j7_context.done` — SHA256 `c703b45bede2ce9baf12c18d72aff00fe41c7a869480af62a1ce4a671f270f91`.

- J8: `/home/dashii/orcd/scratch/feng/model_v2/j8_genome_v2/manifest.json` — SHA256 `8bbe9801ec47c6b9c51ed0691546884c2568e514885daf936d26f6b0db4bbcf7`.

- J9: `/orcd/scratch/orcd/011/dashii/feng/model_v2/j9_smn2_v3/.generations/cdabf4bae9c3b32f42032fe6b282aa32da61391fb2432ddeeb6983808fb04872/manifest.json` — SHA256 `cdabf4bae9c3b32f42032fe6b282aa32da61391fb2432ddeeb6983808fb04872`.

- J10: `/home/dashii/orcd/scratch/feng/model_v2/j10_offtarget_v2/j10_profile.done` — SHA256 `7c3d33bef6872db7b5d1d0b8c7872ed3c3c1e1bb014e6660a40630cdd3867d6b`.

- J11: `/home/dashii/smt-isef/conference/model_v2/j11_report.py` — SHA256 `62db242c7a810559b31cf81aa71cb15b33f0a799e58f65f64c241ba52494e7b6`.

J7 merged parquet recorded SHA256: `660d2b64dfe0ece820ee3e4238ccefe03393d03da59aa6486d126bf3e115cebe` (not rehashed).
J10 aggregate parquet recorded SHA256: `ba1472d398b161ea3e34989989329e58e36027bf0ac5a635e35eb32e957269c9` (not rehashed); JSON SHA256 `3ab6e8787cea3ddda4bd4911a03e1120ae04d7465036516e0551d985a1af76ab` independently matched its completion marker.
J9 `features.parquet`: 233 rows, SHA256 `c3f0367e3b4065d84ad5c4c75f524641307b88ac86907bf90cb138fb3c39ac0a`; read_dataset() verified payload hash, size and row count.
J9 `ranks.parquet`: 132 rows, SHA256 `c84338477fbada6d44ec233c7cdb4f27a961b2515d8f2951ef7ef3507cd508b0`; read_dataset() verified payload hash, size and row count.

J11 output identities are listed in `report/MANIFEST.json`; the predictor identity is in `report/PREDICTION_MANIFEST.json`. This audit is itself locked as an input; it cannot include the resulting snapshot hash without a circular dependency.

## Blocked and evidence boundary

- M6 undefined: no supervised model or efficacy endpoint is available.
- RBP join: BLOCKED pending verified GRCh38 peak bytes, liftover chain hash and mapping QC.
- J5 embeddings: BLOCKED / not_available; no embedding output is authenticated here.
- Supervised training: BLOCKED; no verified negatives. Responders and matched controls do not supply negatives.
- 3UP8 unplaced: needs human confirmation; excluded from the three-ASO ranks.
- ASO chemistry unknown; unmodified-RNA energies do not establish modified-ASO binding or potency.
- SMN coordinates are derived, not GTF-confirmed; paralogs must remain separate.
- July exon-7 gels were previously seen; no prospective claim; excluded from fitting, selection and joins.
- The exon 3 audit is exploratory only; no training, tuning, selection or exon-7 inference.

Descriptive only: sequence-site counts, unmodified-RNA assumptions, accessibility proxies, N-masking reference-model outputs and unlabeled controls. July gels are not prospective validation. No efficacy claims, outcome labels, model fitting, supervised training, or scientific certification are supported.

## Known caveats

- J7 was worker-submitted without prior review; v2 was DO-NOT-ADVANCE. v3 ADVANCE does not erase that history.
- J7 merge 29794069 OOM at 2G; retry 29799681 requested 6G. Tasks 141, 142 and 154 hit transient FileNotFound on J3 shards and were rerun by array 29794068. Merge marker records 11,002,701 rows; scheduler completion and no-writer overlap were not independently verified here.
- J7 review retains a MED concurrency race: shard replacement between verification and consumption. Establish no concurrent shard writer before treating recorded shard hashes as independent proof.
- J9 FASTA identity is bounded: path, size, mtime and FAI SHA, not full FASTA bytes. read_dataset() checks the current immutable generation against the locked manifest; pointer drift fails closed.
- J8 exact-footprint join only: J2 lengths 8/12/18 do not supply 10/14-nt accessibility; missing values stay unavailable.
- J5 v4 gate passed; the historical review remains conditional pending profile/global gate and is not silently relabeled ADVANCE. No J5 embedding output is authenticated.
- J4 historical review/fix history is retained; deep completion is not a new independent scientific review. J2 phase2 total is unavailable.
- v3 review file was absent on initial checks and appeared during the refresh; both explicit v3 ADVANCE verdicts are locked. The 40-minute polling limit was not reached.
- Operational deviations: systemctl returned `Failed to connect to bus: Operation not permitted`; Slurm --test-only connectivity probe did not return and was interrupted (exit 130), with no job submitted. Bare `python` returned `/usr/bin/bash: python: command not found`; subsequent commands used the required interpreter.
- Initial HANDOFF read exposed historical wet-lab outcome text; it was not used. Two initial metadata/listing commands produced excessive output; subsequent reads used bounded summaries. No workbook or sealed artifact was opened.
