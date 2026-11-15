# F1 RESULT — trans_intervention.py
## Files
- conference/model_v2/trans_intervention.py (SCHEMA_VERSION 1.0.0; frozen dataclasses; FaiFastaReader; place / place_at / build_intervention / validate / join_placements; JSON + Parquet; CLI validate|place|schema)
- conference/model_v2/schema/trans_intervention.schema.json (emitted by `schema`; test asserts it equals the generator)
- conference/model_v2/tests/test_trans_intervention.py
## Tests
- `/home/dashii/orcd/pool/smt-isef-env/bin/python -m pytest conference/model_v2/tests/test_trans_intervention.py -q` -> 38 passed, 0 skipped, 0 failed (py3.11.16, pytest 9.1.1, pyarrow 23.0.1).
- Fasta absent (TRANS_INTERVENTION_FASTA=/nonexistent): 30 passed, 8 skipped with reason.
- Bare `python` is not on PATH on the login node; /usr/bin/python3 is 3.6 without pytest. Use the env python above.
## Verified against GRCh38 fasta (tiny .fai slices only)
- chr5:70076584-70076601 (+) == CCAGCATTATGAAAGTGA == revcomp(F18). F14 -> 70076584-70076597, F10 -> 70076584-70076593 (unique in anchor window; target 5' end fixed).
- Donor: 70076575-76 == GT; last exonic base 70076574, first intronic 70076575. Both conventions stored: F18 = +10..+27 (origin last exonic base) / +9..+26 (origin first intronic base).
- Acceptor: 70076519-20 == AG; exon 7 = 70076521-70076574 (54 nt).
- SMN1_exon7 (alternate locus): F18 target at chr5:70952004-70952021 (+), exon 70951941-70951994; exon pos 6 = C (SMN2 = T). Distinct locus_id and site_id; join_placements raises ParalogJoinError.
- 3UP8: status needs_human_confirmation, placement null, candidates enumerated (SMN2 window: 70076581-70076588 only; SMN1 window: 70952001-70952008 only). CLI exit code 3.
## Open questions
1. 3UP8 placement needs human confirmation (single exact candidate per locus, straddles 3 nt upstream of the anchor). No API to record a confirmation yet.
2. Exon-7 start (70076521) and all SMN1 coordinates are DERIVED from sequence (54-nt exon, AG/GT, +875420 offset), not read from the GENCODE v50 GTF; coords_evidence="derived". Confirm against the GTF.
3. ASO chemistry is not in conference/aso_pipeline/sequences.py -> family "unknown", evidence not_available. ASO sequence_evidence defaults to "unlabeled"; say if it should be "observed".
4. Locus search windows (exon -100 / +100 nt) are my choice; ASOs outside them return no_match.
5. conference/model_v2/ has no __init__.py (not owned); `python -m conference.model_v2.trans_intervention` works as a namespace subpackage from repo root; tests insert repo root on sys.path.
## Invariants not satisfied / limits
- No jsonschema library in the env: the schema file is emitted and key/enum-checked, but instances are validated by the module's own strict decoder, not by a JSON-Schema validator.
- Exact matching only (no mismatch-tolerant placement); one exon per locus (distances to that exon's acceptor/donor only).
- Parquet determinism is tested at table level (round-trip equality), not byte-identical files.
- Not committed, not pushed. Nothing outside owned paths touched (plus this result file).
