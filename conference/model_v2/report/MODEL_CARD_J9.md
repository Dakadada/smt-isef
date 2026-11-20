# J9 model card

Inputs: `/orcd/scratch/orcd/011/dashii/feng/model_v2/j9_smn2_v3/.generations/cdabf4bae9c3b32f42032fe6b282aa32da61391fb2432ddeeb6983808fb04872/manifest.json` (SHA256 `cdabf4bae9c3b32f42032fe6b282aa32da61391fb2432ddeeb6983808fb04872`); `/orcd/scratch/orcd/011/dashii/feng/model_v2/j9_smn2_v3/.generations/cdabf4bae9c3b32f42032fe6b282aa32da61391fb2432ddeeb6983808fb04872/ranks.parquet` (SHA256 `c84338477fbada6d44ec233c7cdb4f27a961b2515d8f2951ef7ef3507cd508b0`); `/home/dashii/smt-isef/conference/model_v2/j9_smn2_prediction.py` (SHA256 `c4ae2025fcbafefd0fabb58745804c2cd4c0573ff4ad3d9b7749a94bf163ecaf`).

Evidence class: derived / proxy / model_assumption.

Output: Descriptive M0–M5 per-feature ordering for three placed ASOs and recorded alternatives.

Labels: unlabeled; no outcome labels or verified negatives are read, joined, fitted or selected.

Assumptions: GRCh38 / GENCODE v50; sequence-derived SMN coordinates remain unconfirmed against GTF. Chemistry unknown; any duplex term assumes unmodified RNA. Canonical context and model tissue channels do not establish cell-specific ASO response.

Limitations: No outcome labels, fitted predictor, efficacy test, or placed 3UP8.

Future independent wet-lab comparison must be made against this snapshot. Already-seen July gels are NOT prospective validation.

Review: REVIEWED / ADVANCE; see locked J7J9_V3_REVIEW_RESULT.md. FASTA identity is bounded, not a content hash.
