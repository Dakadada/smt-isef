J7 — implementation prepared; production UNRUN; no submission or commit.
- j7_context_features.py:56: one batch GTF scan builds compact canonical exon JSON; full SHA256 + both J2 manifest hashes; exact padded-exon/complete-set checks and GRCh38/v50 header fence.
- j7_context_features.py:146: derived geometry uses transcript-5′ footprint anchor; donor=last exonic, acceptor=first exonic; downstream positive; true boundary zero retained, missing never zero.
- Geometry: whole-footprint exon/intron/junction-proximal (crossing or ≤10 bp); containing exon/intron lengths; absolute bins 0/1–10/11–50/51–200/201–1000/>1000; ties choose upstream.
- Canonical J2 protein-coding chr1–22,X only (209,762 manifest lines); alternate isoforms not evaluated; unselected genes/regions return not_available. Gene-only model_v2/gtf_genes_v50.gtf is absent; full canonical GTF supplies gene IDs directly.
- j7_context_features.py:186: deterministic background controls match gene/strand/segment/length/region/bins; exclude source overlap; 64 draws then null; evidence=unlabeled, never verified negatives.
- j7_context_features.py:215: RBP BLOCKED: metadata audit found 669 eCLIP sources, all hg19; 114 recorded SHA256; HepG2/K562. No source peak bytes reverified or KD/outcome tables opened.
- Missing: verified GRCh38 peak bytes, source accession/cell/assay/hash mapping, lift chain SHA256 + mapping/drop QC. No join or fabricated zeros; per-row cell/assay/accession/hash/count=null, rbp_missing=true. Manifest audit retained separately.
- Example indexed source only: ENCFF106CGP/FUS/HepG2, hg19, recorded SHA256 51428efbf8e811c9381bd36562599bd89f21d69ad1405a4e728c8fd987ffa75e; not a byte-verification claim.
- No RBP-knockdown == steric-block equivalence claim. No old pre/src imports; feng unchanged.
- j7_context_features.py:255: 512-row streaming sidecar keyed by site_id; one corresponding J3 shard each; atomic Parquet then SHA256/count/provenance .done; locked restart verifies inputs/commit.
- j7_context_features.py:304: merge gate checks exact shard set/schema/count/hash/common provenance, null/bin consistency and disk-backed site-ID uniqueness; atomic merged Parquet/.done; concurrent merges locked.
- tests/test_j7_context_features.py:1: RED→GREEN; final requested gate: 29 passed, 0 skipped, 0 failed (0.55 s). Only requested J7 suite run; broader feng suite not run for this conference-only task.
- Test command: `PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 /home/dashii/orcd/pool/smt-isef-env/bin/python -m pytest conference/model_v2/tests/test_j7_context_features.py -q`
- All three jobs: mit_normal, --requeue, one CPU, 2 GiB; array 0–255%16. bash -n checks passed.
- Submit from repo root ONLY after J3 gate accepts all 256 committed shards. Commands below assume accepted J3 count=11002701 (existing J3 gate target); change only to actual accepted count if different.
`export J7_EXPECTED_ROWS=11002701`
`j7_build=$(sbatch --parsable conference/model_v2/jobs/j7_build_exons.sbatch)`
`j7_array=$(sbatch --parsable --dependency=afterok:${j7_build} conference/model_v2/jobs/j7_array.sbatch)`
`sbatch --dependency=afterok:${j7_array} conference/model_v2/jobs/j7_merge.sbatch`
- Default inputs: ~/orcd/scratch/feng/model_v2/j3_genome_v1; full GTF ~/orcd/scratch/feng/gencode/v50/gencode.v50.primary_assembly.annotation.gtf. No heavy reads performed in this session.
- Default outputs: ~/orcd/scratch/feng/model_v2/j7_context_v1/{canonical_exons.json,shards/}; merged shards/j7_context.parquet. J7_ROOT/J3_ROOT environment overrides supported; use a fresh J7_ROOT for changed inputs.
- Estimate only at 11.0M rows: compact JSON 10–40 MB; context shards 1–4 GB; merged 1–4 GB; reserve 15 GB including temporary merge/SQLite. Runtime/peak RSS/production coverage unmeasured; requested 30m build, 2h array tasks/merge.
- Unfixed limitations: RBP blocked; no production data generated or coverage/size measured; nearest-boundary geometry is canonical-only and anchor-based, not an isoform effect or potency prediction.
