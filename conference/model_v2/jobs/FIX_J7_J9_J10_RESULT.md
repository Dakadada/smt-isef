- ADVANCE code fixes: all 8 findings addressed; production outputs DO-NOT-ADVANCE until rerun/gated. No submissions, commits, wetlab/outcome reads, or feng edits.
- `j7_context_features.py:37,334`: require producer partition + schema/marker agreement; bind reference/universe/parameters and reverify against expected J3 files/exon map at merge.
- `j7_context_features.py:289`: implementation SHA256 participates in shard cache identity and merged common commitment; stale implementations cannot merge.
- `j9_smn2_prediction.py:227`: exact implementation files + dedicated atlas/reference/output roots; reject excluded paths before resolve/stat/open, including symlink targets.
- `j9_smn2_prediction.py:705`: verify J4 marker integrity/run ID/plan, schema provenance and artifact hashes/counts; M5 requires genuine model_evidence=True.
- `j9_smn2_prediction.py:765`: validate placement/endpoint identity; reject duplicate endpoint/tissue keys; require tissues 0,2,4,6; separate endpoints; legacy M5 columns select donor only, never fall back/pool.
- `j9_smn2_prediction.py:1095`: incomplete alternatives or ASO membership yield stable=None/status=not_available and null spread; M2/M3 require all 3 settings.
- `j9_smn2_prediction.py:887,1013`: record/recheck reference path, size, mtime_ns and FAI SHA256 before/after building; synthetic references explicitly identified.
- `j10_offtarget_profile.py:286,337`: reject nonpositive counts/negative expected rows before output I/O; production CLI requires 256 shards; positive synthetic API counts remain supported.
- Tests: initial regressions 34 failed/2 passed before fixes; final focused 150 passed/0 skipped/0 failed. Updated buggy unconditional stability assertion; replaced generic J7 markers and production-dependent J9 atlas fixture with committed synthetic inputs.
- Final full gate: 1239 passed, 1 skipped, 0 failed (19.65s); all 6 sbatch scripts pass bash -n. Initial full-gate legacy-M5 vocabulary failure fixed by explicit donor-only columns.
- Gate command: `PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 /home/dashii/orcd/pool/smt-isef-env/bin/python -m pytest conference/model_v2/tests -q`.
- Invalidated: J7 j7_context_v1 commitments; J9 j9_smn2_v1 M5/stability/reference provenance (rebuild whole dataset); J10 v1 implementation 09fa7cc9 != reviewed 5430cbf3, so historical evidence cannot certify current code; regenerate all v2 partials under current hash.
- Unfixed: none of the 8 findings. Limits: no production reruns/artifact certification; FASTA identity uses bounded stat/index checks, not full-genome hashing; atlas coverage may leave M5 unavailable; RBP remains BLOCKED.
- Exact reruns below are NOT EXECUTED; fresh J7/J9 roots and J10 v2; dependencies gate merges (run from a shell with Slurm access).
```bash
cd /home/dashii/smt-isef
export J7_ROOT="$HOME/orcd/scratch/feng/model_v2/j7_context_v2" J3_ROOT="$HOME/orcd/scratch/feng/model_v2/j3_genome_v1" J7_EXPECTED_ROWS=11002701
j7b=$(sbatch --parsable --partition=mit_preemptable --export=ALL conference/model_v2/jobs/j7_build_exons.sbatch)
j7a=$(sbatch --parsable --partition=mit_preemptable --export=ALL --dependency=afterok:"${j7b%%;*}" conference/model_v2/jobs/j7_array.sbatch)
sbatch --partition=mit_preemptable --export=ALL --dependency=afterok:"${j7a%%;*}" conference/model_v2/jobs/j7_merge.sbatch
j10a=$(sbatch --parsable --partition=mit_preemptable conference/model_v2/jobs/j10_array.sbatch)
sbatch --partition=mit_preemptable --dependency=afterok:"${j10a%%;*}" conference/model_v2/jobs/j10_merge.sbatch
J9_OUT_DIR="$HOME/orcd/scratch/feng/model_v2/j9_smn2_v2" J9_ATLAS_DIR="$HOME/orcd/scratch/feng/model_v2/j4_sentinel_v2" sbatch --partition=mit_preemptable --export=ALL conference/model_v2/jobs/j9_smn2.sbatch
```
