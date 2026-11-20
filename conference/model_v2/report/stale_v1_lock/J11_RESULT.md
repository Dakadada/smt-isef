J11 complete: report snapshot only; all scientific evidence boundaries retained.
conference/model_v2/j11_report.py:100 — bounded input capture; missing/hash/runtime/HEAD drift raises; explicit immutable lock.
conference/model_v2/j11_report.py:238 — four deterministic descriptive figures; threshold repetitions deduplicated without inventing counts.
conference/model_v2/j11_report.py:323 — report, J3/J4/J9/J10 cards, blocked list, reproduction instructions and Slurm template.
conference/model_v2/report/PREDICTION_MANIFEST.json:1 — source/HEAD/recorded J4 model+weight hashes; required future-comparison boundary.
Inputs: 62 files, 1733860 bytes; manifest: 75 verified entries (62 inputs + 13 outputs; manifest excludes itself).
Outputs: 4 figures, 4 model cards, REPORT, INPUT_LOCK, PREDICTION_MANIFEST, MANIFEST, REPRODUCE, generate.sbatch (14 files) + this result.
Tests: 17 passed, 0 failed, 0 skipped; synthetic fixtures only; final run 1.90s.
Test command: `PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 /home/dashii/orcd/pool/smt-isef-env/bin/python -m pytest conference/model_v2/tests/test_j11_report.py -q`
Capture command: `PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 /home/dashii/orcd/pool/smt-isef-env/bin/python conference/model_v2/j11_report.py capture`
Generate command: `PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 /home/dashii/orcd/pool/smt-isef-env/bin/python conference/model_v2/j11_report.py generate`
Verify command: `PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 /home/dashii/orcd/pool/smt-isef-env/bin/python conference/model_v2/j11_report.py verify`
Slurm syntax command: `bash -n conference/model_v2/report/generate.sbatch` (exit 0); supplied, not submitted; all actual inputs fit tiny-read caps.
Snapshot SHA256: def98dcce77c6bfdf51db56f1905bb16fec9cc7c3f655aa3440f45c45086a26c; real-input regeneration preserved every manifest hash.
Nothing requested left ungenerated. Limitations: requested J10 v1 stale (reviewed v2 outside inputs); J2 phase-2 total unavailable; J9 feature count manifest-only; J7 large parquet and external model weights not rehashed (producer-recorded hashes only).
No workbook, sealed artifact, genome shard or weight bytes opened; no wet-lab measurements reproduced; no submissions, commits or edits outside owned paths.
Deviation: initial standing-instruction read of feng/HANDOFF.md exposed historical wet-lab outcome text, contrary to the stricter task boundary; none was used, copied, fitted or compared. The generator reads only locked authorized inputs.
