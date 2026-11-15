# J3 implementation plan

Spec: `../ops/briefs/J3.md`; scientific boundaries: `../SCIENTIFIC_LOCK.md`.
The supplied task authorizes implementation in the named files; no commit or job submission.

Architecture: retain J6 columns and append three independent, nullable channels,
each with status/reason/evidence. Fetch strand-oriented flanks through the existing
indexed FASTA reader. Use a persistent feng-env ViennaRNA worker, with independent
requests for duplexfold and RNAplfold footprint probability. SQLite caches completed
channels by normalized sequence triple within a tool/configuration namespace.
Hash site IDs for deterministic output shards. Validate input completion markers;
publish output parquet then a checksum-bearing atomic completion marker.

Constraints: GRCh38 / GENCODE v50 only; no outcomes or blended score; explicit
unmodified RNA assumption; no heavy login compute; no edits to feng/ or J6/J2.

Review focus: missing/invalid sequences, wrong strand/reference, backend exceptions
and nonfinite outputs, interrupted or stale artifacts, configuration-sensitive cache.

- [x] Write failing tests in `../tests/test_j3_scoring.py`: hand-checkable pairing,
  independent channel failures, normalized cache, retry, strand-oriented context,
  backend footprint indexing, checksum/restart, and sharded union.
- [x] Implement `../j3_scoring.py` and `j3_vienna.py`; run focused tests.
- [x] Add `j3_sentinel.sbatch`, mit_quicktest, ten minutes, --requeue-safe.
- [x] Run model_v2 tests and required feng gate once; inspect owned files and
  write <=30-line `../ops/briefs/J3_RESULT.md` with exact results and submit command.

Rulings: defaults are 80 nt flanks, RNAplfold W=120/L=60 (clipped to context),
37 C ViennaRNA defaults. Cache successful/not-applicable channels only; transient
failures retry. Duplexfold may bind anywhere in target+flank, explicitly labelled;
it does not enforce footprint anchoring or represent modified ASO potency.
J2 output is endpoint-indexed; compute the exact candidate footprint locally rather
than joining approximate J2 u=8/12/18 grids. No J2 outputs are modified.
