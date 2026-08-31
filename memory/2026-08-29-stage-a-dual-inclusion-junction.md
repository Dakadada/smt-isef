# Stage-A dual-inclusion junction debug report

- Status: DONE
- Symptom: a real SMN2 exon-7-spanning fragment aborted targeted quantification after matching both inclusion junction markers.
- Root cause: the targeted counter incorrectly treated the two inclusion markers as mutually exclusive categories, although both support the same included isoform.
- Fix: retain raw per-marker diagnostics, deduplicate either or both inclusion hits into one inclusion fragment, discard and tally only inclusion-plus-skipping contradictions, and compute PSI from fragment-level inclusion and skipping counts.
- Evidence: the read-shaped single-end, contradictory-marker, mixed-PSI, and paired-end regression tests pass; `make smoke-test` passes all 95 tests with Python 3.12.
- Related: targeted progress checkpoints now carry a counting-semantics version; incompatible old partial progress is replayed rather than reused with unreconstructable fragment counts.
