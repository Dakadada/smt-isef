# J8/J10 adversarial rereview
Verdict: ADVANCE — F1–F6 fixed in the reproduced scenarios; no new defect established.
F1 FIXED — j8_feature_join.py:337: recreated 100/1,000/5,000 disjoint-interval miss: 58/41/41 SQLite VM instructions, 0 hits; 216,000-row density test passed.
F2 FIXED — j10_robustness.py:171: original A40/T18/A40 case yields blocked, p_value=None, 50 descriptive decoys.
F3 FIXED — j10_robustness.py:331: missing context/offset yields blocked nulls with upstream reason; two-row CLI regression retains both candidates.
F4 FIXED — j8_feature_join.py:216,240; j10_robustness.py:33: forbidden .done symlink and SQLite -wal/-shm/-journal symlinks raise ValueError before opening; manifest/output symlink tests passed.
F5 FIXED — j10_robustness.py:318,407: incorrect rows or SHA in J3 marker or J9 manifest rejects before feature read; identities include marker hashes.
F6 FIXED — j8_feature_join.py:375; j10_robustness.py:360: missing J8 contract rebuilds; contract-only and damaged J10 checkpoints rebuild, intact checkpoint reused.
Boundary probe: atlas [100,117] vs [99,99]/[99,100]/[117,117]/[117,118]/[118,118] produced 0/1/1/1/0 overlaps.
New findings: none established. Unfixed: none identified.
Requested pytest: 73 passed, 0 skipped, 0 failed (5.68 s).
Limitations: original review supplied repro descriptions, not standalone scripts; cases were recreated and checked alongside regression tests. Synthetic density is not production-wide runtime/RSS proof; sidecar checks do not establish resistance to concurrent symlink replacement.
