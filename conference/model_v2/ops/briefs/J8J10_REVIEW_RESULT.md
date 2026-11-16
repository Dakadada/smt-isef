# J8J10 adversarial review
Verdict: DO-NOT-ADVANCE — 1 HIGH, 5 MED findings; all reproduced with tiny synthetic Python scripts.
Requested pytest: 40 passed, 0 skipped, 0 failed (2.02 s); full feng suite not run for this review.
Scope: source/tests/jobs read-only; only this report added; no submissions or real wet-lab/outcome reads.

HIGH F1 — j8_feature_join.py:324; jobs/j8_full.sbatch:18 — J4 overlap lookup scans the entire earlier chromosome/strand prefix for every J6 row; bounded throughput cannot support the claimed full-genome runtime.
Repro: insert 100/1,000/5,000 disjoint J4 intervals, query beyond their ends: zero hits but 527/5,017/25,017 SQLite VM instructions; EXPLAIN uses only chrom,strand,start<?.
Fix: use an interval/RTree index or a validated maximum-footprint-length lower bound on start; benchmark representative full-density lookups before sizing the full job.

MED F2 — j10_robustness.py:154,175,349 — discovery-selected J6 near-matches receive random-position tails labelled p_value without an observed-position exchangeability/ascertainment model; a stated decoy universe alone does not supply that model.
Repro: flank A*40+T*18+A*40, ASO A*18, start=40,length=18,n=50,seed=1729: p=0.019608; conditioning on J6's <=3 mismatches leaves only six decoys (plus-one tail 1/7, below minimum 20).
Fix: retain descriptive distributions with p_value=None until a prespecified, selection-aware null is justified; overlapping decoys alone are not proof of invalidity.

MED F3 — j10_robustness.py:351 — an upstream failed context fetch aborts the entire report instead of preserving an explicit unavailable null.
Repro: two otherwise valid J3 rows with complementarity_status=ok but target_with_flank=None, footprint_start1_in_context=None: main raises TypeError at None-1; no report exists.
Fix: validate context/offset/channel status per site and emit blocked nulls with the upstream failure reason, retaining the candidate in rank coverage.

MED F4 — j8_feature_join.py:216; j10_robustness.py:33 — derived .done paths are opened before applying safe_path, bypassing the wetlab/terminal boundary through a sidecar symlink.
Repro: allowed.parquet.done -> temporary wetlab/marker.json containing valid window/span: j3_fold_variant returns w120s60; J8 also reads the marker before eventually refusing its hash path.
Fix: resolve and guard every sidecar/catalog/manifest path before opening it; add synthetic forbidden-target symlink tests.

MED F5 — j10_robustness.py:322 — J3 marker parameters are trusted without validating committed rows/hash; hashing the current parquet only records an unverified replacement.
Repro: two-row parquet plus .done {rows:999999,sha256:"wrong",params:{window:120,span:60}} still produces a two-row robustness report.
Fix: verify source commit markers against parquet hashes/counts before reading features; include marker identities in the report contract (and authenticate J9's upstream commit).

MED F6 — j8_feature_join.py:361; j10_robustness.py:329 — incomplete checkpoints can be accepted as completed results.
Repro: remove J8 marker.contract and change context.code_sha256: old parquet is reused; replace a completed J10 JSON with {contract:old_contract}: rerun returns with all result sections absent.
Fix: require an exact, present J8 contract and validate J10 result schema/payload integrity before reuse; rebuild damaged checkpoints.

Coverage gaps: test_j8_feature_join.py:218 checks malformed JSON only; test_j10_robustness.py:137 actually omits the input SHA. Add the above cases and a full-density lookup check.
No additional defect established in strand/midpoint arithmetic, J2 transcript-end indexing, J4 inclusive overlap, Kendall tau-b/sign flips, missing-grid handling, or evaluation_core's supervised refusal.
Unfixed: all findings, by read-only instruction. Limits: synthetic repros only; production-wide cardinality/RSS/runtime not measured; ~37 GB denotes a disk index estimate, not a demonstrated RAM requirement.
