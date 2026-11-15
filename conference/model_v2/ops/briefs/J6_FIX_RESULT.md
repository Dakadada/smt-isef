J6-FIX — five review code findings addressed; production memory target remains unverified.
- `j6_site_universe.py:319`: require pinned GRCh38 primary chromosome length signature before production fetch; reject mismatched assemblies.
- `j6_site_universe.py:390`: synthetic references require explicit opt-in, emit assembly=synthetic and separate site IDs; immutable run modes prohibit mixed reuse; synthetic sentinel acceptance is rejected.
- `j6_site_universe.py:567`: frozen configuration plus atomic run_manifest.json lists every expected shard; changed window/region/inputs, legacy directories, and unlisted shards are rejected without deleting prior results.
- `j6_site_universe.py:457`: SMN evidence uses the full gene universe independently of regional output; missing other paralog yields null instead of false.
- `j6_site_universe.py:378`: cache universe hash, chromosome interval tree, and SMN sequences once per worker; stream unsharded output via iter_sites and Parquet via 512-row batches.
- `j6_site_universe.py:533`: cache run manifest digest/stat once; per-shard verification avoids repeatedly reading/hashing the complete expected-shard list; changed manifests cannot commit.
- `j6_site_universe.py:608`: validate completion marker object/count/hash/expected metadata and Parquet schema; [], null, malformed markers and damaged shards rebuild.
- `j6_site_universe.py:185`: release consumed GTF file cache in bounded 8-MiB chunks with optional POSIX_FADV_DONTNEED; parsing/records remain unchanged.
- `tests/test_j6_site_universe.py:186`: corruption, assembly, restart, paralog, streaming, scaling, worker reuse, cache advice, and manifest-mutation regressions.
Requested gate: **62 passed, 0 skipped, 0 failed**, 3.03 s; exact thread-limited command from the task used.
`jobs/j6_sentinel.sbatch`: bash -n passed; owned-file git diff --check passed; all existing J3 record field names and Arrow types preserved.
`jobs/j6_profile.py`: separate bounded runs with 1,000 / 8,000 genes: setup 0.0113 / 0.0937 s; 20 interval queries 0.0020 / 0.0021 s, 455 / 572 traversal calls.
Each profile also streamed 19,997 dense rows; peak process RSS **118,188 / 124,784 KiB** (115.4 / 121.9 MiB), both below 2 GiB.
Memory cause inference: the known GTF is 4,688,772,094 bytes (4.37 GiB); its sequential scan can populate file cache counted in a Slurm cgroup peak. Prior sentinel had only 6,145 rows in 22 shards.
Unresolved: the reported 4.19-GB production peak was not reproduced/attributed conclusively; sacct failed connecting to slurm003. Production <2 GB requires a fresh sentinel run; POSIX advice is advisory/filesystem-dependent.
Reference limitation: the pinned primary-contig length signature identifies assembly, not FASTA base-content authenticity; immutable FASTA/FAI filesystem identities still apply. enumerate_sites retains its sorted-list compatibility behavior; use iter_sites for streaming.
Prepared sentinel: fresh output .../model_v2/j6_sentinel_fix_v1; fixed region chr5:70000000-71100000, six-anchor acceptance, thread limits, stage RSS and /usr/bin/time -v reporting.
Later submission command: `sbatch conference/model_v2/jobs/j6_sentinel.sbatch` — **not executed**. No commits/submissions made by this session, no sealed-data access, and no edits outside J6 owner files and this requested report.
Concurrent workspace change: HEAD advanced externally to 3c2a3c8 and now contains the tested J6 code; this session did not run git commit. The requested result report remains untracked.
