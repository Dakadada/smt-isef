J5-LIFT: code tested; real-data sentinel DO-NOT-ADVANCE (Slurm submission failed).
`j5_rbp_context.py:151`: accession-keyed lookup is mandatory; missing accession, wrong assay/cell/format/MD5 or assembly conflicts refuse ingestion; no header fallback.
Lookup SHA256: 4dd01bebd4c218c7bad46d01f3300df735a696b1d153ea912100ffa297cf04ca (152 entries: 114 peaks + 38 KD; hg19).
Chain SHA256: 5c0598e500ceb5a78c73086929e8ef993aec309bcafb595139b53d440b125a1d; pinned CLI input at ~/orcd/scratch/feng/liftover/hg19ToHg38.over.chain.gz.
`j5_rbp_context.py:193,335`: uses feng.rbp.liftover read-only; unique ungapped same-size lifts only; derived peak coordinates retain original hg19 BED blocks and chain SHA.
Unmapped/multi-mapped/split/size-changed counts and strand flips are explicit; SQLite lift_rejections preserves rejected source rows/payloads; partial coverage never supplies False misses.
Accounting: input_intervals = outside_region + lifted + rejected_intervals; split/size_changed diagnostics overlap; discarded_event_blocks are lifted blocks withheld with rejected events.
KD inspection ONLY listed ENCFF915LBL.tar.gz members: A3SS/A5SS/MXE/RI/SE.MATS.JunctionCountOnly.txt; no bulk extraction or real event readout.
KD exon blocks use the same lift; failed/incoherent events are withheld atomically; unsupported/empty formats refuse with not_available reason; sample-order direction remains unavailable without documentation.
Cell lines remain separate; PCBP2/K562 stays no_local_assay/missing; KD retains knockdown_not_steric_block and observed raw measurements.
`jobs/j5_sentinel.sbatch:15`: self-contained mit_quicktest regional index + four queries (F18/F14/F10, chr5:69076584-71076584); output ~/orcd/scratch/feng/model_v2/j5_sentinel.
Sentinel scans source files but only maps chain-preimage candidates for the +/-1 Mb window; excluded blocks counted outside_region; nearest distances labelled regional.
Attempted after green J5 tests: `sbatch --parsable conference/model_v2/jobs/j5_sentinel.sbatch`; exit 1: Unable to contact slurm controller (connect failure). Slurm ID: not_available.
`scontrol ping`: Address family '0' not supported; Slurmctld(primary) at slurm003 is DOWN from this session. No sentinel ran; real lift counts: not_available.
`tests/test_j5_rbp_context.py:371`: requested exact thread-limited pytest command: 40 passed, 0 skipped, 0 failed (2.94 s); includes BED boundaries, reverse strand, gap/multi/size, rejection audit, KD, and restart cases.
Required feng gate: 581 passed, 5 skipped, 1 failed (235.23 s): test_cleanroom.py:23/test_implementation_has_no_legacy_reference[source_file42] flags src/ in feng/src/feng/source_auth.py; pre-existing, outside ownership, unchanged.
All three J5 sbatch scripts pass bash -n; owned diff passes git diff --check. No full/index job submitted, no commits, no wet-lab access, no feng source edits.
Limitations: actual KD schemas, regional lift counts, production runtime/RSS remain unverified until Slurm works; GC controls require optional reference/gene bounds; no biological GO claim.
