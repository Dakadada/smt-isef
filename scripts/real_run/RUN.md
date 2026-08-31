# Stage-A real targeted run

This is the default targeted tier. It counts only the run-selected event junction-marker strings directly in raw FASTQ fragments. `build_manifest.py` copies those verified rows byte-for-byte from the master event table into a hashed, run-scoped table; the three-family runnable subset therefore contains only SMN2 exon 7, while DMD remains deferred. It does **not** read or validate a genome FASTA, GTF, STAR index, HISAT2 index, or container. Full alignment remains an opt-in tier with the original all-assets gate.

The two GEO mappings are embedded exactly at GSM level from `docs/dataset-audit.md`. The four verified PRJNA1052769 dCas13d endpoint RNA-seq run mappings are also embedded. The currently runnable subset is `GSE307431_PRJNA1321660,GSE167762_PRJNA705187,PRJNA1052769_dCas13d`. DMD (`PRJEB13980_U7`) is **DEFERRED** pending a verified run-level arm mapping: all 12 ENA runs share the sample title `U7E53_OFFTARGET_ANALYSIS`, so treated and untreated arms cannot be inferred from ENA, and the deposited comparison is vector-control-confounded. Selected families still fail closed on every absent/ambiguous arm, pair set, covariate, URL, MD5, or endpoint count; omitted families are not fetched or validated.

Run from the repository in a networked login/HPC session:

```bash
cd /home/dashii/smt-isef
/usr/bin/python3.12 scripts/real_run/build_manifest.py \
  --out "$SCRATCH/smt-stageA-real" \
  --families GSE307431_PRJNA1321660,GSE167762_PRJNA705187,PRJNA1052769_dCas13d
export RUN_DIR="$SCRATCH/smt-stageA-real"
job_id=$(sbatch --parsable --partition=mit_data_transfer --export=ALL,REPO="$PWD",RUN_DIR="$RUN_DIR" scripts/real_run/run_stageA.sbatch)
scripts/real_run/monitor_requeue.sh "$job_id" "$RUN_DIR"
```

The `mit_data_transfer` partition is required because the job fetches FASTQs from ENA; ordinary compute/preemptable nodes may not have outbound internet. The batch script also declares this partition, while the explicit submit flag makes the requirement visible at submission time.

`--metadata-dir DIR` makes the builder consume saved `DIR/PROJECT.tsv` ENA filereports and is intended for offline testing. On preemption, Slurm signals the job five minutes early; it requests requeue. Downloads use resumable `.partial` files. Completed count shards are checksum-checkpointed, and each accession follows fetch → archive MD5 → marker count → committed count checksum → FASTQ deletion. Re-submit the identical `sbatch` command if site policy terminates without delivering the signal; checkpoints prevent repeated work.

The hard aggregate scratch cap is 300 GiB, comfortably below 700 GB. Processing is one accession at a time and normally retains only one accession's compressed FASTQ plus small count/checkpoint files. Budget roughly 12–36 hours total (archive throughput and FASTQ sizes dominate); a typical live footprint should be tens of GiB and can never pass the 300-GiB guard. The heartbeat updates every 10 minutes.

After completion:

```bash
test -f "$RUN_DIR/STAGE_A_COMPLETE"
column -s $'\t' -t < "$RUN_DIR/work/harmonized/tidy_splicing.tsv" | less -S
awk -F '\t' 'NR==1 || $5=="treated"' "$RUN_DIR/work/harmonized/tidy_splicing.tsv"
```

Interpret `magnitude` as raw inclusion dPSI (treated minus matched control). The SMN2 ASO arms should increase exon-7 inclusion, so their dPSI should be positive. The dCas13d D[333,357] guide targets the ISE-D1 region and should reduce SMN2 exon-7 inclusion, so its dPSI should be negative, as reported by Recinos et al. An unexpected sign is a pipeline-bug alarm, not a result to silently relabel. A zero-marker targeted result is recorded as `missing_junction`, never as a biological zero, and therefore emits no PSI/dPSI label.
