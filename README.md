# Uniform raw-read reanalysis

The headline prospective artifact is the generative splice-switching design engine
in [`src/design_engine/README.md`](src/design_engine/README.md). It ranks bounded
ASO/gRNA/U7 designs using the verified physics adapter and off-target enumerator,
while carrying abstention, uncertainty, therapeutic polarity, Pareto tradeoffs, and
an explicit non-validation banner through every report.

The cross-disease prospective panel and designability map are documented in
[`src/design_panel/README.md`](src/design_panel/README.md). Its target registry is
the path to add newly verified targets during a networked curation session;
coordinate-unverified entries remain explicit required inputs and cannot produce
designs.

This repository contains an offline-testable, resumable implementation of the locked raw-read reanalysis in `docs/architecture.md` and `docs/benchmark-protocol.md`. It resolves ENA runs at execution time, checksum-verifies FASTQs, quantifies the frozen SMN2 exon-7 and DMD exon-53 events, computes study-local PSI/treated-minus-control ΔPSI, and emits one hierarchy-aware tidy table. It does not make a powered generalization or mechanism claim.

## Off-target enumerator contract

The off-target enumerator scans the complete supplied pre-mRNA sequence. Candidate generation is exhaustive only inside the declared `max_edits`, `max_indel_bases`, and inclusive `min_alignment_length`/`max_alignment_length` bounds. Every alignment at every start and length inside those bounds is retained. Expression, localization, RBP support, chemistry score, and their configured thresholds are post-enumeration annotations/ranking signals; they never remove candidates or change recall. Benchmark positives outside a declared bound remain in the denominator and appear in `out_of_bound_misses`.

All candidate, evidence, and positive-site coordinates are **0-based, half-open indices into the strand-oriented transcript sequence**. Transcript sequence is always written 5′→3′, so a minus-strand genomic feature is reverse-oriented when mapped. `genomic_to_transcript_interval` performs checked genomic-to-spliced-transcript conversion and rejects intronic, unmapped, or non-contiguous intervals. Transcript identifiers are version-exact: callers must freeze and supply the same accession version for sequence, evidence, and validation positives; the enumerator never strips or guesses a version. Genomic inputs must first be mapped separately for each transcript isoform using that transcript's exon map.

Known chemistries use the named coarse scoring parameters. Unknown chemistry produces no `estimated_delta_g` by default. Explicit `allow_unknown_chemistry_heuristic=true` opt-in produces an estimate labeled `heuristic-non-calibrated`; downstream consumers must not treat it as chemistry-specific data.

Recall matching requires both candidate endpoints to fall within the declared coordinate tolerance. Empty positive sets have undefined recall. A target-met verdict is issued only when positive-set provenance is nonempty and its completeness and blinding flags are true; otherwise the numerical recall is provisional and `target_met` is undefined. Precision is reported independently and only for a completely labeled prediction universe.

## Smoke test (no network or bioinformatics tools)

```bash
make smoke-test
```

The test uses local `file://` FASTQs, a mock ENA table, exact synthetic junction markers, and blocks HTTP access. No real alignment or scheduler submission occurs.

The cis-prior/nonstructural comparison arena is specified in
`src/baselines/README.md`. Its outer folds are connected components across study,
locus, intervention, chemistry, and biological-sample identities; all dimensions are
asserted disjoint. Labels require a proven matched-control contrast, baselines use a
full regularized representation and real ridge EDF, SpliceAI/Pangolin fail closed
without verified weights plus official golden outputs, and the eight-family CI and
magnitude gate cannot be lowered.

The decisive shared-fold physics-vs-baseline harness and its networked-session
run contract are documented in `src/benchmark/README.md`. It emits deterministic
JSON and text reports and cannot promote the current four-family floor beyond a
descriptive `feasibility-only` verdict.

## Real HPC run

1. Copy `config/pipeline.yaml` to an immutable run configuration. Replace every `REQUIRED_BEFORE_REAL_RUN` provenance hash and reference/index path; freeze SHA-256 hashes for FASTA, GTF, event definitions, container, and indexes. The source hash is the SHA-256 of the sorted per-file SHA-256 manifest (`find src/rawread_reanalysis -type f -name '*.py' -print0 | sort -z | xargs -0 sha256sum | sha256sum`). Inspect the ENA sample fields and finish/freeze the sample rules—unmatched samples deliberately stop computation. For targeted mode, populate exact event `junction_markers`; for full mode, populate `config/event_junctions.GRCh38.tsv` with verified regtools BED coordinates.
2. Set `paths.scratch_dir` to node/local or project scratch and `paths.work_dir`/`checkpoints` to persistent storage. Select `alignment.tier`: `full` with `full_tool` `STAR` or `HISAT2`, or the quota fallback `targeted`. `salmon` and `kallisto` execution is wired, but intentionally stops before PSI unless a frozen event-transcript conversion is supplied; it must not invent junction counts.
3. Resolve and review the frozen manifests before downloads:

   ```bash
   PYTHONPATH=src python3 -m rawread_reanalysis resolve --config run.yaml
   ```

   The resolver calls ENA `filereport?result=read_run`, records accession/URI/archive MD5/size/retrieval time, and filters by `library_strategy`. Confirm endpoint counts: 20, about 12, about 4, and 12. PRJNA1052769 AMPLICON/pooled screen rows must remain excluded.
4. Run stages independently or together. Each completed item has an atomic fsynced checkpoint with config and output hashes; reruns validate rather than trust existing files.

   ```bash
   PYTHONPATH=src python3 -m rawread_reanalysis fetch --config run.yaml
   PYTHONPATH=src python3 -m rawread_reanalysis quantify --config run.yaml
   PYTHONPATH=src python3 -m rawread_reanalysis compute --config run.yaml
   PYTHONPATH=src python3 -m rawread_reanalysis harmonize --config run.yaml
   ```

   Submit CPU work with `sbatch --export=CONFIG=run.yaml,STAGE=all slurm/run_normal.sbatch`, or use `slurm/run_preemptable.sbatch`. Despite cluster GPU availability, raw-read work requests zero GPUs. The preemptable wrapper atomically updates a heartbeat every 10 minutes and preserves run-level checkpoints. A preempted aligner invocation restarts only that incomplete run.
5. Monitor without mutation first, then enable bounded resubmission explicitly:

   ```bash
   python3 slurm/monitor.py --job-id JOBID
   python3 slurm/monitor.py --job-id JOBID --watch --submit --max-retries 3 --interval 60
   ```

   The monitor uses `squeue`/`sacct`, detects a heartbeat stale for 15 minutes, resubmits only retryable jobs, caps retries at three, and quarantines corrupt/checksum-invalid checkpoint manifests. Use a site-specific lowercase partition name if MIT's configured name differs from `PREEMPTABLE`.

## Outputs and conventions

Stages write `manifests/`, scratch FASTQs, per-run `counts/`, `psi/sample_effects.tsv`, and `harmonized/tidy_splicing.tsv`. The latter is keyed by family/sample/intervention/event and carries raw PSI, observed ΔPSI magnitude, primary direction, coverage, batch covariates, three clustering keys, tier/approximation, and confound flags. PSI is cassette-exon inclusion in percent: `100 × (I/2)/(I/2 + S)`, where `I` is the sum of the two inclusion junction counts and `S` the skipping junction count. ΔPSI is treated minus the matched within-study/contrast control mean; data are never rescaled across studies. Primary direction uses `|ΔPSI| < 5` as no-change, with 2/10-point sensitivity labels retained upstream. The ≥20 informative-read flag is per sample; downstream replicate-aware beta-binomial inference remains a separate benchmark-analysis responsibility.

All current studies have `magnitude_eligible=false` because none has matched endpoint occupancy evidence. Observed ΔPSI is retained as a readout, but MAE, rank, calibration, occupancy, and mechanism evaluation must remain disabled. Targeted results carry an approximation label and require concordance against a full-alignment pilot.

PRJEB13980 is never silently treated as clean: every DMD output row carries `VECTOR_CONTROL_CONFOUND`, stating that deposited RNA-seq contrasts rAAV3b-U7-E53 with untreated cells and lacks an RNA-sequenced U7-scramble/vector control. Tissue/donor must also be added as batch covariates after manifest review.

## Footprint and cleanup

The locked envelope is roughly 48 endpoint RNA-seq runs out of 76 archive runs. Plan for **0.5–1 TB persistent** results/references and **2–4 TB scratch** at full tier (the architecture's conservative cluster envelope), 16–32 CPUs and 64–128 GB RAM per shard, and roughly 500–1,500 CPU-hours. The fetch guard requires free space of 2.2 times reported compressed FASTQ bytes (and at least 10 GB). With cleanup enabled, verified FASTQs are deleted only after their per-run count output and checkpoint are atomically committed. Targeted mode can be much smaller but its actual footprint depends on retaining compressed input and must be measured in the pilot.

Pinned expected tool versions live under `tool_lock`; production should run them in an immutable Apptainer image and record the SIF hash. No credentials or secrets are accepted or stored.
