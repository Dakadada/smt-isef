# GENCODE v50 and STAR index

These scripts stage the human GENCODE v50 primary-assembly genome and comprehensive
primary-assembly annotation, then build the Phase 1 STAR index. Large inputs and index
files belong in pool or scratch storage, not HOME.

## Fetch and verify GENCODE

Use the executing account's cluster scratch allocation. `FENG_DATA_ROOT` may be set
explicitly; otherwise the script derives it from `SCRATCH` and fails if neither is
available or writable:

```bash
export FENG_DATA_ROOT="${SCRATCH:?SCRATCH is not configured}/feng/gencode/v50"
bash feng/scripts/fetch_gencode.sh
```

The script resumes `.part` downloads and verifies the official release MD5 before an
archive is atomically promoted. Only after that authoritative check does it compute a
secondary SHA-256 sidecar. Existing archives are rechecked against the pinned MD5 and
are never silently replaced after a mismatch.

Each decompressed file has a provenance sidecar containing the verified archive MD5
and the decompressed-output MD5. Reuse requires both values to match the actual files
and requires format sanity (`>` on the first FASTA line; nine tab-delimited GTF
columns). Missing, stale, or mismatched provenance triggers atomic regeneration.

## Build the STAR index

The default build uses 16 CPUs for two hours (about 32 CPU-hours), requests 60 GB RAM,
and sets `sjdbOverhang=149` for 150-nt reads. Create the relative log directory and
submit from the repository root, replacing both site-specific placeholders:

```bash
mkdir -p feng/logs
export FENG_ENV_PREFIX=/home/dashii/orcd/pool/feng-env
export FENG_DATA_ROOT="${SCRATCH:?SCRATCH is not configured}/feng/gencode/v50"
export FENG_STAR_INDEX_DIR="${SCRATCH}/feng/star/gencode-v50-sjdb149"
export FENG_SJDB_OVERHANG=149
sbatch --export=ALL feng/scripts/star_index.sbatch
```

The job loads `miniforge/25.11.0-0`, activates `FENG_ENV_PREFIX` when supplied, and
resolves STAR and Python from the activated environment or `PATH`. `FENG_STAR_BIN` and
`FENG_PYTHON_BIN` are explicit overrides. The coordinator's current environment prefix
above is a site-specific last resort, not a script default. Storage has no personal-path
fallback: set `SCRATCH`, or explicitly set both data and index locations.

`FENG_STAR_RAM_BYTES` defaults to 55,000,000,000 and can be adjusted within requested
job memory. The checked-in SLURM directives currently use account `mit_general` and
partition `mit_preemptable`, with 16 CPUs, 60 GB RAM, two hours, and requeue enabled.

The job builds in a job-specific staging directory. A completed index is atomically
renamed into place only after the full required STAR file set exists. It then writes
`feng_index_manifest.json` with the decompressed FASTA/GTF MD5s, STAR version, schema
version, and generation parameters (`runMode`, threads, GTF inclusion, `sjdbOverhang`,
and RAM limit). A job skips only when every required file is nonempty and that marker
exactly matches the requested inputs, STAR version, and parameters. Otherwise the final
directory is treated as invalid and rebuilt. Requeued jobs delete only their incomplete
staging directory and restart the non-checkpointable generation safely.
