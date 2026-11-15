# J4-SCALE result (builder + bench + array runner; nothing submitted)

New files only: `conference/model_v2/j4_build_configs.py`, `conference/model_v2/tests/test_j4_build_configs.py`, `conference/model_v2/jobs/j4_scale_bench.{sbatch,py}`, `conference/model_v2/jobs/j4_scale_array.sbatch`. Engine untouched; N-masked rows keep `evidence_class=sequence_ablation_proxy` (engine-set).

**Tests: 31 passed, 0 failed, 0 skipped** — `PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1 /home/dashii/orcd/pool/smt-isef-env/bin/python -m pytest conference/model_v2/tests/test_j4_build_configs.py -q` (0.31 s). Hand-checked toy loci cover +/- strand site orientation, pad-clipped skip, non-contiguous exon indices, single-exon transcript, chrM/scaffold context overrun, grid arithmetic (19/18/17=54), footprint parity with `atlas.tasks()`, engine-dataclass validation, byte-identical rebuild, sha256 chain, CLI.

Conventions: donor = last exonic base, acceptor = first exonic base, transcript orientation (`+`: acceptor=exon_start, donor=exon_end; `-`: swapped), matching `gpu_scores.EXON7_LAST_BASE_1B=70076574` for SMN2 exon 7. Exon recovered as `start1+100 .. end1-100`; `start1<=100` -> explicit `padding_clipped_ambiguous` skip row, never guessed. Terminal non-sites never emitted (exon 1 has no acceptor, exon N no donor; N from contiguous `ex1..exN`, gapped transcript skipped whole).

Interval rule (one absolute genomic tiling interval per config): grid `[site-50, site+50]` 1-based inclusive; window center = site, distance = 50; footprint starts `range(grid_start, grid_end-length+2, stride)` so every footprint lies fully inside the interval. Donor/acceptor of one exon are >50 nt apart in general, so each exon yields up to **two** configs, one per site — never two grids in one config.

Build already run (CPU-only, ~0.3 min, no GPU): `PYTHONPATH=feng/src:. python -m conference.model_v2.j4_build_configs` -> `~/orcd/scratch/feng/model_v2/j4_configs/` with `shards/shard_0000..0318.jsonl`, `j4_configs_manifest.tsv`, `j4_config_shards.tsv`, `j4_build_provenance.json`. 20 395 configs (20 394 exon-site + 1 curated SMN2), 319 shards x 64, 123 skipped loci (all `single_exon_transcript_no_splice_site`; zero clipped in phase 1), loci_sha256 `0320754c…7fc3a`. Curated `SMN2_exon7_issn1_fine`: chr5 +, center/site 70076574 donor, grid 70076554–70076631, stride 1, lengths 8/10/14/18 = 266 footprints.

**Estimated total inferences = 1 121 937** (= footprints + 1 reference per window): 20 394 x 55 = 1 121 670 atlas + 267 curated. Footprints 1 101 542; atlas rows 4 406 168 (= footprints x 4 tissue channels). Per array task: 64 configs, 3 456 footprints, 3 520 inferences.

Exact submit commands (NOT executed):
```
sbatch conference/model_v2/jobs/j4_scale_bench.sbatch
NSHARDS=$(grep -vc '^#' ~/orcd/scratch/feng/model_v2/j4_configs/j4_config_shards.tsv)   # 319
J4_BATCH_SIZE=<winning bench batch size> sbatch --array=0-$((NSHARDS-1))%8 --time=<from bench> \
    conference/model_v2/jobs/j4_scale_array.sbatch
```
Bench: `mit_preemptable`, 1 GPU, 40 min cap, 50 real configs (curated excluded), sweep `--batch-size 4,16,32`, per-sweep private output + `reference_scores` dir (no cross-sweep credit, no `.done` short-circuit), 300 s/sweep budget with partial reporting; prints `inferences/s`, `footprints/s`, peak/reserved GPU GiB, and `model_load_seconds`.

Array runner: no `--array`/`--time` directive by design (both set at submit from bench); hard-fails if `SLURM_ARRAY_TASK_ID` unset. Task = `shard_%04d.jsonl`; pre-pass writes `<atlas>/<config_id>/config.json` + plan TSV whose sha256 equals `config_sha256` in the manifest; one engine invocation per config into `~/orcd/scratch/feng/model_v2/j4_atlas_v1/<config_id>`; idempotent skip requires `shard_0.done` + `ablation_maps/shard_0.parquet` + matching config hash (verified over all four branches offline); requeue trap identical to the PASSED `j4_sentinel.sbatch` and kept installed between configs; a failing config is logged to `failures_task_<id>.tsv` and the shard continues, task exits 1.

Open flag (not fixed, spec said one invocation per config): each invocation reloads the Pangolin ensemble, so startup is paid 20 395 times. Bench reports `model_load_seconds` — if it is a large fraction of ~55-inference-per-config runtime, switch the array to an in-process driver that loads the scorer once per task.
