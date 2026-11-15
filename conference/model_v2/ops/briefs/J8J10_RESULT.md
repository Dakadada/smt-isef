# J8/J10 RESULT — site feature join + robustness engine (claude-opus-5, 2026-11-15)

Status: complete. `40 passed` (J8 22, J10 18). The cut-off Astra draft was largely sound — schemas,
`.done` markers and coordinate conventions all check out against real shards. Fixes: J3 ingest now
projects away `target_with_flank/flank_left/flank_right/duplex_structure` at read time (357 MB ->
~15 MB per shard) and stops calling `to_pylist()` twice per batch; J10 reads the fold window from
the J3 marker `params` instead of hard-coding `w120s60`, refusing unlabelled/undeclared markers;
J10 decoy nulls now really run for complementarity via `j3_scoring.complementarity` (pure
sequence, no backend) behind an orientation guard asserting the rescored observed footprint
reproduces the committed J3 cell; added `jobs/j8_full.sbatch`.

Verified real: J6 52,424 shards / 11,002,701 rows; J3 256 shards / `params.window=120 span=60`;
J2 p1+p2 32+128 shards; J4 `<config_id>/ablation_maps/shard_0.parquet` + `<config_id>/shard_0.done
{ablation}`, 20,395 emitted == 20,395 atlas dirs (no gap); J9 233 rows, all M1..M5 columns present;
`load_loci` on the real manifests 209,760 loci / 0.9 s.

Measured (login node, 1 core, real reads): J3 ingest 14.1 s / 92.7 MB index; J2p1 5.7 s / 14.4 MB;
J2p2 25.5 s / 64.9 MB; J4 config <0.1 s / 0.2 MB; join 1,542 rows/s at 126 B/row. End-to-end
`--bounded`: 48.3 s, peak RSS 1.57 GB, 172 MB index, 344 rows / 2 shards, idempotent rerun 4.5 s
with zero rewrites. J10 bounded: 2.9 s, peak RSS 1.62 GB, 152 KB, idempotent. Full genome
extrapolates to ~4 h, **~37 GB index** (needs `--max-index-gib 48`; the 20 GiB default is too low)
and ~1.4 GB output over 52,424 shards. SCRATCH is 88.4% full, ~119 GB free — fits, not roomy.

Blocked / explicitly degenerate:
- **J2 accessibility join is near-empty genome-wide**: J2 folded u in {8,12,18}, ASO footprints are
  10/14/18 nt, so only F18 matches exactly — 0.01% of a 43,317-row J3 sample (F10 97.24%, F14
  2.74%). Others are explicit `not_available`, no length substitution; surfaced in the output
  manifest. Remedy is a J2 re-run with `ulens=(10,14,18)` — out of scope, not submitted. J3's own
  per-footprint accessibility channel is unaffected and stays the usable source.
- Duplex-energy and accessibility decoy nulls: blocked, both need a folding backend.
- Dinucleotide-shuffle nulls give a distribution but no p-value (randomised Hierholzer draw law is
  non-uniform) — stricter than spec by choice; the shifted-decoy null does emit p-values.
- Supervised metrics blocked by `evaluation_core` (no verified negatives); descriptive only.
- Real J10 output: J3 sample n=300 `unstable` (sign-flip 0.419, 36 comparisons); J9 channel
  agreement `unstable`; M0..M5 all `ok` descriptive; 2/3 probed sites got a shift-null p-value
  (0.0196 at n=50; 0.025 at n=39), the third refused at n=13 < min 20.

Submit (NOT executed): `sbatch conference/model_v2/jobs/j8_bounded.sbatch` ·
`sbatch .../j10_bounded.sbatch` · `sbatch .../j8_full.sbatch` (check `~/orcd/.quota`: needs ~40 GB)
