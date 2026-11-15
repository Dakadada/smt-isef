# ASO Computational Deliverables Implementation Plan

> **For agentic workers:** Execute inline in the existing working tree. Do not
> dispatch subagents, create a worktree, or commit. Track each task with TDD.

**Goal:** Produce candidate-tagged thermodynamic, off-target, and
dinucleotide-scramble deliverables for four ISS-N1 ASOs, plus a future PSI stub.

**Architecture:** A small importable package under `conference/aso_pipeline`
separates sequence definitions, thermodynamics, shuffling, and streamed genomic
overlap. Thin CLIs produce stable TSV/JSON artifacts locally and on Slurm.

**Tech Stack:** Python 3.11, IntaRNA 3.4.1, ViennaRNA/RNAduplex 2.7.2, Bowtie 1,
NCBI BLAST+, GENCODE v50, Slurm.

**Spec:** `conference/docs/2026-09-26-aso-computational-deliverables-design.md`

## Global constraints

- Never access sealed/frozen Arm E material or wet-lab PSI.
- Every modeled output is `[candidate]`; absent values are never fabricated.
- Use RNA:RNA at 37 C, 0.089 M operational TBE salt, and the fixed sensitivity
  grid 0.025/0.050/0.089/0.100/0.150 M.
- Use RNG seed 20260926 and exactly 100 null draws per ASO.
- Exhaustive fixed `-v 3` headline counts are mandatory for F18 and F10.
- Saturated non-headline queries are skipped with a nonnumeric state.
- Run only per-file pytest commands, never the full Feng suite.
- Do not commit.

---

### Task 1: Sequence definitions and thermodynamics

**Files:**
- Create: `conference/aso_pipeline/sequences.py`
- Create: `conference/aso_pipeline/thermo.py`
- Create: `conference/aso_pipeline/item1.py`
- Test: `conference/tests/test_sequences.py`
- Test: `conference/tests/test_thermo.py`

**Interfaces:**
- `derive_target(aso_rna, context_rna) -> TargetWindow`
- `score_duplex(aso, target, intarna, rnaduplex, salts, chemistry) -> rows`
- `run_item1(output_dir) -> list[dict]`

- [ ] Write target-position and real-binary energy tests with hand-derived
  target sequences and validate that missing/zero IntaRNA states stay explicit.
- [ ] Run each test file and observe the expected import failure.
- [ ] Implement normalized immutable ASO records, exact target derivation,
  RNAduplex structure/energy parsing, explicit salt calls, and IntaRNA reuse.
- [ ] Run each test file green and emit raw/sensitivity/ranking tables plus the
  item-1 README and machine-readable provenance.

### Task 2: Dinucleotide null thermodynamics

**Files:**
- Create: `conference/aso_pipeline/shuffle.py`
- Create: `conference/aso_pipeline/item6.py`
- Test: `conference/tests/test_shuffle.py`
- Test: `conference/tests/test_item6.py`

**Interfaces:**
- `dinucleotide_shuffles(sequence, n=100, seed=20260926) -> list[str]`
- `summarize_null(real, values, lower_is_stronger) -> dict`

- [ ] Write deterministic, exact-dinucleotide, length, multiplicity, quantile,
  and direction-aware percentile tests; observe expected import failures.
- [ ] Implement randomized Euler-trail shuffling with replacement, retaining
  draw identity and duplicate multiplicity.
- [ ] Run the two test files green; execute all 400 local thermodynamic null
  draws and write raw and summary tables plus the item-6 README.

### Task 3: Streamed off-target counting and Slurm workflow

**Files:**
- Create: `conference/aso_pipeline/offtarget.py`
- Create: `conference/aso_pipeline/item3.py`
- Create: `conference/scripts/item3_offtargets.sbatch`
- Test: `conference/tests/test_offtarget.py`
- Test: `conference/tests/test_item3.py`

**Interfaces:**
- `expected_random_hits(length, mismatches, genome_size) -> float`
- `classify_query(..., mandatory) -> saturation state`
- `GeneIndex.from_gtf(path)`, `add_hit(...)`, `counts()`
- `run_item3(...) -> compact TSVs`

- [ ] Write literal probability, saturation override, SAM/BLAST coordinate,
  strand-compatibility, and duplicate-query expansion tests; observe failure.
- [ ] Implement preflight manifest generation and streaming parsers/counter.
- [ ] Run both test files green.
- [ ] Write a checkpointed, streamed CPU sbatch that creates Bowtie/BLAST
  indexes once, runs unique nonsaturated queries, and atomically publishes
  compact results. Submit it and record the numeric Slurm job ID.

### Task 4: PSI stub, documentation, and verification

**Files:**
- Create: `conference/aso_pipeline/psi_stub.py`
- Test: `conference/tests/test_psi_stub.py`
- Create: `conference/artifacts/item1_README.md`
- Create: `conference/artifacts/item3_README.md`
- Create: `conference/artifacts/item6_README.md`

**Interfaces:**
- `fit_dg_psi(rows, dg_by_aso)` validates the future row schema then raises
  `PsiDataUnavailable` until explicitly enabled in a later task.

- [ ] Write schema-validation and always-refuse tests; observe import failure.
- [ ] Implement only the typed/refusing stub and run its test file green.
- [ ] Run every conference test file separately and record exact counts.
- [ ] Validate artifact columns, `[candidate]` tags, null counts, checksums,
  software versions, Slurm queue state, and all README assumptions.

