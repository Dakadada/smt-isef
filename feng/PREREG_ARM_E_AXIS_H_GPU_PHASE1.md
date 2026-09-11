# Preregistration seal — Arm E axis H (Pangolin deep-splice scorer), GPU Phase-1 specificity

**STATUS: DRAFT — NOT SEALED.**

Seal fields (placeholders; fill only in the later, separate seal commit):

- Code commit SHA: `<STAMP CODE COMMIT SHA AT SEAL TIME>`
- Authenticated `feng/src` tree SHA-256: `<STAMP AT SEAL TIME>`
- Executed launcher blob SHA-256: `<STAMP AT SEAL TIME>`
- Authenticated recompute job: `<STAMP AT SEAL TIME>`
- Manifest SHA-256: `<STAMP AT SEAL TIME>`
- Seal date: `<STAMP AT SEAL TIME>`

Do NOT treat this document as sealed until re-review returns SEAL-READY and a separate commit changes
the status to SEALED and replaces every placeholder above. A commit cannot attest to its own producing
SHA: this DRAFT therefore ships in the self-contained CODE commit, then a separate later SEAL commit
performs only the DRAFT→SEALED flip and stamps the already-existing CODE commit SHA, authenticated
recompute job, source-tree SHA-256, launcher SHA-256, manifest SHA-256, and seal date.

This preregisters, verbatim from code, the three `candidate` modeling choices of Arm E axis H
**before** any burden value is read. All are additive `candidate` score definitions; they do NOT
touch the frozen Arm E therapeutic-window endpoint. Final freeze is at the month-8 FREEZE tag.

Source of truth: `feng/src/feng/titration/gpu_scores.py`, `feng/src/feng/titration/arm_e_hoff.py`.
Runtime authentication: `feng/src/feng/source_auth.py`, `feng/scripts/arm_e_hoff.sbatch`, and
`feng/scripts/submit_arm_e_hoff.sh`.
Anchor frame: SMN2, GRCh38 GENCODE v50, chr5, + strand.

---

## Frozen anchor constants (`gpu_scores.py:45-53`)

| Constant | Value (1-based) | Meaning |
|---|---|---|
| `CHROM` | chr5 | SMN2 |
| `EXON7_LAST_BASE_1B` | 70076574 | exon 7 last base (+ strand) |
| `EXON7_DONOR_1B` | 70076575 | intron-7 first base (kept for provenance; NOT the readout) |
| `ISS_N1_START_1B` | 70076584 | nusinersen 18-mer target window start |
| `ISS_N1_END_1B` | 70076601 | nusinersen 18-mer target window end |
| `PANGOLIN_CONTEXT` | 5000 | flank each side fed to Pangolin |
| `DEFAULT_DISTANCE` (`d`) | 50 | Pangolin `-d`; output covers ±d around center |
| `DEFAULT_READOUT_1B` | = `EXON7_LAST_BASE_1B` (**70076574**) | readout position — Pangolin labels 5′ **donor usage at the LAST EXONIC base**, not the first intronic base (review round-1 finding 1) |

Context window fed to Pangolin: `window_bounds_1b` = `[center−(context+d), center+(context+d)]`,
length `2*(context+d)+1`, so `compute_score` returns exactly `2*d+1` positions centered on `center`.

---

## Choice 1 — N-ablation footprint proxy (`footprint_mask`, `gpu_scores.py:79`)

**Named honestly (review round-1 finding 2): this is an N-ABLATION SENSITIVITY PROXY, not a validated
steric-occlusion prediction.** Replacing the ASO footprint with `'N'` presents Pangolin with
missing/unknown sequence at those positions; the resulting splice-score change is a proxy for how much
the model's splice call *depends on* that footprint's sequence. It is NOT a mechanistic claim that an
occupied base behaves as an all-zero input — an ASO-bound base is physically still present. Reports must
not upgrade this proxy to a steric-occlusion prediction without an occupancy-aware model or validation.

**Definition.** Replace footprint bases `[footprint_start_1b, footprint_end_1b]` (1-based inclusive) in
the reference context window with `'N'` (one-hot all-zero column, `IN_MAP[0]`). Length-preserving (no
splice-coordinate reindexing), fail-closed (footprint outside window raises, never silently masks
nothing, `:92-100`), and the **only** perturbation applied.

**Signal = ref−alt Pangolin delta** at the readout position (reference window vs its N-ablated copy).

---

## Choice 2 — Exon-7 readout position (`readout_index`, `gpu_scores.py:105`; `DEFAULT_READOUT_1B`)

**Definition (frozen).** The induced splice signal is read at the exon-7 5′ donor, which Pangolin
labels at the **last exonic base** = `EXON7_LAST_BASE_1B = 70076574` (review round-1 finding 1 — the
prior 70076575 was an off-by-one; Pangolin indexes donor usage at the last exonic, not first intronic,
base). The readout is fixed a priori to this position; it is **not** selected to maximize signal.
`readout_index` maps a genomic position to the `2*d+1` output array index and **bounds-checks** ±d of
window center (`:107-109`). The auxiliary `pangolin_max_gain` and its offset are **diagnostic-only** and
may NOT replace this primary donor endpoint after inspection.

---

## Choice 3 — Per-mer off-target burden statistic (`pangolin_off_target_burden`, `gpu_scores.py:261`)

**Definition (frozen).** For a given ASO (mer) and its set of near-complement off-target sites:

- Each site is authenticated before scoring: strand ∈ {+,−}; the declared `target_5to3` must match
  the strand-oriented reference (`:304-311`); and the ASO must bind by **hybridization** (reverse-
  complement), rejected only if a *strictly* closer **direct-ASO** match exists (`rc_distance >
  direct_distance` → reject; ties proceed as hybridization — the palindrome/near-symmetric case,
  `:312-321`).
- Each site's footprint is N-ablated (Choice 1) and scored; the site contribution is the **directional**
  `max(0, max_gain) + max(0, −max_loss)` (review round-1 finding 3 — the prior `|max_gain|+|max_loss|`
  double-counted same-sign extrema, e.g. counted a positive `max_loss` as a loss). Only a genuine
  positive splice-gain excursion plus a genuine negative splice-loss excursion contribute.
- **Burden = the SUM of the directional per-site contribution over all supplied off-target sites.** This
  is a **TOTAL** burden (conflates site prevalence with per-site severity), NOT a mean site effect; it
  varies by mer with the number of near-complement sites, by design. `n_sites` is a canonical-manifest
  count (equality-enforced against manifest rows), **not** a denominator.

**Honesty guards (frozen).**
- **No post-hoc threshold** on per-site deltas — every scanned site contributes; no cutoff is applied.
- **Authenticated empty states are distinct.** An authenticated scan artifact with `complete=True` and
  exact-integer `off_target_hits=0` proves a genuine complete-zero site set. `run_sweep` emits a terminal
  completed H-off artifact with `pangolin_off_target_incomplete=False`, `n_sites=0`, and burden `0.0`.
  An authenticated below-minimum-length scan has `complete=False`, an empty manifest, and the fixed
  `length_below_minimum_informative_length` reason; `run_sweep` emits a terminal diagnostic with
  `pangolin_off_target_incomplete=True`, `n_sites=0`, and burden `None`. The low-level scorer by itself
  still returns incomplete/`None` for an unauthenticated empty site list. Other GPU-env, worker/fetch,
  duplicate/missing-result, or missing/non-finite-value failures also return incomplete/`None`; they do
  not manufacture a zero.
- **Checkpoint resume is strict (review round-1 finding 4).** Fresh-state initialization occurs only
  when `load_checkpoint` returns `None`; it initializes committed `site_index=0` and burden `0.0`. A
  signal before the first completed batch may therefore write a legitimate index-0 checkpoint. Once
  any checkpoint file exists, `site_index` must be present with `type(site_index) is int` (including at
  zero), `input_digest` must match, and `pangolin_off_target_burden` must be present with exact type
  `float`, finite, and ≥0. Missing fields, booleans, strings, floats for `site_index`, and all coercions
  are refused; resumed state never defaults to zero.
- `n_sites` is the manifest site count (equality-enforced), NOT a denominator; burden is a **sum**.
  The per-mer site set is fixed by the frozen near-complement scan (max-edit ≤ 5) before the burden is read.
- **Pangolin-ONLY.** Burden uses no IntaRNA/SPLAM/SpliceAI signal — no ensemble stacking (project
  do-not). The 12-model aggregation inside is Pangolin's own internal model/head ensemble, not cross-scorer stacking.
- **Code + weight + version provenance (review round-1 finding 5).** The authenticated input digest
  pins the Pangolin package at executable-import-candidate content granularity:
  `pangolin_package_sha256` hashes the sorted `(package-relative path, SHA-256(file content))` pairs for
  every `.py`, `.pyc`, `.pyo`, `.pyd`, and platform extension-module suffix under the installed package;
  for `.pyc`/`.pyo`, “file content” is source-normalized bytecode with exactly the volatile 16-byte
  CPython 3.7+ header removed, while `.py` and extension candidates use full-content hashes. Thus an
  unchanged-source bytecode recompile does not change identity, but a changed bytecode body does.
  Provenance records the complete file list and count. The reviewed hash is committed as
  `PINNED_PANGOLIN_PACKAGE_SHA256` inside the independently authenticated `feng/src` tree and is the sole
  trust root: a fresh live measurement cannot authorize itself. The digest also pins the Pangolin package
  version string plus a SHA-256 of each of the 12 Pangolin weight files, so an artifact's identity changes
  if package code, executable cache/extension candidates, version, or weights change; missing package
  code or weights fail closed at scoring time.
- **Runtime-authenticated committed source.** The documented submit wrapper streams
  `git show <CODE commit>:feng/scripts/arm_e_hoff.sbatch` directly to `sbatch` and passes the raw-byte
  SHA-256 of that committed launcher blob. At task start the launcher hashes its executing script and
  exits nonzero before source staging or scoring on a mismatch. Each task materializes `feng/src` with
  `git archive <CODE commit>`, derives a manifest from that commit's Git object store, and compares the
  materialized bytes and complete file set with every committed blob. Permission bits are not treated
  as an authentication mechanism.
- **Point-of-use verification.** A stdlib-only bootstrap re-hashes every materialized `feng/src` entry
  immediately before importing the parent H-off module. Every Pangolin GPU worker is launched through
  that same bootstrap and independently repeats the complete-tree check in its own process before it
  imports `feng.titration.gpu_scores`. The parent provenance probe and `_worker_main()` entry each re-hash
  the installed Pangolin package's complete executable-candidate set and compare it directly with the
  committed pin before importing NumPy, Torch, or Pangolin. Verification additionally requires the
  Pangolin package and importable submodules to resolve through `SourceFileLoader`; sourceless, extension,
  absent, or other non-source loaders are refused. Any missing, changed, extra, malformed, or
  unauthenticated executable entry exits nonzero without scoring, including on resume/re-entry.
- **Digest-bound identity versus execution provenance.** `source_tree_sha256` is SHA-256 over sorted
  `(repo-relative path, SHA-256(raw committed blob bytes))` pairs. `input_digest` covers the visible
  `inputs` object, including this complete-tree identity, `executed_git_commit`,
  the verified-live `pangolin_package_sha256`, `pangolin_package_pinned_sha256`, and executable-candidate
  manifest, `DEFAULT_DISTANCE=50`, `PANGOLIN_CONTEXT=5000`, and the biological/model inputs listed above.
  It does **not** cover run identifiers: `execution_provenance` separately records `run_id`, Slurm array
  job id, array task id, and `launcher_blob_sha256`, identifying the particular execution while allowing
  identical inputs and code to reproduce an identical `input_digest`.

**Threat model.** `input_digest` authenticates the full `feng/src` tree plus the installed Pangolin
package identity (executable import-candidate content hashes), verified against its committed pin at task
start; this defends against code drift and at-rest substitution. It does **not** defend against an attacker
with live write access to the
compute account during execution (TOCTOU between hash and import, or a hostile interpreter startup /
pre-bootstrap import via inherited `PYTHONPATH`, `.pth`, or `sitecustomize`); those attacks are out of
scope for a single-user preregistration integrity claim.
- `STATE = "candidate"` is stamped on every output AND enforced by the artifact validator.

---

## Reading rule (freeze discipline)

The earlier 15 H-off outputs, including job 22343195, are SUPERSEDED and do not support this seal or
the registered endpoint. Path A requires a clean recompute from the
committed CODE snapshot described above. At the first start of the new timestamped run, each task holds
a per-length lock and moves any prior output, `.done`, and checkpoint into
`hoff/superseded/<run-id>/LXX/` before scoring. A destination with no live source is accepted as a
completed move after interruption; only after all moves finish is the preparation marker published by
atomic rename. The lock is retained through scoring, so duplicate launches cannot race or re-supersede
current-run output. Only outputs passing launcher, parent, and worker cryptographic checks count.

Sequence: (1) commit the complete code plus this DRAFT; (2) submit the launcher directly from that
commit and run the authenticated recompute without reading any burden value or ranking; (3) validate
structural completion and provenance; (4) in a separate commit, change this DRAFT to SEALED and stamp
the CODE commit SHA, source-tree SHA-256, launcher SHA-256, recompute job, manifest SHA-256, and date. The per-mer
ranking is NOT inspected before the month-9 Phase-9 / FREEZE step. No definition is refit after the
burden is read.

---

## Open items carried to FREEZE
- GPL-3.0 Pangolin vendoring flag (from commit `f0b0d48`) — resolve before any release use.
- Any non-blocking notes from the structural review are folded in above before seal.
