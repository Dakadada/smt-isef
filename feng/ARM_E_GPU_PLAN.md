# Arm E — GPU-intensive enhancement plan (feasibility, pre-results)

**Status: FEASIBILITY / PLAN ONLY. Nothing implemented. Written 2026-09-06, deliberately
BEFORE examining Arm E scan results (freeze-safe: no plan fitted to the outcome).**

Scope decided with user (2026-09-06): plan BOTH GPU axes, phased — specificity first,
MD second; feasibility depth (go/no-go), not a spec; both axes are **additive** `candidate`
score defs preregistered at the month-8 FREEZE tag (do NOT touch the existing frozen Arm E
therapeutic-window endpoint or existing scores); budget may go big (H100/H200 preemptable).

Reframe that motivates this: on-target axes lean *longer-acts-harder* (splice saturation
ambiguous, per commit `f0b0d48`), so the live "shorter-mer-is-better" hypothesis now rests
almost entirely on the **specificity axis**. GPU intensity is worth spending there.

---

## GPU inventory (verified on this ORCD/Slurm allocation, 2026-09-06)

| Partition | GPU types | Walltime | Notes |
|---|---|---|---|
| `mit_normal_gpu` | H100×4, **H200×8**, L40S×4 | 6 h | **QOS cap 2 GPU/user** |
| `mit_preemptable` | A100 (2–8), H100 (2–8), **H200×8 (31 nodes)**, L40S, A40, L4, RTX-Pro-6000 | 2 day | effectively unlimited GPU **iff** checkpoint+requeue |

=> "Think bigger" is real: H100/H200 are plentiful on `mit_preemptable`. Heavy sweeps belong
there, using the existing SIGTERM@60 checkpoint+requeue discipline
(`feng/scripts/templates/gpu_preemptable.sbatch`). Respect shared-quota courtesy.

## Already built (key finding — specificity is a scale-up, not greenfield)

- `feng/src/feng/titration/gpu_scores.py` = Pangolin axis H:
  - **H-on** (on-target potency): ASO footprint modeled as steric N-masking of its target
    window; Pangolin exon-7 donor usage gain = predicted inclusion shift. Answers the
    saturation question.
  - **H-off** (off-target splicing burden): same delta-splice prediction at each
    near-complement off-target site, aggregated to a per-mer burden. Per-site function is
    built + wired.
- GPU env + template + checkpoint infra exist (`feng-gpu-env`, cuda12.9/cudnn,
  `gpu_preemptable.sbatch`, `feng.gpu.checkpoint`, `check_gpu_env` gating that returns
  `value=None` with a reason rather than fabricating a number).
- Pending job `22054047` IS the H-off genome sweep — currently throttled to 1×L40S, `%2`.
- **MD axis: not built at all.** True greenfield.

---

## Phase 1 — Specificity (off-target burden sweep) · **GO**

- **Runs on our GPUs:** yes. Re-target the existing H-off sweep from 1×L40S to
  H100/H200 preemptable with wider array parallelism. Infra already requeue-safe.
- **Cost:** cheap. Pangolin = one forward pass per ~10 kb window, batched thousands/sec on
  H100. Near-complement sites (max-edit ≤5) across the ladder × transcriptome ≈ low-hundreds
  of GPU-hours thorough; trivial on h200:8 preemptable.
- **Freeze fit:** additive `candidate` score def. **The masking scheme + exon-7 readout
  position are flagged `candidate` modeling choices in the code — they MUST pass a structural
  review + prereg-seal BEFORE the sweep result counts.** That review is the real gate, not
  compute.
- **Serves:** the live falsifiable claim. "shorter-better" rests here.
- **Risk:** low compute. Scientific risk = defining the per-mer off-target burden statistic
  honestly (no post-hoc threshold picking) and avoiding the ensemble trap (do NOT stack
  Pangolin + IntaRNA + SPLAM without a defined causal role + locked ablation — project do-not).

**Phase-1 steps (when authorized):**
1. Structural (adversarial) review of the masking scheme + exon-7 readout choices → prereg-seal.
2. Prereg the per-mer off-target burden statistic (definition frozen before the sweep is read).
3. Re-run the H-off genome sweep at H100/H200 preemptable scale, wider parallelism, checkpointed.
4. Validate sentinels/manifest before accepting; do not read burden result until prereg is sealed.

## Phase 2 — MD (duplex free-energy / kinetics) · **GO-WITH-HEAVY-CAVEATS, defer**

- **Runs on our GPUs:** yes hardware-wise (H100/H200), but the **stack is not built** — needs
  OpenMM + force fields, incl. modified-nucleotide params (nusinersen is 2'-MOE / phosphoro-
  thioate, nonstandard), in a new env.
- **Cost:** heavy. Explicit-solvent duplex MD ~50–100k atoms × ~15 lengths × replicas ×
  100s ns = hundreds-to-low-thousands of GPU-hours. Feasible on preemptable, not free.
- **Risk (the real problem):** tiny ladder n + MD convergence = high overfit /
  unfalsifiable-noise risk; modified-nt force fields are shaky; marginal payoff over IntaRNA's
  thermodynamic point estimate unless framed as a preregistered free-energy ablation.
- **Recommendation:** do NOT start before Phase 1 lands and reads out. If pursued, scope as
  enhanced-sampling ΔG of duplex melting (umbrella / MBAR), NOT brute MD, with a locked
  ablation vs IntaRNA.

---

## Bottom line

Bigger-GPU Arm E is worth it — but the lever is **Phase 1 specificity, a compute-scale-up of
already-built code, not new science**. Do Phase 1 first (review → prereg-seal masking/readout
→ rerun H-off at H100/H200 scale → prereg the burden statistic). Phase 2 MD is an optional
deep follow-up: heaviest cost, weakest small-n payoff — defer past the Phase-1 readout. Both
additive, preregistered at FREEZE, endpoint untouched.

## Constraints carried (do not violate)

- Additive `candidate` defs only; no change to the frozen Arm E endpoint; prereg at FREEZE.
- Do NOT read sealed/terminal or currently-running Arm E results to shape this plan.
- Clean-room; honesty architecture; no chemistry-pooling; no ensemble without locked ablation.
- Resolve the GPL-3.0 Pangolin vendoring flag (noted in `f0b0d48`) before any release use.
