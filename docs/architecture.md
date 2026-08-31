# Locked Task-1 architecture

**Lock date:** 2026-08-28. Amendments are append-only and frozen before affected outcomes open. Nothing under `pre/` is part of this prospective lock.

## Scope and claim boundary

The guaranteed deliverable is an **SMN2-centered case study**, a public benchmark assembled for sequence-programmed **trans steric-blockade splicing interventions**, feasibility/negative-result reporting, and immutable timestamped prospective designs. A multi-locus benchmark activates only after a parallel verification retrieves and checksum-verifies raw reads for **>=3 independent study families across >=2 loci**. Candidate upgrade domains are SMN2 plus DMD/U7 and dCas13d; STMN2/UNC13A remain leads unless their intervention reads qualify. Novelty and prior-art claims are limited to the locked analysis in `docs/prior-art-defense.md`.

`docs/dataset-audit.md` establishes the current floor: **0 byte-verified families**, about **4 archive-anchored candidates**, **2 loci (SMN2, DMD)**, about **5 unique endpoint interventions**, and zero exact-sequence-and-FASTQ-byte-verified interventions. ENA metadata confirms `fastq_ftp` URLs for all 76 archive runs across the four candidates, but no FASTQ was downloaded and 0 runs are byte/MD5-verified. No general predictor, biological rule, mechanism, therapeutic benefit, patentability, or FTO claim is currently permitted. The benchmark is the durable contribution; biological claims require future independent validation. Novelty/IP language is governed by corrected `docs/prior-art-defense.md`.

## Evaluated quantity and frozen model

The primary feasibility endpoint is **direction only**: the sign of cross-fitted predicted treated-minus-control delta-PSI for one preregistered event per intervention, evaluated by directional balanced accuracy and directional concordance. The model still emits a treated-minus-control delta-PSI curve,

\[
\widehat{\Delta\Psi}_{i}=100\{g_{-s(i)}(P_{i,T})-g_{-s(i)}(P_{i,C})\}.
\]

where `s(i)` is the held-out study family and the monotone bounded link `g_-s:[0,1]->[0,1]` is fitted only on outer-training families and frozen. Output is clipped only to `[-100,100]`. The one primary model is an **equilibrium full-transcript model plus a constrained PSI link**. Its exact class map is

\[
P_{i,c}=(1-\theta_{i,c})P_{U,i,c}+\theta_{i,c}P_{B,i,c};\quad
\widehat{\Delta\Psi}_i=100[g(P_{i,T})-g(P_{i,C})].
\]

The primary direction must be unchanged over the frozen theta grid. Delta-PSI **magnitude** evaluation (including MAE, rank correlation, and calibration slope) is a conditional secondary endpoint that activates only for interventions with documented matched occupancy evidence under the contract below; none of the current candidates qualifies. The co-transcriptional model is secondary. eCLIP/probing is only an ablation/feature variant of either model, never a third voter. No majority vote is used.

The PSI link is a monotone logistic spline with at most **3 interior knots**, global intercept/slope, and chemistry/study random intercepts with zero-centered regularizing priors. Locus-ID, oligo-ID, study-by-chemistry, and unrestricted categorical interactions are forbidden. Effective degrees of freedom and priors are frozen per fold. Matched-capacity nonphysical baselines receive the same metadata, coefficients/knots (within one degree of freedom), nesting, and tuning budget.

## Occupancy (`theta`) contract

The one accepted external proxy is the **fraction of target RNA molecules occupied at endpoint sampling**, measured in the matched cell/system, dose, chemistry/effector, and time by quantitative target-site protection or pull-down/CLIP with input and spike-in normalization. Its public raw assay or preregistered matched companion assay must provide accession, molecule identity, calibration curve, and uncertainty. Nominal dose, predicted affinity, endpoint PSI, and unrelated-cell binding are not occupancy measurements.

Every manifest contains proxy accession, assay, time, dose, matched system, estimate, uncertainty, and transfer flag. If absent/unmatched, output is downgraded to **theta-conditional effects on `{0,.25,.50,.75,1}` plus robust direction** (unchanged over the grid). No magnitude, calibration, correlation, effector-mechanism, or occupancy claim is allowed for that intervention.

| Claimed parameter | Observable evidence | Prior/constraint | If missing |
|---|---|---|---|
| observed delta-PSI | treated/control junction counts and replicates | beta-binomial; frozen event | effect with uncertainty if coverage passes |
| `theta` | matched protection or pull-down/CLIP above | `[0,1]`, measurement error | no point value; grid only |
| `P_U`,`P_B` | sequence/structure inputs, frozen state definitions | physical bounds/partition model | score only; no observed-mechanism claim |
| PSI link | training-family trans outcomes | monotone, <=3 knots, regularized groups | no magnitude transfer outside support |
| effector effect | repeated interventions at comparable sites/loci | zero-centered hierarchy | no mechanism claim |
| calibration/generalization | eligible cross-fitted interventions/loci | frozen model/grouping | only above protocol floors |

## Pillars

- **P3a:** strand-correct full-pre-mRNA off-target enumeration across mismatch/bulge/indel and energy bands.
- **P3b (flagship):** immutable timestamped SMN2 designs with sequence, chemistry, coordinates, theta assumption/grid, predictions, controls, off-target dossier, hashes, and future assay. Multi-locus designs wait for the upgrade gate.
- **P1 effector-class variance: CONDITIONAL-ON-DATA.** SMN2-only data cannot separate class from study/site/chemistry. Activate only under the multi-locus upgrade with repeated interventions in >=2 effector classes at >=2 loci.
- **P4 span robustness/flip pairs: CONDITIONAL-ON-DATA.** SMN2-only data cannot identify it. Activate only under the upgrade with matched concordant/flip pairs at >=2 loci; events never substitute for pairs.
- **P5:** compute-only minigene tooling; public minigene results retain a domain label.
- **P2 C2B34:** shelved under the archived causal-timing objection.

## Baseline contract

SpliceAI/Pangolin act on the **unchanged reference sequence as frozen feature extractors**, not direct trans predictors. Report untreated splice propensity and a target-region importance vector from a frozen in-silico occlusion transform: replace each occupied base one at a time with each alternative, average signed event-score changes, and aggregate with the frozen footprint kernel. This is a sensitivity feature, not an assertion that substitution equals binding.

Matched-capacity baselines are: untreated propensity; propensity plus target-region importance; affinity/complementarity, accessibility, splice-element distance/orientation, chemistry, dose, length, and k-mers; and their capacity-matched stack. Physics has value only if it beats the best outer-fold-selected matched-capacity baseline. Ablations remove bound-state physics, occupancy, footprint, co-transcriptional terms, or eCLIP/probing.

## Locked compute and infrastructure

MIT normal allocation is **2 GPUs; no H200s**. Favor frozen embeddings, then parameter-efficient fine-tuning. Training from scratch is disallowed unless a locked pilot shows material benefit and fits the envelope. Larger jobs use **PREEMPTABLE** with checkpoints every **5-15 minutes** and a queue/requeue monitor.

| Stage | CPU/RAM | GPU-hours | Persistent/scratch | Approach |
|---|---:|---:|---:|---|
| read reanalysis | 16-32 CPU, 64-128 GB; roughly 500-1,500 CPU-hours | 0-20 | 0.5-1 TB persistent / 2-4 TB scratch | download, checksum/manifest, align, QC, and splice-junction quantify; screen 76 archive runs and process the ~48 RNA-seq endpoint runs (20 + ~12 + ~4 + 12) |
| cis model | 16 CPU, 64 GB | 100-300/2 GPUs | 100 GB / 0.5 TB | frozen first; PEFT only if pilot wins |
| partition functions | 32 CPU, 128 GB | 0-80 | 100 GB / 1 TB | shard transcript/intervention/theta |
| physics outer folds | 16 CPU, 64-128 GB | 100-250/2 GPUs | 200 GB / 1 TB | fixed seeds; fold/shard jobs |

Checkpoints go to a temporary path, are fsynced and atomically renamed, with code/container/data hashes, fold, shard, RNG/optimizer state, completed IDs, and output checksums. Resume validates hashes and a deterministic probe batch. Idempotent shards are preemption-safe. The monitor polls jobs/heartbeats, requeues stale/preempted shards, caps retries at 3, quarantines partial/checksum-invalid output, and logs failures.

Raw-read work is predominantly CPU alignment, not GPU work. Run each family as an independently checkpointed download/alignment/QC/junction shard on PREEMPTABLE capacity, persist manifests and completed-run outputs, and requeue only incomplete family shards. The 76-run archive envelope comprises PRJNA1321660 (20 RNA-seq), PRJNA705187 (28 total, about 12 endpoint RNA-seq), PRJNA1052769 (16 total, about 4 individual-guide endpoint RNA-seq), and PRJEB13980 (12 RNA-seq). If full genomic alignment cannot fit the allocation or scratch quota, use a locked reduced-processing path based on splice-aware pseudoalignment and/or targeted junction counting for the preregistered events, retain the approximation label, and report concordance against a full-alignment pilot family.

If full-transcript partitioning exceeds **30 minutes or 32 GB RAM per intervention-state** in the pilot, fall back to 2-kb event windows plus annotated competing elements; ensemble 0.5/1/2-kb windows, report sensitivity, and abstain when direction changes. Label this an approximation.

## Reproducibility

Raw objects receive accession, retrieval time, URI, byte size, and SHA-256. References, containers, commands, sample sheets, eligibility, theta evidence, and predictions are versioned. Outcomes remain inaccessible until predictions/hashes deposit. Negative results persist; amendments apply only to later untouched data.
