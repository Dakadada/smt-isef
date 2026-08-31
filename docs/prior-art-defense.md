# Prior-art defense: compute-native design of steric splice modulators

**Purpose and claim boundary.** This document states what this project can and cannot fairly claim. It is written for an ISEF Scientific Review Committee (SRC), a skeptical scientific reviewer, or a patent examiner. The project is not presented as the invention of splice-switching oligonucleotides, programmable RNA occlusion, STMN2 or UNC13A rescue, RNA-accessibility prediction, target-site discovery, or sequence-based splicing prediction. Its proposed durable contribution is a rigorously evaluated *intervention-prediction framework* centered on a public, chemistry-aware benchmark of trans-acting steric interventions. Cis/trans separation, genuinely held-out tests, prospective immutable designs, and calibrated uncertainty under domain shift are required safeguards, not separate inventions. Any occupancy/footprint representation is a computational hypothesis whose biological utility would remain contingent on future independent experimental confirmation.

## 1. The closest prior art already establishes the biological premise

### STMN2: ASOs, programmable steric replacement, and functional rescue are prior art

Baughn et al. showed that loss of nuclear TDP-43 exposes a cryptic splice acceptor/polyadenylation program in human **STMN2**. They corrected processing with splice-modulating ASOs, and ASO restoration of stathmin-2 rescued axonal regeneration after axotomy in TDP-43-deficient human motor neurons. They also showed that a nuclease-dead RfxCas13d effector (dCasRx), positioned at the relevant pre-mRNA region, can substitute for the steric-blocking function of TDP-43 and restore normal STMN2 processing. This is direct prior art against any broad claim to “using an ASO on STMN2,” “using dCas13 to block a splice element,” or “discovering that steric occupancy can correct STMN2.” The conservative reading is important: the paper supports functional axonal-regeneration rescue for the ASO intervention and molecular splicing correction for dCasRx; this document does **not** assert that a separate dCasRx axonal-regeneration experiment was performed. [Baughn et al., 2023, *Science* 379:1140–1149, PMID 36927019, DOI 10.1126/science.abq5622; verified via PubMed/PMC and the journal metadata.]

### dCas13d splicing screens: dense proximal/distal mapping and SMN2 comparison are prior art

Recinos et al. performed a CRISPR-dCas13d screen of proximal and distal splicing-regulatory elements in an SMN2 reporter, using a 1,937-guide library and fluorescence sorting. The study explicitly framed dCas13d/gRNA action as steric hindrance, validated individual guides, and compared the platform with ASOs directed at known SMN2 regulatory elements, including ISS-N1, the target of nusinersen. It therefore anticipates broad claims to screening dCas13 guide positions, discovering distal regulatory sites, or comparing dCas13 steric modulation with an ASO at SMN2. Raw screen/RNA-seq data are deposited under **PRJNA1052769**, which makes this study a key benchmark source rather than an unoccupied niche. [Recinos et al., 2024, *Nature Communications* 15:3839, PMID 38714659, PMCID PMC11076525, DOI 10.1038/s41467-024-47140-8; accession verified in the article's Data Availability statement.]

### UNC13A: splice-switching rescue and candidate sequences are prior art

Keuss et al. showed that TDP-43 loss drives UNC13A cryptic splicing and synaptic dysfunction in human neuron models. Their M18-17 and M21-4 splice-switching ASOs restored UNC13A protein to different degrees and rescued synaptic transmission. Targeting UNC13A, suppressing its cryptic exon, or proposing an UNC13A splice-switching ASO is thus not novel. [Keuss et al., 2024 preprint, PMID 38979232, PMCID PMC11230273, DOI 10.1101/2024.06.20.599684; verified via PubMed/PMC.]

At least one associated filing is identifiable: **US20250051775A1, “Therapeutics for the treatment of neurodegenerative disorders,”** issued as **US12522826B2 on 2026-01-13**, with priority **2020-12-09**; its inventors are **Pietro Fratta, Anna Leigh Brown, Oscar Wilkins, Matthew Keuss, Gemma Ward, and Andrew Hill**, and its assignees are **UCL Business Ltd and US HHS**. Its exact claim scope and overlap with the paper's particular ASOs still require a professional patent-family review and counsel before any freedom-to-operate or patentability conclusion. Earlier UNC13A/TDP-43 filings also exist—for example **WO2022216759A1, “Compositions and methods for treating TDP-43 proteinopathy,”** which discloses UNC13A-targeting ASOs and is assigned to **Maze Therapeutics**; its US family member is **US20250011773A1**. The PCT record is listed as ceased, while national rights require a separate check. It is not attributed here to Keuss et al.; it is cited only to show a crowded claim landscape. Claim scope, national rights, and FTO implications remain **to verify with counsel**.

### Nusinersen/ISS-N1: the canonical steric splice-switching precedent

Singh and colleagues identified ISS-N1 near the 5′ end of SMN2 intron 7 as a potent intronic splicing silencer; ASO masking of the element promotes exon 7 inclusion. Subsequent Krainer-lab work optimized antisense modulation in this region, and the 18-mer 2′-O-methoxyethyl/phosphorothioate drug nusinersen targets intron-7 nucleotides 10–27. Thus, “find a silencer, bind it with a steric ASO, and restore exon inclusion” is a mature therapeutic design pattern, not this project's invention. [Singh et al., 2006, *Molecular and Cellular Biology* 26(4):1333–1346, PMID **16449646**, DOI **10.1128/MCB.26.4.1333-1346.2006**. The ASO-10-27 primary citation is Hua et al., 2008, *American Journal of Human Genetics* 82, PMID **18371932**, DOI **10.1016/j.ajhg.2008.01.014** (not the 2007 paper).]

## 2. Sequence models are strong prior art—but mostly answer a different causal question

SpliceAI predicts splice junctions and variant-induced changes from long primary-sequence context. [Jaganathan et al., 2019, *Cell* 176:535–548.e24, PMID 30661751, DOI 10.1016/j.cell.2018.12.015; verified via PubMed.] Pangolin predicts tissue-aware splice-site usage from DNA sequence; Pangolin reports improved splice-site prediction over SpliceAI on its benchmark. [Zeng & Li, 2022, *Genome Biology* 23:103, PMID 35449021, DOI 10.1186/s13059-022-02664-4; verified via PubMed.]

MMSplice uses modular neural components to predict variant effects on splicing. [Cheng et al., 2019, *Genome Biology* 20:48, PMID 30823901, DOI 10.1186/s13059-019-1653-z; verified via PubMed.] Massively parallel splicing assays (MPSAs)—including Vex-seq, MaPSy, MFASS, and saturation assays—have both trained and benchmarked splicing-effect predictors. A useful synthesis benchmark found that Pangolin and SpliceAI generally agreed best with saturation-MPSA measurements, while performance varied by exon, region, and task. [Smith & Kitzman, *Genome Biology* 24:294 (2023), PMID **38129864**, PMCID **PMC10734170**, DOI **10.1186/s13059-023-03144-z**.]

These methods are indispensable baselines. They are not, however, sufficient evidence that a model predicts the outcome of placing a chemically specified, finite-size trans-acting object on an unchanged pre-mRNA. A nucleotide substitution changes the cis sequence and can simultaneously alter splice-site strength, motifs, pairing, structure, and transcription-linked features. An ASO or dCas13 RNP instead occupies a region with chemistry-dependent affinity, kinetics, mismatch tolerance, and an excluded-volume footprint. Conflating those estimands is a testable modeling choice, not a biological identity.

## 3. Chemistry and mismatch tolerance prevent a sequence-only notion of “the target”

The Singh-lab study associated with **GSE307431** compares three ASOs aimed at ISS-N1—F18OMe (2′-O-methyl/phosphorothioate), F18MOE (2′-O-methoxyethyl/phosphorothioate), and F20PMO (morpholino)—and reports chemistry-dependent transcriptome perturbation. In particular, F18MOE induced POLR2H exon-2 skipping through imperfect pairing and steric blockade, with tolerance involving terminal mismatches and internal wobble pairs; increasing complementarity strengthened effects in validation constructs. This directly warns against treating sequence identity, predicted duplex energy, or a single “off-target count” as chemistry-invariant. [Ottesen et al., *NAR Molecular Medicine* 3(1):ugag002 (2026), PMID **41551726**, PMCID **PMC12805893**, DOI **10.1093/narmme/ugag002**; GEO **GSE307431**.]

Consequently, an honest intervention benchmark must encode at least oligonucleotide chemistry/backbone, guide or ASO length, concentration or expression regime when available, mismatch pattern, target accessibility, and effector footprint. Missing fields must remain missingness indicators—not be silently imputed as if all interventions were equivalent.

## 4. Explicit non-novelty statement

This project does **not** claim novelty for:

1. using an ASO on STMN2;
2. using dCas13/dCasRx to alter splicing;
3. targeting UNC13A or its TDP-43-dependent cryptic exon;
4. predicting RNA accessibility;
5. identifying an ASO-binding site; or
6. building another cis-splicing foundation model.

Nor does it claim that affinity, accessibility, motif content, distance to a splice site, or a pre-trained splice model is individually new. Each is a necessary comparator.

## 5. The defensible novelty is the cross-study benchmark

Novelty here is conditional: each point below is a contribution only if implemented exactly, released with auditable provenance, and supported by the stated held-out results. A proposal or architecture diagram alone does not establish it.

### 5.1 A cross-study, chemistry-aware public **STERIC-INTERVENTION** benchmark

The benchmark unit is an intervention, not a variant: `(pre-mRNA context, locus, cell/model system, effector modality, chemistry, targeting sequence, occupied interval, dose/expression, outcome, replicate)`. Public raw reads from eligible ASO and dCas13 steric-occlusion studies are reprocessed through one versioned workflow: adapter/QC policy, alignment or guide-count procedure, event definition, PSI or other outcome estimator, replicate aggregation, and uncertainty. Published summary values are retained for reconciliation, not mixed uncritically with recomputed endpoints.

The benchmark must publish inclusion criteria, exclusions, data licenses, accession-to-sample maps, genome/transcript versions, checksums, workflow/container versions, and a per-study audit trail. Cross-study generalization is evaluated with grouped splits by study, locus, chemistry, and effector—not random rows that leak neighboring guides or shared constructs. PRJNA1052769 is an anchor dataset; GSE307431 is relevant chemistry/off-target evidence. Baughn and Keuss data are included only where raw public data and metadata support uniform reprocessing. The search protocol and search manifest must also be published so that “no prior art found” is supported by a reproducible method rather than assertion.

The novelty claim is therefore not “we gathered data.” It is: **no identified prior work has yet been shown here to provide a public, uniformly reprocessed, cross-study benchmark whose labels are outcomes of steric trans interventions and whose schema preserves effector chemistry.** A formal systematic search could still defeat this claim, and the claim must be updated if overlapping work is found.

## 6. Required design and validation safeguards (not novelty claims)

### 6.1 Strict separation of cis-mutation learning from trans-occlusion learning

Two datasets and two tasks are kept explicit:

- **Cis task:** predict the effect of changing the RNA/DNA sequence.
- **Trans task:** predict the effect of occupying a region while the underlying pre-mRNA sequence is unchanged.

Cis-trained outputs from SpliceAI, Pangolin, MMSplice, or MPSA-trained models may be frozen features or baselines. They must not be relabeled as measured trans outcomes. Training, calibration, and test partitions for trans prediction use only steric-intervention labels. A pre-registered ablation tests whether cis-model features help, hurt, or merely proxy distance/motif density. This separation is required dataset hygiene because it prevents an apparently large sample size built from causally different perturbations; it is not an invention.

### 6.2 A constrained, falsifiable occupancy/footprint hypothesis

Before testing, the project must operationally define “occupancy,” “productive state,” “footprint,” and “steric effect,” including their measurable inputs, units or scales, and boundary conditions. It must constrain and freeze the footprint kernel before testing rather than tune its interpretation after observing held-out outcomes. The proposed model would predict an occupancy-weighted productive-state estimand rather than “binding” in the abstract. For candidate intervention (i), define an occupancy distribution (p_i(x)) over target positions/states, conditioned on chemistry, complementarity, accessibility, and dose/expression. Define a productive steric state (q_i(x)): the probability that occupancy at (x) blocks or permits the relevant splice-regulatory assembly. The predicted intervention effect is an aggregation such as

\[
\widehat{\Delta\Psi}_i = f\!\left(\sum_x p_i(x)\,q_i(x),\; F_i,\; C_i\right),
\]

where (F_i) is an excluded-volume footprint kernel and (C_i) is sequence/cellular context. The footprint is larger than the nominal hybridized bases for an RNP and modality-specific for an ASO; it masks or perturbs nearby protein/RNA interactions according to distance from the occupied interval. This is an explicit, ablatable representation—not a claim that the simplified kernel is a complete molecular simulation.

This computational hypothesis merits attention only if the pre-constrained representation improves **held-out intervention prediction** over strong nonphysical baselines under identical grouped splits. Required matched-capacity nonphysical baselines include: SpliceAI; Pangolin; MMSplice where applicable; duplex affinity/complementarity; RNA accessibility; motif disruption/overlap; distance and orientation to splice sites; chemistry and length; study/locus covariates; and flexible combinations of these features. Report continuous metrics (for example held-out MAE and rank correlation), classification metrics only when thresholds are predeclared, calibration, confidence intervals from grouping-aware bootstrap, and paired tests across held-out groups. The decisive ablations remove occupancy weighting, footprint, chemistry, or productive-state terms one at a time. Even computational improvement would not establish the proposed mechanism or a biological design rule without future independent experimental confirmation.

A gain on random row splits is not evidence. A gain confined to one locus, one study, or one modality is not a general design rule. The claim should be narrowed to the domains in which grouped, leakage-resistant results reproduce.

### 6.3 Immutable, timestamped **prospective** designs

After model selection and before any new outcome is observed, the project exports a signed, content-addressed design packet. It contains model and code commit hashes; container/data-manifest hashes; exact transcript/genome build; candidate sequences and coordinates; chemistry/effector assumptions; ranked predicted effects; uncertainty intervals; inclusion/exclusion rules; and a declaration that the file is frozen. The packet is timestamped in a durable third-party or public repository. Any later revision is a new version, never an overwrite.

Candidates used during training, feature engineering, threshold choice, or retrospective demonstration are not prospective. A future laboratory can test the frozen set blind to subsequent tuning. This is an anti-hindsight evidentiary safeguard, not scientific novelty. Until independent experimental confirmation, the output is a computational prospective prediction and an auditable experimental wager—not an experimentally validated molecule or biological design rule.

### 6.4 Quantified out-of-domain uncertainty

The model must state when a proposed `(chemistry, locus, effector)` combination is unlike training data. Evaluation therefore includes leave-one-chemistry-out, leave-one-locus-out, leave-one-study-out, and leave-one-effector-out tests where sample support permits. Predictions carry both outcome uncertainty and a domain-distance/OOD score. Calibration is reported separately in-domain and OOD; intervals should widen or the model should abstain when combinations are unsupported.

Useful implementations may include grouped conformal intervals, hierarchical variance components, ensembles, and explicit nearest-domain distances, but method choice is secondary to coverage testing. These are standard ML validation methods and are required here, not claimed as novelty. “Uncertainty” is not defensible if 90% intervals fail to approach 90% empirical coverage on held-out groups. For a truly unseen effector, the honest output may be “unsupported extrapolation,” not a ranked candidate.

## 7. What would falsify or materially narrow the claim

The central claim is weakened or defeated if:

- a prior public benchmark already provides uniformly reprocessed, chemistry-aware trans-occlusion outcomes with comparable grouped evaluation;
- the occupancy/footprint model does not outperform matched-capacity nonphysical baselines on held-out interventions;
- improvement vanishes under study-, locus-, chemistry-, or effector-held-out splits;
- results depend on leakage among overlapping guides, shared reporters, replicates, or near-identical sequences;
- uncertainty intervals are miscalibrated under domain shift;
- raw-read reprocessing cannot reproduce published directional effects; or
- the “prospective” set was inspected, revised, or selected after outcome data became available.

Negative results remain scientifically useful: the benchmark can show that current public data are insufficient to identify an occupancy/footprint advantage, or that simple accessibility/affinity/distance baselines are already adequate. The project should report that result rather than retroactively redefining success.

## 8. Biological ceiling without new wet-lab validation

No analysis of public data can establish that a newly proposed molecule is better in a disease-relevant organism, safe, deliverable, durable, or clinically useful. With no new wet lab, the strongest current deliverable would be: **a reproducible model that, if the pre-specified evaluation succeeds, performs better on specified held-out public interventions and generates an immutable, timestamped set of prospective predictions for later testing.** It would remain computational-only evidence, not an in-house-validated ASO/guide. Any biological claim about this project's candidates or proposed mechanism must await future independent experimental confirmation. Functional rescue established by Baughn, Keuss, nusinersen, or others validates their interventions and biological mechanisms; it does not validate this project's new candidates.

## 9. Patentability and freedom-to-operate limits

A computational representation of occupancy and footprint, by itself, is vulnerable as intellectual property: it may be viewed as an abstract model, an obvious aggregation of known biophysical features, or insufficiently enabled without demonstrated utility. Publication may also become prior art against later claims, while it does not create freedom to operate around existing composition-of-matter or method-of-treatment claims.

A computational representation alone is a weak IP position. Any potentially defensible IP position—still requiring counsel and a full prior-art/FTO search—would need a **specific non-obvious algorithmic step, enabling evidence, and a specifically disclosed candidate set**, all of which remain pending. Relevant evidence would need to include unexpected held-out computational improvement, a pre-specified chemistry/effector regime, exact candidates, and future independent experimental confirmation against appropriate comparators before any biological performance claim. Candidate composition, target interval, chemistry, dose, delivery, indication, and method claims may each encounter different prior art, including the STMN2, UNC13A, SMN2, and dCas13 work above. The identified UNC13A filings make broad target-level claims especially risky. No patentability or non-infringement conclusion is made here.

## 10. Conditional claim language

If the required results are obtained, a defensible summary is:

> We will release a public, uniformly reprocessed benchmark of chemically and mechanistically annotated steric splice interventions and its search protocol and manifest. We will keep cis mutations separate from trans occlusion and test a pre-constrained occupancy/footprint hypothesis against matched-capacity nonphysical baselines under grouped held-out evaluation. If those computational tests succeed, we will freeze and timestamp prospective candidates with calibrated OOD uncertainty. Any biological interpretation and candidate performance claim will await future independent experimental confirmation.

The following language is not supported:

> We invented STMN2/UNC13A splice rescue, invented dCas13 splice modulation, discovered the best ASO, or proved a clinically superior therapy.

## 11. Citation-verification ledger

“Verified” means that the identifier and the claim used above were checked against the named primary or authoritative record during preparation; it is not a replication of the experiment.

| Citation or record | Status | Verification basis / limitation |
|---|---|---|
| Baughn et al., 2023, PMID 36927019, DOI 10.1126/science.abq5622 | **Verified** | PubMed-linked PMC full text and Science metadata; distinguished ASO functional rescue from dCasRx molecular correction. |
| Recinos et al., 2024, PMID 38714659, DOI 10.1038/s41467-024-47140-8, PRJNA1052769 | **Verified** | PubMed, PMC full text, Nature Communications record, and article Data Availability statement. |
| Keuss et al., 2024, PMID 38979232, DOI 10.1101/2024.06.20.599684 | **Verified** | PubMed and PMC preprint record/full text. It is cited as a preprint, not silently upgraded to a journal article. |
| US20250051775A1; issued as US12522826B2 (2026-01-13); priority 2020-12-09; inventors Fratta, Brown, Wilkins, Keuss, Ward, Hill; assignees UCL Business Ltd + US HHS | **Resolved bibliographic record; claim scope/FTO needs counsel** | Issuance, priority, complete inventor list, and assignees resolved; exact claim overlap and FTO remain to verify with counsel. |
| WO2022216759A1; Maze Therapeutics; US family member US20250011773A1 | **Resolved bibliographic record; national rights need separate check** | PCT record listed ceased; national rights, claim construction, and FTO remain to verify with counsel. |
| Singh et al., 2006 ISS-N1 discovery, *Molecular and Cellular Biology* 26(4):1333–1346, PMID 16449646, DOI 10.1128/MCB.26.4.1333-1346.2006 | **Resolved** | Primary ISS-N1 citation. |
| Singh, Howell, Androphy & Singh, 2017, DOI 10.1038/gt.2017.34 | **Verified** | PMC/Gene Therapy record; used as a historical review of ISS-N1 and nusinersen development. |
| Hua et al., 2008, *American Journal of Human Genetics* 82, PMID 18371932, DOI 10.1016/j.ajhg.2008.01.014 | **Resolved** | Primary ASO-10-27 citation; not the 2007 paper. |
| Jaganathan et al., 2019, PMID 30661751, DOI 10.1016/j.cell.2018.12.015 | **Verified** | PubMed/Cell bibliographic record. |
| Zeng & Li, 2022, PMID 35449021, DOI 10.1186/s13059-022-02664-4 | **Verified** | PubMed/Genome Biology record. |
| Cheng et al., 2019 MMSplice, PMID 30823901, DOI 10.1186/s13059-019-1653-z | **Verified** | PubMed/Genome Biology record. |
| Smith & Kitzman, *Genome Biology* 24:294 (2023), PMID 38129864, PMCID PMC10734170, DOI 10.1186/s13059-023-03144-z | **Resolved** | MPSA predictor benchmark. |
| Ottesen et al., *NAR Molecular Medicine* 3(1):ugag002 (2026), PMID 41551726, PMCID PMC12805893, DOI 10.1093/narmme/ugag002; GSE307431 | **Resolved** | F18MOE/POLR2H study and associated GEO accession. |

All citation records specifically resolved for this review are reflected above. Patent claim scope, national rights, patentability, and FTO remain matters for separate counsel review.
