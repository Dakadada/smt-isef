# Public steric-intervention dataset catalog

**Snapshot:** 2026-08-28. Reconciled to `docs/dataset-audit.md`. “Archive-anchored” means a repository/paper declares raw reads. ENA metadata now confirms FASTQ URLs for all four candidates, but no FASTQ bytes or checksums and no exact intervention sequences have been verified.

## Locked scope and counts

The guaranteed deliverable is an **SMN2-centered case study**, public trans steric-blockade feasibility/negative-result benchmark, and immutable timestamped designs. The multi-locus upgrade (SMN2 -> +DMD/U7/dCas13d) activates only after raw reads are retrieved/checksum-verified for **>=3 independent families across >=2 loci**.

Current floor: **4 metadata-retrievable archive-anchored candidates but 0 byte/MD5-verified families; 2 loci (SMN2, DMD); about 5 unique endpoint interventions; 0 exact-sequence-and-FASTQ-byte-verified interventions.** Candidates are GSE307431, GSE167762/PRJNA705187, PRJNA1052769 individual-guide RNA-seq, and PRJEB13980. This cannot support a trainable general adapter, nested tuning, or terminal study.

## Eligibility and units

A trans label requires a sequence-programmed non-cleaving steric effector, matched control, quantitative endpoint splice reads, exact intervention/chemistry, dose/time/system, and sample map. Strict counting requires every relevant raw run, repository checksum, and manifest validation. Targeted reads retain a domain flag. RNase-H gapmers, active Cas13, DNA editing, small molecules, RBP knockdown, uncontrolled before/after clinical data, and guide-enrichment without splice-product reads are not primary labels.

Hierarchy: **study family -> biological sample -> oligo/intervention -> event**. Aliases and campaigns sharing controls/samples/lab/constructs form one family. Chemistry variants are interventions. Events never create intervention replication.

## Archive-anchored candidates

| Family | Locus/intervention | Evidence and unresolved fields | Locked role |
|---|---|---|---|
| **GSE307431 / PRJNA1321660** | SMN2 ISS-N1; F18OMe, F18MOE, F20PMO (3 implementations) | ENA metadata returns 20/20 RNA-seq runs with `fastq_ftp`. No FASTQ downloaded or byte/MD5 verified; exact strings remain unverified. | frozen SMN2 case-study anchor; no outcome-driven tuning |
| **GSE167762 / SRP308402 / PRJNA705187** | SMN2 ISS-N1 ASO1; VPA arm is not pure steric | ENA metadata returns 28/28 runs with `fastq_ftp` across RNA-Seq, ChIP-Seq, and OTHER; the about 12-run RNA-seq subset is the endpoint. No FASTQ downloaded or byte/MD5 verified; sequence and dose/time remain unresolved. Repeats nusinersen-like molecule/site; no new unique intervention. | SMN2 replication/feasibility; POINT-seq/ChIP auxiliary |
| **PRJNA1052769** | SMN2 dCas13d D[333,357], one individual-guide intervention | ENA metadata returns 16/16 runs with `fastq_ftp`; RNA-Seq is the about four-run individual-guide endpoint and AMPLICON is the pooled screen. No FASTQ downloaded or byte/MD5 verified; spacer remains unverified. The 1,937-guide enrichment screen is not endpoint labels. | conditional dCas13d upgrade evidence |
| **PRJEB13980** | DMD exon 53 U7snRNA-E53 | ENA metadata returns 12/12 RNA-seq runs with `fastq_ftp`. No FASTQ downloaded or byte/MD5 verified; cassette and arms/replicates remain unresolved; no RNA-sequenced U7-scramble and vector confounding. | conditional second-locus/U7 upgrade, sensitivity flag |

## Contextual/ineligible

| Record | Disposition |
|---|---|
| GSE221900 clinical nusinersen CSF | uncontrolled/confounded longitudinal context only |
| PMC8395096 Ottesen 2021 | supplementary tables, no verified raw accession; not GSE307431 |
| Holgersen 2021 | 81 SBOs/105 experiments, but proprietary/conditional raw data |
| STMN2 Baughn; UNC13A Keuss | biological precedent; no verified public raw intervention endpoint reads |
| DMD PMO beyond PRJEB13980 | no additional qualifying public endpoint series verified |
| GSE202178, GSE235673, cis MPSAs | ChIP, DNA editing, or cis perturbation; excluded |

GTEx/sQTL, ClinVar, MPSA, GSE121569, PRJEB42763/GSE126543, and eCLIP/probing supply cis/context/frozen features only. They never increase trans counts. eCLIP/probing is an ablation feature, not a model vote.

## Occupancy requirement

For every intervention seek the one proxy in `docs/architecture.md`: matched endpoint-time fraction of target RNA occupied via quantitative target-site protection or pull-down/CLIP with input/spike-in normalization. Record accession, assay, molecule, system, dose, time, estimate, uncertainty, transfer. None of the four candidates currently has this documented. Direction-only evaluation (sign of treated-minus-control delta-PSI, using directional balanced accuracy/concordance) is therefore the primary feasibility endpoint. Restrict each prediction to theta grid `{0,.25,.50,.75,1}` and require a robust direction; magnitude is a conditional secondary endpoint, so exclude MAE, Spearman, calibration, occupancy, and mechanism claims until matched theta is documented.

## Verification and activation

In a networked session on 2026-08-28, ENA REST `filereport` with `result=read_run` was queried for metadata only; **no FASTQ was downloaded**. It returned PRJNA1321660: 20 runs, all 20 with `fastq_ftp`, RNA-Seq; PRJNA705187: 28 runs, all with `fastq_ftp`, RNA-Seq + ChIP-Seq + OTHER (the RNA-seq subset is the endpoint); PRJNA1052769: 16 runs, all with `fastq_ftp`, RNA-Seq (individual-guide endpoint) + AMPLICON (pooled screen); and PRJEB13980: 12 runs, all with `fastq_ftp`, RNA-Seq. Thus archive retrievability is confirmed because FASTQ URLs exist, but **0 runs are byte/MD5-verified**, and exact intervention sequences remain to transcribe. This preserves the existing strict byte-verification and sequence caveats; see `docs/dataset-audit.md` section 7.

Export run table; retrieve every endpoint FASTQ; match repository MD5 and record local SHA-256; resolve controls/replicates; independently verify sequence/chemistry, dose, time, system, strand, event, license; fingerprint overlap; publish signed manifest. One fetched run does not verify a family.

Then apply `docs/benchmark-protocol.md`: N=2 paired feasibility; N=3 across 2 loci mechanical LOSO feasibility; N=4 across 2 loci/5 interventions descriptive stress test only. Higher floors govern correlation, calibration, cluster intervals, and terminal claims.

Search GEO/SRA/ENA/BioStudies using steric/splice terms plus SMN2, DMD, STMN2, UNC13A. Version queries, dates, results, exclusions, aliases, and data statements. New data activate the upgrade only through byte verification.

## Novelty and ceiling

The corrected novelty/FTO source is `docs/prior-art-defense.md`. The defensible contribution is a uniformly reprocessed, chemistry-aware, leakage-resistant public benchmark and search/manifest protocol. Known locus/effector work is prior art. Prospective packets are computational wagers; molecule, mechanism, rescue, safety, therapeutic, patentability, and superiority claims require future independent validation (and qualified counsel for IP/FTO).
