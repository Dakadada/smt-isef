# Full-pre-mRNA off-target enumerator (P3a)

This module enumerates hybridization candidates at **every offset of every complete
expressed pre-mRNA sequence supplied**. It uses deterministic Levenshtein dynamic
programming over target lengths around the oligo length, recovering mismatches,
RNA/oligo bulges, indels, and configurable edit distances above four. It does not
use the retired signed-PSI polarity rule: the estimand is recall of validated
binding/splicing events. Effect sign is `unknown` unless overlapping eCLIP or
perturbation evidence supplies a sign; conflicting signed evidence also yields
`unknown`.

Chemistry-specific coarse enumeration priors are represented in `chemistry.py` for
2'-MOE, 2'-OMe, PMO, and LNA. They are transparent ranking approximations rather
than validated duplex thermodynamics. Any other chemistry uses a conservative
fallback, is penalized, and emits `chemistry_parameters_known=false`. Scores are
recall-oriented ranking scores; `confidence` is a bounded heuristic and is not
claimed as clinical/statistical calibration. The benchmark reports recall and
reports precision only when the prediction universe is completely labeled.

## Offline test

From the repository root (no downloads or alignment):

```bash
PYTHONPATH=src python3.12 -m unittest tests.test_offtarget -v
```

## Later networked GENCODE run

In a separate network-enabled, provenance-recorded session, download and checksum
a pinned GENCODE genome FASTA and matching comprehensive GTF. Construct one
strand-correct JSONL record per expressed transcript with `transcript_id`,
`gene_id`, the **full pre-mRNA genomic span including all introns and exons**,
`expression_tpm`, and `nuclear_fraction`. Do not provide spliced cDNA. Intersect a
pinned eCLIP BED/metadata release with those transcript coordinates and convert it
to optional evidence JSONL. Then run:

```bash
PYTHONPATH=src python3 -m offtarget \
  --config config/offtarget.example.json \
  --transcripts run/gencode_full_premrna.jsonl \
  --evidence run/eclip_evidence.jsonl > run/offtargets.tsv
```

Record input hashes and tool/interpreter versions. For large transcriptomes, shard
by transcript without prefiltering sequence offsets, merge all shard outputs, and
rerank globally. The exhaustive reference implementation favors recall over speed.

## Held-out control

POLR2H exon 2 / F18MOE (GSE307431) is a held-out known-positive check. Historically
this control forced expansion beyond a four-mismatch budget. It must remain absent
from enumerator logic and tuning fixtures: add it only to the blinded evaluation
manifest after candidate generation, with a target of at least 90% held-out recall.
