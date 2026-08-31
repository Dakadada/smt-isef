"""JSON/FASTA command line entry point with deterministic TSV output."""

import argparse
import csv
import json
import sys

from .enumerator import enumerate_offtargets
from .models import EnumeratorConfig, Evidence, Transcript


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--transcripts", required=True,
                        help="JSONL records with complete pre-mRNA sequence and annotations")
    parser.add_argument("--evidence", help="optional JSONL eCLIP/perturbation records")
    args = parser.parse_args()
    with open(args.config, encoding="utf-8") as handle:
        raw = json.load(handle)
    with open(args.transcripts, encoding="utf-8") as handle:
        transcripts = [Transcript(**json.loads(line)) for line in handle if line.strip()]
    evidence = []
    if args.evidence:
        with open(args.evidence, encoding="utf-8") as handle:
            evidence = [Evidence(**json.loads(line)) for line in handle if line.strip()]
    config = EnumeratorConfig(**raw.get("enumerator", {}))
    candidates = enumerate_offtargets(raw["effector_sequence"], transcripts, config, evidence)
    fields = list(CandidateFieldNames)
    writer = csv.DictWriter(sys.stdout, fieldnames=fields, dialect="excel-tab")
    writer.writeheader()
    for c in candidates:
        writer.writerow({name: getattr(c, name) for name in fields})


CandidateFieldNames = (
    "rank", "transcript_id", "gene_id", "start", "end", "target_sequence",
    "mismatches", "insertions", "deletions", "edit_distance", "cigar",
    "chemistry", "chemistry_parameters_known", "estimated_delta_g",
    "expression_tpm", "nuclear_fraction", "expression_localization_weight",
    "effect_sign", "recall_score", "confidence",
)


if __name__ == "__main__":
    main()
