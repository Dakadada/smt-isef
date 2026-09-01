"""Committed production contexts; no offline fixtures are permitted here."""
from __future__ import annotations

import json
from pathlib import Path

from design_engine import ObjectiveDirection, TargetDefinition, TherapeuticObjective
from physics_adapter import EventDefinition, Interval, TranscriptContext

from .models import RealTranscriptContext
from .panel import bind_real_transcript_context


DEFAULT_SMN2_CONTEXT = Path(__file__).resolve().parents[2] / "config" / "smn2_issn1_context.json"


def load_smn2_context(path: Path = DEFAULT_SMN2_CONTEXT) -> RealTranscriptContext:
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    required = {"schema_version", "target_id", "locus", "assembly", "chrom", "strand",
                "transcript_version", "sequence", "genomic_start_1based_inclusive",
                "genomic_end_1based_inclusive", "sequence_provenance"}
    if set(raw) != required or raw["schema_version"] != "real-transcript-context/1":
        raise ValueError("unsupported or incomplete SMN2 real-context schema")
    sequence = raw["sequence"].upper()
    start = int(raw["genomic_start_1based_inclusive"])
    end = int(raw["genomic_end_1based_inclusive"])
    if (raw["target_id"] != "SMN2_exon7_inclusion" or raw["locus"] != "SMN2"
            or raw["assembly"] != "GRCh38" or raw["chrom"] != "chr5"
            or raw["strand"] != "+" or raw["transcript_version"] != "NM_017411.4"
            or end - start + 1 != len(sequence)):
        raise ValueError("SMN2 real-context identity or coordinate mismatch")
    objective = TherapeuticObjective(raw["target_id"], raw["locus"],
                                     ObjectiveDirection.INCREASE_INCLUSION)
    definition = TargetDefinition(raw["target_id"], raw["locus"], objective,
        TranscriptContext(raw["transcript_version"], sequence),
        EventDefinition(raw["target_id"], (), (), (Interval(0, len(sequence)),)),
        0, len(sequence))
    mapping = tuple((start + offset, offset) for offset in range(len(sequence)))
    return bind_real_transcript_context(definition=definition, assembly=raw["assembly"],
        chrom=raw["chrom"], strand=raw["strand"],
        transcript_version=raw["transcript_version"],
        sequence_provenance=raw["sequence_provenance"],
        genomic_to_transcript_mapping=mapping)


def default_real_contexts() -> dict[str, RealTranscriptContext]:
    return {"SMN2_exon7_inclusion": load_smn2_context()}
