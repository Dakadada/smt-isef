"""Deterministic reporter-context compression into design-engine coordinates."""
from __future__ import annotations

from dataclasses import dataclass

from design_engine import TargetDefinition
from physics_adapter import EventDefinition, Interval

from .models import (CLAIM_BOUNDARY, CompressedReporterContext, ContextDomain,
                     ReporterSourceContext, ReporterTranscriptContext,
                     RetainedBlock, _event_intervals, _hash)


@dataclass(frozen=True)
class _Derived:
    mapping: tuple[tuple[int, int], ...]
    target: TargetDefinition
    manifest_sha256: str


def _validate_blocks(source: ReporterSourceContext,
                     blocks: tuple[RetainedBlock, ...]) -> None:
    if not blocks:
        raise ValueError("at least one retained reporter block is required")
    if len({block.label for block in blocks}) != len(blocks):
        raise ValueError("retained block labels must be unique")
    if tuple(sorted(blocks, key=lambda block: block.source_interval.start)) != blocks:
        raise ValueError("retained blocks must be supplied in source-coordinate order")
    previous_end = -1
    for block in blocks:
        interval = block.source_interval
        if interval.end > len(source.sequence):
            raise ValueError("retained block lies outside the reporter source sequence")
        if interval.start < previous_end:
            raise ValueError("retained reporter blocks must not overlap")
        previous_end = interval.end


def _remap(interval: Interval, blocks: tuple[RetainedBlock, ...]) -> Interval:
    """Map one interval only when a single retained block covers it exactly."""
    compressed_offset = 0
    for block in blocks:
        retained = block.source_interval
        if retained.start <= interval.start and interval.end <= retained.end:
            return Interval(compressed_offset + interval.start - retained.start,
                            compressed_offset + interval.end - retained.start)
        compressed_offset += retained.end - retained.start
    raise ValueError("event/design interval is not wholly retained in one reporter block")


def _derive(source: ReporterSourceContext,
            blocks: tuple[RetainedBlock, ...]) -> _Derived:
    source.__post_init__()
    _validate_blocks(source, blocks)
    mapping: list[tuple[int, int]] = []
    sequence_parts: list[str] = []
    compressed_offset = 0
    for block in blocks:
        interval = block.source_interval
        sequence_parts.append(source.sequence[interval.start:interval.end])
        mapping.extend((source_position, compressed_offset + source_position - interval.start)
                       for source_position in range(interval.start, interval.end))
        compressed_offset += interval.end - interval.start
    sequence = "".join(sequence_parts)
    artificial_junction_offsets = tuple(
        sum(block.source_interval.end - block.source_interval.start for block in blocks[:index])
        for index in range(1, len(blocks))
        if blocks[index - 1].source_interval.end != blocks[index].source_interval.start
    )
    event = EventDefinition(
        source.event.event_id,
        tuple(_remap(item, blocks) for item in source.event.u1_contacts),
        tuple(_remap(item, blocks) for item in source.event.branchpoint_ppt),
        tuple(_remap(item, blocks) for item in source.event.silencers),
    )
    design = _remap(source.design_interval, blocks)
    manifest = {
        "schema_version": "minigene-context-compression-input/1",
        "domain_label": ContextDomain.MINIGENE_REPORTER,
        "claim_boundary": CLAIM_BOUNDARY,
        "source": source,
        "retained_blocks": blocks,
    }
    manifest_sha256 = _hash(manifest)
    context = ReporterTranscriptContext(
        transcript_id=f"{source.reporter_id}::compressed", sequence=sequence,
        domain_label=ContextDomain.MINIGENE_REPORTER,
        source_sequence_sha256=source.sequence_sha256,
        source_sequence_provenance=source.sequence_provenance,
        retained_blocks=tuple((block.label, block.source_interval.start, block.source_interval.end)
                              for block in blocks),
        source_to_context_mapping=tuple(mapping),
        artificial_junction_offsets=artificial_junction_offsets,
        compression_manifest_sha256=manifest_sha256,
        claim_boundary=CLAIM_BOUNDARY,
    )
    target = TargetDefinition(source.event.event_id, source.locus, source.objective,
                              context, event, design.start, design.end)
    return _Derived(tuple(mapping), target, manifest_sha256)


def compress_reporter_context(
        source: ReporterSourceContext,
        retained_blocks: tuple[RetainedBlock, ...]) -> CompressedReporterContext:
    """Compress caller-supplied reporter blocks and remap them for the engine.

    This function evaluates no public result and makes no claim that the
    compressed construct preserves endogenous or experimental behaviour.
    """
    if not isinstance(source, ReporterSourceContext):
        raise TypeError("source must be a ReporterSourceContext")
    blocks = tuple(retained_blocks)
    if any(not isinstance(block, RetainedBlock) for block in blocks):
        raise TypeError("retained blocks must be typed RetainedBlock values")
    derived = _derive(source, blocks)
    return CompressedReporterContext(
        "minigene-context-compression/1", ContextDomain.MINIGENE_REPORTER,
        CLAIM_BOUNDARY, source, blocks, derived.mapping, derived.target,
        derived.manifest_sha256,
    )
