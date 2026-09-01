# Minigene reporter-context tooling

This package deterministically concatenates explicitly retained regions of a
caller-supplied reporter sequence, remaps event/design intervals into the
compressed coordinates, and returns a `TargetDefinition` accepted by the design
engine. Every source, wrapper, and engine `TranscriptContext` retains the frozen
`minigene_reporter` domain label. Coordinates are zero-based, half-open reporter
construct coordinates—not genomic coordinates.

The input sequence, its SHA-256 digest, provenance text, retained blocks, and
coordinate mapping are validated fail-closed. Event and design intervals must be
wholly contained in one retained block; this prevents a design window or modeled
element from silently crossing an artificial compression junction.

Use `design_compressed_context`, rather than calling the generic engine directly,
to produce a `ReporterDesignReport`. The wrapper repeats the reporter domain and
claim boundary at top level, validates the complete compression receipt against
the nested engine manifest, and rejects every generated target footprint or
dCas13 shadow that touches an artificial junction. Its JSON and text renderers
therefore cannot silently serialize the result as an endogenous design.

This is compute-only representation tooling. It does not download or assert any
public reporter data, evaluate a published result, establish preservation of
splicing behavior, or validate a biological or therapeutic claim.
