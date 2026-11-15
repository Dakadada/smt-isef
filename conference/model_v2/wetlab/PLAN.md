# Exon-3 audit implementation plan

Spec: ../ops/briefs/WETLAB_exon3.md. This task explicitly authorizes implementation.
Design: read only the exon-3 workbook with openpyxl in formula and cached modes.
Group raw bands by gel/lane, independently compute conditional signal fractions,
and compare both formula results and caches without modifying the source.
Preserve a cell inventory and column-level evidence schema. No sequence or
coordinate inference; no computational joins. Plot individual positional values
with gel separation and no statistics. Outputs are ignored.

- [x] Write synthetic tests for missing/zero/invalid bands, duplicate bands,
  formula references, stale caches, metadata, positions and output evidence.
- [x] Run tests red, implement parser, run focused tests green (20 passed).
- [x] Run source audit, inspect outputs and unchanged SHA256.
- [x] Run the feng gate once; document exact counts and limitations.

Execution rulings: scratch is read-only, so use ignored owned `outputs/` and
document the regeneration command. Full gate: 581 passed, 5 skipped, 1 failed;
the cleanroom substring check flags unchanged `feng/src/feng/source_auth.py`.
Leave that unrelated code untouched; repository-wide acceptance remains blocked.

Review focus: absent bands cannot become zeros; invalid volumes cannot produce
scores; raw grouping must be independent of quantification formula references;
cached formulas may be missing or stale; no biological/technical replicate
classification or primer-derived size may be invented.
