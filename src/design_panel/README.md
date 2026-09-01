# Prospective cross-disease design panel

This package runs the stable `design_engine.design_for` API only for panel targets
whose coordinates resolve exactly against the hash-pinned registry **and** whose
committed transcript context matches the registry's
`sequence_sha256`. Targets missing either proof remain pending and never reach the
design engine. The package emits one frozen-clock JSON report and a human-readable
designability table. These are prospective predictions, not experimentally validated
results. The default panel now loads the real, sequence-bound SMN2 ISS-N1 context from
`config/smn2_issn1_context.json`; DMD and every other target remain pending. SMN2 is
the only locus intended for real benchmark support; every other target is exploratory.

The checked-in panel is [`target_panel.json`](target_panel.json). In a networked
curation session, add a real target there **only after** adding its independently
verified, transcript-version-exact GRCh38 event to
[`../../config/event_junctions.GRCh38.tsv`](../../config/event_junctions.GRCh38.tsv)
and committing the exact transcript-context sequence digest in its `sequence_sha256`
column. Supply that real `TargetDefinition`, retrieval provenance, and genomic-to-
transcript mapping through `bind_real_transcript_context`. Change `coordinate_status`
to `verified`, replace the `REQUIRED-INPUT` provenance with the registry digest-bound
value, and set `event_id` plus the transcript-version-exact identity fields. Until all
those inputs exist, keep `event_id: null`: the loader rejects unverified entries that
claim an event, while verified rows with an empty sequence digest fail into pending
without calling the design engine.

```python
from datetime import datetime, timezone
from pathlib import Path
from design_panel import deposit, run_panel

report = run_panel(clock=lambda: datetime(2026, 1, 2, tzinfo=timezone.utc))
json_path, text_path = deposit(report, Path("artifacts/design-panels"))
```

`deposit` uses the frozen timestamp in both filenames and refuses to overwrite an
existing report.
