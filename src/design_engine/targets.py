"""Load target labels from the frozen junction config and provide offline fixtures."""
from __future__ import annotations

import csv
from pathlib import Path

from offtarget import Transcript
from physics_adapter import EventDefinition, Interval, TranscriptContext

from .models import ObjectiveDirection, TargetDefinition, TherapeuticObjective


DEFAULT_EVENT_CONFIG = Path(__file__).resolve().parents[2] / "config" / "event_junctions.GRCh38.tsv"


def load_target_polarities(path: Path = DEFAULT_EVENT_CONFIG) -> dict[str, str]:
    """Compatibility view of the typed, fixed built-in objectives."""
    polarities = {"SMN2_exon7_inclusion": ObjectiveDirection.INCREASE_INCLUSION,
                  "DMD_exon53_inclusion": ObjectiveDirection.DECREASE_INCLUSION}
    with path.open(newline="", encoding="utf-8") as handle:
        event_ids = {row["event"] for row in csv.DictReader(handle, delimiter="\t")}
    missing = set(polarities) - event_ids
    if missing:
        raise ValueError(f"therapeutic targets missing from event config: {sorted(missing)}")
    return {event_id: polarities[event_id].value for event_id in sorted(event_ids) if event_id in polarities}


def offline_target(target_id: str, path: Path = DEFAULT_EVENT_CONFIG) -> TargetDefinition:
    """Small deterministic context for CI; callers supply real transcript contexts in production."""
    polarity = load_target_polarities(path).get(target_id)
    if polarity is None:
        raise KeyError(f"unknown built-in target: {target_id}")
    locus = "SMN2" if target_id.startswith("SMN2") else "DMD"
    # Both regulatory kinds are present so therapeutic polarity selects opposite designs.
    sequence = ("ACGT" * 30)[:120]
    context = TranscriptContext(f"{locus}_OFFLINE_FIXTURE", sequence)
    event = EventDefinition(target_id, (Interval(28, 38), Interval(72, 82)),
                            (Interval(44, 50),), (Interval(54, 66),))
    # Same synthetic motif content guarantees that the real enumerator contributes
    # a non-zero burden axis during offline integration tests.
    decoy = Transcript("OFFTARGET_FIXTURE", "DECOY", sequence, 20.0, 0.8)
    objective = TherapeuticObjective(target_id, locus, ObjectiveDirection(polarity))
    return TargetDefinition(target_id, locus, objective, context, event, 20, 90, (decoy,))
