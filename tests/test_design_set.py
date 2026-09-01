from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
from types import SimpleNamespace
from tempfile import TemporaryDirectory
import unittest
from unittest import mock

from design_engine import PROSPECTIVE_BANNER, PROSPECTIVE_STATUS, design_for
from design_panel.panel import committed_real_contexts
from design_set import (build_design_set, deposit, prediction_axis_tiers,
                        render_json, render_text, superior_candidate_ids)
from design_set.builder import _reference_placement


FIXED_TIME = datetime(2026, 9, 1, 12, 0, 0, tzinfo=timezone.utc)


def scored(identifier: str, confidence: float, burden: float, robustness: float):
    candidate = SimpleNamespace(candidate_id=identifier)
    return SimpleNamespace(candidate=candidate, directional_confidence=confidence,
                           off_target_burden=burden, span_robustness=robustness,
                           abstained=False, goal_alignment="supports_goal")


class ProspectiveDesignSetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.report = build_design_set(clock=lambda: FIXED_TIME)

    def test_frozen_clock_reproduces_identical_artifact_bytes(self):
        first = build_design_set(clock=lambda: FIXED_TIME)
        second = build_design_set(clock=lambda: FIXED_TIME)
        self.assertEqual(render_json(first).encode(), render_json(second).encode())
        self.assertEqual(render_text(first).encode(), render_text(second).encode())

    def test_superior_set_includes_dominator_and_excludes_dominated_candidate(self):
        reference = scored("reference", .8, .4, .7)
        dominating = scored("dominating", .9, .3, .7)
        dominated = scored("dominated", .7, .5, .6)
        tied = scored("tied", .8, .4, .7)
        candidates = (reference, dominating, dominated, tied)
        self.assertEqual(superior_candidate_ids(candidates, reference), ("dominating",))
        tiers = prediction_axis_tiers(candidates)
        self.assertEqual(tiers["dominating"], 1)
        self.assertGreater(tiers["dominated"], tiers["dominating"])

    def test_exact_banner_status_reach_artifact_and_every_candidate(self):
        payload = json.loads(render_json(self.report))
        self.assertEqual(payload["banner"], PROSPECTIVE_BANNER)
        self.assertEqual(payload["prediction_status"], PROSPECTIVE_STATUS)
        self.assertTrue(payload["candidates"])
        for row in payload["candidates"]:
            self.assertEqual(row["banner"], PROSPECTIVE_BANNER)
            self.assertEqual(row["prediction_status"], PROSPECTIVE_STATUS)
        self.assertIn("NOT experimentally validated; NOT efficacy/potency",
                      payload["claim_boundary"])

    def test_renderer_refuses_non_prospective_artifact(self):
        forged = build_design_set(clock=lambda: FIXED_TIME)
        object.__setattr__(forged, "prediction_status", "validated")
        with self.assertRaisesRegex(ValueError, "prospective"):
            render_json(forged)

    def test_deposit_refuses_same_timestamp_overwrite(self):
        with TemporaryDirectory() as directory:
            first = deposit(self.report, Path(directory))
            self.assertTrue(all(path.is_file() for path in first))
            with self.assertRaises(FileExistsError):
                deposit(self.report, Path(directory))

    def test_builder_sources_authenticated_committed_context_seam(self):
        contexts = committed_real_contexts()
        with mock.patch("design_set.builder.committed_real_contexts",
                        return_value=contexts) as source:
            report = build_design_set(clock=lambda: FIXED_TIME)
        source.assert_called_once_with()
        self.assertEqual(report.provenance["committed_context_source"],
                         "design_panel.panel.committed_real_contexts")
        self.assertEqual(report.provenance["binding_sequence_sha256"],
                         contexts["SMN2_exon7_inclusion"].binding_sequence_sha256)

    def test_reference_resolves_to_committed_nusinersen_class_18mer(self):
        reference = self.report.reference
        self.assertEqual(reference.coordinate_status, "resolved-full-18mer")
        self.assertEqual((reference.start, reference.end), (0, 18))
        self.assertEqual(reference.chemistry, "2MOE")
        self.assertEqual(reference.candidate_id,
                         "SMN2_exon7_inclusion:ASO:2MOE:0-18")

    def test_undersized_window_withholds_reference_coordinates_and_comparison(self):
        context = committed_real_contexts()["SMN2_exon7_inclusion"]
        engine_report = design_for(context.definition, clock=lambda: FIXED_TIME)
        short_context = SimpleNamespace(
            definition=SimpleNamespace(
                target_id=context.definition.target_id, locus=context.definition.locus,
                context=SimpleNamespace(sequence=context.definition.context.sequence[:17])),
            assembly=context.assembly, chrom=context.chrom, strand=context.strand,
            transcript_version=context.transcript_version,
            genomic_to_transcript_mapping=context.genomic_to_transcript_mapping[:17])
        placement, reference = _reference_placement(engine_report, short_context)
        self.assertEqual(placement.coordinate_status,
                         "reference-position-within-committed-window")
        self.assertIsNone(placement.start)
        self.assertIsNone(placement.end)
        self.assertIsNone(placement.candidate_id)
        self.assertIsNone(reference)
        self.assertIn("cannot verify the expected ISS-N1 18-mer", placement.limitation)
        self.assertEqual(superior_candidate_ids(engine_report.candidates, reference), ())

    def test_reference_resolution_fails_closed_on_mapping_or_sequence_mismatch(self):
        context = committed_real_contexts()["SMN2_exon7_inclusion"]
        engine_report = design_for(context.definition, clock=lambda: FIXED_TIME)
        for mapping, sequence in (
                (tuple((genomic + 1, offset) for genomic, offset
                       in context.genomic_to_transcript_mapping),
                 context.definition.context.sequence),
                (context.genomic_to_transcript_mapping,
                 "A" + context.definition.context.sequence[1:])):
            forged = SimpleNamespace(
                definition=SimpleNamespace(target_id=context.definition.target_id,
                    locus=context.definition.locus,
                    context=SimpleNamespace(sequence=sequence)),
                assembly=context.assembly, chrom=context.chrom, strand=context.strand,
                transcript_version=context.transcript_version,
                genomic_to_transcript_mapping=mapping)
            placement, reference = _reference_placement(engine_report, forged)
            self.assertEqual(placement.coordinate_status,
                             "reference-position-within-committed-window")
            self.assertIsNone(placement.start)
            self.assertIsNone(placement.end)
            self.assertIsNone(reference)

    def test_current_tied_prediction_axes_make_superior_set_honestly_empty(self):
        context = committed_real_contexts()["SMN2_exon7_inclusion"]
        self.assertEqual(self.report.superior_candidate_ids, ())
        self.assertTrue(all(row.pareto and row.rank == 1 for row in self.report.candidates))
        self.assertEqual({(row.directional_confidence, row.off_target_burden,
                           row.span_robustness) for row in self.report.candidates},
                         {(0.9, 0.0, 1.0)})
        self.assertEqual(len(context.definition.offtarget_transcripts), 0)
        self.assertEqual(self.report.provenance["off_target_corpus_size"], 0)
        self.assertNotIn("works better", render_text(self.report).lower())

    def test_renderer_rejects_forged_superiority_annotation(self):
        forged = build_design_set(clock=lambda: FIXED_TIME)
        candidate = next(row for row in forged.candidates if not row.is_reference)
        object.__setattr__(candidate, "off_target_burden", 999.0)
        object.__setattr__(candidate, "predicted_superior_to_reference", True)
        object.__setattr__(forged, "superior_candidate_ids", (candidate.candidate_id,))
        with self.assertRaisesRegex(
                ValueError, "recomputed three-axis weak dominance against reference"):
            render_json(forged)

    def test_empty_off_target_corpus_disclosure_is_in_both_rendered_outputs(self):
        payload = json.loads(render_json(self.report))
        provenance = payload["provenance"]
        self.assertEqual(provenance["off_target_corpus_size"], 0)
        self.assertIs(provenance["off_target_burden_evaluated"], False)
        disclosure = provenance["off_target_disclosure"]
        self.assertIn("NOT evaluated", disclosure)
        self.assertIn("non-discriminating", disclosure)
        self.assertIn("must NOT be read as 'no better molecule exists.'", disclosure)
        text = render_text(self.report)
        self.assertIn(f"Off-target disclosure: {disclosure}", text)


if __name__ == "__main__":
    unittest.main()
