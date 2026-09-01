import hashlib
import json
import unittest
from dataclasses import replace
from datetime import datetime, timezone
from unittest.mock import patch

from design_engine import (DesignEngine, GenerationConfig, ObjectiveDirection,
                           TherapeuticObjective)
from design_engine.engine import _hash as _engine_hash
from minigene_context import (CLAIM_BOUNDARY, ContextDomain, ReporterSourceContext,
                              ReporterTranscriptContext, RetainedBlock,
                              compress_reporter_context, design_compressed_context,
                              render_json, render_text)
from physics_adapter import EventDefinition, Interval


class MinigeneContextTests(unittest.TestCase):
    def source(self, **changes):
        sequence = ("ACGT" * 30)[:120]
        objective = TherapeuticObjective(
            "fixture_reporter_event", "FIXTURE", ObjectiveDirection.INCREASE_INCLUSION,
            override_provenance="test-only objective; no biological claim")
        values = dict(
            reporter_id="fixture-reporter-v1", locus="FIXTURE", sequence=sequence,
            sequence_sha256=hashlib.sha256(sequence.encode("ascii")).hexdigest(),
            sequence_provenance="caller-supplied test fixture; no public-data assertion",
            objective=objective,
            event=EventDefinition("fixture_reporter_event", (Interval(62, 66),),
                                  (Interval(72, 77),), (Interval(82, 88),)),
            design_interval=Interval(55, 95),
        )
        values.update(changes)
        return ReporterSourceContext(**values)

    def blocks(self):
        return (RetainedBlock("upstream", Interval(0, 30)),
                RetainedBlock("event-and-design", Interval(50, 110)))

    def test_compresses_sequence_coordinates_and_retains_domain_end_to_end(self):
        result = compress_reporter_context(self.source(), self.blocks())
        self.assertEqual(result.domain_label, ContextDomain.MINIGENE_REPORTER)
        self.assertIsInstance(result.target.context, ReporterTranscriptContext)
        self.assertEqual(result.target.context.domain_label, ContextDomain.MINIGENE_REPORTER)
        self.assertEqual(result.target.context.sequence,
                         self.source().sequence[:30] + self.source().sequence[50:110])
        self.assertEqual(result.target.design_start, 35)
        self.assertEqual(result.target.design_end, 75)
        self.assertEqual(result.target.event.u1_contacts, (Interval(42, 46),))
        self.assertEqual(result.source_to_context_mapping[30], (50, 30))
        payload = result.to_dict()
        self.assertEqual(payload["domain_label"], "minigene_reporter")
        self.assertEqual(payload["target"]["context"]["domain_label"], "minigene_reporter")
        self.assertEqual(payload["claim_boundary"], CLAIM_BOUNDARY)

    def test_engine_manifest_retains_reporter_domain_label(self):
        result = compress_reporter_context(self.source(), self.blocks())
        report = DesignEngine(
            generation_config=GenerationConfig(max_candidates=8),
            clock=lambda: datetime(2026, 8, 31, tzinfo=timezone.utc),
        ).design_for(result.target)
        context = report.provenance["run_manifest"]["target_definition"]["context"]
        self.assertEqual(context["domain_label"], "minigene_reporter")
        self.assertEqual(context["source_coordinate_convention"],
                         "reporter-construct-zero-based-half-open")
        self.assertEqual(context["source_sequence_sha256"], self.source().sequence_sha256)
        self.assertEqual(context["source_sequence_provenance"], self.source().sequence_provenance)
        self.assertEqual(context["retained_blocks"],
                         (("upstream", 0, 30), ("event-and-design", 50, 110)))
        self.assertEqual(context["source_to_context_mapping"][30], (50, 30))
        self.assertEqual(context["artificial_junction_offsets"], (30,))
        self.assertEqual(context["compression_manifest_sha256"], result.input_manifest_sha256)
        self.assertEqual(context["claim_boundary"], CLAIM_BOUNDARY)
        self.assertTrue(report.candidates)

    def test_reporter_wrapper_and_renderers_surface_round_trip_receipt(self):
        compressed = compress_reporter_context(self.source(), self.blocks())
        wrapped = design_compressed_context(compressed, DesignEngine(
            clock=lambda: datetime(2026, 8, 31, tzinfo=timezone.utc)))
        payload = json.loads(render_json(wrapped))
        self.assertEqual(payload["domain_label"], "minigene_reporter")
        self.assertEqual(payload["claim_boundary"], CLAIM_BOUNDARY)
        receipt = payload["compression_receipt"]
        self.assertEqual(receipt["source_to_context_mapping"][30], [50, 30])
        self.assertEqual(receipt["target"]["context"]["source_sequence_provenance"],
                         self.source().sequence_provenance)
        self.assertEqual(receipt["target"]["context"]["artificial_junction_offsets"], [30])
        nested = payload["design_report"]["provenance"]["run_manifest"]
        self.assertEqual(nested["target_definition"]["context"]["claim_boundary"], CLAIM_BOUNDARY)
        text = render_text(wrapped)
        self.assertIn("DOMAIN: minigene_reporter", text)
        self.assertIn("not an endogenous-genome representation", text)
        self.assertIn(self.source().sequence_provenance, text)

    def test_forged_or_omitted_engine_metadata_is_rejected_even_if_rehashed(self):
        compressed = compress_reporter_context(self.source(), self.blocks())
        wrapped = design_compressed_context(compressed)
        base = wrapped.design_report
        provenance = base.to_dict()["provenance"]
        del provenance["run_manifest"]["target_definition"]["context"]["claim_boundary"]
        provenance["run_manifest_sha256"] = _engine_hash(provenance["run_manifest"])
        forged = replace(base, provenance=provenance)
        with self.assertRaisesRegex(ValueError, "target definition differs"):
            replace(wrapped, design_report=forged)

    def test_rehashed_event_and_design_interval_substitution_is_rejected(self):
        wrapped = design_compressed_context(
            compress_reporter_context(self.source(), self.blocks()))
        for field, value in (("design_start", 36), ("event", {
                "event_id": "fixture_reporter_event", "u1_contacts": [[43, 46]],
                "branchpoint_ppt": [[52, 57]], "silencers": [[62, 68]]})):
            base = wrapped.design_report
            provenance = base.to_dict()["provenance"]
            provenance["run_manifest"]["target_definition"][field] = value
            provenance["run_manifest_sha256"] = _engine_hash(provenance["run_manifest"])
            forged = replace(base, provenance=provenance)
            with self.assertRaisesRegex(ValueError, "target definition differs"):
                replace(wrapped, design_report=forged)

    def test_rehashed_candidate_score_substitution_fails_deterministic_rerun(self):
        wrapped = design_compressed_context(
            compress_reporter_context(self.source(), self.blocks()))
        base = wrapped.design_report
        changed = replace(base.candidates[0],
                          directional_confidence=base.candidates[0].directional_confidence + 0.01)
        candidates = (changed,) + base.candidates[1:]
        provenance = base.to_dict()["provenance"]
        provenance["run_manifest"]["scored_candidates"] = [
            base.to_dict()["candidates"][index] if index else {
                **base.to_dict()["candidates"][0],
                "directional_confidence": changed.directional_confidence}
            for index in range(len(candidates))]
        provenance["run_manifest_sha256"] = _engine_hash(provenance["run_manifest"])
        forged = replace(base, candidates=candidates, provenance=provenance)
        with self.assertRaisesRegex(ValueError, "deterministic approved-engine output"):
            replace(wrapped, design_report=forged)

    def test_rehashed_grna_shadow_change_is_rejected(self):
        compressed = compress_reporter_context(self.source(), self.blocks())
        with self.assertRaisesRegex(ValueError, "configuration is not approved"):
            design_compressed_context(compressed, DesignEngine(
                generation_config=GenerationConfig(grna_shadow=0)))
        wrapped = design_compressed_context(compressed)
        base = wrapped.design_report
        provenance = base.to_dict()["provenance"]
        provenance["run_manifest"]["generation_config"]["grna_shadow"] = 0
        provenance["run_manifest_sha256"] = _engine_hash(provenance["run_manifest"])
        forged = replace(base, provenance=provenance)
        with self.assertRaisesRegex(ValueError, "deterministic approved-engine output"):
            replace(wrapped, design_report=forged)

    def test_construction_and_render_ignore_replaced_mutable_engine_seams(self):
        compressed = compress_reporter_context(self.source(), self.blocks())
        with patch("minigene_context.design.DesignEngine", side_effect=AssertionError("seam")), \
             patch("minigene_context.design.generate_candidates", side_effect=AssertionError("seam")), \
             patch("minigene_context.design._assert_approved_engine", side_effect=AssertionError("seam")), \
             patch("minigene_context.design._rerun_approved", side_effect=AssertionError("seam")):
            wrapped = design_compressed_context(compressed)
            self.assertIn("minigene_reporter", render_json(wrapped))

    def test_source_sequence_digest_and_provenance_fail_closed(self):
        with self.assertRaisesRegex(ValueError, "digest mismatch"):
            self.source(sequence_sha256="0" * 64)
        with self.assertRaisesRegex(ValueError, "provenance"):
            self.source(sequence_provenance="")
        with self.assertRaisesRegex(ValueError, "uppercase DNA"):
            self.source(sequence="ACGX", sequence_sha256=hashlib.sha256(b"ACGX").hexdigest())

    def test_block_order_overlap_bounds_and_identity_fail_closed(self):
        source = self.source()
        with self.assertRaisesRegex(ValueError, "source-coordinate order"):
            compress_reporter_context(source, tuple(reversed(self.blocks())))
        with self.assertRaisesRegex(ValueError, "must not overlap"):
            compress_reporter_context(source, (RetainedBlock("a", Interval(0, 70)),
                                                RetainedBlock("b", Interval(60, 110))))
        with self.assertRaisesRegex(ValueError, "outside"):
            compress_reporter_context(source, (RetainedBlock("a", Interval(0, 121)),))
        with self.assertRaisesRegex(ValueError, "labels must be unique"):
            compress_reporter_context(source, (RetainedBlock("same", Interval(0, 30)),
                                                RetainedBlock("same", Interval(50, 110))))

    def test_deleted_or_junction_crossing_model_intervals_are_rejected(self):
        source = self.source(event=EventDefinition(
            "fixture_reporter_event", (Interval(25, 55),), (), ()))
        with self.assertRaisesRegex(ValueError, "not wholly retained"):
            compress_reporter_context(source, self.blocks())

    def test_candidate_footprint_or_dcas_shadow_touching_join_fails_closed(self):
        footprint_source = self.source(design_interval=Interval(10, 30))
        footprint = compress_reporter_context(footprint_source, self.blocks())
        with self.assertRaisesRegex(ValueError, "candidate footprint touches"):
            design_compressed_context(footprint)

        shadow_source = self.source(design_interval=Interval(54, 100))
        shadow = compress_reporter_context(shadow_source, self.blocks())
        with self.assertRaisesRegex(ValueError, "dCas13 shadow touches"):
            design_compressed_context(shadow)
        source = self.source(design_interval=Interval(25, 60))
        with self.assertRaisesRegex(ValueError, "not wholly retained"):
            compress_reporter_context(source, self.blocks())

    def test_derived_mapping_target_and_manifest_cannot_be_substituted(self):
        result = compress_reporter_context(self.source(), self.blocks())
        with self.assertRaisesRegex(ValueError, "coordinate mapping"):
            replace(result, source_to_context_mapping=result.source_to_context_mapping[:-1])
        with self.assertRaisesRegex(ValueError, "deterministic compression output"):
            replace(result, target=replace(result.target, design_start=result.target.design_start + 1))
        with self.assertRaisesRegex(ValueError, "manifest digest"):
            replace(result, input_manifest_sha256="0" * 64)
        with self.assertRaisesRegex(ValueError, "domain"):
            replace(result, domain_label="endogenous_transcript")


if __name__ == "__main__":
    unittest.main()
