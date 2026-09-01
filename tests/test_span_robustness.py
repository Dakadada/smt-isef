import unittest
from dataclasses import replace
import span_robustness.method as method
from span_robustness import PredictedDirection,evaluate_sliding_windows

class SpanRobustnessTests(unittest.TestCase):
    def test_fixture_is_non_tautological_and_uses_actual_adapter(self):
        report=evaluate_sliding_windows()
        directions={x.direction for x in report.windows}
        self.assertIn(PredictedDirection.INCREASES_INCLUSION,directions)
        self.assertIn(PredictedDirection.DECREASES_INCLUSION,directions)
        self.assertIn(PredictedDirection.ABSTAIN,directions)
        provenance=dict(report.engine_provenance)
        self.assertIn("EffectorPhysicsAdapter",provenance["engine_identity"])
        self.assertEqual(provenance["backend_version"],"design-engine-deterministic-stub/1")
        self.assertEqual(provenance["outcomes_consumed"],"false")
        report.to_dict()

    def test_exact_grid_and_explicit_event_geometry_bound(self):
        report=evaluate_sliding_windows(window_size=3,stride=1)
        self.assertEqual(report.accessible_elements,((2,5),)); self.assertEqual(report.silencers,((9,12),))
        self.assertEqual([(x.start,x.end) for x in report.windows],[(i,i+3) for i in range(14)])

    def test_arbitrary_ids_accessions_and_leads_rejected(self):
        for value in ("fixture","STMN2","unc13-a","NM_001 accession","PRJNA1052769",""):
            with self.subTest(value=value),self.assertRaisesRegex(ValueError,"allowlist"):
                evaluate_sliding_windows(fixture_id=value)

    def test_no_arbitrary_sequence_argument_exists(self):
        with self.assertRaises(TypeError): evaluate_sliding_windows(sequence="ACGT")

    def test_stride_is_overlapping_or_contiguous(self):
        with self.assertRaisesRegex(ValueError,"stride"): evaluate_sliding_windows(window_size=2,stride=3)
        evaluate_sliding_windows(window_size=2,stride=2)

    def test_forged_direction_and_summary_rejected(self):
        report=evaluate_sliding_windows()
        opposite=(PredictedDirection.INCREASES_INCLUSION if report.windows[0].direction is not PredictedDirection.INCREASES_INCLUSION else PredictedDirection.DECREASES_INCLUSION)
        forged=replace(report.windows[0],direction=opposite)
        with self.assertRaisesRegex(ValueError,"approved-engine output"): replace(report,windows=(forged,)+report.windows[1:])
        with self.assertRaisesRegex(ValueError,"summary"): replace(report,evaluable_window_count=999)

    def test_fixture_context_tampering_rejected(self):
        report=evaluate_sliding_windows()
        for field,value in (("fixture_id","OTHER"),("sequence","T"+report.sequence[1:]),
                            ("sequence_sha256","0"*64),("accessible_elements",((0,1),)),
                            ("silencers",((0,1),))):
            with self.subTest(field=field),self.assertRaisesRegex(ValueError,"allowlisted"):
                replace(report,**{field:value})

    def test_provenance_code_config_tampering_rejected(self):
        report=evaluate_sliding_windows(); provenance=dict(report.engine_provenance)
        for key in ("model_config_json","model_protocol_sha256","backend_version",
                    "adapter_source_sha256","package_source_manifest_sha256"):
            changed=dict(provenance); changed[key]="tampered"
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,"provenance"):
                replace(report,engine_provenance=tuple(sorted(changed.items())))

    def test_actual_predictor_reassignment_cannot_authorize_forgery(self):
        original_report=evaluate_sliding_windows(); original=method._predict_direction
        try:
            method._predict_direction=lambda *args: PredictedDirection.INCREASES_INCLUSION
            rebuilt=evaluate_sliding_windows()
            self.assertEqual(rebuilt.windows,original_report.windows)
            original_report.to_dict()
            all_increase=tuple(replace(row,direction=PredictedDirection.INCREASES_INCLUSION)
                               for row in original_report.windows)
            with self.assertRaisesRegex(ValueError,"approved-engine output"):
                replace(original_report,windows=all_increase)
        finally:
            method._predict_direction=original

    def test_claims_are_fixture_only_and_leads_unevaluated(self):
        report=evaluate_sliding_windows()
        self.assertIn("only the frozen non-locus method fixture",report.claim_boundary)
        self.assertEqual([(x.locus,x.role) for x in report.leads],
            [("STMN2","lead-only-not-evaluated"),("UNC13A","lead-only-not-evaluated")])

if __name__=="__main__": unittest.main()
