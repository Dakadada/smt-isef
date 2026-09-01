import json
import tempfile
import unittest
from copy import deepcopy
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import patch

from effector_characterization import (CLAIM_BOUNDARY, P1_STATUS,
    build_characterization, deposit, derive_characterization, render_json, render_text)


FIXED_TIME = datetime(2026, 8, 31, 13, 0, tzinfo=timezone.utc)


class EffectorCharacterizationTests(unittest.TestCase):
    def setUp(self):
        self.report = build_characterization(clock=lambda: FIXED_TIME)

    def test_honesty_boundary_and_p1_gate_are_structural_in_both_outputs(self):
        payload = json.loads(render_json(self.report))
        self.assertEqual(payload["claim_boundary"], CLAIM_BOUNDARY)
        self.assertEqual(payload["p1_status"], P1_STATUS)
        self.assertFalse(payload["provenance"]["p1_activation_met"])
        text = render_text(self.report)
        for phrase in ("ENGINE-INTERNAL", "NOT biological effector-class effects",
                       "P1 remains gated OFF"):
            self.assertIn(phrase, text)

    def test_actual_committed_tiling_and_component_scores(self):
        self.assertEqual(self.report.context_length, 21)
        self.assertEqual(len(self.report.candidates), 10)
        self.assertEqual({row.effector_class for row in self.report.candidates},
                         {"ASO", "U7"})
        observed = {(row.effector_class, row.chemistry, row.length, row.engine_score)
                    for row in self.report.candidates}
        for length in (18, 20):
            self.assertIn(("ASO", "PMO", length, .90), observed)
            self.assertIn(("ASO", "2MOE", length, .88), observed)
            self.assertIn(("ASO", "OMe", length, .86), observed)
            self.assertIn(("ASO", "LNA", length, .80), observed)
            self.assertIn(("U7", "U7 snRNA cassette", length, .67), observed)
        self.assertTrue(all(row.length_penalty == 0 for row in self.report.candidates))

    def test_class_distributions_include_absent_dcas_without_fabricating_values(self):
        summaries = {row.effector_class: row for row in self.report.class_summaries}
        self.assertEqual(summaries["ASO"].candidate_count, 8)
        self.assertEqual((summaries["ASO"].minimum, summaries["ASO"].maximum,
                          summaries["ASO"].mean), (.8, .9, .86))
        self.assertEqual(summaries["U7"].candidate_count, 2)
        dcas = summaries["dCas13_gRNA"]
        self.assertEqual(dcas.candidate_count, 0)
        self.assertEqual(dcas.observed_scores, ())
        self.assertIsNone(dcas.mean)
        self.assertIn("21 nt is shorter", dcas.limitation)
        self.assertIn("4-nt gRNA shadow", dcas.limitation)

    def test_configured_formula_records_base_deltas_penalty_and_shadow(self):
        self.assertEqual(dict(self.report.class_bases),
                         {"ASO": .82, "U7": .67, "dCas13_gRNA": .58})
        self.assertEqual(dict(self.report.aso_chemistry_deltas),
                         {"PMO": .08, "2MOE": .06, "OMe": .04, "LNA": -.02})
        self.assertEqual(self.report.length_penalty_threshold, 20)
        self.assertEqual(self.report.length_penalty_per_nt, .01)
        self.assertEqual(self.report.grna_shadow_nt, 4)
        self.assertIn("changes the physics footprint, not the manufacturability scalar",
                      render_text(self.report))
        self.assertEqual(self.report.provenance["tiling_start_positions"], (0,))
        self.assertIn("not independent positions or biological replicates",
                      self.report.provenance["tiling_disclosure"])

    def test_provenance_binds_authenticated_context_and_actual_engine_run(self):
        provenance = self.report.provenance
        for key in ("engine_report_sha256", "engine_run_manifest_sha256",
                    "manufacturability_function_source_sha256",
                    "score_rule_probe_sha256", "design_engine_module_sha256",
                    "approved_derivation_source_sha256",
                    "approved_derivation_module_sha256",
                    "binding_sequence_sha256",
                    "mapping_sha256", "characterization_manifest_sha256"):
            self.assertRegex(provenance[key], r"^[0-9a-f]{64}$")
        self.assertEqual(provenance["candidate_count"], len(self.report.candidates))

    def test_renderer_rejects_forged_candidate_and_self_consistent_summary(self):
        row = self.report.candidates[0]
        object.__setattr__(row, "class_base", 0.81)
        object.__setattr__(row, "engine_score", 0.89)
        summary = next(x for x in self.report.class_summaries
                       if x.effector_class == row.effector_class)
        object.__setattr__(summary, "observed_scores", tuple(
            0.89 if value == 0.9 else value for value in summary.observed_scores))
        object.__setattr__(summary, "maximum", 0.89)
        object.__setattr__(summary, "mean",
                           round(sum(summary.observed_scores) / len(summary.observed_scores), 6))
        with self.assertRaisesRegex(ValueError, "deterministic authenticated SMN2"):
            render_json(self.report)

    def test_renderer_rejects_forged_absent_class_claim(self):
        dcas = next(row for row in self.report.class_summaries
                    if row.effector_class == "dCas13_gRNA")
        object.__setattr__(dcas, "availability", "observed-in-engine-output")
        object.__setattr__(dcas, "candidate_count", 1)
        object.__setattr__(dcas, "observed_scores", (.58,))
        object.__setattr__(dcas, "minimum", .58)
        object.__setattr__(dcas, "maximum", .58)
        object.__setattr__(dcas, "mean", .58)
        with self.assertRaisesRegex(ValueError, "deterministic authenticated SMN2"):
            render_text(self.report)

    def test_deposit_is_immutable_and_deterministic(self):
        second = build_characterization(clock=lambda: FIXED_TIME)
        self.assertEqual(render_json(self.report), render_json(second))
        with tempfile.TemporaryDirectory() as directory:
            paths = deposit(self.report, Path(directory))
            self.assertTrue(all(path.exists() for path in paths))
            self.assertTrue(all("T130000Z" in path.name for path in paths))
            with self.assertRaises(FileExistsError):
                deposit(self.report, Path(directory))

    def test_public_builder_derivation_substitution_cannot_authorize_construction(self):
        forged = deepcopy(derive_characterization(FIXED_TIME))
        forged["candidates"][0]["class_base"] = .81
        forged["candidates"][0]["engine_score"] = .87
        distribution = next(row for row in forged["distributions"]
            if row["effector_class"] == "ASO" and row["chemistry"] == "2MOE"
            and row["length"] == 18)
        distribution.update(observed_scores=(.87,), minimum=.87, maximum=.87, mean=.87)
        summary = next(row for row in forged["class_summaries"]
                       if row["effector_class"] == "ASO")
        values = list(summary["observed_scores"])
        values[values.index(.88)] = .87
        summary["observed_scores"] = tuple(sorted(values))
        summary["minimum"] = min(values)
        summary["maximum"] = max(values)
        summary["mean"] = round(sum(values) / len(values), 6)

        with patch("effector_characterization.builder._APPROVED_DERIVATION",
                   return_value=forged), self.assertRaisesRegex(
                       ValueError, "deterministic authenticated SMN2"):
            build_characterization(clock=lambda: FIXED_TIME)

    def test_runtime_engine_context_and_validation_seam_patches_do_not_substitute_truth(self):
        baseline = render_json(self.report)

        def forbidden(*args, **kwargs):
            raise AssertionError("mutable substituted seam was called")

        with patch("effector_characterization.builder.derive_characterization", forbidden), \
             patch("effector_characterization.approved.design_for", forbidden), \
             patch("effector_characterization.approved.committed_real_contexts", forbidden), \
             patch("design_engine.engine.design_for", forbidden), \
             patch("design_panel.panel.committed_real_contexts", forbidden), \
             patch("effector_characterization.models._APPROVED_DERIVATION", forbidden):
            self.assertEqual(render_json(self.report), baseline)

    def test_substituted_seams_cannot_make_forged_report_render_or_deposit(self):
        row = self.report.candidates[0]
        object.__setattr__(row, "class_base", .81)
        object.__setattr__(row, "engine_score", .87)

        def self_consistent(*args, **kwargs):
            return self.report.to_dict()

        patches = (
            patch("effector_characterization.builder.derive_characterization", self_consistent),
            patch("effector_characterization.approved.design_for", self_consistent),
            patch("effector_characterization.approved.committed_real_contexts", self_consistent),
            patch("effector_characterization.models._APPROVED_DERIVATION", self_consistent),
        )
        with patches[0], patches[1], patches[2], patches[3]:
            with self.assertRaisesRegex(ValueError, "deterministic authenticated SMN2"):
                render_json(self.report)
            with tempfile.TemporaryDirectory() as directory:
                with self.assertRaisesRegex(ValueError,
                                            "deterministic authenticated SMN2"):
                    deposit(self.report, Path(directory))
                self.assertEqual(list(Path(directory).iterdir()), [])

    def test_naive_clock_fails_closed(self):
        with self.assertRaisesRegex(ValueError, "timezone-aware"):
            build_characterization(clock=lambda: datetime(2026, 8, 31, 13, 0))


if __name__ == "__main__":
    unittest.main()
