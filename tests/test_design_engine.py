import json
import unittest
from dataclasses import replace
from datetime import datetime, timezone
from itertools import permutations

from design_engine import (DesignEngine, DesignReport, GeneratedCandidate,
    GenerationConfig, MagnitudeClaim, ObjectiveDirection, PROSPECTIVE_BANNER,
    PROSPECTIVE_STATUS, ScoredCandidate, TargetDefinition, TherapeuticObjective,
    generate_candidates, offline_target, pareto_front, pareto_tiers, render_json,
    render_text)
from offtarget import EnumeratorConfig, Transcript, enumerate_offtargets
from physics_adapter import (EventDefinition, FoldingResult, Interval, Probability,
                             TranscriptContext, expected_receipt)


FIXED_TIME = datetime(2026, 1, 2, 3, 4, 5, tzinfo=timezone.utc)
HASH_A = "a" * 64
HASH_B = "b" * 64


class AbstainingFolder:
    backend_version = "always-abstain-test/1"

    def productive_probability(self, request):
        probability = Probability(.5, .2, .8)
        return FoldingResult(probability, expected_receipt(request, self.backend_version))


class StatefulFolder:
    backend_version = "stateful-test"

    def __init__(self, shift):
        self.shift = shift

    def productive_probability(self, request):
        value = .5 + (self.shift if request.bound_complex else 0)
        return FoldingResult(Probability(value, value, value),
                             expected_receipt(request, self.backend_version))


def enumerator_a(*args, **kwargs):
    return enumerate_offtargets(*args, **kwargs)


def enumerator_b(*args, **kwargs):
    return enumerate_offtargets(*args, **kwargs)


def scored(name, efficacy, burden, manufacture, robustness, *, abstained=False,
           opposing=False, claim=None):
    candidate = GeneratedCandidate(name, "ASO", "PMO", 0, 18, "A" * 18, "T" * 18)
    alignment = "abstained" if abstained else ("opposes_goal" if opposing else "supports_goal")
    return ScoredCandidate(candidate, "increases_inclusion",
        ObjectiveDirection.INCREASE_INCLUSION, alignment, efficacy, claim, HASH_A, HASH_B,
        abstained, (), burden, 1, "unknown", manufacture, robustness)


class DesignEngineTests(unittest.TestCase):
    def engine(self, **kwargs):
        return DesignEngine(clock=lambda: FIXED_TIME,
            generation_config=GenerationConfig(max_candidates=12), **kwargs)

    def test_generation_is_deterministic_bounded_and_covers_classes_and_chemistries(self):
        target = offline_target("SMN2_exon7_inclusion")
        config = GenerationConfig(stride=7, max_candidates=80)
        first = generate_candidates(target, config)
        self.assertEqual(first, generate_candidates(target, config))
        self.assertLessEqual(len(first), 80)
        self.assertEqual({x.effector_class for x in first}, {"ASO", "dCas13_gRNA", "U7"})
        self.assertTrue({"2MOE", "OMe", "PMO", "LNA"}.issubset({x.chemistry for x in first}))

    def test_typed_therapeutic_objective_fixes_opposite_signs_and_dmd(self):
        smn = self.engine().design_for("SMN2_exon7_inclusion")
        dmd = self.engine().design_for("DMD_exon53_inclusion")
        self.assertEqual(smn.objective.direction, ObjectiveDirection.INCREASE_INCLUSION)
        self.assertEqual(dmd.objective.direction, ObjectiveDirection.DECREASE_INCLUSION)
        self.assertTrue(all(x.therapeutic_goal is ObjectiveDirection.INCREASE_INCLUSION for x in smn.candidates))
        self.assertTrue(all(x.therapeutic_goal is ObjectiveDirection.DECREASE_INCLUSION for x in dmd.candidates))
        with self.assertRaises(ValueError):
            TherapeuticObjective("DMD_exon53_inclusion", "DMD",
                                 ObjectiveDirection.INCREASE_INCLUSION)

    def test_custom_target_requires_override_and_event_locus_consistency(self):
        with self.assertRaises(ValueError):
            TherapeuticObjective("CUSTOM", "GENE", ObjectiveDirection.INCREASE_INCLUSION)
        objective = TherapeuticObjective("CUSTOM", "GENE", ObjectiveDirection.INCREASE_INCLUSION,
                                         override_provenance="IRB protocol X / curator Y")
        context = TranscriptContext("TX", "A" * 100)
        event = EventDefinition("CUSTOM", (Interval(10, 20),))
        TargetDefinition("CUSTOM", "GENE", objective, context, event, 5, 50)
        with self.assertRaises(ValueError):
            TargetDefinition("WRONG", "GENE", objective, context, event, 5, 50)

    def test_all_abstain_has_no_ranks_or_rank_one_and_renders_diagnostics(self):
        report = self.engine(folder=AbstainingFolder()).design_for("SMN2_exon7_inclusion")
        self.assertTrue(all(x.abstained and x.rank is None for x in report.candidates))
        self.assertEqual(report.pareto_candidate_ids, ())
        text = render_text(report)
        self.assertIn("no rank-1 candidate exists", text)
        self.assertIn("UNRANKED DIAGNOSTICS", text)

    def test_pareto_is_permutation_invariant_and_exact_ties_share_tier(self):
        a = scored("a", .9, .5, .7, .7)
        b = scored("b", .8, .2, .8, .8)
        tie = scored("tie", .9, .5, .7, .7)
        dominated = scored("dominated", .7, .6, .6, .6)
        expected = {"a", "b", "tie"}
        for ordering in permutations((a, b, tie, dominated)):
            self.assertEqual(set(pareto_front(ordering)), expected)
            tiers = pareto_tiers(ordering)
            self.assertEqual(tiers["a"], tiers["tie"])
            self.assertEqual(tiers["a"], 1)
            self.assertGreater(tiers["dominated"], 1)

    def test_pareto_rejects_duplicate_ids_and_handles_full_dominance(self):
        best = scored("best", 1, 0, 1, 1)
        weak = scored("weak", .1, 9, .1, .1)
        self.assertEqual(pareto_front((weak, best)), ("best",))
        with self.assertRaises(ValueError):
            pareto_front((best, replace(best)))

    def test_increased_offtarget_burden_worsens_dominance(self):
        low = scored("low", .8, .1, .8, .8)
        high = scored("high", .8, .9, .8, .8)
        self.assertEqual(pareto_front((high, low)), ("low",))
        self.assertGreater(pareto_tiers((high, low))["high"], pareto_tiers((high, low))["low"])

    def test_banner_status_enforced_on_construction_replace_and_render(self):
        report = self.engine().design_for("SMN2_exon7_inclusion")
        for banner in ("", "validated", None):
            with self.assertRaises(ValueError):
                replace(report, banner=banner)
        with self.assertRaises(ValueError):
            replace(report, prediction_status="validated")
        with self.assertRaises(ValueError):
            replace(report, provenance={**report.provenance, "prediction_status": "validated"})
        object.__setattr__(report, "banner", "")
        with self.assertRaises(ValueError):
            render_text(report)
        with self.assertRaises(ValueError):
            render_json(report)

    def test_renderers_emit_canonical_banner_independently(self):
        report = self.engine().design_for("SMN2_exon7_inclusion")
        self.assertTrue(render_text(report).startswith(PROSPECTIVE_BANNER))
        payload = json.loads(render_json(report))
        self.assertEqual(payload["banner"], PROSPECTIVE_BANNER)
        self.assertEqual(payload["prediction_status"], PROSPECTIVE_STATUS)

    def test_magnitude_is_absent_by_default_and_claim_hashes_are_schema_gated(self):
        report = self.engine().design_for("SMN2_exon7_inclusion")
        self.assertTrue(all(x.magnitude_claim is None for x in report.candidates))
        valid = MagnitudeClaim(12.0, "delta-PSI", "assay:plate-7", HASH_A, HASH_B)
        scored("audited", .8, .1, .8, .8, claim=valid)
        with self.assertRaises(ValueError):
            scored("mismatch", .8, .1, .8, .8,
                   claim=replace(valid, input_hash="c" * 64))

    def test_manifest_binds_target_sequence_event_and_transcriptome(self):
        base = offline_target("SMN2_exon7_inclusion")
        digest = self.engine().design_for(base).provenance["run_manifest_sha256"]
        changed_sequence = replace(base, context=replace(base.context,
            sequence="T" + base.context.sequence[1:]))
        changed_event = replace(base, event=replace(base.event,
            silencers=(Interval(55, 66),)))
        changed_corpus = replace(base, offtarget_transcripts=(replace(
            base.offtarget_transcripts[0], transcript_id="OTHER_TX"),))
        for changed in (changed_sequence, changed_event, changed_corpus):
            self.assertNotEqual(digest,
                self.engine().design_for(changed).provenance["run_manifest_sha256"])

    def test_render_recomputes_manifest_digest_and_nested_provenance_is_immutable(self):
        report = self.engine().design_for("SMN2_exon7_inclusion")
        with self.assertRaises(TypeError):
            report.provenance["run_manifest"]["backend"]["qualified_name"] = "forged"
        forged = report.to_dict()["provenance"]
        forged["run_manifest"]["target_definition"]["locus"] = "FORGED"
        forged_report = replace(report, provenance=forged)
        with self.assertRaisesRegex(ValueError, "run manifest content"):
            render_json(forged_report)

    def test_manifest_binds_backend_enumerator_and_configs(self):
        first = self.engine(folder=StatefulFolder(.1), enumerator=enumerator_a).design_for(
            "SMN2_exon7_inclusion")
        backend_changed = self.engine(folder=StatefulFolder(.2), enumerator=enumerator_a).design_for(
            "SMN2_exon7_inclusion")
        enumerator_changed = self.engine(folder=StatefulFolder(.1), enumerator=enumerator_b).design_for(
            "SMN2_exon7_inclusion")
        config_changed = self.engine(folder=StatefulFolder(.1), enumerator=enumerator_a,
            offtarget_config=EnumeratorConfig(max_edits=3, max_indel_bases=1)).design_for(
                "SMN2_exon7_inclusion")
        digests = {x.provenance["run_manifest_sha256"] for x in
                   (first, backend_changed, enumerator_changed, config_changed)}
        self.assertEqual(len(digests), 4)

    def test_full_report_is_deterministic_with_frozen_clock(self):
        first = render_json(self.engine().design_for("DMD_exon53_inclusion"))
        second = render_json(self.engine().design_for("DMD_exon53_inclusion"))
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
