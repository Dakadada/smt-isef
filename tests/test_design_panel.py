import csv
import hashlib
import json
import tempfile
import unittest
from contextlib import contextmanager
from dataclasses import asdict, replace
from datetime import datetime, timezone
from pathlib import Path
from unittest.mock import patch

from design_engine import (DesignEngine, GenerationConfig, ObjectiveDirection,
                           TargetDefinition, TherapeuticObjective, offline_target)
from design_engine.engine import _hash as _engine_hash, _implementation_identity
from design_panel import (PENDING_CONTEXT_STATUS, DesignabilityState,
    PENDING_UNCOMMITTED_CONTEXT_STATUS, PanelReport, bind_real_transcript_context,
    deposit, load_panel, registry_provenance_digest, render_json, render_text, run_panel)
from physics_adapter import (EventDefinition, FoldingResult, Interval, Probability,
                             TranscriptContext, expected_receipt)
import design_panel.panel as panel_module
from design_panel.panel import _hash, _pending, _summary

FIXED_TIME = datetime(2026, 2, 3, 4, 5, 6, tzinfo=timezone.utc)


class NetworkedTestFolder:
    backend_version = "networked-test-model/1"
    def productive_probability(self, request):
        base = .52 - (.02 if request.max_bp_span is not None else 0)
        if request.bound_complex:
            target = request.bound_complex.transcript_coordinates
            base += (.22 if any(target.overlaps(x) for x in request.required_silencer_free) else 0)
            base -= (.22 if any(target.overlaps(x) for x in request.required_accessible) else 0)
        return FoldingResult(Probability(base, base-.025, base+.025),
                             expected_receipt(request, self.backend_version))


class DesignPanelTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.registry_path = Path(self.temp.name) / "registry.tsv"
        with Path("config/event_junctions.GRCh38.tsv").open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        for row in rows:
            row["sequence_sha256"] = self.context(row["event"]).binding_sequence_sha256
        with self.registry_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
            writer.writeheader(); writer.writerows(rows)
        registry_hash = hashlib.sha256(self.registry_path.read_bytes()).hexdigest()
        config = json.loads(Path("src/design_panel/target_panel.json").read_text())
        by_event = {row["event"]: row for row in rows}
        for target in config["targets"]:
            if target["coordinate_status"] != "verified": continue
            row = by_event[target["event_id"]]
            context = self.context(row["event"])
            if row["sequence_sha256"]:
                target["mapping_sha256"] = context.mapping_sha256
            target["coordinate_provenance"] = (
                f"event_junctions.GRCh38.tsv:sha256:{registry_hash}:{target['event_id']}")
            target["provenance_digest"] = registry_provenance_digest(
                row, target["transcript_version"])
        self.config_path = Path(self.temp.name) / "panel.json"
        self.config_path.write_text(json.dumps(config))

    def tearDown(self):
        self.temp.cleanup()

    def run_test_panel(self, **kwargs):
        registry_hash = hashlib.sha256(self.registry_path.read_bytes()).hexdigest()
        panel_hash = hashlib.sha256(self.config_path.read_bytes()).hexdigest()
        design_api = kwargs["design_api"]
        engine_api = f"{design_api.__module__}.{getattr(design_api, '__qualname__', design_api.__class__.__qualname__)}"
        backend_identity = _implementation_identity(NetworkedTestFolder())["identity_sha256"]
        with patch("design_panel.panel.DEFAULT_COORDINATE_REGISTRY", self.registry_path), \
             patch("design_panel.panel.EXPECTED_REGISTRY_SHA256", registry_hash), \
             patch("design_panel.panel.DEFAULT_PANEL_CONFIG", self.config_path), \
             patch("design_panel.panel.EXPECTED_PANEL_SHA256", panel_hash), \
             patch("design_panel.panel.APPROVED_ENGINE_API", engine_api), \
             patch("design_panel.panel.APPROVED_BACKEND_IDENTITY_SHA256", backend_identity), \
             patch("design_panel.panel.APPROVED_DESIGN_ENGINE", design_api), \
             patch("design_panel.panel.committed_real_contexts",
                   lambda: {"SMN2_exon7_inclusion": self.context(),
                            "DMD_exon53_inclusion": self.context("DMD_exon53_inclusion")}):
            return run_panel(config_path=self.config_path, registry_path=self.registry_path,
                             _test_only_allow_unapproved_engine=True,
                             **kwargs)

    @contextmanager
    def frozen_test_anchor(self, report):
        engine_api = report.provenance["panel_manifest"]["engine_api"]
        backend_identity = (_implementation_identity(NetworkedTestFolder())["identity_sha256"]
                            if not report.designs else
                            report.designs[0].provenance["folding_backend_identity_sha256"])
        with patch("design_panel.panel.DEFAULT_COORDINATE_REGISTRY", self.registry_path), \
             patch("design_panel.panel.EXPECTED_REGISTRY_SHA256", hashlib.sha256(self.registry_path.read_bytes()).hexdigest()), \
             patch("design_panel.panel.DEFAULT_PANEL_CONFIG", self.config_path), \
             patch("design_panel.panel.EXPECTED_PANEL_SHA256", hashlib.sha256(self.config_path.read_bytes()).hexdigest()), \
             patch("design_panel.panel.APPROVED_ENGINE_API", engine_api), \
             patch("design_panel.panel.APPROVED_BACKEND_IDENTITY_SHA256", backend_identity), \
             patch("design_panel.panel.APPROVED_DESIGN_ENGINE", self.api()), \
             patch("design_panel.panel.committed_real_contexts",
                   lambda: {"SMN2_exon7_inclusion": self.context(),
                            "DMD_exon53_inclusion": self.context("DMD_exon53_inclusion")}):
            yield

    def rebind_candidates(self, report, candidates, pareto_candidate_ids=None):
        pareto_ids = (report.pareto_candidate_ids if pareto_candidate_ids is None
                      else pareto_candidate_ids)
        provenance = report.to_dict()["provenance"]
        provenance["run_manifest"]["scored_candidates"] = [asdict(x) for x in candidates]
        provenance["run_manifest"]["pareto_candidate_ids"] = list(pareto_ids)
        provenance["run_manifest_sha256"] = _engine_hash(provenance["run_manifest"])
        return replace(report, candidates=tuple(candidates), pareto_candidate_ids=tuple(pareto_ids),
                       provenance=provenance)

    def load_test_panel(self, path):
        registry_hash = hashlib.sha256(self.registry_path.read_bytes()).hexdigest()
        with patch("design_panel.panel.DEFAULT_COORDINATE_REGISTRY", self.registry_path), \
             patch("design_panel.panel.EXPECTED_REGISTRY_SHA256", registry_hash):
            return load_panel(path, self.registry_path, _test_only_allow_unfrozen_config=True)

    def definition(self, target_id="SMN2_exon7_inclusion"):
        locus, transcript, direction = (("SMN2", "NM_017411.4", ObjectiveDirection.INCREASE_INCLUSION)
            if target_id.startswith("SMN2") else ("DMD", "NM_004006.3", ObjectiveDirection.DECREASE_INCLUSION))
        sequence = ("GATTACAACCTGAGTCCGTA" * 7)[:130]
        objective = TherapeuticObjective(target_id, locus, direction)
        event = EventDefinition(target_id, (Interval(28, 38), Interval(72, 82)),
                                (Interval(44, 50),), (Interval(54, 66),))
        return TargetDefinition(target_id, locus, objective, TranscriptContext(transcript, sequence), event, 20, 90)

    def context(self, target_id="SMN2_exon7_inclusion"):
        definition = self.definition(target_id)
        return bind_real_transcript_context(definition=definition, assembly="GRCh38",
            chrom="chr5" if target_id.startswith("SMN2") else "chrX",
            strand="+" if target_id.startswith("SMN2") else "-",
            transcript_version=definition.context.transcript_id, sequence_provenance="NCBI retrieval receipt sha256:real",
            genomic_to_transcript_mapping=tuple((100 + i if target_id.startswith("SMN2") else 300 - i, i)
                                                for i in range(len(definition.context.sequence))))

    def api(self, mutate=None):
        def call(target, **options):
            report = DesignEngine(folder=NetworkedTestFolder(), clock=options["clock"],
                generation_config=GenerationConfig(max_candidates=20)).design_for(target)
            return mutate(report) if mutate else report
        return call

    def test_default_emits_real_smn2_design_and_leaves_dmd_pending(self):
        report = run_panel(clock=lambda: FIXED_TIME)
        self.assertEqual([x.target_id for x in report.designs], ["SMN2_exon7_inclusion"])
        definition = report.designs[0].provenance["run_manifest"]["target_definition"]
        self.assertEqual(definition["context"]["transcript_id"], "NM_017411.4")
        self.assertEqual(definition["context"]["sequence"], "CCAGCATTATGAAAGTGAATC")
        self.assertNotIn("OFFLINE_FIXTURE", str(definition))
        self.assertIn("DMD_exon53_inclusion", {x.target_id for x in report.pending_targets
                      if x.status == PENDING_UNCOMMITTED_CONTEXT_STATUS})

    def test_real_definition_not_string_is_designed_and_manifest_is_reproducible(self):
        seen = []
        def api(target, **options):
            seen.append(target); return self.api()(target, **options)
        contexts = {"SMN2_exon7_inclusion": self.context()}
        first = self.run_test_panel(clock=lambda: FIXED_TIME, design_api=api, target_contexts=contexts)
        second = self.run_test_panel(clock=lambda: FIXED_TIME, design_api=api, target_contexts=contexts)
        self.assertTrue(all(isinstance(x, TargetDefinition) for x in seen))
        with self.frozen_test_anchor(first):
            self.assertEqual(render_json(first), render_json(second))
        self.assertEqual(first.provenance["panel_manifest"]["coordinate_registry_sha256"],
                         hashlib.sha256(self.registry_path.read_bytes()).hexdigest())
        bound = first.provenance["panel_manifest"]["bound_target_definitions"][0]
        self.assertEqual(bound["genomic_to_transcript_mapping"][0], (100, 0))
        self.assertEqual(bound["genomic_to_transcript_mapping"][-1], (229, 129))

    def test_renamed_or_fixture_context_cannot_bypass_committed_sequence_digest(self):
        fixture = offline_target("SMN2_exon7_inclusion")
        fixture = replace(fixture, context=replace(fixture.context, transcript_id="NM_017411.4"))
        asserted = bind_real_transcript_context(definition=fixture, assembly="GRCh38",
            chrom="chr5", strand="+",
            transcript_version="NM_017411.4", sequence_provenance="renamed synthetic input",
            genomic_to_transcript_mapping=tuple((1 + i, i) for i in range(len(fixture.context.sequence))))
        with self.assertRaisesRegex(ValueError, "committed registry digest"):
            self.run_test_panel(clock=lambda: FIXED_TIME,
                target_contexts={"SMN2_exon7_inclusion": asserted}, design_api=self.api())

    def test_self_asserted_verified_absent_from_registry_is_rejected(self):
        config = json.loads(Path("src/design_panel/target_panel.json").read_text())
        config["targets"][0]["target_id"] = config["targets"][0]["event_id"] = "invented"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"panel.json"; path.write_text(json.dumps(config))
            with self.assertRaises(ValueError): load_panel(path)

    def test_registry_override_and_tampering_fail_closed_without_test_flag(self):
        with self.assertRaisesRegex(ValueError, "frozen project registry"):
            load_panel(self.config_path, self.registry_path)
        original = Path("config/event_junctions.GRCh38.tsv").read_text()
        tampered = Path(self.temp.name) / "tampered.tsv"
        tampered.write_text(original + "\n")
        with self.assertRaisesRegex(ValueError, "frozen project registry"):
            load_panel(registry_path=tampered)
        with patch("design_panel.panel.EXPECTED_REGISTRY_SHA256", "0" * 64):
            with self.assertRaisesRegex(ValueError, "registry digest mismatch"):
                load_panel()

    def test_identity_polarity_and_clock_mismatches_are_rejected(self):
        contexts = {"SMN2_exon7_inclusion": self.context()}
        wrong_time = lambda r: replace(r, provenance={**r.provenance, "generated_at": "2020-01-01T00:00:00+00:00"})
        with self.assertRaisesRegex(ValueError, "frozen panel clock"):
            self.run_test_panel(clock=lambda: FIXED_TIME, target_contexts=contexts, design_api=self.api(wrong_time))
        wrong_identity = lambda r: replace(r, target_id="DMD_exon53_inclusion")
        with self.assertRaisesRegex(ValueError, "identity/polarity"):
            self.run_test_panel(clock=lambda: FIXED_TIME, target_contexts=contexts, design_api=self.api(wrong_identity))
        dmd_context = {"DMD_exon53_inclusion": self.context("DMD_exon53_inclusion")}
        polarity_mutation = lambda r: object.__setattr__(r.objective, "direction", ObjectiveDirection.INCREASE_INCLUSION) or r
        with self.assertRaises(ValueError):
            self.run_test_panel(clock=lambda: FIXED_TIME, target_contexts=dmd_context, design_api=self.api(polarity_mutation))

    def test_panel_rejects_nested_design_with_forged_manifest_content(self):
        def forge(report):
            provenance = report.to_dict()["provenance"]
            provenance["run_manifest"]["target_definition"]["locus"] = "FORGED"
            return replace(report, provenance=provenance)
        with self.assertRaisesRegex(ValueError, "run manifest content"):
            self.run_test_panel(clock=lambda: FIXED_TIME,
                target_contexts={"SMN2_exon7_inclusion": self.context()},
                design_api=self.api(forge))

    def test_report_must_be_bound_to_exact_supplied_sequence(self):
        committed = self.context()
        other_definition = replace(committed.definition, context=replace(
            committed.definition.context,
            sequence="C" + committed.definition.context.sequence[1:]))
        substituted = self.api()(other_definition, clock=lambda: FIXED_TIME)
        with self.assertRaisesRegex(ValueError, "bound to supplied real context"):
            self.run_test_panel(clock=lambda: FIXED_TIME,
                target_contexts={"SMN2_exon7_inclusion": committed},
                design_api=lambda *a, **k: substituted)

    def test_panel_rejects_injected_design_and_incoherent_designability(self):
        context = {"SMN2_exon7_inclusion": self.context()}
        panel = self.run_test_panel(clock=lambda: FIXED_TIME,
                                    target_contexts=context, design_api=self.api())
        with self.frozen_test_anchor(panel):
            with self.assertRaisesRegex(ValueError, "report manifest list"):
                replace(panel, designs=panel.designs + (panel.designs[0],))
        pending_index = next(i for i, row in enumerate(panel.designability_map)
                             if row.state is DesignabilityState.PENDING)
        rows = list(panel.designability_map)
        rows[pending_index] = replace(rows[pending_index], state=DesignabilityState.DESIGNABLE)
        with self.frozen_test_anchor(panel):
            with self.assertRaisesRegex(ValueError, "pending targets cannot be marked designable"):
                replace(panel, designability_map=tuple(rows))

    def test_tampered_panel_binding_and_wrong_assembly_fail_closed(self):
        config = json.loads(self.config_path.read_text())
        config["targets"][0]["provenance_digest"] = "0" * 64
        bad_digest = Path(self.temp.name) / "bad-digest.json"
        bad_digest.write_text(json.dumps(config))
        with self.assertRaisesRegex(ValueError, "not exactly resolved"):
            self.load_test_panel(bad_digest)

        config["assembly"] = "GRCh37"
        bad_assembly_config = Path(self.temp.name) / "bad-assembly.json"
        bad_assembly_config.write_text(json.dumps(config))
        with self.assertRaisesRegex(ValueError, "schema/assembly"):
            self.load_test_panel(bad_assembly_config)

        wrong_assembly = replace(self.context(), assembly="GRCh37")
        with self.assertRaisesRegex(ValueError, "identity mismatch"):
            self.run_test_panel(clock=lambda: FIXED_TIME,
                target_contexts={"SMN2_exon7_inclusion": wrong_assembly},
                design_api=self.api())

    def test_honesty_is_immutable_and_revalidated_by_render_and_deposit(self):
        report = run_panel(clock=lambda: FIXED_TIME)
        with self.assertRaises(TypeError): report.provenance["prediction_status"] = "validated"
        object.__setattr__(report, "banner", "VALIDATED")
        with self.assertRaises(ValueError): render_json(report)
        with self.assertRaises(ValueError): render_text(report)
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaises(ValueError): deposit(report, Path(directory))

    def test_designability_thresholds_and_tied_best_candidates(self):
        context = {"SMN2_exon7_inclusion": self.context()}
        baseline = self.api()(self.definition(), clock=lambda: FIXED_TIME)
        candidates = list(baseline.candidates)
        supported = [x for x in candidates if x.rank is not None]
        if len(supported) >= 2:
            tier = supported[0].rank
            candidates[candidates.index(supported[1])] = replace(supported[1], rank=tier)
            tied_report = self.rebind_candidates(baseline, candidates)
            panel = self.run_test_panel(clock=lambda: FIXED_TIME, target_contexts=context, design_api=lambda *a, **k: tied_report)
            row = next(x for x in panel.designability_map if x.target_id.startswith("SMN2"))
            self.assertGreaterEqual(len(row.best_candidate_ids), 2)
        abstained = tuple(replace(x, abstained=True, rank=None, pareto=False,
                                  pareto_tier=None, goal_alignment="abstained") for x in baseline.candidates)
        no_go = self.rebind_candidates(baseline, abstained, ())
        panel = self.run_test_panel(clock=lambda: FIXED_TIME, target_contexts=context, design_api=lambda *a, **k: no_go)
        row = next(x for x in panel.designability_map if x.target_id.startswith("SMN2"))
        self.assertEqual(row.state, DesignabilityState.NON_DESIGNABLE)
        self.assertEqual(row.best_candidate_ids, ())

        # Each trust threshold independently prevents a positive designability claim.
        support = next(x for x in baseline.candidates if x.rank is not None)
        cases = {
            "mostly-abstained": tuple(
                x if i == 0 else replace(x, abstained=True, rank=None, pareto=False,
                                         pareto_tier=None, goal_alignment="abstained")
                for i, x in enumerate(baseline.candidates)),
            "high-burden": tuple(replace(x, off_target_burden=2.0) for x in baseline.candidates),
            "no-goal-support": tuple(replace(x, rank=None, pareto=False, pareto_tier=None,
                                              goal_alignment="opposes_goal")
                                     for x in baseline.candidates),
        }
        for name, rows in cases.items():
            pareto_ids = tuple(x.candidate.candidate_id for x in rows if x.pareto)
            case_report = self.rebind_candidates(baseline, rows, pareto_ids)
            with self.subTest(name=name):
                panel = self.run_test_panel(clock=lambda: FIXED_TIME, target_contexts=context,
                                  design_api=lambda *a, _report=case_report, **k: _report)
                row = next(x for x in panel.designability_map if x.target_id.startswith("SMN2"))
                self.assertEqual(row.state, DesignabilityState.NON_DESIGNABLE)
                self.assertEqual(row.best_candidate_ids, ())

    def test_deposit_and_render_include_honesty_goal_and_state(self):
        report = run_panel(clock=lambda: FIXED_TIME)
        text = render_text(report)
        self.assertIn("NOT EXPERIMENTALLY VALIDATED", text)
        self.assertIn("goal | state", text)
        with tempfile.TemporaryDirectory() as directory:
            deposit(report, Path(directory))
            with self.assertRaises(FileExistsError): deposit(report, Path(directory))

    def test_production_rejects_unapproved_engine_and_options(self):
        with self.assertRaisesRegex(ValueError, "approved immutable engine manifest"):
            run_panel(clock=lambda: FIXED_TIME, design_api=self.api())
        with self.assertRaisesRegex(ValueError, "approved immutable engine manifest"):
            run_panel(clock=lambda: FIXED_TIME, design_options={"generation_config": GenerationConfig()})

    def test_production_rejects_custom_or_tampered_panel_config(self):
        with self.assertRaisesRegex(ValueError, "frozen project panel"):
            load_panel(self.config_path)
        with patch("design_panel.panel.DEFAULT_PANEL_CONFIG", self.config_path), \
             patch("design_panel.panel.EXPECTED_PANEL_SHA256", "0" * 64):
            with self.assertRaisesRegex(ValueError, "panel digest mismatch"):
                load_panel(self.config_path)

    def test_mapping_requires_independent_digest_identity_and_full_coverage(self):
        definition = self.definition()
        with self.assertRaisesRegex(ValueError, "exactly and directionally cover"):
            bind_real_transcript_context(definition=definition, assembly="GRCh38", chrom="chr5",
                strand="+", transcript_version="NM_017411.4", sequence_provenance="receipt",
                genomic_to_transcript_mapping=((0, 0),))
        context = replace(self.context(), mapping_sha256="0" * 64)
        with self.assertRaisesRegex(ValueError, "context hash mismatch"):
            self.run_test_panel(clock=lambda: FIXED_TIME,
                target_contexts={"SMN2_exon7_inclusion": context}, design_api=self.api())

    def test_every_designability_field_is_recomputed(self):
        panel = self.run_test_panel(clock=lambda: FIXED_TIME,
            target_contexts={"SMN2_exon7_inclusion": self.context()}, design_api=self.api())
        rows = list(panel.designability_map)
        rows[0] = replace(rows[0], locus="FORGED", best_candidate_ids=("forged",),
                          pareto_candidate_count=999)
        with self.frozen_test_anchor(panel):
            with self.assertRaisesRegex(ValueError, "deterministic panel outputs"):
                replace(panel, designability_map=tuple(rows))

    def test_report_manifest_entries_are_target_bound_and_positional(self):
        panel = self.run_test_panel(clock=lambda: FIXED_TIME,
            target_contexts={"SMN2_exon7_inclusion": self.context()}, design_api=self.api())
        with self.frozen_test_anchor(panel):
            provenance = panel.to_dict()["provenance"]
        provenance["panel_manifest"]["report_manifests"][0]["target_id"] = "DMD_exon53_inclusion"
        digest = _hash(provenance["panel_manifest"])
        provenance["panel_manifest_sha256"] = digest
        with self.frozen_test_anchor(panel):
            with self.assertRaisesRegex(ValueError, "report manifest list"):
                replace(panel, provenance=provenance, panel_manifest_sha256=digest)

    def test_forged_bound_ids_cannot_upgrade_pending_target_to_designable(self):
        pending_panel = self.run_test_panel(clock=lambda: FIXED_TIME,
            target_contexts={}, design_api=self.api())
        designed_panel = self.run_test_panel(clock=lambda: FIXED_TIME,
            target_contexts={"SMN2_exon7_inclusion": self.context()}, design_api=self.api())
        with self.frozen_test_anchor(pending_panel):
            provenance = pending_panel.to_dict()["provenance"]
            designed_provenance = designed_panel.to_dict()["provenance"]
        provenance["panel_manifest"]["bound_target_definitions"] = (
            designed_provenance["panel_manifest"]["bound_target_definitions"])
        digest = _hash(provenance["panel_manifest"])
        provenance["panel_manifest_sha256"] = digest
        with self.frozen_test_anchor(pending_panel):
            metadata_only_forge = replace(pending_panel, provenance=provenance,
                                          panel_manifest_sha256=digest)
        smn2 = next(row for row in metadata_only_forge.designability_map
                    if row.target_id == "SMN2_exon7_inclusion")
        self.assertIs(smn2.state, DesignabilityState.PENDING)
        rows = tuple(replace(row, state=DesignabilityState.DESIGNABLE,
                             status="forged designable claim")
                     if row.target_id == "SMN2_exon7_inclusion" else row
                     for row in metadata_only_forge.designability_map)
        with self.frozen_test_anchor(metadata_only_forge):
            with self.assertRaisesRegex(ValueError, "pending targets cannot be marked designable"):
                replace(metadata_only_forge, designability_map=rows)

    def test_direct_design_injection_requires_frozen_sequence_evidence(self):
        base = DesignEngine(clock=lambda: FIXED_TIME).design_for(
            offline_target("SMN2_exon7_inclusion"))
        candidates = tuple(replace(row, abstained=False, goal_alignment="supports_goal",
            rank=1, pareto=True, pareto_tier=1, off_target_burden=0.0,
            directional_confidence=1.0) for row in base.candidates)
        pareto_ids = tuple(row.candidate.candidate_id for row in candidates)
        forged = self.rebind_candidates(base, candidates, pareto_ids)
        targets = load_panel()
        pending = []
        rows = []
        for target in targets:
            if target.target_id == forged.target_id:
                rows.append(_summary(target, forged))
            else:
                status = (PENDING_UNCOMMITTED_CONTEXT_STATUS if not target.sequence_sha256
                          else PENDING_CONTEXT_STATUS)
                item, row = _pending(target, status)
                pending.append(item)
                rows.append(row)
        manifest = {
            "panel_version": panel_module.PANEL_VERSION,
            "engine_api": panel_module.APPROVED_ENGINE_API,
            "options": {},
            "ordered_targets": [{**asdict(target), "expected_goal":
                panel_module.EXPECTED_GOALS.get(
                    target.event_id or target.target_id,
                    ObjectiveDirection.INCREASE_INCLUSION).value}
                for target in targets],
            "bound_target_definitions": [],
            "coordinate_registry_sha256": panel_module.EXPECTED_REGISTRY_SHA256,
            "report_manifests": [{"target_id": forged.target_id,
                "run_manifest_sha256": forged.provenance["run_manifest_sha256"]}],
            "generated_at": FIXED_TIME.isoformat(),
        }
        digest = _hash(manifest)
        provenance = {"panel_manifest": manifest, "panel_manifest_sha256": digest,
            "prediction_status": panel_module.PANEL_PREDICTION_STATUS,
            "coordinate_policy": "forged"}
        with self.assertRaisesRegex(
                ValueError, "sequence does not match frozen committed digest"):
            PanelReport("prospective-design-panel-report/2", panel_module.PANEL_BANNER,
                panel_module.PANEL_PREDICTION_STATUS, FIXED_TIME.isoformat(), (forged,), tuple(pending),
                tuple(rows), provenance, digest)

    def test_render_revalidates_nested_backend_and_identity(self):
        panel = run_panel(clock=lambda: FIXED_TIME)
        design = panel.designs[0]

        def replace_design(forged):
            provenance = panel.to_dict()["provenance"]
            provenance["panel_manifest"]["report_manifests"][0]["run_manifest_sha256"] = (
                forged.provenance["run_manifest_sha256"])
            digest = _hash(provenance["panel_manifest"])
            provenance["panel_manifest_sha256"] = digest
            return replace(panel, designs=(forged,), provenance=provenance,
                           panel_manifest_sha256=digest)

        provenance = design.to_dict()["provenance"]
        provenance["run_manifest"]["backend"]["identity_sha256"] = "0" * 64
        provenance["run_manifest_sha256"] = _engine_hash(provenance["run_manifest"])
        backend_forge = replace(design, provenance=provenance)
        with self.assertRaisesRegex(ValueError, "lacks approved backend receipt"):
            replace_design(backend_forge)

        identity_forge = replace(design, locus="FORGED")
        with self.assertRaisesRegex(ValueError, "identity/polarity mismatch"):
            replace_design(identity_forge)

    def test_panel_manifest_policy_metadata_is_revalidated(self):
        panel = run_panel(clock=lambda: FIXED_TIME)
        mutations = {
            "panel_version": "forged-panel-version",
            "engine_api": "forged.engine",
            "options": {"forged": True},
            "coordinate_registry_sha256": "0" * 64,
        }
        for field, value in mutations.items():
            with self.subTest(field=field):
                provenance = panel.to_dict()["provenance"]
                provenance["panel_manifest"][field] = value
                digest = _hash(provenance["panel_manifest"])
                provenance["panel_manifest_sha256"] = digest
                with self.assertRaisesRegex(
                        ValueError, "manifest metadata does not match frozen policy"):
                    replace(panel, provenance=provenance, panel_manifest_sha256=digest)

    def test_panel_and_nested_timestamps_must_match_canonical_utc(self):
        panel = run_panel(clock=lambda: FIXED_TIME)
        with self.assertRaisesRegex(ValueError, "timestamps must be identical"):
            replace(panel, generated_at="2026-02-03T05:05:06+00:00")
        with self.assertRaisesRegex(ValueError, "canonical UTC"):
            replace(panel, generated_at="2026-02-03T04:05:06Z")

    def test_pending_records_are_exactly_reconstructed(self):
        panel = run_panel(clock=lambda: FIXED_TIME)
        pending = list(panel.pending_targets)
        pending[0] = replace(pending[0], locus="FORGED", status="forged")
        with self.assertRaisesRegex(ValueError, "pending target records"):
            replace(panel, pending_targets=tuple(pending))

    def test_forged_self_validating_panel_manifest_cannot_replace_frozen_truth(self):
        panel = run_panel(clock=lambda: FIXED_TIME)
        provenance = panel.to_dict()["provenance"]
        dmd = next(x for x in provenance["panel_manifest"]["ordered_targets"]
                   if x["target_id"] == "DMD_exon53_inclusion")
        dmd["locus"] = "FORGED"
        pending = tuple(replace(x, locus="FORGED") if x.target_id == dmd["target_id"] else x
                        for x in panel.pending_targets)
        rows = tuple(replace(x, locus="FORGED") if x.target_id == dmd["target_id"] else x
                     for x in panel.designability_map)
        digest = _hash(provenance["panel_manifest"])
        provenance["panel_manifest_sha256"] = digest
        with self.assertRaisesRegex(ValueError, "frozen panel"):
            replace(panel, provenance=provenance, designability_map=rows,
                    pending_targets=pending, panel_manifest_sha256=digest)

    def test_dmd_design_cannot_be_injected_across_deposit_boundary(self):
        panel = run_panel(clock=lambda: FIXED_TIME)
        test_panel = self.run_test_panel(clock=lambda: FIXED_TIME, design_api=self.api(),
            target_contexts={"SMN2_exon7_inclusion": self.context(),
                             "DMD_exon53_inclusion": self.context("DMD_exon53_inclusion")})
        dmd_design = next(x for x in test_panel.designs if x.target_id == "DMD_exon53_inclusion")
        dmd_row = next(x for x in test_panel.designability_map if x.target_id == dmd_design.target_id)
        provenance = panel.to_dict()["provenance"]
        provenance["panel_manifest"]["report_manifests"].append({
            "target_id": dmd_design.target_id,
            "run_manifest_sha256": dmd_design.provenance["run_manifest_sha256"]})
        digest = _hash(provenance["panel_manifest"])
        provenance["panel_manifest_sha256"] = digest
        rows = tuple(dmd_row if x.target_id == dmd_design.target_id else x
                     for x in panel.designability_map)
        pending = tuple(x for x in panel.pending_targets if x.target_id != dmd_design.target_id)
        with self.assertRaisesRegex(ValueError, "committed sequence digest"):
            replace(panel, designs=panel.designs + (dmd_design,), pending_targets=pending,
                    designability_map=rows, provenance=provenance, panel_manifest_sha256=digest)

    def test_engine_metadata_spoof_does_not_bypass_callable_identity(self):
        def fake(*args, **kwargs):
            raise AssertionError("spoofed engine must never be invoked")
        fake.__module__ = "design_engine.engine"
        fake.__qualname__ = "design_for"
        with self.assertRaisesRegex(ValueError, "approved immutable engine manifest"):
            run_panel(clock=lambda: FIXED_TIME, design_api=fake)

    def test_designs_cannot_be_reordered_with_matching_manifest_rehash(self):
        panel = self.run_test_panel(clock=lambda: FIXED_TIME, design_api=self.api(),
            target_contexts={"SMN2_exon7_inclusion": self.context(),
                             "DMD_exon53_inclusion": self.context("DMD_exon53_inclusion")})
        with self.frozen_test_anchor(panel):
            provenance = panel.to_dict()["provenance"]
        provenance["panel_manifest"]["report_manifests"].reverse()
        digest = _hash(provenance["panel_manifest"])
        provenance["panel_manifest_sha256"] = digest
        with self.frozen_test_anchor(panel):
            with self.assertRaisesRegex(ValueError, "frozen target order"):
                replace(panel, designs=tuple(reversed(panel.designs)), provenance=provenance,
                        panel_manifest_sha256=digest)

    def test_fabricated_candidate_scores_fail_deterministic_rederivation(self):
        # Teeth for the render-time re-derivation (P0 close): a hand-built report whose
        # nested design is internally consistent but carries fabricated candidate scores
        # cannot reproduce the approved engine's deterministic output, so it is rejected.
        context = {"SMN2_exon7_inclusion": self.context()}
        panel = self.run_test_panel(clock=lambda: FIXED_TIME, target_contexts=context,
                                    design_api=self.api())
        # Baseline: the genuine panel re-derives byte-for-byte and validates.
        with self.frozen_test_anchor(panel):
            panel.validate_for_rendering()
        design = next(x for x in panel.designs if x.target_id.startswith("SMN2"))
        bumped = list(design.candidates)
        bumped[0] = replace(bumped[0],
                            directional_confidence=min(0.999, bumped[0].directional_confidence + 0.1))
        forged_design = self.rebind_candidates(design, bumped)
        # The forged nested report is internally self-consistent (its own digest recomputed).
        forged_design.validate_for_rendering()
        with self.frozen_test_anchor(panel):
            provenance = panel.to_dict()["provenance"]
        provenance["panel_manifest"]["report_manifests"][0]["run_manifest_sha256"] = (
            forged_design.provenance["run_manifest_sha256"])
        digest = _hash(provenance["panel_manifest"])
        provenance["panel_manifest_sha256"] = digest
        # ... yet it is not what the approved engine deterministically produces.
        with self.frozen_test_anchor(panel):
            with self.assertRaisesRegex(ValueError, "deterministic engine output"):
                replace(panel, designs=(forged_design,), provenance=provenance,
                        panel_manifest_sha256=digest)

    def test_frozen_inputs_are_each_read_once_for_hashing_and_parsing(self):
        original = Path.read_bytes
        counts = {}
        def tracked(path):
            resolved = path.resolve()
            counts[resolved] = counts.get(resolved, 0) + 1
            return original(path)
        with patch.object(Path, "read_bytes", tracked):
            load_panel()
        self.assertEqual(counts[Path("config/event_junctions.GRCh38.tsv").resolve()], 1)
        self.assertEqual(counts[Path("src/design_panel/target_panel.json").resolve()], 1)


if __name__ == "__main__": unittest.main()
