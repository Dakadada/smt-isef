import unittest

from physics_adapter import (AppliedConstraintType, Baseline,
    CapacityConstrainedLink, CapacityError, ConstraintReceipt,
    CoordinateConvention, DuplexOrientation, EffectorPhysicsAdapter,
    EffectorSpec, EventDefinition, FROZEN_THETA_GRID, FoldingResult,
    IntermolecularComplexSpec, IntermolecularPair, Interval,
    MagnitudeEvidence, ModelClass, ModelConfig, OccupancySource, Probability,
    TableFolder, TranscriptContext, prediction_input_hash, protocol_hash,
    robust_interval_direction)


def complex_spec(target=Interval(20, 25), aso_id="aso"):
    return IntermolecularComplexSpec(
        target, CoordinateConvention.ZERO_BASED_HALF_OPEN, "UUUUU", "tx", aso_id,
        DuplexOrientation.ANTIPARALLEL, target.start, "2'-MOE phosphorothioate",
        tuple(IntermolecularPair(target.start + i, 4 - i) for i in range(5)))


def prediction_inputs():
    context = TranscriptContext("tx", "A" * 40)
    effector = EffectorSpec("aso", Interval(20, 25), complex_spec())
    event = EventDefinition("event", (Interval(2, 5),), (Interval(10, 12),),
                            (Interval(15, 17),))
    return context, effector, event


def fixture(values=None, folder_class=TableFolder, evidence=None):
    values = values or {
      (ModelClass.EQUILIBRIUM_FULL_TRANSCRIPT.value, False): Probability(.2, .18, .22),
      (ModelClass.EQUILIBRIUM_FULL_TRANSCRIPT.value, True): Probability(.8, .75, .85),
      (ModelClass.COTRANSCRIPTIONAL.value, False): Probability(.3, .25, .35),
      (ModelClass.COTRANSCRIPTIONAL.value, True): Probability(.6, .55, .65),
      (ModelClass.ECLIP_PROBING_ABLATION.value, False): Probability(.1, .08, .12),
      (ModelClass.ECLIP_PROBING_ABLATION.value, True): Probability(.7, .65, .75),
    }
    folder = folder_class(values)
    context, effector, event = prediction_inputs()
    result = EffectorPhysicsAdapter(folder).predict(
        context, effector, event, magnitude_evidence=evidence)
    return folder, result


def valid_rows(studies=4, repetitions=6):
    rows = []
    for study_index in range(studies):
        for chemistry in ("moe", "pmo"):
            for repetition in range(repetitions):
                label = (repetition + study_index) % 2
                rows.append({"p": .2 + .6 * label, "psi": label,
                             "study_family": f"study-{study_index}",
                             "chemistry": chemistry, "fold": study_index % 2})
    return rows


class ForcedUnpairedFolder(TableFolder):
    """Contract-negative backend: reports a region was forced unpaired."""

    def productive_probability(self, request):
        result = super().productive_probability(request)
        if request.bound_complex is None:
            return result
        receipt = result.receipt
        return FoldingResult(result.probability, ConstraintReceipt(
            AppliedConstraintType.FORCED_UNPAIRED, receipt.coordinate_convention,
            receipt.transcript_sequence_hash, receipt.aso_sequence_hash,
            receipt.strand_hash, receipt.pairing_hash, receipt.applied_orientation,
            receipt.applied_register, receipt.chemistry_hash,
            receipt.backend_version, False))


class PhysicsAdapterTest(unittest.TestCase):
    def test_bound_request_is_a_complete_intermolecular_complex(self):
        folder, result = fixture()
        bound = [r for r in folder.requests if r.bound_complex is not None]
        self.assertEqual(len(bound), 3)
        spec = bound[0].bound_complex
        self.assertEqual(spec.transcript_coordinates, Interval(20, 25))
        self.assertEqual(spec.aso_sequence, "UUUUU")
        self.assertEqual(spec.orientation, DuplexOrientation.ANTIPARALLEL)
        self.assertEqual(len(spec.pairing_constraints), 5)
        self.assertTrue(all(Interval(20, 25) in r.unavailable_to_factors for r in bound))
        state = result.classes[ModelClass.EQUILIBRIUM_FULL_TRANSCRIPT].bound_state
        self.assertEqual(state.target_paired_to_aso, (Interval(20, 25),))

    def test_forced_unpaired_backend_is_rejected_fail_closed(self):
        with self.assertRaisesRegex(ValueError, "receipt does not match"):
            fixture(folder_class=ForcedUnpairedFolder)

    def test_direction_and_theta_grid_are_correct(self):
        _, result = fixture()
        primary = result.classes[ModelClass.EQUILIBRIUM_FULL_TRANSCRIPT]
        self.assertEqual(result.theta_grid, FROZEN_THETA_GRID)
        self.assertAlmostEqual(primary.theta_probabilities[.5].estimate, .5)
        self.assertEqual(result.primary_direction, 1)
        self.assertFalse(result.magnitude_claim_allowed)

    def test_theta_sign_flip_abstains(self):
        self.assertIsNone(robust_interval_direction({0: (0, 0), .25: (.1, .2),
                                                     .5: (-.2, -.1)}))

    def test_uncertainty_driven_abstention_uses_intervals(self):
        values = {(m.value, False): Probability(.4, .2, .6) for m in ModelClass}
        values.update({(m.value, True): Probability(.7, .5, .9) for m in ModelClass})
        _, result = fixture(values)
        primary = result.classes[ModelClass.EQUILIBRIUM_FULL_TRANSCRIPT]
        self.assertGreater(primary.theta_deltas[1], 0)
        self.assertIsNone(result.primary_direction)

    def test_concordance_is_diagnostic_and_ablation_never_votes(self):
        values = {(m.value, b): Probability(v, v, v) for m, b, v in [
          (ModelClass.EQUILIBRIUM_FULL_TRANSCRIPT, False, .8), (ModelClass.EQUILIBRIUM_FULL_TRANSCRIPT, True, .2),
          (ModelClass.COTRANSCRIPTIONAL, False, .2), (ModelClass.COTRANSCRIPTIONAL, True, .8),
          (ModelClass.ECLIP_PROBING_ABLATION, False, .9), (ModelClass.ECLIP_PROBING_ABLATION, True, .1)]}
        _, result = fixture(values)
        self.assertEqual(result.primary_direction, -1)
        self.assertFalse(result.directional_concordance)

    def test_label_pure_two_row_studies_are_rejected(self):
        rows = []
        for i in range(12):
            label = i % 2
            for chemistry in ("moe", "pmo"):
                rows.append({"p": .5, "psi": label, "study_family": f"s{i}",
                             "chemistry": chemistry, "fold": i % 2})
        with self.assertRaisesRegex(CapacityError, "label-pure"):
            CapacityConstrainedLink().validate_training_design(rows)

    def test_sparse_joint_study_chemistry_identities_are_rejected(self):
        rows = valid_rows()
        for row in rows:
            row["chemistry"] = "c-" + row["study_family"]
        with self.assertRaisesRegex(CapacityError, "crossed"):
            CapacityConstrainedLink().validate_training_design(rows)

    def test_parameter_ratio_and_group_disjoint_cross_fitting_are_enforced(self):
        with self.assertRaisesRegex(CapacityError, "degree of freedom"):
            CapacityConstrainedLink().validate_training_design(valid_rows(repetitions=2))
        rows = valid_rows()
        rows[0]["fold"] = 1
        with self.assertRaisesRegex(CapacityError, "disjoint"):
            CapacityConstrainedLink().validate_training_design(rows)

    def test_evidence_bound_magnitude_authorization(self):
        context, effector, event = prediction_inputs()
        config = ModelConfig()
        evidence = MagnitudeEvidence("assay:plate-7", MagnitudeEvidence.MATCHING_CRITERION,
            OccupancySource.MATCHED_EXPERIMENTAL,
            prediction_input_hash(context, effector, event), protocol_hash(config))
        _, result = fixture(evidence=evidence)
        self.assertTrue(result.magnitude_claim_allowed)
        wrong = MagnitudeEvidence("assay:plate-7", MagnitudeEvidence.MATCHING_CRITERION,
            OccupancySource.MATCHED_EXPERIMENTAL, "0" * 64, protocol_hash(config))
        with self.assertRaisesRegex(ValueError, "does not match prediction inputs"):
            fixture(evidence=wrong)
        with self.assertRaises(TypeError):
            EffectorPhysicsAdapter(TableFolder({})).predict(
                context, effector, event, matched_theta_evidence=True)

    def test_fitted_link_is_monotone_and_capacity_constrained(self):
        fitted = CapacityConstrainedLink().fit(valid_rows(), iterations=100)
        predictions = [fitted.predict(x) for x in (0, .25, .5, .75, 1)]
        self.assertEqual(predictions, sorted(predictions))
        budgets = [CapacityConstrainedLink().baseline_capacity(b) for b in Baseline]
        self.assertEqual(len({(x["global_df"], x["max_group_levels"],
                               x["min_observations_per_degree_of_freedom"],
                               x["ridge_precision"]) for x in budgets}), 1)

    def test_invalid_complex_pair_register_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "antiparallel"):
            IntermolecularComplexSpec(Interval(20, 25),
                CoordinateConvention.ZERO_BASED_HALF_OPEN, "UUUUU", "tx", "aso",
                DuplexOrientation.ANTIPARALLEL, 20, "PMO",
                tuple(IntermolecularPair(20 + i, i) for i in range(5)))


if __name__ == "__main__":
    unittest.main()
