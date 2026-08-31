import unittest

from offtarget import (DeclaredBounds, Evidence, EnumeratorConfig, GenomicExon,
                       PositiveSetMetadata, PositiveSite, Transcript, benchmark,
                       enumerate_offtargets, genomic_to_transcript_interval)
from offtarget.alignment import reverse_complement


class EnumeratorTests(unittest.TestCase):
    OLIGO = "AACCGGTTAACCGGTT"
    PERFECT = reverse_complement(OLIGO)

    def scan(self, site, **kwargs):
        transcript = Transcript("tx", "gene", "TTT" + site + "AAA", 20.0, 0.8)
        config = EnumeratorConfig(max_edits=8, max_indel_bases=2, **kwargs)
        return enumerate_offtargets(self.OLIGO, [transcript], config)

    def test_recovers_bulge_indel_and_more_than_four_mismatches(self):
        mismatch = list(self.PERFECT)
        for index in (0, 2, 4, 6, 8):
            mismatch[index] = "A" if mismatch[index] != "A" else "C"
        cases = [
            (self.PERFECT[:7] + "A" + self.PERFECT[7:], "insertions"),
            (self.PERFECT[:7] + self.PERFECT[8:], "deletions"),
            ("".join(mismatch), "mismatches"),
        ]
        for site, feature in cases:
            with self.subTest(feature=feature):
                hits = self.scan(site)
                local = [h for h in hits if h.start == 3]
                self.assertTrue(local)
                self.assertGreater(getattr(local[0], feature), 0)
        self.assertGreater([h for h in self.scan("".join(mismatch)) if h.start == 3][0].mismatches, 4)

    def test_sign_requires_signed_rbp_evidence(self):
        tx = Transcript("tx", "gene", self.PERFECT, 10, 1)
        unsigned = Evidence("tx", 0, len(self.PERFECT), source="motif", effect_sign="enhancer")
        hit = enumerate_offtargets(self.OLIGO, [tx], evidence=[unsigned])[0]
        self.assertEqual(hit.effect_sign, "unknown")
        signed = Evidence("tx", 0, len(self.PERFECT), source="eCLIP", rbp="RBP1",
                          effect_sign="silencer")
        self.assertEqual(enumerate_offtargets(self.OLIGO, [tx], evidence=[signed])[0].effect_sign,
                         "silencer")

    def test_expression_and_nuclear_localization_raise_weight_and_rank(self):
        transcripts = [
            Transcript("low", "g1", self.PERFECT, 1, 0.1),
            Transcript("high", "g2", self.PERFECT, 50, 0.9),
        ]
        hits = [h for h in enumerate_offtargets(self.OLIGO, transcripts) if h.start == 0]
        by_id = {h.transcript_id: h for h in hits}
        self.assertGreater(by_id["high"].expression_localization_weight,
                           by_id["low"].expression_localization_weight)
        self.assertGreater(by_id["high"].recall_score, by_id["low"].recall_score)
        self.assertLess(by_id["high"].rank, by_id["low"].rank)

    def test_unknown_chemistry_is_flagged(self):
        hit = self.scan(self.PERFECT, chemistry="experimental-X")[0]
        self.assertFalse(hit.chemistry_parameters_known)
        self.assertIsNone(hit.estimated_delta_g)
        self.assertEqual(hit.energy_estimate_kind, "none")

    def test_unknown_chemistry_heuristic_requires_opt_in_and_is_labeled(self):
        hit = self.scan(self.PERFECT, chemistry="experimental-X",
                        allow_unknown_chemistry_heuristic=True)[0]
        self.assertIsNotNone(hit.estimated_delta_g)
        self.assertEqual(hit.energy_estimate_kind, "heuristic-non-calibrated")

    def test_same_start_competing_lengths_both_survive(self):
        hits = self.scan(self.PERFECT)
        lengths = {hit.end - hit.start for hit in hits if hit.start == 3}
        self.assertIn(len(self.PERFECT), lengths)
        self.assertIn(len(self.PERFECT) + 1, lengths)

    def test_expression_and_score_filters_are_annotations_not_recall_gates(self):
        tx = Transcript("tx", "gene", "TTT" + self.PERFECT + "AAA", 0.0, 0.8)
        permissive = enumerate_offtargets(self.OLIGO, [tx], EnumeratorConfig())
        annotated = enumerate_offtargets(
            self.OLIGO, [tx], EnumeratorConfig(min_expression_tpm=1000, min_score=1.0))
        self.assertEqual([(x.start, x.end, x.cigar) for x in permissive],
                         [(x.start, x.end, x.cigar) for x in annotated])
        positive = [PositiveSite("tx", 3, 3 + len(self.PERFECT))]
        self.assertEqual(benchmark(permissive, positive).recall,
                         benchmark(annotated, positive).recall)
        self.assertTrue(any(not x.passes_expression_annotation for x in annotated))
        self.assertTrue(any(not x.passes_score_annotation for x in annotated))

    def test_out_of_bound_positive_is_an_explicit_miss(self):
        bounds = DeclaredBounds(2, 1, 15, 17)
        positive = PositiveSite("tx", 3, 21, required_edits=4,
                                required_indel_bases=2)
        result = benchmark([], [positive], declared_bounds=bounds)
        self.assertEqual(result.missed, (positive,))
        self.assertEqual(result.out_of_bound_misses, (positive,))

    def test_empty_positive_set_has_undefined_recall(self):
        result = benchmark([], [])
        self.assertIsNone(result.recall)
        self.assertIsNone(result.target_met)
        self.assertEqual(result.recall_status, "undefined-empty-positive-set")

    def test_incomplete_positive_set_is_provisional(self):
        positive = PositiveSite("tx", 0, len(self.PERFECT))
        result = benchmark([], [positive], positive_set_metadata=PositiveSetMetadata(
            provenance="frozen-set-v1", complete=False, blinded=True))
        self.assertEqual(result.recall, 0.0)
        self.assertIsNone(result.target_met)
        self.assertEqual(result.recall_status, "provisional-incomplete-metadata")

    def test_endpoint_matching_rejects_loose_overlap(self):
        hits = self.scan(self.PERFECT)
        # This interval overlaps many candidates but matches neither endpoint pair.
        result = benchmark(hits, [PositiveSite("tx", 4, 30)])
        self.assertEqual(result.recovered, 0)

    def test_conflicting_signed_evidence_is_unknown(self):
        tx = Transcript("tx", "gene", self.PERFECT, 10, 1)
        evidence = [Evidence("tx", 0, 16, effect_sign="enhancer"),
                    Evidence("tx", 0, 16, effect_sign="silencer")]
        self.assertEqual(enumerate_offtargets(self.OLIGO, [tx], evidence=evidence)[0].effect_sign,
                         "unknown")

    def test_perturbation_only_signed_evidence_is_accepted(self):
        tx = Transcript("tx", "gene", self.PERFECT, 10, 1)
        evidence = Evidence("tx", 0, 16, source="knockdown",
                            effect_sign="enhancer", perturbation_supported=True)
        self.assertEqual(enumerate_offtargets(self.OLIGO, [tx], evidence=[evidence])[0].effect_sign,
                         "enhancer")

    def test_default_sign_is_unknown(self):
        self.assertEqual(self.scan(self.PERFECT)[0].effect_sign, "unknown")

    def test_minus_strand_coordinate_conversion(self):
        exons = [GenomicExon(100, 110), GenomicExon(200, 220)]
        self.assertEqual(genomic_to_transcript_interval(214, 219, "-", exons), (1, 6))

    def test_mock_blinded_positive_recall_and_separate_precision(self):
        hits = self.scan(self.PERFECT)
        positives = [PositiveSite("tx", 3, 3 + len(self.PERFECT)),
                     PositiveSite("missing", 1, 10)]
        result = benchmark(hits, positives)
        self.assertEqual(result.recall, 0.5)
        self.assertIsNone(result.target_met)
        self.assertIsNone(result.precision)
        labeled = benchmark(hits, positives, precision_labels_complete=True)
        self.assertIsNotNone(labeled.precision)

    def test_complete_blinded_metadata_can_issue_verdict(self):
        hits = self.scan(self.PERFECT)
        result = benchmark(hits, [PositiveSite("tx", 3, 19)],
                           positive_set_metadata=PositiveSetMetadata(
                               "sha256:example", complete=True, blinded=True))
        self.assertTrue(result.target_met)
        self.assertEqual(result.recall_status, "validated")


if __name__ == "__main__":
    unittest.main()
