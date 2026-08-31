"""A5 tests: transcriptome-wide off-target sieve.

Deliberately self-contained -- no test here depends on the ~3 GB downloaded
GENCODE genome or the built ``transcripts.jsonl`` index being present on
disk. Structural tests build small synthetic ``Transcript``/``IntronFlank``
objects directly (same dataclasses :mod:`smt.sieve.transcriptome` builds from
the real GTF+FASTA, just populated by hand with a known, hand-checkable
structure), so the tests are fast and reproducible without a multi-GB
dependency.
"""

from __future__ import annotations

import pytest

from smt.sieve.polarity import (
    PROXIMITY_NT,
    Category,
    classify_hit,
    score_polarity,
)
from smt.sieve.scan import (
    Candidate,
    default_seed_k,
    hamming,
    score_candidate,
    seed_scan_transcript,
    shard_of,
    tiling_seeds,
)
from smt.sieve.transcriptome import Exon, IntronFlank, Transcript, revcomp

# ---------------------------------------------------------------------------
# a small synthetic "+" strand transcript with a known structure
#
#   exon1  g 1000-1049 (50nt)  t   1- 50   5'UTR g 1000-1009
#   intron1 (short, "whole")   g 1050-1069 (20nt)
#   exon2  g 1100-1149 (50nt)  t  51-100
#   intron2 (long, donor+acceptor flanks of 100nt each)
#   exon3  g 1200-1259 (60nt)  t 101-160   3'UTR g 1250-1259
# ---------------------------------------------------------------------------


def _make_transcript() -> Transcript:
    exons = (
        Exon(1, 1000, 1049, 1, 50),
        Exon(2, 1100, 1149, 51, 100),
        Exon(3, 1200, 1259, 101, 160),
    )
    mrna_seq = "ACGT" * 40  # 160 nt, matches exon t-space total exactly
    assert len(mrna_seq) == 160

    whole_flank = IntronFlank(1, "whole", 1050, 1069, "A" * 20, 20)
    donor_flank = IntronFlank(2, "donor", 1150, 1249, "C" * 100, 300)
    acceptor_flank = IntronFlank(2, "acceptor", 1150, 1249, "G" * 100, 300)

    return Transcript(
        transcript_id="TEST1",
        gene_id="GENE1",
        gene_name="TESTGENE",
        chrom="chrTEST",
        strand="+",
        biotype="protein_coding",
        exons=exons,
        five_utr=((1000, 1009),),
        three_utr=((1250, 1259),),
        mrna_seq=mrna_seq,
        flanks=(whole_flank, donor_flank, acceptor_flank),
    )


def _make_minus_transcript() -> Transcript:
    """Same shape, minus strand, to exercise the strand-flip coordinate math."""
    # On a minus-strand gene, transcript 5'->3' order is DECREASING genomic
    # coordinate, so exon "1" (first in transcript order) is the HIGHEST
    # genomic interval.
    exons = (
        Exon(1, 1200, 1249, 1, 50),
        Exon(2, 1100, 1149, 51, 100),
        Exon(3, 1000, 1059, 101, 160),
    )
    mrna_seq = "ACGT" * 40
    return Transcript(
        transcript_id="TEST1_MINUS",
        gene_id="GENE1M",
        gene_name="TESTGENE_MINUS",
        chrom="chrTEST",
        strand="-",
        biotype="protein_coding",
        exons=exons,
        five_utr=(),
        three_utr=(),
        mrna_seq=mrna_seq,
        flanks=(),
    )


# ---------------------------------------------------------------------------
# exon/intron/UTR classification correctness
# ---------------------------------------------------------------------------


class TestClassification:
    def test_five_utr(self):
        t = _make_transcript()
        # t_pos=5 -> genomic 1000+4=1004, inside the annotated 5'UTR (1000-1009).
        cat, dist = classify_hit(t, "mrna", 5, 5)
        assert cat is Category.FIVE_UTR

    def test_three_utr(self):
        t = _make_transcript()
        # t_pos=155 -> exon3 offset 155-101=54 -> genomic 1200+54=1254, in 3'UTR.
        cat, dist = classify_hit(t, "mrna", 155, 155)
        assert cat is Category.THREE_UTR

    def test_exonic_proximal_near_internal_splice_site(self):
        t = _make_transcript()
        # exon1's downstream edge (t_end=50) IS a real 5'ss (not poly-A, since
        # a later exon exists). t_pos=45 is 5 nt from it.
        cat, dist = classify_hit(t, "mrna", 45, 45)
        assert cat is Category.EXONIC_PROXIMAL
        assert dist == 5

    def test_exonic_core_far_from_any_splice_site(self):
        t = _make_transcript()
        # exon3 is the LAST exon: only its upstream edge (t_start=101) is a
        # real splice site (downstream is poly-A, excluded). t_pos=140 is
        # 39 nt from the one real splice site -> beyond PROXIMITY_NT (30).
        cat, dist = classify_hit(t, "mrna", 140, 140)
        assert cat is Category.EXONIC_CORE
        assert dist == 39
        assert dist > PROXIMITY_NT

    def test_intronic_proximal_near_donor_splice_site(self):
        t = _make_transcript()
        cat, dist = classify_hit(t, "intron:2:donor", 5, 5)
        assert cat is Category.INTRONIC_PROXIMAL
        assert dist == 5

    def test_intronic_core_deep_in_intron(self):
        t = _make_transcript()
        cat, dist = classify_hit(t, "intron:2:donor", 80, 80)
        assert cat is Category.INTRONIC_CORE
        assert dist == 80

    def test_intronic_acceptor_side_distance_from_far_end(self):
        t = _make_transcript()
        # acceptor flank: splice site sits at the flank's LAST base (local
        # pos 100 of a 100-nt flank), so local_pos=95 is 6 nt away
        # (n - (pos+span-1) + 1 = 100 - 95 + 1 = 6).
        cat, dist = classify_hit(t, "intron:2:acceptor", 95, 95)
        assert cat is Category.INTRONIC_PROXIMAL
        assert dist == 6

    def test_minus_strand_mrna_to_genomic_mapping(self):
        t = _make_minus_transcript()
        # t_pos=1 (5' end of mRNA) must map to the HIGH genomic coordinate
        # (exon1's g_end=1249) on a minus-strand transcript.
        assert t.mrna_pos_to_genomic(1) == 1249
        # t_pos=50 (end of exon1 in transcript space) -> exon1's g_start.
        assert t.mrna_pos_to_genomic(50) == 1200
        # t_pos=160 (very 3' end) -> exon3's LOWEST genomic coordinate.
        assert t.mrna_pos_to_genomic(160) == 1000


# ---------------------------------------------------------------------------
# stage 1: seed stage recovers a planted perfect match
# ---------------------------------------------------------------------------


class TestSeedStage:
    def test_recovers_planted_perfect_match_in_mrna(self):
        # Non-periodic on purpose: a repeating oligo (e.g. ACGTACGT...) makes
        # its own reverse complement partially self-similar at shifted
        # registers, which plants spurious near-matches against any filler
        # and defeats the point of this test.
        oligo = "GTCAGCTTAGGCATGA"  # 16-mer, arbitrary, non-periodic
        target = revcomp(oligo)
        # Build a transcript whose mRNA is target planted inside filler.
        filler_before = "TTTTTTTTTT"  # 10 nt
        filler_after = "GGGGGGGGGG"  # 10 nt
        mrna = filler_before + target + filler_after
        exons = (Exon(1, 5000, 5000 + len(mrna) - 1, 1, len(mrna)),)
        t = Transcript(
            "T2", "G2", "GENE2", "chrTEST", "+", "protein_coding",
            exons, (), (), mrna, (),
        )
        cands = seed_scan_transcript(t, oligo)
        assert len(cands) == 1
        c = cands[0]
        assert c.mismatches == 0
        assert c.unit == "mrna"
        # Planted at 0-based offset len(filler_before) -> 1-based local_start.
        assert c.local_start == len(filler_before) + 1
        assert c.local_end == len(filler_before) + len(target)

    def test_no_false_positive_on_unrelated_sequence(self):
        oligo = "ACGTACGTACGTACGT"
        mrna = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"  # poly-T, unrelated
        exons = (Exon(1, 1, len(mrna), 1, len(mrna)),)
        t = Transcript("T3", "G3", "GENE3", "chrTEST", "+", "protein_coding",
                        exons, (), (), mrna, ())
        assert seed_scan_transcript(t, oligo) == []

    def test_tiling_seeds_cover_whole_target_with_default_k(self):
        oligo = "TCACTTTCATAATGCTGG"  # nusinersen, 18-mer
        target = revcomp(oligo)
        k = default_seed_k(len(oligo))
        assert 7 <= k <= 10
        seeds = tiling_seeds(target, k)
        covered = set()
        for off, s in seeds:
            covered.update(range(off, off + len(s)))
        assert covered == set(range(len(target)))

    def test_pigeonhole_guarantee_for_single_mismatch(self):
        """A target with exactly 1 mismatch must still be found by one of the
        two non-overlapping tiles (the mismatch can only spoil one tile)."""
        oligo = "TCACTTTCATAATGCTGG"
        target = list(revcomp(oligo))
        target[3] = "A" if target[3] != "A" else "C"  # corrupt one base
        target = "".join(target)
        mrna = "N" * 5 + target + "N" * 5
        exons = (Exon(1, 1, len(mrna), 1, len(mrna)),)
        t = Transcript("T4", "G4", "GENE4", "chrTEST", "+", "protein_coding",
                        exons, (), (), mrna, ())
        cands = seed_scan_transcript(t, oligo, max_mismatches=3)
        assert len(cands) == 1
        assert cands[0].mismatches == 1


# ---------------------------------------------------------------------------
# ddG_net: reuse of existing thermo primitives, sign conventions
# ---------------------------------------------------------------------------


class TestDdGNet:
    def test_ddG_net_is_sum_of_hybrid_and_unfold(self):
        """score_candidate must not reimplement any energy -- it only calls
        smt.thermo.hybridization and smt.thermo.accessibility and adds."""
        oligo = "TCACTTTCATAATGCTGG"
        target = revcomp(oligo)
        filler = "AAAAAAAAAAAAAAAAAAAA"  # neutral, foldable filler (not 'N')
        mrna = filler + target + filler
        exons = (Exon(1, 1, len(mrna), 1, len(mrna)),)
        t = Transcript("T5", "G5", "GENE5", "chrTEST", "+", "protein_coding",
                        exons, (), (), mrna, ())
        c = Candidate("T5", "GENE5", "chrTEST", "+", "mrna",
                      len(filler) + 1, len(filler) + len(target), 0,
                      len(filler) + 1, len(filler) + len(target))
        row = score_candidate(t, c, oligo)
        assert row is not None
        dg_hybrid = float(row["dG_hybrid"])
        dg_unfold = float(row["dG_unfold"])
        ddg_net = float(row["ddG_net"])
        assert ddg_net == pytest.approx(dg_hybrid + dg_unfold, abs=1e-3)
        # A perfect 18-mer duplex must be strongly favorable (negative).
        assert dg_hybrid < -10.0
        # Opening cost is a -RT*ln(probability) and is >= 0 by construction.
        assert dg_unfold >= 0.0


# ---------------------------------------------------------------------------
# polarity: sign convention and the exonic/intronic flip
# ---------------------------------------------------------------------------


class TestPolaritySignConvention:
    def test_intronic_hit_is_positive_inclusion_promoting(self):
        call = score_polarity(Category.INTRONIC_PROXIMAL, distance_nt=10, ddG_net=-20.0)
        assert call.sign == 1
        assert call.signed_effect > 0

    def test_exonic_hit_is_negative_skipping_promoting(self):
        call = score_polarity(Category.EXONIC_PROXIMAL, distance_nt=10, ddG_net=-20.0)
        assert call.sign == -1
        assert call.signed_effect < 0

    def test_utr_hit_has_no_polarity(self):
        call = score_polarity(Category.FIVE_UTR, distance_nt=0, ddG_net=-20.0)
        assert call.sign == 0
        assert call.signed_effect == 0.0

    def test_unfavorable_binding_gives_zero_effect_regardless_of_category(self):
        call = score_polarity(Category.INTRONIC_PROXIMAL, distance_nt=5, ddG_net=+5.0)
        assert call.signed_effect == 0.0

    def test_same_ddG_net_same_distance_flips_sign_exonic_vs_intronic(self):
        """The core claim of this module: identical binding energetics at an
        identical distance from the splice site predict OPPOSITE splicing
        outcomes depending only on which side of the splice site the site is
        on. This is the exonic-enhancer / intronic-silencer flip the whole
        sieve exists to make explicit.
        """
        exonic = score_polarity(Category.EXONIC_PROXIMAL, distance_nt=15, ddG_net=-18.0)
        intronic = score_polarity(Category.INTRONIC_PROXIMAL, distance_nt=15, ddG_net=-18.0)
        assert exonic.signed_effect < 0
        assert intronic.signed_effect > 0
        assert exonic.signed_effect == pytest.approx(-intronic.signed_effect)

    def test_magnitude_decays_with_distance(self):
        near = score_polarity(Category.INTRONIC_PROXIMAL, distance_nt=2, ddG_net=-20.0)
        far = score_polarity(Category.INTRONIC_CORE, distance_nt=100, ddG_net=-20.0)
        assert near.signed_effect > far.signed_effect > 0


# ---------------------------------------------------------------------------
# sharding: exact partition, copied audit from length_sweep's discipline
# ---------------------------------------------------------------------------


class TestSharding:
    def test_partition_union_and_disjoint(self):
        items = [f"transcript_{i}" for i in range(1237)]  # odd number on purpose
        nshards = 37
        shards = [shard_of(items, s, nshards) for s in range(nshards)]

        union = set()
        total_len = 0
        for s in shards:
            union.update(s)
            total_len += len(s)

        assert total_len == len(items)  # no double-counting -> disjoint
        assert union == set(items)  # union == full set

    def test_out_of_range_shard_rejected(self):
        items = list(range(10))
        with pytest.raises(ValueError):
            shard_of(items, 10, 10)
        with pytest.raises(ValueError):
            shard_of(items, -1, 10)

    def test_single_shard_is_identity(self):
        items = list(range(50))
        assert shard_of(items, 0, 1) == items


# ---------------------------------------------------------------------------
# small helper functions
# ---------------------------------------------------------------------------


def test_hamming_distance():
    assert hamming("ACGT", "ACGT") == 0
    assert hamming("ACGT", "AGGT") == 1
    assert hamming("ACGT", "TGCA") == 4
    with pytest.raises(ValueError):
        hamming("ACG", "ACGT")


def test_revcomp():
    assert revcomp("ACGT") == "ACGT"
    assert revcomp("TCACTTTCATAATGCTGG") == "CCAGCATTATGAAAGTGA"
