"""Tests for R3: paralog-aware SMN2 exon 7 quantification.

The load-bearing facts here are NOT things this code decides -- they are
properties of GRCh38 and of the annotation STAR was indexed against. So each
one is re-derived from that source file inside the test rather than restated
from ``smt.rnaseq.targets``. A test that only repeats the module's own
constants would have passed just as happily on the truncated UNC13A exon
table (see docs/HANDOFF.md), which is exactly the failure mode to avoid.

Three facts carry the whole design and each gets its own test:
  * SMN exon 6 is IDENTICAL between the paralogs -> a 6->8 skipping read
    carries no paralog information -> the loci must be pooled;
  * SMN exon 7 differs at exactly ONE base, c.840 -> that base, and only
    that base, can assign an inclusion read;
  * POLR2H reproduces rMATS's dPSI -> the counter is measuring what rMATS
    measured, so the SMN2 number it produces means something.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from smt.rnaseq import targets as T
from smt.rnaseq.junction_psi import (
    SAMTOOLS, Cassette, _parse, base_census, count_junctions, form_lengths,
    pooled_counts, psi, read_length, JunctionCounts,
)
from smt.rnaseq.report import CONTROL_TOL, wilson_upper
from smt.seq.reference import Interval

GENOME = Path("/home/dashii/orcd/pool/smt-isef-data/transcriptome/GRCh38.primary_assembly.genome.fa")
GTF = Path("/home/dashii/orcd/pool/smt-isef-data/transcriptome/gencode.v44.basic.annotation.gtf")

needs_genome = pytest.mark.skipif(
    not (GENOME.exists() and SAMTOOLS.exists()),
    reason="GRCh38 fasta / rna-env samtools not present (cluster only)",
)
needs_bams = pytest.mark.skipif(
    not Path(T.bam_path(T.F18MOE[0])).exists() or not SAMTOOLS.exists(),
    reason="STAR BAMs not present (cluster only)",
)
needs_gtf = pytest.mark.skipif(not GTF.exists(), reason="gencode GTF not present (cluster only)")


def faidx(region: str) -> str:
    out = subprocess.run([str(SAMTOOLS), "faidx", str(GENOME), region],
                         capture_output=True, text=True, check=True).stdout
    return "".join(out.splitlines()[1:]).upper()


# ---------------------------------------------------------------- geometry


def test_cigar_parse_blocks_and_gaps():
    """A spliced CIGAR splits into the right reference blocks and N gaps."""
    blocks, gaps = _parse(100, "30M170N70M")
    assert blocks == [(100, 129, 0), (300, 369, 30)]
    assert gaps == [(300 - 170, 299)]


def test_cigar_parse_soft_clip_shifts_query_not_reference():
    blocks, gaps = _parse(100, "5S20M")
    assert blocks == [(100, 119, 5)]
    assert gaps == []


def test_cassette_gaps_are_the_introns():
    cas = Cassette("X", "chrZ", Interval(100, 200), Interval(300, 350), Interval(500, 600))
    assert cas.inc_up_gap == (201, 299)
    assert cas.inc_down_gap == (351, 499)
    assert cas.skip_gap == (201, 499)


def test_form_lengths_match_rmats_own_output():
    """rMATS reported 183/99 for the 84 nt POLR2H cassette at 100 nt reads.

    And 198/99 for POLR2H's 693 nt cassette in the same file, which is what
    pins the ``min(2*(L-1), exon+L-1)`` form rather than just ``exon+L-1``.
    """
    assert form_lengths(T.POLR2H_EXON3, 100) == (183, 99)
    big = Cassette("X", "chr3", Interval(1, 10), Interval(100, 792), Interval(1000, 1010))
    assert len(big.cassette) == 693
    assert form_lengths(big, 100) == (198, 99)


def test_psi_is_none_when_there_are_no_reads():
    """No reads is not the same as no inclusion, so it must not return 0 or 1."""
    assert psi(JunctionCounts(0, 0, 0, __import__("collections").Counter()),
               T.POLR2H_EXON3, 100) is None


def test_wilson_upper_bounds_a_zero_count():
    """Zero skipping out of n is a bound, not a proof of zero."""
    assert wilson_upper(0, 0) == 1.0
    assert 0.0 < wilson_upper(0, 118) < 0.05
    assert wilson_upper(0, 500) < wilson_upper(0, 118)


# ------------------------------------------------- annotation ground truth


@pytest.fixture(scope="session")
def gencode_smn_exons():
    """Re-read the SMN exons from the GTF STAR was actually indexed against."""
    wanted = {"ENST00000380743": "SMN2", "ENST00000380707": "SMN1"}
    found: dict[str, list[tuple[int, int]]] = {"SMN1": [], "SMN2": []}
    with GTF.open() as fh:
        for line in fh:
            if "\texon\t" not in line:
                continue
            for tid, gene in wanted.items():
                if tid in line:
                    f = line.split("\t")
                    found[gene].append((int(f[3]), int(f[4])))
    return {g: sorted(v) for g, v in found.items()}


@needs_gtf
def test_smn_cassettes_match_the_indexed_annotation(gencode_smn_exons):
    """Both pinned cassettes are exons 7/8/9 of the 9-exon canonical form.

    GENCODE exon 8 == the literature's exon 7. The nine-exon count and the
    54 nt cassette are what make that mapping checkable rather than assumed.
    """
    for gene, cas in (("SMN2", T.SMN2_EXON7), ("SMN1", T.SMN1_EXON7)):
        exons = gencode_smn_exons[gene]
        assert len(exons) == 9, f"{gene}: expected 9 exons, got {len(exons)}"
        assert exons[6] == (cas.upstream.start, cas.upstream.end)
        assert exons[7] == (cas.cassette.start, cas.cassette.end)
        assert exons[8] == (cas.downstream.start, cas.downstream.end)
        assert len(cas.cassette) == 54


# -------------------------------------------------- the paralog facts


@needs_genome
def test_exon6_carries_no_paralog_information():
    """SMN1 and SMN2 exon 6 are identical.

    This is WHY the two loci must be pooled: a 6->8 skipping read overlaps
    only exon 6 and the first bases of exon 8, so nothing in it says which
    paralog it came from. If this test ever fails, per-read assignment of
    the exclusion isoform becomes possible and the pooling can be dropped.
    """
    a = faidx(f"chr5:{T.SMN1_EXON7.upstream.start}-{T.SMN1_EXON7.upstream.end}")
    b = faidx(f"chr5:{T.SMN2_EXON7.upstream.start}-{T.SMN2_EXON7.upstream.end}")
    assert len(a) == len(b) == 111
    assert a == b


@needs_genome
def test_exon7_differs_at_exactly_one_base_and_it_is_c840():
    """One difference, at position 6, SMN1 C / SMN2 T. That is c.840C>T."""
    a = faidx(f"chr5:{T.SMN1_EXON7.cassette.start}-{T.SMN1_EXON7.cassette.end}")
    b = faidx(f"chr5:{T.SMN2_EXON7.cassette.start}-{T.SMN2_EXON7.cassette.end}")
    assert len(a) == len(b) == 54
    diffs = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
    assert diffs == [5], f"expected one difference at exon position 6, got {diffs}"
    assert (a[5], b[5]) == ("C", "T")
    assert T.C840_SMN1 == T.SMN1_EXON7.cassette.start + 5
    assert T.C840_SMN2 == T.SMN2_EXON7.cassette.start + 5


@needs_genome
def test_c840_positions_read_back_the_expected_bases():
    assert faidx(f"chr5:{T.C840_SMN1}-{T.C840_SMN1}") == "C"
    assert faidx(f"chr5:{T.C840_SMN2}-{T.C840_SMN2}") == "T"


# --------------------------------------------------------- against the data


@needs_bams
def test_smn1_is_absent_from_this_dataset():
    """0 C at c.840 is the premise that makes pooling legitimate.

    Measured on the three ScrMOE samples (the untreated-splicing arm, where
    any SMN1 would be easiest to see), both loci.
    """
    total = 0
    for acc in T.SCRMOE:
        bam = Path(T.bam_path(acc))
        for cas, pos in ((T.SMN2_EXON7, T.C840_SMN2), (T.SMN1_EXON7, T.C840_SMN1)):
            census = base_census(bam, cas.chrom, pos, cas.cassette)
            assert census["C"] == 0, f"{acc}: {census['C']} SMN1 reads -- pooling invalid"
            total += census["T"]
    assert total > 100, f"only {total} SMN2 reads at c.840 -- too thin to conclude"


@needs_bams
def test_skipping_reads_land_at_smn1_and_inclusion_at_smn2():
    """The tie-break this design exists to defeat, asserted as a fact.

    In a SMN1-null line every read is SMN2, yet the aligner puts the
    inclusion junctions at the SMN2 locus and the paralog-uninformative
    6->8 skipping junction at SMN1. Counting either locus alone would give
    a fabricated answer: SMN2 alone reads as PSI 1.0 in EVERY sample.
    """
    bam = Path(T.bam_path(T.SCRMOE[0]))
    at_smn2 = count_junctions(bam, T.SMN2_EXON7)
    at_smn1 = count_junctions(bam, T.SMN1_EXON7)
    assert at_smn2.ijc > 50 and at_smn2.sjc == 0
    assert at_smn1.sjc > 10
    pooled = pooled_counts(bam, T.SMN_EXON7_LOCI)
    assert pooled.sjc == at_smn2.sjc + at_smn1.sjc > 0


@needs_bams
def test_polr2h_control_reproduces_rmats():
    """The gate on everything else: same event, same normalization, same answer.

    Counts run above rMATS's because multimappers are included and no
    anchor filter is applied, so the tolerance is on the PSI, not the counts.
    """
    means = []
    for accs in (T.F18MOE, T.SCRMOE):
        vals = []
        for acc in accs:
            bam = Path(T.bam_path(acc))
            counts = count_junctions(bam, T.POLR2H_EXON3)
            vals.append(psi(counts, T.POLR2H_EXON3, read_length(bam, T.POLR2H_EXON3.region)))
        means.append(sum(vals) / len(vals))
    dpsi = means[0] - means[1]
    assert abs(dpsi - T.POLR2H_RMATS_DPSI) < CONTROL_TOL, (
        f"POLR2H dPSI {dpsi:+.3f} vs rMATS {T.POLR2H_RMATS_DPSI:+.3f}"
    )


@needs_bams
def test_exact_gap_matching_does_not_miss_a_shifted_junction():
    """An exact-coordinate matcher returns a SILENT zero if the aligner put
    the junction one base off -- and a silent zero here reads as "no
    skipping", i.e. as a perfect result. So assert that the three counted
    junctions really are the dominant ones in the region."""
    bam = Path(T.bam_path(T.SCRMOE[0]))
    counts = count_junctions(bam, T.SMN2_EXON7)
    top = [gap for gap, _ in counts.gaps.most_common(4)]
    assert T.SMN2_EXON7.inc_up_gap in top
    assert T.SMN2_EXON7.inc_down_gap in top


@needs_bams
def test_smn2_exon7_inclusion_rises_under_the_iss_n1_aso():
    """The actual R3 deliverable, asserted as a direction and a floor.

    F18MOE is nusinersen's chemistry targeting ISS-N1, so exon 7 inclusion
    must RISE. Bounds are deliberately loose -- this asserts the finding,
    not the third decimal place.
    """
    means = []
    for accs in (T.F18MOE, T.SCRMOE):
        vals = []
        for acc in accs:
            bam = Path(T.bam_path(acc))
            counts = pooled_counts(bam, T.SMN_EXON7_LOCI)
            vals.append(psi(counts, T.SMN2_EXON7, read_length(bam, T.SMN2_EXON7.region)))
        means.append(sum(vals) / len(vals))
    f18, scr = means
    assert 0.5 < scr < 0.8, f"ScrMOE baseline PSI {scr:.3f} outside the observed band"
    assert f18 > 0.95, f"F18MOE PSI {f18:.3f} -- expected near-complete inclusion"
    assert f18 - scr > 0.25, f"dPSI {f18 - scr:+.3f} -- expected a large rise"
