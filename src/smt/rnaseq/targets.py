"""Pinned GRCh38 coordinates for the events counted by ``smt.rnaseq``.

Every coordinate here was read out of the annotation the BAMs were actually
aligned against -- ``gencode.v44.basic.annotation.gtf``, the same file passed
to ``STAR --sjdbGTFfile`` in ``slurm/star_index.sbatch`` -- and every
paralog-discriminating base was then confirmed against
``GRCh38.primary_assembly.genome.fa`` itself. Nothing here is from memory.

SMN EXON NUMBERING. The SMN literature numbers exons 1, 2a, 2b, 3, 4, 5, 6,
7, 8; GENCODE numbers the same nine exons 1..9. So literature exon 7 -- the
54 nt exon this project is about -- is GENCODE ``exon_number 8`` of the
9-exon canonical transcripts ENST00000380743.9 (SMN2) and ENST00000380707.9
(SMN1). The 54 nt length is the check that the mapping is right.
"""

from __future__ import annotations

from smt.rnaseq.junction_psi import Cassette
from smt.seq.reference import Interval

#: SMN2 exon 7 cassette. GENCODE exons 7/8/9 of ENST00000380743.9 == the
#: literature's exons 6/7/8.
SMN2_EXON7 = Cassette(
    gene="SMN2", chrom="chr5",
    upstream=Interval(70070641, 70070751),    # exon 6, 111 nt
    cassette=Interval(70076521, 70076574),    # exon 7,  54 nt
    downstream=Interval(70077019, 70077595),  # exon 8, 577 nt
)

#: The SMN1 copy of the same event, ENST00000380707.9. Counted only to be
#: POOLED with SMN2 -- see the ``smt.rnaseq.junction_psi`` module docstring
#: for why the exclusion junction cannot be assigned to a paralog at all.
SMN1_EXON7 = Cassette(
    gene="SMN1", chrom="chr5",
    upstream=Interval(70946066, 70946176),    # exon 6, 111 nt
    cassette=Interval(70951941, 70951994),    # exon 7,  54 nt
    downstream=Interval(70952439, 70953015),  # exon 8, 577 nt
)

SMN_EXON7_LOCI = (SMN2_EXON7, SMN1_EXON7)

#: c.840C>T, exon 7 position 6. SMN1 carries C, SMN2 carries T. Verified by
#: faidx on both loci, and by diffing the two 54 nt exon 7 sequences: they
#: differ at exactly ONE position, this one.
C840_SMN2 = 70076526
C840_SMN1 = 70951946

#: POLR2H exon 3 (GENCODE-canonical) == exon 2 in five transcripts; see
#: docs/HANDOFF.md. This is the positive control: rMATS calls this event at
#: dPSI -0.920, FDR 0, so any counter that cannot reproduce it is wrong.
#: Taken from rMATS's own SE.MATS.JC.txt row, exonStart_0base converted +1.
POLR2H_EXON3 = Cassette(
    gene="POLR2H", chrom="chr3",
    upstream=Interval(184362873, 184363565),
    cassette=Interval(184364966, 184365049),
    downstream=Interval(184365133, 184365226),
)
POLR2H_RMATS_DPSI = -0.920

#: GSE307431 group membership, from slurm/rmats b1.txt / b2.txt -- the same
#: split that produced the -0.920. b1 is F18MOE (nusinersen's chemistry),
#: b2 is ScrMOE (its scrambled control), so dPSI is reported b1 - b2 to
#: match rMATS's IncLevelDifference sign convention.
F18MOE = ("SRR35298700", "SRR35298705", "SRR35298714")
SCRMOE = ("SRR35298699", "SRR35298704", "SRR35298715")

BAM_DIR = "/home/dashii/orcd/scratch/smt-rna/bam"


def bam_path(accession: str) -> str:
    return f"{BAM_DIR}/{accession}.Aligned.sortedByCoord.out.bam"
