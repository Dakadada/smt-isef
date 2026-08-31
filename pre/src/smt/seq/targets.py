"""A7: pinned reference sequences for the TDP-43-dependent cryptic exon
targets STMN2 and UNC13A.

Mirrors ``smt.seq.reference``'s SMN2 pinning pattern (accession + pinned
length assertion + literature-anchored feature coordinates), generalized to
two more genes. Kept as a separate module rather than folded into
``reference.py``/``construct.py`` because those two modules are hard-wired to
NG_008728/SMN2 (``construct._build`` calls ``load_smn2()`` directly); this
module is additive and does not touch either.

THE SIGN FLIP -- why this module exists
----------------------------------------
``smt.thermo.constrained.a_5ss`` was built for SMN2 exon 7: an effector binds
an intronic silencer element (ISS-N1 / the GC core) and the therapeutic goal
is to FREE the exon 7 5' splice site so U1 snRNP can bind and the exon is
INCLUDED. Here the goal is the opposite: an effector should OCCLUDE a cryptic
splice site so it is NOT recognised and the cryptic exon is EXCLUDED. Same
machinery (see ``smt.thermo.occlusion``), opposite sign, different gene. This
module supplies the "different gene" half: real, literature-anchored
coordinates for the two loci, with every coordinate's provenance and
confidence stated explicitly.

No RefSeqGene (NG_) record exists for either STMN2 or UNC13A -- checked
directly (``esearch db=nuccore term=<GENE>[gene] AND Homo sapiens[orgn] AND
NG_[accn]`` returns 0 hits for both, 2026-08-22). The pinned accession is
therefore the RefSeq chromosome record (NC_000008.11 / NC_000019.10, GRCh38.p14)
sliced to the gene body via NCBI eutils ``efetch`` with an explicit
``seq_start``/``seq_stop``/``strand``. This is the same kind of object
NG_008728 is for SMN2 (a curated genomic slice with an annotated feature
table) -- it just lacks the LRG curation layer, so provenance is documented
here instead of inherited from an LRG record.

Strand handling (the reason UNC13A needs no special-casing downstream)
------------------------------------------------------------------------
STMN2 is transcribed on the plus genomic strand; UNC13A on the minus strand
(confirmed via NCBI Gene esummary: UNC13A's ``chrstart`` (17,688,353) >
``chrstop`` (17,601,335)). Requesting UNC13A's slice with ``strand=2`` makes
NCBI's efetch return the REVERSE COMPLEMENT of that genomic range with every
annotation feature already remapped to increasing local coordinates in the
gene's own 5'->3' direction. All coordinates below for UNC13A are therefore
already in "reference space" in exactly the same sense NG_008728's
coordinates are for SMN2: 1-based, increasing 5'->3', no minus-strand
bookkeeping required anywhere else in the codebase. This was verified, not
just assumed: the reference allele NCBI dbSNP reports on the forward genomic
strand for both pinned SNPs (see below), complemented, matches the actual
base found at the derived local position in the fetched slice.

Verification standard applied to every literature coordinate below
----------------------------------------------------------------------
A coordinate is marked CONFIDENCE: HIGH only when, in addition to being
quoted from the primary paper, it was independently cross-checked against
the actual pinned sequence (canonical splice-site dinucleotide present at the
implied boundary, computed feature length matching the paper's stated
length, or both). Anything not cleared to that standard is marked lower and
explained. Nothing here is invented; where a claim could not be pinned to a
specific base pair, it is stated as approximate/derived, not silently
rounded to fit.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from smt.seq.construct import Construct
from smt.seq.reference import DATA, Interval

# ---------------------------------------------------------------------------
# STMN2 -- chr8, NC_000008.11, PLUS strand.
# Gene body chr8:79,611,116-79,666,157 (GRCh38.p14, NCBI Gene ID 11075,
# esummary db=gene id=11075, RefSeq annotation release RS_2025_08, fetched
# live 2026-08-22). This span is stable across every annotation release back
# to RS_2023_03 for this assembly (checked in the same esummary response).
# ---------------------------------------------------------------------------

STMN2_ACCESSION = "NC_000008.11"
STMN2_CHR_START = 79_611_116  # 1-based inclusive, GRCh38, plus strand
STMN2_CHR_END = 79_666_157
STMN2_LENGTH = 55042
STMN2_GB = DATA / "STMN2_NC_000008.11_79611116_79666157.gb"

# Local (gene 5'->3') 1-based inclusive exon coordinates, canonical RefSeq
# transcript NM_007029.4 (5 exons). Read directly out of this slice's own
# GenBank feature table (Bio.SeqIO), not hand-transcribed from elsewhere.
STMN2_EXONS: dict[str, tuple[int, int]] = {
    "1": (2, 99),
    "2": (25687, 25782),
    "3": (30263, 30435),
    "4": (43756, 43947),
    "5": (53700, 55042),
}

# Cryptic exon "2a", the intron-1-embedded poison exon whose inclusion
# truncates STMN2 and triggers NMD. Melamed et al. 2019 (Nat Neurosci 22:180,
# PMC6348009), quoted directly from the article text via a full-text fetch:
# "exon 2a (hg38; chr8:79,616,822-79,617,048)", stated length 227 nt.
# Local = chr - STMN2_CHR_START + 1 (plus strand, no flip):
#   chr8:79,616,822 -> local 5707   chr8:79,617,048 -> local 5933
#
# CONFIDENCE: HIGH. Cross-checked directly against the pinned sequence, not
# just quoted:
#   (a) length 5933-5707+1 = 227, matches the paper's stated length exactly;
#   (b) the 2 nt immediately 5' of local 5707 are "AG" -- the canonical 3'
#       splice-site acceptor dinucleotide, exactly where a spliced-in exon's
#       upstream intron must end;
#   (c) the paper states an AUUAAA polyadenylation signal sits "24
#       nucleotides upstream of the polyadenylation site"; an ATTAAA hexamer
#       is found at local 5910-5915, 18-23 nt upstream of the exon's 3'
#       boundary (5933) depending on which end of the hexamer the paper
#       counts from -- consistent to within its own rounding, not an exact
#       coordinate match but strong corroboration, not merely a name match.
STMN2_CRYPTIC_EXON = Interval(5707, 5933)
STMN2_CRYPTIC_EXON_SHA256 = (
    "d8c89c73082be82dd8b2e16df99f4fc3a760823b3855ef9ddd3a3e97dc4deeea"
)

# 3' splice site (acceptor) readout window: last 6 nt of intron 1 (ending in
# the AG dinucleotide confirmed above, local 5705-5706) + first 3 nt of the
# cryptic exon. Same 9-nt style as smt.seq.construct.SS5_EXON7, mirrored for
# an acceptor rather than a donor because the STMN2 cryptic exon's regulated
# decision is 3'ss recognition (it terminates by intrinsic polyadenylation,
# not by a downstream 5'ss) -- there is no analogous donor to model.
STMN2_CRYPTIC_3SS = Interval(5701, 5709)

# Candidate TDP-43 GU-repeat binding site (Melamed et al. describe "three
# GUGUGU hexamers ... spanning 30 nucleotides downstream of the 3' splice
# site" and show TDP-43 binding there sterically blocks recognition of this
# same acceptor -- i.e. TDP-43's own wild-type mechanism IS an occlusion
# effector on this exact site, which is why an ASO mimicking it is the
# published therapeutic strategy, Baughn et al. 2023, Science 379:1140).
# NOT independently pinned to an exact base range here -- derived by scanning
# the already-pinned STMN2_CRYPTIC_EXON sequence for GT-repeat runs and
# locating three GTGT/GTGTGT clusters at exon-local offsets 53, 63 and 69.
# That is further downstream than the paper's "30 nt" language suggests.
# FLAGGED EXPLICITLY: UNVERIFIED-PRECISE. Treat as an approximate region of
# interest, not a literature-pinned coordinate; do not report its exact
# boundary as sourced.
STMN2_TDP43_GU_SITE = Interval(5759, 5786)  # ASSUMPTION: approximate, see above

# ---------------------------------------------------------------------------
# UNC13A -- chr19, NC_000019.10, MINUS strand.
# Gene body chr19:17,601,335-17,688,353 (GRCh38.p14, NCBI Gene ID 23025,
# esummary db=gene id=23025, fetched live 2026-08-22). Minus strand:
# chrstart (17,688,353) > chrstop (17,601,335) in NCBI's own record.
# ---------------------------------------------------------------------------

UNC13A_ACCESSION = "NC_000019.10"
UNC13A_CHR_START = 17_601_335  # smaller genomic bound, as passed to efetch
UNC13A_CHR_END = 17_688_353  # larger genomic bound == transcription start (minus strand gene)
UNC13A_STRAND = 2  # efetch strand=2: returns revcomp, already 5'->3' local
UNC13A_LENGTH = 87019
UNC13A_GB = DATA / "UNC13A_NC_000019.10_17601335_17688353_revcomp.gb"

# Local (gene 5'->3') 1-based inclusive exon coordinates for the canonical /
# reference transcript NM_001080421.3, read directly out of this slice's own
# GenBank feature table (the ``mRNA`` feature whose ``/transcript_id`` is
# NM_001080421.3; its 44 join() parts are already in ascending local
# coordinates because the slice was fetched with strand=2). The count, 44,
# matches NCBI Gene's own genomicinfo.exoncount=44 for gene ID 23025 exactly
# -- an independent cross-check that this is the transcript the literature's
# "exon 20/21" numbering refers to.
#
# HISTORY -- why this used to be two entries. Only exons "20" and "21" were
# pinned originally, because A7 (cryptic-exon retarget) needs nothing but the
# pair that brackets the cryptic exon. That was correct for A7 and wrong for
# A6: ``smt.minigene.compress.compress_gene`` spans first-to-last exon it is
# given, so whole-gene compression silently produced a 743 nt compression of
# the 1.5 kb window 45394..46961 rather than of the 87019 nt locus. The two
# original values are reproduced here unchanged and still agree with the
# GenBank table to the base -- the pins were right, the table was short.
UNC13A_EXONS: dict[str, tuple[int, int]] = {
    "1": (1, 176),
    "2": (12313, 12342),
    "3": (13598, 13697),
    "4": (15859, 15976),
    "5": (18678, 18801),
    "6": (20164, 20237),
    "7": (21650, 21704),
    "8": (24787, 24822),
    "9": (30085, 30292),
    "10": (31956, 32471),
    "11": (32972, 33080),
    "12": (35677, 35723),
    "13": (38767, 38845),
    "14": (39010, 39015),
    "15": (39371, 39442),
    "16": (39704, 39923),
    "17": (40862, 41089),
    "18": (42243, 42384),
    "19": (42511, 42680),
    "20": (45394, 45509),
    "21": (46798, 46961),
    "22": (47693, 47843),
    "23": (48446, 48514),
    "24": (48829, 48918),
    "25": (49137, 49271),
    "26": (52197, 52330),
    "27": (55161, 55246),
    "28": (55446, 55572),
    "29": (57604, 57700),
    "30": (58066, 58209),
    "31": (59031, 59114),
    "32": (60414, 60491),
    "33": (60757, 60845),
    "34": (61569, 61721),
    "35": (63402, 63525),
    "36": (64807, 64812),
    "37": (66484, 66522),
    "38": (67632, 67661),
    "39": (69392, 69448),
    "40": (69853, 69933),
    "41": (70505, 70652),
    "42": (76499, 76591),
    "43": (78255, 78414),
    "44": (82000, 87018),
}

#: 5'->3' exon order, as key strings. Explicit rather than ``.keys()`` so
#: nothing downstream depends on dict insertion order surviving an edit.
UNC13A_EXON_ORDER: tuple[str, ...] = tuple(str(i) for i in range(1, 45))

#: The pair bracketing the cryptic exon. A7 uses only these two; naming them
#: keeps that lane independent of the full table above.
UNC13A_EXON_20 = Interval(45394, 45509)
UNC13A_EXON_21 = Interval(46798, 46961)

UNC13A_CANONICAL_TRANSCRIPT = "NM_001080421.3"

# Two alternative cryptic-exon isoforms, both between canonical exons 20 and
# 21, sharing one 5' donor and differing only in which of two 3' acceptors is
# used. Ma et al. 2022 (Nat Neurosci, PMID 35410343) and Brown et al. 2022
# (Nature 603:131, PMC8891019); PMC8891019 quoted directly via full-text
# fetch: "128-bp or a 178-bp cryptic exon between the canonical exons 20 and
# 21 (hg38; chr19:17642414-17642541 (128bp); chr19:17642414-17642591 (178bp))";
# donor "near chr19:17642414"; acceptors "near chr19:17642541" (128bp form)
# and "near chr19:17642591" (178bp form).
# Local = UNC13A_CHR_END - chr + 1 (minus strand: local increases as genomic
# coordinate decreases):
#   chr19:17,642,414 -> local 45940 (shared donor-side boundary)
#   chr19:17,642,541 -> local 45813 (128bp acceptor)
#   chr19:17,642,591 -> local 45763 (178bp acceptor)
#
# CONFIDENCE: HIGH, same independent-verification standard as STMN2:
#   (a) both computed lengths match the paper's stated lengths exactly
#       (45940-45813+1 = 128; 45940-45763+1 = 178);
#   (b) both candidate acceptors end in the canonical "AG" dinucleotide in
#       the pinned (already strand-corrected) sequence, and the shared donor
#       boundary (local 45940) is immediately followed by the canonical
#       "GT" -- all three checked directly against the fetched sequence;
#   (c) the pair sits entirely inside intron 20 as independently defined by
#       NM_001080421.3's own exon table (exon 20 ends 45509, exon 21 starts
#       46798), matching the paper's "between exons 20 and 21" framing using
#       the same exon-count transcript NCBI reports for this gene.
UNC13A_CRYPTIC_EXON_128 = Interval(45813, 45940)
UNC13A_CRYPTIC_EXON_178 = Interval(45763, 45940)
UNC13A_CRYPTIC_EXON_128_SHA256 = (
    "e92bdccd84ae9f2d0e9bdc24372fe4602a49df20c590111cde71d13647757fa9"
)
UNC13A_CRYPTIC_EXON_178_SHA256 = (
    "1293b34734cfc9ad4a9439ea80854bfaec235b8e2614f7a05ac0a6bb105744b5"
)

# Splice site readout windows, 9 nt, same style as SS5_EXON7 / STMN2_CRYPTIC_3SS.
UNC13A_CRYPTIC_DONOR = Interval(45938, 45946)  # exon last3 + intron first6; GT confirmed at 45941-45942
UNC13A_CRYPTIC_ACCEPTOR_128 = Interval(45807, 45815)  # intron last6 (AG at 45811-45812) + exon first3
UNC13A_CRYPTIC_ACCEPTOR_178 = Interval(45757, 45765)  # intron last6 (AG at 45761-45762) + exon first3

# ---------------------------------------------------------------------------
# UNC13A risk SNPs rs12608932 / rs12973192.
# Positions from dbSNP (NCBI Variation Services API,
# https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/<id>, GRCh38 SPDI
# placement on NC_000019.10, fetched live 2026-08-22).
# CONFIDENCE on POSITION and REFERENCE ALLELE: HIGH -- independently
# verified, not assumed: the SPDI-reported forward-strand reference base,
# complemented for UNC13A's minus strand, was checked against the actual
# base at the derived local position in the pinned sequence and matches
# exactly for both SNPs.
# CONFIDENCE on which allele is "risk" (vs. reference/major): MEDIUM --
# carried from the literature's own A>C / C>G notation (Ma et al. 2022,
# Brown et al. 2022) rather than independently re-derived from a GWAS
# summary statistics table here; flagged rather than silently asserted.
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SNP:
    rsid: str
    chrom_pos: int  # 1-based, GRCh38 forward (plus) strand, NC_000019.10
    local_pos: int  # 1-based, this module's UNC13A local coordinates
    ref_forward: str  # dbSNP reference allele, forward genomic strand
    alt_forward: str  # literature-cited alternate allele, forward strand
    ref_local: str  # complement of ref_forward (UNC13A is minus strand)
    alt_local: str
    source: str
    confidence: str


# rs12608932: dbSNP SPDI position 17,641,879 (0-based) -> 1-based 17,641,880.
# Local = 17,688,353 - 17,641,880 + 1 = 46,474. Falls in intron 20, between
# the cryptic exon's shared donor (45,940) and canonical exon 21 (46,798) --
# i.e. in the intron but NOT inside either cryptic exon isoform, matching the
# literature's framing ("located in the same intron that harbours CE").
RS12608932 = SNP(
    rsid="rs12608932",
    chrom_pos=17_641_880,
    local_pos=46474,
    ref_forward="A",
    alt_forward="C",
    ref_local="T",
    alt_local="G",
    source=(
        "dbSNP refsnp/12608932 (GRCh38, NC_000019.10); "
        "Brown et al. 2022 Nature 603:131 (PMC8891019); Ma et al. 2022"
    ),
    confidence="position/ref allele HIGH (dbSNP + direct sequence cross-check); risk-allele direction MEDIUM",
)

# rs12973192: dbSNP SPDI position 17,642,429 (0-based) -> 1-based 17,642,430.
# Local = 17,688,353 - 17,642,430 + 1 = 45,924. Falls INSIDE both cryptic
# exon isoforms (45,813-45,940 and 45,763-45,940 both contain 45,924),
# matching the literature's explicit statement that this SNP sits within
# the cryptic exon itself.
RS12973192 = SNP(
    rsid="rs12973192",
    chrom_pos=17_642_430,
    local_pos=45924,
    ref_forward="C",
    alt_forward="G",
    ref_local="G",
    alt_local="C",
    source=(
        "dbSNP refsnp/12973192 (GRCh38, NC_000019.10); "
        "Brown et al. 2022 Nature 603:131 (PMC8891019); Ma et al. 2022"
    ),
    confidence="position/ref allele HIGH (dbSNP + direct sequence cross-check); risk-allele direction MEDIUM",
)


# ---------------------------------------------------------------------------
# sequence loading (mirrors smt.seq.reference.load_smn2)
# ---------------------------------------------------------------------------


def _load(gb_path: Path, accession: str, pinned_length: int) -> str:
    if not gb_path.exists():
        raise FileNotFoundError(
            f"{gb_path} missing. Fetch with NCBI efetch db=nuccore id={accession} "
            "(see this module's docstring for the exact seq_start/seq_stop/strand)."
        )
    from Bio import SeqIO

    rec = SeqIO.read(gb_path, "genbank")
    seq = str(rec.seq).upper()
    if len(seq) != pinned_length:
        raise ValueError(
            f"{accession} slice length {len(seq)} != pinned {pinned_length}; "
            "reference changed, re-verify all coordinates before proceeding."
        )
    return seq


def load_stmn2() -> str:
    """Return the pinned STMN2 gene-body slice, uppercase DNA."""
    return _load(STMN2_GB, STMN2_ACCESSION, STMN2_LENGTH)


def load_unc13a() -> str:
    """Return the pinned UNC13A gene-body slice, uppercase DNA (already
    reverse-complemented to 5'->3' gene orientation; see module docstring)."""
    return _load(UNC13A_GB, UNC13A_ACCESSION, UNC13A_LENGTH)


def slice_seq(seq: str, iv: Interval) -> str:
    """Extract a 1-based inclusive interval (same convention as slice_ref)."""
    return seq[iv.start - 1 : iv.end]


def verify_checksums() -> None:
    """Re-derive the pinned cryptic-exon checksums from disk and assert.

    Cheap extra guard against a coordinate slip: if STMN2_CRYPTIC_EXON,
    UNC13A_CRYPTIC_EXON_128 or _178 above were ever edited by one base
    without updating the checksum, this raises immediately instead of
    silently shifting every downstream footprint score by a nucleotide.
    """
    stmn2 = load_stmn2()
    got = hashlib.sha256(slice_seq(stmn2, STMN2_CRYPTIC_EXON).encode()).hexdigest()
    if got != STMN2_CRYPTIC_EXON_SHA256:
        raise ValueError(f"STMN2_CRYPTIC_EXON checksum mismatch: {got} != {STMN2_CRYPTIC_EXON_SHA256}")

    unc13a = load_unc13a()
    got128 = hashlib.sha256(slice_seq(unc13a, UNC13A_CRYPTIC_EXON_128).encode()).hexdigest()
    if got128 != UNC13A_CRYPTIC_EXON_128_SHA256:
        raise ValueError(f"UNC13A_CRYPTIC_EXON_128 checksum mismatch: {got128} != {UNC13A_CRYPTIC_EXON_128_SHA256}")
    got178 = hashlib.sha256(slice_seq(unc13a, UNC13A_CRYPTIC_EXON_178).encode()).hexdigest()
    if got178 != UNC13A_CRYPTIC_EXON_178_SHA256:
        raise ValueError(f"UNC13A_CRYPTIC_EXON_178 checksum mismatch: {got178} != {UNC13A_CRYPTIC_EXON_178_SHA256}")


# ---------------------------------------------------------------------------
# construct builders -- compact folding windows around each cryptic exon.
#
# Same purpose as smt.seq.construct.compact_window for SMN2's exon 7 (a
# window comfortably under GLOBAL_PF_SOFT_LIMIT for the exact constrained
# partition function), kept here rather than added to construct.py because
# construct._build() is hard-wired to load_smn2()/SMN2_EXONS and generalizing
# it was out of scope for this addition; smt.seq.construct.Construct itself
# (the dataclass and its coordinate-mapping methods) is fully generic and is
# reused unchanged.
# ---------------------------------------------------------------------------


def stmn2_cryptic_window(
    margin_up: int = 400, margin_down: int = 200, seq: Optional[str] = None
) -> Construct:
    """Compact window around the STMN2 cryptic exon 2a, for global folding.

    Default span 827 nt (400 nt of intron-1 context upstream of the acceptor,
    the 227 nt cryptic exon, 200 nt further into intron 1 downstream) -- well
    under GLOBAL_PF_SOFT_LIMIT (2000 nt), same design point as SMN2's 798 nt
    compact_window.
    """
    ref = seq if seq is not None else load_stmn2()
    win = Interval(STMN2_CRYPTIC_EXON.start - margin_up, STMN2_CRYPTIC_EXON.end + margin_down)
    sub = slice_seq(ref, win)
    features = {
        "STMN2_CRYPTIC_EXON": STMN2_CRYPTIC_EXON,
        "STMN2_CRYPTIC_3SS": STMN2_CRYPTIC_3SS,
        "STMN2_TDP43_GU_SITE": STMN2_TDP43_GU_SITE,
    }
    return Construct(f"STMN2_cryptic_u{margin_up}_d{margin_down}", win, (), sub, features)


def unc13a_cryptic_window(
    margin_up: int = 400, margin_down: int = 200, seq: Optional[str] = None
) -> Construct:
    """Compact window spanning both UNC13A cryptic exon isoforms and their
    shared donor, for global folding. Default span 784 nt."""
    ref = seq if seq is not None else load_unc13a()
    win = Interval(
        UNC13A_CRYPTIC_EXON_178.start - margin_up, UNC13A_CRYPTIC_DONOR.end + margin_down
    )
    sub = slice_seq(ref, win)
    features = {
        "UNC13A_CRYPTIC_EXON_128": UNC13A_CRYPTIC_EXON_128,
        "UNC13A_CRYPTIC_EXON_178": UNC13A_CRYPTIC_EXON_178,
        "UNC13A_CRYPTIC_DONOR": UNC13A_CRYPTIC_DONOR,
        "UNC13A_CRYPTIC_ACCEPTOR_128": UNC13A_CRYPTIC_ACCEPTOR_128,
        "UNC13A_CRYPTIC_ACCEPTOR_178": UNC13A_CRYPTIC_ACCEPTOR_178,
    }
    return Construct(f"UNC13A_cryptic_u{margin_up}_d{margin_down}", win, (), sub, features)
