"""Tests for the Arm E StericRuler length-titration ladder builder.

The ladder is generated symmetrically around the true nusinersen 18-mer at its
confirmed ISS-N1 anchor (chr5:70076584-70076601, + strand, GRCh38 v50). These
tests cross-check the builder against:

* the fail-closed nusinersen identity at L=18,
* independently recomputed symmetric-trim coordinates for every even L in 2..30,
* an INDEPENDENT reverse-complement oracle (hand-rolled here, not imported), and
* the published +/- flank string
  chr5:70076575-70076610 = GTAAGTCTGCCAGCATTATGAAAGTGAATCTTACTT
  which spans every rung's window, so the real-fasta fetch is validated
  end-to-end against a known sequence.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from feng.bench.scores import NUSINERSEN_5P_TO_3P, NUSINERSEN_RC, OLIGO_LENGTH
from feng.honesty.states import ALLOWED_RESULT_STATES
from feng.titration import ladder as ladder_mod
from feng.titration.ladder import build_ladder


# --- known-truth constants used only as test oracles -------------------------

GENOME_FASTA = Path(
    "/home/dashii/orcd/scratch/feng/gencode/v50/GRCh38.primary_assembly.genome.fa"
)

# Published + strand flank around ISS-N1 (plan / phase2 recon README).
FLANK_START_1BASED = 70076575
FLANK_SEQUENCE = "GTAAGTCTGCCAGCATTATGAAAGTGAATCTTACTT"

ANCHOR_START_1BASED = 70076584
ANCHOR_END_1BASED = 70076601

_needs_fasta = pytest.mark.skipif(
    not GENOME_FASTA.is_file() or GENOME_FASTA.stat().st_size == 0,
    reason=f"reference fasta not available: {GENOME_FASTA}",
)


def _independent_rc(seq: str) -> str:
    """A hand-rolled reverse complement, deliberately not the module's."""
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement[base] for base in reversed(seq))


def _expected_window(length: int) -> tuple[int, int]:
    """Symmetric trim/extend coords (1-based inclusive), recomputed here."""
    half = (OLIGO_LENGTH - length) // 2
    return ANCHOR_START_1BASED + half, ANCHOR_END_1BASED - half


def _flank_target(start_1based: int, end_1based: int) -> str:
    """Slice the published + strand flank for a 1-based inclusive window."""
    lo = start_1based - FLANK_START_1BASED
    hi = end_1based - FLANK_START_1BASED + 1
    assert 0 <= lo < hi <= len(FLANK_SEQUENCE)
    return FLANK_SEQUENCE[lo:hi]


# --- fasta-free unit tests (coordinate/API contract) -------------------------


def test_flank_reproduces_nusinersen_target():
    # Sanity on the test's own oracle: the L=18 flank slice IS the RC target.
    assert _flank_target(ANCHOR_START_1BASED, ANCHOR_END_1BASED) == NUSINERSEN_RC


def test_independent_rc_oracle_matches_nusinersen():
    assert _independent_rc(NUSINERSEN_RC) == NUSINERSEN_5P_TO_3P


def test_no_aso_literal_hardcoded_in_module():
    source = Path(ladder_mod.__file__).read_text(encoding="utf-8")
    forbidden = [
        NUSINERSEN_5P_TO_3P,      # the 18-mer ASO
        NUSINERSEN_RC,            # its target
        "CACTTTCATAATGCTG",       # L=16 ASO
        "GCCAGCATTATGAAAGTGAA",   # L=20 target
    ]
    hits = [lit for lit in forbidden if lit in source]
    assert not hits, f"ladder.py hardcodes sequence literals: {hits}"


def test_rejects_odd_bounds():
    with pytest.raises(ValueError):
        build_ladder(genome_fasta=GENOME_FASTA, min_len=3, max_len=30)
    with pytest.raises(ValueError):
        build_ladder(genome_fasta=GENOME_FASTA, min_len=2, max_len=29)


def test_rejects_range_excluding_18():
    # Without L=18 the nusinersen anchor invariant cannot be checked.
    with pytest.raises(ValueError):
        build_ladder(genome_fasta=GENOME_FASTA, min_len=2, max_len=16)
    with pytest.raises(ValueError):
        build_ladder(genome_fasta=GENOME_FASTA, min_len=20, max_len=30)


# --- fasta-backed tests (the real end-to-end build) --------------------------


@_needs_fasta
def test_builds_fifteen_even_rungs():
    rows = build_ladder(genome_fasta=GENOME_FASTA)
    assert len(rows) == 15
    assert [r["length"] for r in rows] == list(range(2, 31, 2))


@_needs_fasta
def test_l18_reproduces_nusinersen_exactly():
    rows = build_ladder(genome_fasta=GENOME_FASTA)
    row18 = next(r for r in rows if r["length"] == 18)
    assert row18["aso_5to3"] == NUSINERSEN_5P_TO_3P
    assert row18["target_5to3"] == NUSINERSEN_RC
    assert row18["start"] == ANCHOR_START_1BASED
    assert row18["end"] == ANCHOR_END_1BASED


@_needs_fasta
def test_symmetric_coordinates_every_rung():
    rows = build_ladder(genome_fasta=GENOME_FASTA)
    for row in rows:
        exp_start, exp_end = _expected_window(row["length"])
        assert row["start"] == exp_start
        assert row["end"] == exp_end
        assert row["end"] - row["start"] + 1 == row["length"]
        assert len(row["target_5to3"]) == row["length"]
        assert len(row["aso_5to3"]) == row["length"]


@_needs_fasta
def test_aso_is_reverse_complement_of_target():
    rows = build_ladder(genome_fasta=GENOME_FASTA)
    for row in rows:
        assert row["aso_5to3"] == _independent_rc(row["target_5to3"])


@_needs_fasta
def test_target_matches_published_flank_every_rung():
    rows = build_ladder(genome_fasta=GENOME_FASTA)
    for row in rows:
        assert row["target_5to3"] == _flank_target(row["start"], row["end"])


@_needs_fasta
def test_spot_values_l16_and_l20():
    rows = build_ladder(genome_fasta=GENOME_FASTA)
    row16 = next(r for r in rows if r["length"] == 16)
    row20 = next(r for r in rows if r["length"] == 20)
    # Derived from the flank, not copied from the plan.
    exp16_target = _flank_target(*_expected_window(16))
    exp20_target = _flank_target(*_expected_window(20))
    assert row16["target_5to3"] == exp16_target
    assert row16["aso_5to3"] == _independent_rc(exp16_target)
    assert row20["target_5to3"] == exp20_target
    assert row20["aso_5to3"] == _independent_rc(exp20_target)


@_needs_fasta
def test_schema_and_honesty_state():
    rows = build_ladder(genome_fasta=GENOME_FASTA)
    expected_keys = {
        "length",
        "chrom",
        "start",
        "end",
        "strand",
        "target_5to3",
        "aso_5to3",
        "state",
    }
    for row in rows:
        assert set(row) == expected_keys
        assert row["chrom"] == "chr5"
        assert row["strand"] == "+"
        assert row["state"] in ALLOWED_RESULT_STATES
        assert row["state"] == "candidate"
        assert isinstance(row["length"], int)
        assert isinstance(row["start"], int)
        assert isinstance(row["end"], int)
