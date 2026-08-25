"""Tests for S5's pinned inputs: the U1 duplex arm and the risdiplam ligand.

Nothing in S5 may be typed from memory, so these tests re-derive every claim
from the pinned records rather than restating the module's constants.

The load-bearing check is ``test_u1_arm_pairs_perfectly_with_the_consensus``.
``U1_ARM = (3, 11)`` is the one window choice in this module that could be
wrong in a way that still looks plausible -- an off-by-one would still produce
a 9-mer and still "pair" partially. Requiring 9/9 against the independently
known 5'ss consensus pins the register from outside the module.

Covers:
* the U1 record's length is asserted, so a silently updated RefSeq entry fails
  loudly instead of shifting the arm;
* the arm pairs 9/9 with the consensus CAGGUAAGU -- which is what fixes the
  register -- and 7/9 with SMN2;
* SMN2's two mismatches land at exactly the positions that are SMN2-specific,
  and the second of them is re-derived to be reference coordinate 32052, the
  bulged -1 adenosine risdiplam stacks on. This is the structural premise of
  the whole arm and it is DERIVED, not asserted;
* ``pairs_with`` really anneals ANTIPARALLEL -- a parallel implementation is
  caught by a fixture built to pair one way and not the other;
* the ligand validates against PubChem, and a deliberately corrupted ligand is
  REJECTED. A validator that cannot fail is not a validator.
"""

from __future__ import annotations

import pytest

from smt.seq.construct import SS5_EXON7
from smt.seq.reference import load_smn2, slice_ref
from smt.structure import s5_inputs
from smt.structure.s5_inputs import (
    U1_ARM,
    U1_LENGTH,
    LigandError,
    load_u1,
    pairs_with,
    u1_duplex_arm,
    validate_ligand,
)

#: The mammalian 5' splice site consensus. Independent of this project's code
#: and of the U1 record -- that independence is the point.
CONSENSUS_5SS = "CAGGUAAGU"

pytestmark = pytest.mark.skipif(
    not s5_inputs.U1_GB.exists(), reason="RNU1-1 record not on disk"
)


def test_u1_record_length_is_pinned():
    assert len(load_u1()) == U1_LENGTH


def test_u1_is_rna_not_dna():
    seq = load_u1()
    assert "T" not in seq
    assert set(seq) <= set("ACGU")


def test_u1_arm_pairs_perfectly_with_the_consensus():
    """Fixes the register from OUTSIDE the module. An off-by-one fails here."""
    arm = u1_duplex_arm()
    assert len(arm) == U1_ARM[1] - U1_ARM[0] + 1 == len(CONSENSUS_5SS)
    assert all(pairs_with(arm, CONSENSUS_5SS)), "U1 nt 3-11 must read CAGGUAAGU 9/9"


def test_neighbouring_registers_pair_worse_than_the_pinned_one():
    """The pinned window must be a strict optimum, not merely acceptable."""
    seq = load_u1()
    lo = U1_ARM[0]
    best = sum(pairs_with(seq[lo - 1 : lo - 1 + 9], CONSENSUS_5SS))
    assert best == 9
    for shift in (-2, -1, 1, 2):
        s = lo - 1 + shift
        if s < 0 or s + 9 > len(seq):
            continue
        assert sum(pairs_with(seq[s : s + 9], CONSENSUS_5SS)) < best


def test_smn2_five_ss_mismatches_the_duplex_at_two_positions():
    five_ss = slice_ref(load_smn2(), SS5_EXON7).replace("T", "U")
    paired = pairs_with(u1_duplex_arm(), five_ss)

    assert len(five_ss) == 9
    assert sum(paired) == 7
    unpaired = [i for i, ok in enumerate(paired) if not ok]
    assert unpaired == [0, 2]


def test_the_bulged_minus_one_a_is_reference_32052_and_is_unpaired():
    """The structural premise of S5, re-derived from coordinates and sequence."""
    smn2 = load_smn2()
    five_ss = slice_ref(smn2, SS5_EXON7).replace("T", "U")
    paired = pairs_with(u1_duplex_arm(), five_ss)

    # Third position of the 9-mer, by construction of the interval.
    ref_pos = SS5_EXON7.start + 2
    assert ref_pos == 32052
    assert five_ss[2] == "A", "ref 32052 must be the adenosine"
    assert not paired[2], "the -1 A must be UNPAIRED in the U1 duplex (bulged)"

    # And it is the last base of exon 7, i.e. genuinely position -1.
    from smt.seq.reference import exon
    assert exon("7").end == ref_pos


def test_pairs_with_is_antiparallel_not_parallel():
    """A parallel implementation scores this fixture 0/4; antiparallel gives 4/4."""
    arm = "GGGA"
    site = "UCCC"          # reverse(arm) = AGGG, complement -> UCCC
    assert all(pairs_with(arm, site))
    assert sum(pairs_with(arm, arm[::-1])) < 4


def test_pairs_with_rejects_length_mismatch():
    with pytest.raises(ValueError):
        pairs_with("ACGU", "ACG")


def test_gu_wobble_is_optional_and_actually_toggles():
    assert pairs_with("G", "U", allow_wobble=True) == [True]
    assert pairs_with("G", "U", allow_wobble=False) == [False]


def test_risdiplam_validates_against_pubchem():
    rdkit = pytest.importorskip("rdkit")
    info = validate_ligand()
    assert info["formula"] == "C22H23N7O"
    assert info["inchikey"] == "ASKZRYGFUPSJPN-UHFFFAOYSA-N"
    assert info["heavy_atoms"] == 30


def test_a_corrupted_ligand_is_rejected(monkeypatch):
    """A validator that cannot fail is not a validator."""
    pytest.importorskip("rdkit")
    # Benzene: parses cleanly, wrong molecule.
    monkeypatch.setattr(s5_inputs, "RISDIPLAM_SMILES", "c1ccccc1")
    with pytest.raises(LigandError):
        validate_ligand()


def test_unparseable_ligand_is_rejected(monkeypatch):
    pytest.importorskip("rdkit")
    monkeypatch.setattr(s5_inputs, "RISDIPLAM_SMILES", "not a smiles")
    with pytest.raises(LigandError):
        validate_ligand()
