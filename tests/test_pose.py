"""Tests for S5 pose grading. A grader that cannot fail is not a gate.

The whole point of ``smt.structure.pose`` is to detect a failure mode that
base-pair F1 cannot see: correct fold, wrong ligand placement. So the
load-bearing test here is not that a deposited structure scores well against
itself -- it is that a DELIBERATELY DISPLACED ligand is REJECTED. Without
that, the module could be returning PASS unconditionally and every other test
would still be green.

Covers:
* self-grading a deposited complex is perfect (contact F1 1.000, RMSD 0.000),
  which pins the metric's zero point;
* a ligand translated 10 A FAILS, and fails on the placement criterion
  specifically -- the test that proves the gate can gate;
* a ligand nudged 0.5 A still PASSES, so the threshold is not simply rejecting
  everything that is not bit-identical;
* PER-COPY SCORING IS LIVE. 2GDI holds two TPP copies; pooling them doubles
  the contact count and is exactly the 6DCL bug that made S3's gate first read
  a false MARGINAL 0.400. The test asserts pooled == 2x per-copy, so if
  someone reverts to pooling it fails loudly;
* solvent and ions are not mistaken for the ligand;
* structures that cannot be compared residue-by-residue are REFUSED rather
  than silently mis-indexed, and a ligand-free structure is refused too.
"""

from __future__ import annotations

from pathlib import Path

import pytest

gemmi = pytest.importorskip("gemmi")

from smt.structure.pose import (
    PASS_LIGAND_RMSD,
    CompositionMismatch,
    grade_pose,
    ligand_contacts,
    ligand_residues,
    ligand_rmsd,
)

PDB = Path("/home/dashii/orcd/pool/smt-isef-data/pdb")
TPP = PDB / "2GDI.cif"      # TPP riboswitch, TWO ligand copies
SAM = PDB / "3GX5.cif"      # SAM riboswitch, ONE ligand copy

pytestmark = pytest.mark.skipif(
    not (TPP.exists() and SAM.exists()),
    reason="deposited RNA:ligand test complexes not on disk",
)


def _displace_ligand(src: Path, dst: Path, shift: float) -> Path:
    """Copy ``src`` with every ligand heavy atom translated along +x."""
    st = gemmi.read_structure(str(src))
    st.setup_entities()
    names = {r.name for r in ligand_residues(src)}
    for chain in st[0]:
        for res in chain:
            if res.name in names:
                for atom in res:
                    atom.pos = gemmi.Position(
                        atom.pos.x + shift, atom.pos.y, atom.pos.z)
    st.setup_entities()
    dst.write_text(st.make_mmcif_document().as_string())
    return dst


def test_self_grade_is_perfect():
    g = grade_pose(SAM, SAM, bp_f1=1.0)
    assert g.contact_f1 == pytest.approx(1.0)
    assert g.ligand_rmsd == pytest.approx(0.0, abs=1e-6)
    assert g.passed
    assert g.verdict == "PASS"


def test_a_displaced_ligand_is_REJECTED(tmp_path):
    """The test that proves the gate can gate."""
    moved = _displace_ligand(SAM, tmp_path / "moved.cif", shift=10.0)
    g = grade_pose(SAM, moved, bp_f1=1.0)

    assert not g.passed
    assert g.ligand_rmsd > PASS_LIGAND_RMSD
    assert "placement" in g.verdict


def test_a_slightly_nudged_ligand_still_passes(tmp_path):
    """The threshold must not be rejecting everything non-identical."""
    moved = _displace_ligand(SAM, tmp_path / "nudge.cif", shift=0.5)
    g = grade_pose(SAM, moved, bp_f1=1.0)

    assert g.ligand_rmsd == pytest.approx(0.5, abs=0.05)
    assert g.passed


def test_rmsd_scales_with_the_displacement(tmp_path):
    for shift in (1.0, 3.0, 7.0):
        moved = _displace_ligand(SAM, tmp_path / f"s{shift}.cif", shift)
        rmsd, n = ligand_rmsd(SAM, moved)
        assert rmsd == pytest.approx(shift, abs=0.05)
        assert n > 10


def test_per_copy_scoring_is_live_and_pooling_would_double_count():
    """2GDI holds two TPP. Pooling is the 6DCL bug; this asserts it is not done."""
    assert len(ligand_residues(TPP)) == 2

    pooled = ligand_contacts(TPP)
    copy0 = ligand_contacts(TPP, copy=0)
    copy1 = ligand_contacts(TPP, copy=1)

    assert len(copy0) > 0 and len(copy1) > 0
    assert not (copy0 & copy1), "the two copies must touch disjoint residues"
    assert len(pooled) == len(copy0) + len(copy1)

    # And the graded result uses ONE copy, not the pool.
    g = grade_pose(TPP, TPP, bp_f1=1.0)
    assert g.n_copies == (2, 2)
    assert g.n_ref_contacts == len(copy0)
    assert g.contact_f1 == pytest.approx(1.0), "pooling would cap this at 0.5"


def test_solvent_and_ions_are_not_mistaken_for_the_ligand():
    for path in (TPP, SAM):
        names = {r.name.strip().upper() for r in ligand_residues(path)}
        assert names
        assert not (names & {"HOH", "WAT", "MG", "NA", "K", "CL", "SO4"})


def test_composition_mismatch_is_refused():
    with pytest.raises(CompositionMismatch):
        grade_pose(TPP, SAM, bp_f1=1.0)


def test_ligand_free_structure_is_refused(tmp_path):
    """Grading a pose with no pose is meaningless and must raise."""
    st = gemmi.read_structure(str(SAM))
    st.setup_entities()
    names = {r.name for r in ligand_residues(SAM)}
    for chain in st[0]:
        for i in range(len(chain) - 1, -1, -1):
            if chain[i].name in names:
                del chain[i]
    st.setup_entities()
    stripped = tmp_path / "noligand.cif"
    stripped.write_text(st.make_mmcif_document().as_string())

    assert ligand_residues(stripped) == []
    with pytest.raises(CompositionMismatch):
        grade_pose(SAM, stripped, bp_f1=1.0)


def test_all_three_criteria_are_required():
    """A perfect pose with a broken fold must still FAIL."""
    g = grade_pose(SAM, SAM, bp_f1=0.10)
    assert g.contact_f1 == pytest.approx(1.0)
    assert not g.passed
    assert "fold" in g.verdict


def test_different_ligands_are_never_paired(tmp_path):
    """7QSH carries a uridine AND its real ligand. Comparing one against the
    other is meaningless, not merely noisy, so it must be refused."""
    st = gemmi.read_structure(str(SAM))
    st.setup_entities()
    names = {r.name for r in ligand_residues(SAM)}
    for chain in st[0]:
        for res in chain:
            if res.name in names:
                res.name = "XXX"
    st.setup_entities()
    renamed = tmp_path / "renamed.cif"
    renamed.write_text(st.make_mmcif_document().as_string())

    assert [r.name for r in ligand_residues(renamed)] == ["XXX"]
    with pytest.raises(CompositionMismatch):
        grade_pose(SAM, renamed, bp_f1=1.0)


def test_crystallisation_additives_are_not_ligands():
    """1J8G is held together by spermine; it must not read as a complex."""
    spm = PDB / "1J8G.cif"
    if not spm.exists():
        pytest.skip("1J8G not on disk")
    assert ligand_residues(spm) == []
