"""S5 pose grading: is the LIGAND in the right place, not just the fold?

WHY THIS MODULE HAD TO EXIST. ``smt.structure.gate`` grades base pairs, and
base-pair F1 grades the RNA FOLD. A predicted complex can score F1 = 1.000
with the ligand sitting in entirely the wrong pocket, or in solvent. Reading
an S5 prediction on base-pair F1 alone would therefore be a gate that cannot
detect the failure S5 is specifically about.

THREE METRICS, ALL REQUIRED, ALL DECLARED BEFORE ANY PREDICTION IS READ:

  base-pair F1  >= 0.70   the fold          (inherits gate.PASS_F1)
  contact F1    >= 0.70   the pose          (which residues the ligand touches)
  ligand RMSD   <= 2.0 A  the placement     (standard docking success criterion)

A prediction passes only if ALL THREE hold. Any one of them alone is
game-able: contacts without RMSD tolerates a flipped ligand in the right
pocket; RMSD without contacts is undefined when the pocket itself moved.

TWO HONEST LIMITS, stated here rather than appended to the report later.

1. LIGAND RMSD IS AN UPPER BOUND, not a symmetry-corrected value. Atoms are
   matched by NAME. A molecule with a symmetric substituent (risdiplam has a
   spiro-cyclopropyl and two methyls) can be chemically identical under an
   atom relabelling that this metric scores as displaced. So a LOW RMSD is
   trustworthy and a HIGH one is not conclusive on its own -- which is the
   safe direction for a gate, since it can only make passing harder.
2. LIGAND COPIES ARE SCORED PER COPY, NEVER POOLED. An asymmetric unit often
   holds several copies of the same ligand (2GDI holds two TPP). Pooling their
   contacts caps recall at 1/n_copies and reads as a bad pose when the pose is
   fine -- this is exactly the 6DCL failure that made S3's gate first report a
   false MARGINAL 0.400. Every copy pairing is scored and the BEST is reported,
   with the copy count printed so pooling can never hide.

3. Residues are matched by POSITION IN POLYMER ORDER across chains, not by
   author numbering, because predicted and deposited files rarely agree on
   chain names or numbering. ``check_composition`` asserts the two structures
   have the same per-chain residue counts and refuses to grade if they do not,
   rather than silently comparing misaligned indices.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import gemmi

from smt.structure.gate import PASS_F1

#: Heavy-atom contact distance, angstroms. 4.5 A is the usual cutoff for
#: "in contact" in structure-based ligand analysis.
CONTACT_CUTOFF = 4.5

#: Declared before grading, all three required.
PASS_CONTACT_F1 = 0.70
PASS_LIGAND_RMSD = 2.0

#: Residue names that are solvent/ions, never the ligand of interest.
_NOT_LIGAND = {"HOH", "WAT", "DOD", "MG", "NA", "K", "CL", "ZN", "CA",
               "MN", "SO4", "PO4", "GOL", "EDO", "IOD", "BR"}


class CompositionMismatch(RuntimeError):
    """Two structures cannot be compared residue-by-residue. Stop, do not guess."""


@dataclass(frozen=True)
class PoseGrade:
    bp_f1: float
    contact_f1: float
    ligand_rmsd: float
    n_ref_contacts: int
    n_pred_contacts: int
    n_ligand_atoms: int
    n_copies: tuple[int, int] = (1, 1)

    @property
    def passed(self) -> bool:
        return (self.bp_f1 >= PASS_F1
                and self.contact_f1 >= PASS_CONTACT_F1
                and self.ligand_rmsd <= PASS_LIGAND_RMSD)

    @property
    def verdict(self) -> str:
        if self.passed:
            return "PASS"
        failed = []
        if self.bp_f1 < PASS_F1:
            failed.append(f"fold (bp F1 {self.bp_f1:.3f})")
        if self.contact_f1 < PASS_CONTACT_F1:
            failed.append(f"pose (contact F1 {self.contact_f1:.3f})")
        if self.ligand_rmsd > PASS_LIGAND_RMSD:
            failed.append(f"placement (RMSD {self.ligand_rmsd:.2f} A)")
        return "FAIL on " + ", ".join(failed)


def _polymer_residues(st: gemmi.Structure) -> list:
    """All polymer residues, in chain order then sequence order."""
    out = []
    for chain in st[0]:
        for res in chain.get_polymer():
            out.append(res)
    return out


def _chain_counts(st: gemmi.Structure) -> list[int]:
    return [len(chain.get_polymer()) for chain in st[0]
            if len(chain.get_polymer()) > 0]


def ligand_residues(path: Path, exclude: set[str] | None = None) -> list:
    """Non-polymer, non-solvent residues -- the ligand(s)."""
    exclude = _NOT_LIGAND if exclude is None else exclude
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    out = []
    for chain in st[0]:
        poly = {r.seqid.num for r in chain.get_polymer()}
        for res in chain:
            name = res.name.strip().upper()
            if res.seqid.num in poly or name in exclude:
                continue
            if not res.het_flag == "H" and len(res) < 5:
                continue
            out.append(res)
    return out


def ligand_contacts(path: Path, cutoff: float = CONTACT_CUTOFF,
                    copy: int | None = None) -> set[int]:
    """Polymer-residue indices within ``cutoff`` of a ligand copy's heavy atoms.

    ``copy`` selects ONE ligand copy by index; ``None`` pools every copy, which
    is correct only when there is exactly one. See the module docstring on why
    pooling multiple copies is a bug, not a convenience.

    Indices are positions in polymer order across chains -- see the module
    docstring on why author numbering is not used.
    """
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    residues = _polymer_residues(st)
    ligands = ligand_residues(path)
    if copy is not None:
        ligands = [ligands[copy]]

    lig_atoms = [a.pos for res in ligands for a in res
                 if a.element != gemmi.Element("H")]
    if not lig_atoms:
        return set()

    contacts: set[int] = set()
    for idx, res in enumerate(residues):
        for atom in res:
            if atom.element == gemmi.Element("H"):
                continue
            if any(atom.pos.dist(p) <= cutoff for p in lig_atoms):
                contacts.add(idx)
                break
    return contacts


def check_composition(deposited: Path, predicted: Path) -> None:
    """Refuse to grade structures whose polymers do not correspond."""
    a = _chain_counts(gemmi.read_structure(str(deposited)))
    b = _chain_counts(gemmi.read_structure(str(predicted)))
    if a != b:
        raise CompositionMismatch(
            f"per-chain polymer residue counts differ: {a} vs {b}. "
            "Residue indices would not correspond; grading refused."
        )


def ligand_rmsd(deposited: Path, predicted: Path) -> tuple[float, int]:
    """Ligand heavy-atom RMSD after superposing on polymer C1' atoms.

    Returns ``(rmsd, n_atoms_matched)``. Upper bound only -- atoms are matched
    by name, so a symmetry-equivalent relabelling reads as displacement.
    """
    dep = gemmi.read_structure(str(deposited))
    pred = gemmi.read_structure(str(predicted))
    dep.setup_entities()
    pred.setup_entities()

    dep_res, pred_res = _polymer_residues(dep), _polymer_residues(pred)
    if len(dep_res) != len(pred_res):
        raise CompositionMismatch(
            f"polymer lengths differ: {len(dep_res)} vs {len(pred_res)}")

    ref_pos, mov_pos = [], []
    for a, b in zip(dep_res, pred_res):
        pa, pb = a.find_atom("C1'", "*"), b.find_atom("C1'", "*")
        if pa is not None and pb is not None:
            ref_pos.append(pa.pos)
            mov_pos.append(pb.pos)
    if len(ref_pos) < 3:
        raise CompositionMismatch("fewer than 3 C1' atoms to superpose on")

    sup = gemmi.superpose_positions(ref_pos, mov_pos)

    def atoms_of(res):
        return {a.name: a.pos for a in res if a.element != gemmi.Element("H")}

    dep_copies = [atoms_of(r) for r in ligand_residues(deposited)]
    pred_copies = [atoms_of(r) for r in ligand_residues(predicted)]
    if not dep_copies or not pred_copies:
        return float("inf"), 0

    # PER COPY, never pooled. Best pairing wins; see the module docstring.
    best, best_n = float("inf"), 0
    for dep_lig in dep_copies:
        for pred_lig in pred_copies:
            shared = sorted(set(dep_lig) & set(pred_lig))
            if not shared:
                continue
            total = 0.0
            for name in shared:
                moved = gemmi.Position(sup.transform.apply(pred_lig[name]))
                total += dep_lig[name].dist(moved) ** 2
            rmsd = (total / len(shared)) ** 0.5
            if rmsd < best:
                best, best_n = rmsd, len(shared)
    return best, best_n


def grade_pose(deposited: Path, predicted: Path, bp_f1: float) -> PoseGrade:
    """All three metrics. ``bp_f1`` comes from ``smt.structure.gate.grade``."""
    check_composition(deposited, predicted)

    n_dep = len(ligand_residues(deposited))
    n_pred = len(ligand_residues(predicted))
    if n_dep == 0 or n_pred == 0:
        raise CompositionMismatch(
            f"no ligand found (deposited {n_dep}, predicted {n_pred}); "
            "pose grading is meaningless without one")

    # PER COPY, never pooled -- the 6DCL lesson. Best pairing wins.
    best = None
    for i in range(n_dep):
        ref = ligand_contacts(deposited, copy=i)
        for j in range(n_pred):
            pred = ligand_contacts(predicted, copy=j)
            tp = len(ref & pred)
            prec = tp / len(pred) if pred else 0.0
            rec = tp / len(ref) if ref else 0.0
            f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
            if best is None or f1 > best[0]:
                best = (f1, len(ref), len(pred))

    f1, n_ref, n_pred_c = best
    rmsd, n_atoms = ligand_rmsd(deposited, predicted)
    return PoseGrade(bp_f1=bp_f1, contact_f1=f1, ligand_rmsd=rmsd,
                     n_ref_contacts=n_ref, n_pred_contacts=n_pred_c,
                     n_ligand_atoms=n_atoms, n_copies=(n_dep, n_pred))
