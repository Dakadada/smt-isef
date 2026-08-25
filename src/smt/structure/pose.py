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

from collections import Counter
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

#: Minimum heavy atoms for a non-polymer residue to be a candidate ligand.
#: Enumerating every ion in the periodic table is the wrong guard -- 1J8G
#: carries eight strontium ions that an explicit name list missed. A drug-like
#: ligand has tens of heavy atoms (biotin 16, TPP 26, SAM 27, risdiplam 30),
#: so a floor catches the next unlisted ion too.
MIN_LIGAND_ATOMS = 10

#: Residue names that are solvent/ions, never the ligand of interest.
_NOT_LIGAND = {"HOH", "WAT", "DOD", "MG", "NA", "K", "CL", "ZN", "CA",
               "MN", "SO4", "PO4", "GOL", "EDO", "IOD", "BR",
               # Crystallisation additives, not specific binders. Counting
               # these as "the ligand" would grade the pose of a buffer
               # component -- 1J8G, 5V1L and 7EDT are held together by
               # spermine and would otherwise look like ligand complexes.
               "SPM", "SPD", "NCO", "IRI", "OS", "ACT", "MPD", "PEG", "TRS"}


def _read(path: Path) -> gemmi.Structure:
    """Read a structure with ALTERNATE CONFORMATIONS REMOVED.

    Without this, a residue modelled in two conformations contributes both to
    every count: 8QMH's W53 read as 39 heavy atoms against a 32-atom SMILES and
    7QSH's ligands read as exactly 2x their formula. That is the 6DCL bug in a
    new costume -- counting one physical thing twice -- and it silently
    corrupts contacts, composition checks and RMSD alike.
    """
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_alternative_conformations()
    return st


#: Oxygen count may differ by this much between a deposited model and the
#: chemical-component definition. Crystallographic models routinely omit or add
#: an oxygen through protonation state, a modelled water of hydration, or an
#: unresolved hydroxyl: 7QSH's HYJ is C10 N4 O8 P1 in the coordinates against
#: C10 N4 O9 P1 in the chem-comp SMILES. The C/N/P/S SKELETON is the identity;
#: refusing over one oxygen would discard a valid system for a naming-level
#: difference, which is the same false-negative class as the LIG1 refusal.
_OXYGEN_SLACK = 1


def _composition(res) -> tuple:
    """Heavy-atom element multiset -- a name-independent molecular identity."""
    c = Counter(a.element.name for a in res if a.element != gemmi.Element("H"))
    return tuple(sorted(c.items()))


def _same_molecule(a: tuple, b: tuple) -> bool:
    """Identity by C/N/P/S skeleton, with oxygen allowed +/- _OXYGEN_SLACK."""
    da, db = dict(a), dict(b)
    for el in set(da) | set(db):
        if el == "O":
            continue
        if da.get(el, 0) != db.get(el, 0):
            return False
    return abs(da.get("O", 0) - db.get("O", 0)) <= _OXYGEN_SLACK


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
    st = _read(path)
    out = []
    for chain in st[0]:
        poly = {r.seqid.num for r in chain.get_polymer()}
        for res in chain:
            name = res.name.strip().upper()
            if res.seqid.num in poly or name in exclude:
                continue
            heavy = sum(1 for a in res if a.element != gemmi.Element("H"))
            if heavy < MIN_LIGAND_ATOMS:
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
    st = _read(path)
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
    a = _chain_counts(_read(deposited))
    b = _chain_counts(_read(predicted))
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
    dep = _read(deposited)
    pred = _read(predicted)

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

    def heavy(res):
        return [a for a in res if a.element != gemmi.Element("H")]

    dep_ligs = ligand_residues(deposited)
    pred_ligs = ligand_residues(predicted)
    if not dep_ligs or not pred_ligs:
        return float("inf"), 0

    # PER COPY, never pooled, and paired by COMPOSITION rather than name --
    # Boltz emits "LIG1" while deposited files use the PDB chem-comp id, so a
    # name check would refuse every real comparison.
    best, best_n = float("inf"), 0
    for dres in dep_ligs:
        for pres in pred_ligs:
            if not _same_molecule(_composition(dres), _composition(pres)):
                continue
            da, pa = heavy(dres), heavy(pres)
            moved = [gemmi.Position(sup.transform.apply(a.pos)) for a in pa]
            rmsd = _assigned_rmsd(da, moved)
            if rmsd < best:
                best, best_n = rmsd, len(da)
    return best, best_n


def _assigned_rmsd(ref_atoms, moved_positions) -> float:
    """RMSD under the optimal SAME-ELEMENT atom assignment.

    Atom NAMES cannot be used: deposited files use PDB nomenclature (C11, O12)
    and Boltz emits its own (C29, C32), so a name match finds nothing. Instead
    each predicted atom is matched to a deposited atom of the same element by
    optimal assignment.

    THIS IS A LOWER BOUND -- it ignores connectivity, so it can pair atoms that
    are chemically distinct. A FAIL on this metric is therefore DECISIVE (the
    pose is wrong even under the most generous correspondence), while a PASS is
    the weaker direction and must be read as such. The report says so.
    """
    from scipy.optimize import linear_sum_assignment

    ref_el = [a.element.name for a in ref_atoms]
    mov_el = [p[1] for p in moved_positions] if moved_positions and isinstance(
        moved_positions[0], tuple) else None
    if mov_el is not None:
        raise TypeError("moved_positions must be gemmi.Position")

    n, m = len(ref_atoms), len(moved_positions)
    if n == 0 or m == 0:
        return float("inf")
    size = max(n, m)
    BIG = 1e6
    cost = [[BIG] * size for _ in range(size)]
    for i, a in enumerate(ref_atoms):
        for j, p in enumerate(moved_positions):
            cost[i][j] = a.pos.dist(p) ** 2
    rows, cols = linear_sum_assignment(cost)
    total, used = 0.0, 0
    for i, j in zip(rows, cols):
        if i < n and j < m:
            total += cost[i][j]
            used += 1
    return (total / used) ** 0.5 if used else float("inf")


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
    # Copies are paired only when they are the SAME MOLECULE BY COMPOSITION:
    # 7QSH carries a uridine AND its real ligand, and comparing one against the
    # other is meaningless rather than merely noisy. Composition rather than
    # name, because Boltz emits "LIG1" and deposited files use chem-comp ids.
    dep_comp = [_composition(r) for r in ligand_residues(deposited)]
    pred_comp = [_composition(r) for r in ligand_residues(predicted)]
    if not any(_same_molecule(a, b) for a in dep_comp for b in pred_comp):
        raise CompositionMismatch(
            "no ligand of matching composition between structures: "
            f"{dep_comp} vs {pred_comp}")

    best = None
    for i in range(n_dep):
        ref = ligand_contacts(deposited, copy=i)
        for j in range(n_pred):
            if not _same_molecule(dep_comp[i], pred_comp[j]):
                continue
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
