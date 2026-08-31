"""Chain identification and guide/target assignment for Cas13 ternary complexes.

Scope note: this module is pure structural bookkeeping for a *non-enzymatic*
steric-blockade model. No turnover chemistry is modelled anywhere in this
package; see ``tests/test_no_catalysis.py``.

Units and conventions
---------------------
* All distances are Angstroms (A).
* ``seqid`` is the author residue number as deposited in the mmCIF.
* ``register`` is the duplex-register coordinate defined by
  :func:`duplex_register`: 0 is the 5'-most base-paired target nucleotide and
  the coordinate advances by +1 per *resolved* target nucleotide toward the 3'
  end.  Nucleotides 5' of the hybridized span therefore carry negative
  register values.  Because it counts resolved nucleotides, a stretch of
  missing density does not advance the coordinate; every such break is
  recorded in :attr:`ComplexRecord.flags`.

Nothing here is hardcoded to a particular chain letter: polymer kind comes
from residue composition, and the guide/target assignment comes from
base-pairing geometry plus the direct-repeat hairpin signature.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Sequence

import gemmi
import numpy as np

# --- Physical constants -----------------------------------------------------

# Bondi van der Waals radii (A).  Source: A. Bondi, J. Phys. Chem. 68, 441
# (1964).  Mg from S. Alvarez, Dalton Trans. 42, 8617 (2013), which extends
# Bondi's set to the metals.  These are the standard radii used by freesasa
# and by essentially every SASA implementation; they are element-level
# approximations that ignore hybridisation and local charge.
VDW_RADII_A: dict[str, float] = {
    "H": 1.20,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "P": 1.80,
    "S": 1.80,
    "MG": 1.73,
}
# Fallback for any element not tabulated above.  Carbon-like; APPROXIMATION.
DEFAULT_VDW_RADIUS_A: float = 1.70

# Maximum donor-acceptor N...N distance accepted as a Watson-Crick pair (A).
# Canonical WC N1(purine)...N3(pyrimidine) hydrogen bonds are 2.8-3.0 A
# (Saenger, "Principles of Nucleic Acid Structure", 1984).  3.4 A is a
# deliberately loose cutoff so that 3.1-3.5 A cryo-EM models still register
# their pairs; it is an APPROXIMATION tuned for tolerance, not precision.
WC_PAIR_CUTOFF_A: float = 3.4

# O3'(i)-P(i+1) covalent bond length is ~1.6 A.  Anything beyond 2.5 A is a
# genuine chain break (unmodelled density), regardless of residue numbering.
BACKBONE_BREAK_CUTOFF_A: float = 2.5

# Minimum separation along a chain before two residues are allowed to be
# called an intramolecular pair; excludes trivial i,i+1 geometry.
_MIN_HAIRPIN_SEPARATION: int = 3

# Minimum consecutive base pairs required to call a stack a "duplex ladder".
_MIN_LADDER_LENGTH: int = 4

# Heavy-atom contact distance used only as a corroborating (never decisive)
# guide/target discriminator (A).
_CONTACT_CUTOFF_A: float = 4.5

_RNA_RESIDUES = {"A", "C", "G", "U", "I", "N"}
_DNA_RESIDUES = {"DA", "DC", "DG", "DT", "DU", "DI", "DN"}

_PURINE_WC_ATOM = "N1"
_PYRIMIDINE_WC_ATOM = "N3"
_PURINES = {"A", "G", "DA", "DG", "I", "DI"}
_PYRIMIDINES = {"C", "U", "T", "DC", "DT", "DU"}


# --- Records ----------------------------------------------------------------


@dataclass(frozen=True)
class ChainRecord:
    """One deposited polymer chain, described independently of its chain ID."""

    chain_id: str
    entity_description: str
    polymer_kind: str  # 'protein' | 'rna' | 'dna' | 'other'
    n_residues: int
    seqid_first: int
    seqid_last: int
    sequence: str
    numbering_gaps: tuple[tuple[int, int], ...]
    backbone_breaks: tuple[tuple[int, int], ...]
    role: str  # 'protein' | 'guide' | 'target' | 'guide+target fusion' | 'ambiguous' | 'other'
    role_evidence: str


@dataclass(frozen=True)
class Segment:
    """A contiguous stretch of one chain, identified by explicit residue numbers.

    A ``Segment`` is how this package refers to "the guide" or "the target"
    even when both live on a single covalently fused chain (8FTI).
    """

    chain_id: str
    seqids: tuple[int, ...]

    def __len__(self) -> int:
        return len(self.seqids)

    @property
    def seqid_range(self) -> tuple[int, int]:
        return (min(self.seqids), max(self.seqids))


@dataclass(frozen=True)
class ComplexRecord:
    """Full classification of one Cas13 ternary-complex coordinate file."""

    pdb_id: str
    path: str
    title: str
    resolution_A: float | None
    chains: tuple[ChainRecord, ...]
    protein_chain_ids: tuple[str, ...]
    guide: Segment | None
    target: Segment | None
    #: (target_seqid, guide_seqid) pairs of the guide:target duplex ladder.
    duplex_pairs: tuple[tuple[int, int], ...]
    fused: bool
    flags: tuple[str, ...] = field(default_factory=tuple)

    @property
    def detected_pair_seqids(self) -> tuple[int, ...]:
        """Target residue numbers with a geometrically detected base pair."""
        return tuple(sorted({t for t, _ in self.duplex_pairs}))

    @property
    def hybridized_target_seqids(self) -> tuple[int, ...]:
        """Contiguous target footprint of the guide:target duplex.

        Spans the resolved target nucleotides from the first to the last
        ladder pair inclusive.  The span is filled in because a duplex is
        physically contiguous: at 3.1-3.5 A a handful of interior pairs fall
        outside :data:`WC_PAIR_CUTOFF_A` and drop out of
        :attr:`detected_pair_seqids`, but they are still hybridized.
        """
        detected = self.detected_pair_seqids
        if not detected or self.target is None:
            return ()
        order = {s: i for i, s in enumerate(self.target.seqids)}
        lo, hi = order[detected[0]], order[detected[-1]]
        return tuple(s for s in self.target.seqids if lo <= order[s] <= hi)

    @property
    def ok(self) -> bool:
        """True when both a guide and a target strand were confidently assigned."""
        return self.guide is not None and self.target is not None


# --- Loading ----------------------------------------------------------------


def load_structure(path: str | Path) -> gemmi.Structure:
    """Read an mmCIF and normalise it for geometry work.

    Removes hydrogens (cryo-EM models at 3.1-3.5 A have none anyway),
    alternative conformations, ligands and waters.  Entities are set up so
    that ``chain.get_polymer()`` is populated.
    """
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.remove_ligands_and_waters()
    if len(st) == 0:
        raise ValueError(f"{path}: no models in file")
    return st


def _entity_descriptions(block: gemmi.cif.Block) -> dict[str, str]:
    """Map entity id -> ``_entity.pdbx_description`` (the depositor's own label)."""
    out: dict[str, str] = {}
    for row in block.find("_entity.", ["id", "pdbx_description"]):
        out[gemmi.cif.as_string(row[0])] = gemmi.cif.as_string(row[1])
    return out


def _entity_description(
    st: gemmi.Structure, chain: gemmi.Chain, by_id: dict[str, str]
) -> str:
    names = {s.subchain_id() for s in chain.subchains() if len(s)}
    for ent in st.entities:
        if names & set(ent.subchains):
            return by_id.get(ent.name, ent.name or "")
    return ""


def _polymer_kind(residues: Sequence[gemmi.Residue]) -> str:
    """Infer polymer kind from residue composition alone."""
    rna = dna = aa = 0
    for r in residues:
        name = r.name.strip().upper()
        if name in _DNA_RESIDUES:
            dna += 1
        elif name in _RNA_RESIDUES:
            rna += 1
        else:
            info = gemmi.find_tabulated_residue(name)
            if info is not None and info.is_amino_acid():
                aa += 1
            elif info is not None and info.is_nucleic_acid():
                rna += 1
    total = max(rna + dna + aa, 1)
    if aa / total > 0.7:
        return "protein"
    if dna > rna:
        return "dna"
    if rna / total > 0.7:
        return "rna"
    return "other"


def _numbering_gaps(residues: Sequence[gemmi.Residue]) -> tuple[tuple[int, int], ...]:
    """Residue-number jumps.  A jump across 0 is a numbering convention, not a gap."""
    out = []
    for a, b in zip(residues, residues[1:]):
        na, nb = a.seqid.num, b.seqid.num
        if nb == na + 1:
            continue
        if na == -1 and nb == 1:  # deposited numbering skips zero
            continue
        out.append((na, nb))
    return tuple(out)


def _backbone_breaks(residues: Sequence[gemmi.Residue]) -> tuple[tuple[int, int], ...]:
    """True unmodelled breaks, detected from O3'(i)-P(i+1) / C(i)-N(i+1) distance."""
    out = []
    for a, b in zip(residues, residues[1:]):
        pa = a.find_atom("O3'", "*") or a.find_atom("C", "*")
        pb = b.find_atom("P", "*") or b.find_atom("N", "*")
        if pa is None or pb is None:
            out.append((a.seqid.num, b.seqid.num))
            continue
        if pa.pos.dist(pb.pos) > BACKBONE_BREAK_CUTOFF_A:
            out.append((a.seqid.num, b.seqid.num))
    return tuple(out)


# --- Base pairing -----------------------------------------------------------


def _wc_atom(residue: gemmi.Residue) -> gemmi.Atom | None:
    name = residue.name.strip().upper()
    if name in _PURINES:
        return residue.find_atom(_PURINE_WC_ATOM, "*")
    if name in _PYRIMIDINES:
        return residue.find_atom(_PYRIMIDINE_WC_ATOM, "*")
    return None


def watson_crick_pairs(
    a: Sequence[gemmi.Residue],
    b: Sequence[gemmi.Residue],
    same_chain: bool = False,
) -> list[tuple[int, int]]:
    """Detect base pairs geometrically from N1(purine)...N3(pyrimidine) distance.

    Geometry, not sequence: this works even where a deposited model contains
    non-canonical or mismatched pairs.  Returns ``(seqid_a, seqid_b)`` with
    distances below :data:`WC_PAIR_CUTOFF_A`.
    """
    pos_a = [(r.seqid.num, _wc_atom(r)) for r in a]
    pos_b = [(r.seqid.num, _wc_atom(r)) for r in b]
    pairs: list[tuple[int, int]] = []
    for na, xa in pos_a:
        if xa is None:
            continue
        for nb, xb in pos_b:
            if xb is None:
                continue
            if same_chain and nb - na < _MIN_HAIRPIN_SEPARATION:
                continue
            if xa.pos.dist(xb.pos) <= WC_PAIR_CUTOFF_A:
                pairs.append((na, nb))
    return pairs


def _index_map(residues: Sequence[gemmi.Residue]) -> dict[int, int]:
    return {r.seqid.num: i for i, r in enumerate(residues)}


def longest_ladder(
    pairs: Iterable[tuple[int, int]],
    idx_a: dict[int, int],
    idx_b: dict[int, int],
    max_slip: int = 2,
) -> list[tuple[int, int]]:
    """Longest antiparallel stack of consecutive base pairs.

    Walks pairs sorted by strand-A index, accepting a step when strand A
    advances and strand B retreats, tolerating up to ``max_slip`` skipped
    positions on either strand (bulges and unmodelled residues).  Returns the
    longest such run, which for a Cas13 ternary complex is the guide:target
    spacer duplex.
    """
    uniq = sorted({p for p in pairs}, key=lambda p: (idx_a[p[0]], -idx_b[p[1]]))
    if not uniq:
        return []
    best: list[tuple[int, int]] = []
    # Dynamic programming over the sorted pair list.
    length = [1] * len(uniq)
    prev = [-1] * len(uniq)
    for j in range(len(uniq)):
        aj, bj = idx_a[uniq[j][0]], idx_b[uniq[j][1]]
        for i in range(j):
            ai, bi = idx_a[uniq[i][0]], idx_b[uniq[i][1]]
            da, db = aj - ai, bi - bj
            if 1 <= da <= 1 + max_slip and 1 <= db <= 1 + max_slip:
                if length[i] + 1 > length[j]:
                    length[j] = length[i] + 1
                    prev[j] = i
    end = int(np.argmax(length))
    chain: list[tuple[int, int]] = []
    while end != -1:
        chain.append(uniq[end])
        end = prev[end]
    best = list(reversed(chain))
    return best if len(best) >= _MIN_LADDER_LENGTH else []


# --- Guide / target assignment ---------------------------------------------


def _protein_contact_fraction(
    residues: Sequence[gemmi.Residue], protein_coords: np.ndarray
) -> float:
    """Fraction of residues with any heavy atom within 4.5 A of the protein."""
    if len(protein_coords) == 0 or not residues:
        return 0.0
    from scipy.spatial import cKDTree

    tree = cKDTree(protein_coords)
    hits = 0
    for r in residues:
        xyz = np.array([[a.pos.x, a.pos.y, a.pos.z] for a in r], dtype=float)
        if len(xyz) and tree.query_ball_point(xyz, _CONTACT_CUTOFF_A, return_length=True).sum():
            hits += 1
    return hits / len(residues)


def _hairpin_count(residues: Sequence[gemmi.Residue], exclude: set[int]) -> int:
    """Intramolecular base pairs outside the spacer duplex.

    The crRNA direct repeat folds into a hairpin that is clamped by the
    protein; the target strand has no such element.  This is the primary
    guide/target discriminator and it is structural, not sequence-based.
    """
    subset = [r for r in residues if r.seqid.num not in exclude]
    return len(watson_crick_pairs(subset, subset, same_chain=True))


def classify_structure(path: str | Path) -> ComplexRecord:
    """Identify protein / guide / target from residue composition and geometry.

    Handles three deposited architectures without special-casing chain IDs:

    1. protein + separate crRNA chain + separate target chain;
    2. protein + a single covalently fused ``crRNA-linker-target`` chain, which
       is split at the apex of the intramolecular spacer duplex;
    3. anything else, which is flagged rather than guessed at.
    """
    path = Path(path)
    st = load_structure(path)
    pdb_id = st.name or path.stem
    title = ""
    doc = gemmi.cif.read(str(path))
    block = doc.sole_block()
    raw_title = block.find_value("_struct.title")
    if raw_title:
        title = gemmi.cif.as_string(raw_title)
    entity_names = _entity_descriptions(block)

    flags: list[str] = []
    model = st[0]

    polymers: list[tuple[gemmi.Chain, list[gemmi.Residue], str]] = []
    for ch in model:
        pol = ch.get_polymer()
        residues = list(pol)
        if not residues:
            continue
        polymers.append((ch, residues, _polymer_kind(residues)))

    protein_coords = np.array(
        [
            [a.pos.x, a.pos.y, a.pos.z]
            for _, residues, kind in polymers
            if kind == "protein"
            for r in residues
            for a in r
        ],
        dtype=float,
    ).reshape(-1, 3)

    rna_chains = [(ch, residues) for ch, residues, kind in polymers if kind == "rna"]

    guide: Segment | None = None
    target: Segment | None = None
    duplex_pairs: tuple[tuple[int, int], ...] = ()
    fused = False
    evidence: dict[str, str] = {}

    if len(rna_chains) >= 2:
        guide, target, duplex_pairs, notes = _assign_two_chain(
            rna_chains, protein_coords
        )
        flags.extend(notes)
        if guide is not None:
            evidence[guide.chain_id] = "guide"
        if target is not None:
            evidence[target.chain_id] = "target"
    elif len(rna_chains) == 1:
        ch, residues = rna_chains[0]
        guide, target, duplex_pairs, notes = _assign_fused(ch, residues, protein_coords)
        fused = guide is not None and target is not None and guide.chain_id == target.chain_id
        flags.extend(notes)
        if fused:
            evidence[ch.name] = "guide+target fusion"
    else:
        flags.append("no RNA polymer chain found; cannot assign guide or target")

    chain_records: list[ChainRecord] = []
    for ch, residues, kind in polymers:
        if kind == "protein":
            role, why = "protein", "residue composition is >70% amino acid"
        elif ch.name in evidence:
            role = evidence[ch.name]
            why = _ROLE_EVIDENCE[role]
        elif kind in ("rna", "dna"):
            role, why = "ambiguous", "nucleic acid chain not assigned to the spacer duplex"
            flags.append(f"chain {ch.name}: nucleic acid with no guide/target assignment")
        else:
            role, why = "other", "unrecognised polymer"
        seq = gemmi.one_letter_code([r.name for r in residues])
        breaks = _backbone_breaks(residues)
        if breaks:
            flags.append(
                f"chain {ch.name}: unmodelled backbone break(s) at "
                + ", ".join(f"{a}->{b}" for a, b in breaks)
            )
        chain_records.append(
            ChainRecord(
                chain_id=ch.name,
                entity_description=_entity_description(st, ch, entity_names),
                polymer_kind=kind,
                n_residues=len(residues),
                seqid_first=residues[0].seqid.num,
                seqid_last=residues[-1].seqid.num,
                sequence=seq.upper(),
                numbering_gaps=_numbering_gaps(residues),
                backbone_breaks=breaks,
                role=role,
                role_evidence=why,
            )
        )

    return ComplexRecord(
        pdb_id=pdb_id.upper(),
        path=str(path),
        title=title,
        resolution_A=st.resolution or None,
        chains=tuple(chain_records),
        protein_chain_ids=tuple(c.chain_id for c in chain_records if c.role == "protein"),
        guide=guide,
        target=target,
        duplex_pairs=duplex_pairs,
        fused=fused,
        flags=tuple(flags),
    )


_ROLE_EVIDENCE = {
    "guide": "carries the direct-repeat hairpin and pairs with the target via its spacer",
    "target": "pairs with the guide spacer and carries no direct-repeat hairpin",
    "guide+target fusion": "single chain whose two arms form the spacer duplex intramolecularly",
}


def _assign_two_chain(
    rna_chains: list[tuple[gemmi.Chain, list[gemmi.Residue]]],
    protein_coords: np.ndarray,
) -> tuple[Segment | None, Segment | None, tuple[tuple[int, int], ...], list[str]]:
    notes: list[str] = []
    best: tuple[int, tuple, tuple] | None = None
    for i in range(len(rna_chains)):
        for j in range(i + 1, len(rna_chains)):
            (ci, ri), (cj, rj) = rna_chains[i], rna_chains[j]
            pairs = watson_crick_pairs(ri, rj)
            ladder = longest_ladder(pairs, _index_map(ri), _index_map(rj))
            if best is None or len(ladder) > best[0]:
                best = (len(ladder), (ci, ri), (cj, rj)), ladder  # type: ignore[assignment]
    if best is None:
        return None, None, (), ["fewer than two RNA chains to pair"]
    (nlad, (ca, ra), (cb, rb)), ladder = best  # type: ignore[misc]
    if nlad < _MIN_LADDER_LENGTH:
        return None, None, (), ["no guide:target duplex of >=4 bp found between RNA chains"]

    paired_a = {p[0] for p in ladder}
    paired_b = {p[1] for p in ladder}
    hp_a = _hairpin_count(ra, paired_a)
    hp_b = _hairpin_count(rb, paired_b)
    con_a = _protein_contact_fraction(ra, protein_coords)
    con_b = _protein_contact_fraction(rb, protein_coords)

    if hp_a != hp_b:
        guide_is_a = hp_a > hp_b
        why = f"direct-repeat hairpin pairs {hp_a} vs {hp_b}"
    elif abs(con_a - con_b) > 0.10:
        guide_is_a = con_a > con_b
        why = f"protein-contact fraction {con_a:.2f} vs {con_b:.2f} (hairpin test tied)"
        notes.append(f"chains {ca.name}/{cb.name}: guide assigned on protein contact, not hairpin")
    else:
        notes.append(
            f"AMBIGUOUS guide/target for chains {ca.name}/{cb.name}: "
            f"hairpin {hp_a} vs {hp_b}, contact {con_a:.2f} vs {con_b:.2f}"
        )
        return None, None, (), notes

    if guide_is_a:
        guide_ch, guide_res, tgt_ch, tgt_res = ca, ra, cb, rb
        pairs_out = tuple((t, g) for g, t in ladder)
    else:
        guide_ch, guide_res, tgt_ch, tgt_res = cb, rb, ca, ra
        pairs_out = tuple((t, g) for t, g in ladder)
    notes.append(
        f"guide={guide_ch.name} target={tgt_ch.name} via {why}; "
        f"spacer duplex {len(ladder)} bp"
    )
    return (
        Segment(guide_ch.name, tuple(r.seqid.num for r in guide_res)),
        Segment(tgt_ch.name, tuple(r.seqid.num for r in tgt_res)),
        pairs_out,
        notes,
    )


def _assign_fused(
    ch: gemmi.Chain,
    residues: list[gemmi.Residue],
    protein_coords: np.ndarray,
) -> tuple[Segment | None, Segment | None, tuple[tuple[int, int], ...], list[str]]:
    """Split a covalently fused crRNA-linker-target chain into its two arms.

    The spacer duplex is formed intramolecularly, so the chain folds back on
    itself.  The longest intramolecular ladder identifies the two arms; the
    chain is cut at the apex between them.  The arm carrying the remaining
    direct-repeat hairpin is the guide.
    """
    notes: list[str] = []
    idx = _index_map(residues)
    pairs = watson_crick_pairs(residues, residues, same_chain=True)
    ladder = longest_ladder(pairs, idx, idx)
    if len(ladder) < _MIN_LADDER_LENGTH:
        return None, None, (), [
            f"chain {ch.name}: single RNA chain with no intramolecular duplex; "
            "cannot separate guide from target"
        ]
    five_side = [p[0] for p in ladder]
    three_side = [p[1] for p in ladder]
    apex_lo, apex_hi = max(five_side), min(three_side)
    cut = (idx[apex_lo] + idx[apex_hi]) / 2.0
    arm5 = [r for r in residues if idx[r.seqid.num] <= cut]
    arm3 = [r for r in residues if idx[r.seqid.num] > cut]
    notes.append(
        f"chain {ch.name} is a fused crRNA-linker-target: split after resolved "
        f"index {math.floor(cut)} (seqid {apex_lo} | {apex_hi}); "
        f"arms {len(arm5)} nt and {len(arm3)} nt, spacer duplex {len(ladder)} bp"
    )

    hp5 = _hairpin_count(arm5, set(five_side))
    hp3 = _hairpin_count(arm3, set(three_side))
    con5 = _protein_contact_fraction(arm5, protein_coords)
    con3 = _protein_contact_fraction(arm3, protein_coords)
    if hp5 != hp3:
        guide_is_5 = hp5 > hp3
        why = f"direct-repeat hairpin pairs {hp5} vs {hp3}"
    elif abs(con5 - con3) > 0.10:
        guide_is_5 = con5 > con3
        why = f"protein-contact fraction {con5:.2f} vs {con3:.2f} (hairpin test tied)"
        notes.append(f"chain {ch.name}: guide arm assigned on protein contact, not hairpin")
    else:
        notes.append(
            f"AMBIGUOUS guide/target arms in fused chain {ch.name}: "
            f"hairpin {hp5} vs {hp3}, contact {con5:.2f} vs {con3:.2f}"
        )
        return None, None, (), notes

    guide_res, tgt_res = (arm5, arm3) if guide_is_5 else (arm3, arm5)
    pairs_out = (
        tuple((g, t) for g, t in ladder) if not guide_is_5 else tuple((t, g) for g, t in ladder)
    )
    notes.append(f"guide arm = {'5' if guide_is_5 else '3'}' arm via {why}")
    return (
        Segment(ch.name, tuple(r.seqid.num for r in guide_res)),
        Segment(ch.name, tuple(r.seqid.num for r in tgt_res)),
        pairs_out,
        notes,
    )


# --- Register ---------------------------------------------------------------


def duplex_register(record: ComplexRecord) -> dict[int, int]:
    """Map target ``seqid`` -> duplex-register coordinate.

    Register 0 is the 5'-most base-paired target nucleotide; the coordinate
    advances +1 per resolved target nucleotide toward the 3' end (so it is a
    resolved-nucleotide index, not an arithmetic seqid offset, and unmodelled
    residues do not advance it).  This puts every structure on one axis:
    0..N-1 spans the hybridized footprint, negatives are 5' flank, values
    beyond the footprint are 3' flank.
    """
    if record.target is None or not record.duplex_pairs:
        return {}
    seqids = list(record.target.seqids)
    order = {s: i for i, s in enumerate(seqids)}
    paired = [s for s in record.hybridized_target_seqids if s in order]
    if not paired:
        return {}
    anchor = order[min(paired, key=lambda s: order[s])]
    return {s: order[s] - anchor for s in seqids}


def describe(record: ComplexRecord) -> str:
    """Human-readable classification report for one structure."""
    lines = [
        f"{record.pdb_id}  {record.title}",
        f"  resolution: {record.resolution_A} A   fused guide/target: {record.fused}",
    ]
    for c in record.chains:
        lines.append(
            f"  chain {c.chain_id:<3} {c.polymer_kind:<8} n={c.n_residues:<4} "
            f"seqid {c.seqid_first}..{c.seqid_last}  role={c.role}"
        )
        lines.append(f"        entity: {c.entity_description}")
        lines.append(f"        why:    {c.role_evidence}")
    if record.guide:
        lines.append(f"  guide  segment: chain {record.guide.chain_id} {record.guide.seqid_range} ({len(record.guide)} nt)")
    if record.target:
        lines.append(f"  target segment: chain {record.target.chain_id} {record.target.seqid_range} ({len(record.target)} nt)")
    span = record.hybridized_target_seqids
    if span:
        lines.append(
            f"  spacer duplex: {len(record.duplex_pairs)} pairs detected, "
            f"contiguous footprint = {len(span)} nt over target seqid {span[0]}..{span[-1]}"
        )
    for f in record.flags:
        lines.append(f"  ! {f}")
    return "\n".join(lines)


def classify_all(pdb_dir: str | Path, pattern: str = "*.cif") -> list[ComplexRecord]:
    """Classify every coordinate file in a directory (sorted by filename)."""
    return [classify_structure(p) for p in sorted(Path(pdb_dir).glob(pattern))]
