"""Per-nucleotide steric shadow of an RNA-bound effector.

Scope
-----
This package models NON-ENZYMATIC steric blockade only.  The effector is
assumed catalytically dead: it occupies space and nothing else.  No turnover
chemistry, and no RNase-H-style degradation, is modelled or modellable here;
that is deliberately out of scope (see ``tests/test_no_catalysis.py``).

What a "steric shadow" is
-------------------------
``s(i)`` is the fraction of nucleotide *i*'s solvent-accessible surface that
the effector removes:

    s(i) = 1 - SASA_bound(i) / SASA_free(i)

Both terms are rolling-probe solvent-accessible surface areas (Lee & Richards,
J. Mol. Biol. 55, 379 (1971); Shrake & Rupley, J. Mol. Biol. 79, 351 (1973))
computed over the *same* target atoms, so the ratio is like-for-like.
``SASA_free`` uses the target strand alone; ``SASA_bound`` uses the target
plus whatever the occluding body is.

The probe radius is the whole point.  A water probe (1.4 A) only reports
direct contact.  The physiological competitor for a splice site is U1 snRNP,
a ~250 kDa particle whose own radius is tens of Angstroms, and a probe that
large is excluded from a wide skirt of space around the effector.  The shadow
therefore extends well past the hybridized footprint, and asymmetrically.

Units and conventions
---------------------
* Distances Angstroms (A), areas A^2, masses kDa.
* ``seqid`` is the deposited author residue number; ``register`` is the
  duplex-register coordinate from
  :func:`smt.steric.structures.duplex_register`.
* ``s(i)`` is dimensionless in [0, 1].
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Literal, Sequence

import gemmi
import numpy as np

from .structures import (
    VDW_RADII_A,
    DEFAULT_VDW_RADIUS_A,
    ComplexRecord,
    classify_structure,
    duplex_register,
    load_structure,
)

# --- Probe radii ------------------------------------------------------------

# Water probe.  Lee & Richards (1971) used 1.4 A for the radius of water and
# every SASA convention since has kept it.
WATER_PROBE_RADIUS_A: float = 1.4

# Mean partial specific volume of globular protein, expressed as a density.
# 1.35 g/cm^3 is the large-protein asymptote reported by Fischer, Polikarpov &
# Craievich, Protein Sci. 13, 2825 (2004); small proteins run nearer 1.41.
# APPROXIMATION: it ignores hydration shell and any RNA content, and U1 snRNP
# is roughly one third RNA (denser, ~1.7 g/cm^3), so this slightly
# over-estimates the volume and hence the radius of a real snRNP.
PROTEIN_DENSITY_G_PER_CM3: float = 1.35

# Avogadro's number (CODATA 2018, exact by definition of the mole).
AVOGADRO: float = 6.02214076e23

# Mass of the U1 snRNP particle in kDa.  U1 is ~1.7 MDa in some assemblies but
# the core 17S/12S U1 particle (U1 snRNA + seven Sm proteins + U1-70K, U1-A,
# U1-C) is ~250 kDa.  APPROXIMATION: this is the figure the project brief
# specifies and it is quoted at one significant figure.
U1_SNRNP_MASS_KDA: float = 250.0

# Axial rise per base pair in A-form RNA duplex, A.  Saenger, "Principles of
# Nucleic Acid Structure" (1984), Table 9-3: A-RNA rise = 2.81 A.  Used only
# to convert a spatial extent into an equivalent nucleotide count, and flagged
# as an approximation wherever it is applied.
A_FORM_RISE_PER_NT_A: float = 2.81


def probe_radius_from_mass(
    mass_kda: float, density_g_per_cm3: float = PROTEIN_DENSITY_G_PER_CM3
) -> float:
    """Radius (A) of an equivalent-volume sphere for a particle of ``mass_kda``.

    V = M / (rho * N_A), then R = (3V / 4pi)^(1/3), with 1 cm^3 = 1e24 A^3.

    APPROXIMATION on two counts: it assumes the particle is a compact sphere
    (U1 snRNP is not - it is a lopsided particle with a protruding stem-loop),
    and it assumes a single uniform density.  Treat the result as an upper
    bound on the effective steric probe, not a measurement.
    """
    if mass_kda <= 0:
        raise ValueError("mass_kda must be positive")
    volume_A3 = (mass_kda * 1000.0) / (density_g_per_cm3 * AVOGADRO) * 1e24
    return float((3.0 * volume_A3 / (4.0 * np.pi)) ** (1.0 / 3.0))


#: ~41.9 A for 250 kDa at 1.35 g/cm^3.  See :func:`probe_radius_from_mass`.
U1_SNRNP_PROBE_RADIUS_A: float = probe_radius_from_mass(U1_SNRNP_MASS_KDA)

#: Ladder spanning a water probe up to a U1-snRNP-sized probe (A).
DEFAULT_PROBE_RADII_A: tuple[float, ...] = (
    WATER_PROBE_RADIUS_A,
    3.0,
    6.0,
    12.0,
    20.0,
    round(U1_SNRNP_PROBE_RADIUS_A, 1),
)

#: Shrake-Rupley sphere sampling density.  200 points gives ~1-2% relative
#: error on a per-atom area, which is far below the differences this module
#: reports.  Raise it for publication figures.
DEFAULT_SPHERE_POINTS: int = 200

State = Literal["free", "duplex", "rnp"]


# --- Geometry core ----------------------------------------------------------


def _fibonacci_sphere(n: int) -> np.ndarray:
    """``n`` near-uniform unit vectors via the golden-angle spiral."""
    i = np.arange(n, dtype=np.float64) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    theta = np.pi * (1.0 + 5.0**0.5) * i
    return np.stack(
        [np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1
    )


def sasa_for_subset(
    coords: np.ndarray,
    radii: np.ndarray,
    subset: np.ndarray,
    probe_radius_A: float,
    n_points: int = DEFAULT_SPHERE_POINTS,
) -> np.ndarray:
    """Shrake-Rupley SASA (A^2) for ``subset`` atoms occluded by all of ``coords``.

    A test point on the expanded sphere of radius ``r_i + probe`` is accessible
    when it lies outside every other expanded sphere.  Because the surface
    sampled is the locus of *probe centres*, a large probe is automatically
    excluded from any crevice narrower than itself - which is exactly the
    physics wanted here.

    Only ``subset`` atoms are evaluated, so the cost scales with the size of
    the target strand rather than the whole complex.  ``freesasa`` would
    compute the full complex and does not expose a per-atom subset, so this
    direct implementation is used instead; it also removes any doubt about
    whether the library clamps large probe radii.
    """
    coords = np.ascontiguousarray(coords, dtype=np.float64)
    radii = np.ascontiguousarray(radii, dtype=np.float64)
    expanded = radii + probe_radius_A
    unit = _fibonacci_sphere(n_points)
    out = np.zeros(len(subset), dtype=np.float64)

    for k, i in enumerate(subset):
        Ri = expanded[i]
        delta = coords - coords[i]
        dist = np.sqrt(np.einsum("ij,ij->i", delta, delta))
        # Only spheres that actually reach atom i's expanded sphere matter.
        cand = np.flatnonzero((dist < Ri + expanded) & (dist > 1e-9))
        pts = coords[i] + unit * Ri
        alive = np.ones(n_points, dtype=bool)
        if cand.size:
            # Nearest occluders first: they kill most test points immediately,
            # so the surviving point set shrinks fast.
            cand = cand[np.argsort(dist[cand], kind="stable")]
            for start in range(0, cand.size, 64):
                idx = np.flatnonzero(alive)
                if idx.size == 0:
                    break
                blk = cand[start : start + 64]
                d2 = np.einsum(
                    "pbc,pbc->pb",
                    pts[idx][:, None, :] - coords[blk][None, :, :],
                    pts[idx][:, None, :] - coords[blk][None, :, :],
                )
                alive[idx[(d2 < expanded[blk][None, :] ** 2).any(axis=1)]] = False
        out[k] = 4.0 * np.pi * Ri * Ri * alive.sum() / n_points
    return out


# --- Atom bookkeeping -------------------------------------------------------


@dataclass(frozen=True)
class AtomSet:
    """Heavy atoms of one complex, tagged by role and residue number."""

    coords: np.ndarray  # (N, 3) float64, Angstroms
    radii: np.ndarray  # (N,) float64, van der Waals radii in Angstroms
    role: np.ndarray  # (N,) str, one of 'protein' / 'guide' / 'target'
    seqid: np.ndarray  # (N,) int, author residue number

    def where(self, role: str) -> np.ndarray:
        return np.flatnonzero(self.role == role)


def build_atom_set(record: ComplexRecord) -> AtomSet:
    """Collect heavy atoms and label each as protein, guide or target.

    Works for the fused ``crRNA-linker-target`` architecture too, because the
    guide and target are addressed by residue number within a chain rather
    than by chain ID.
    """
    if record.guide is None or record.target is None:
        raise ValueError(
            f"{record.pdb_id}: guide and/or target unassigned; cannot build atom set. "
            f"flags={record.flags}"
        )
    st = load_structure(record.path)
    guide_key = {(record.guide.chain_id, s) for s in record.guide.seqids}
    target_key = {(record.target.chain_id, s) for s in record.target.seqids}

    coords: list[list[float]] = []
    radii: list[float] = []
    roles: list[str] = []
    seqids: list[int] = []
    for ch in st[0]:
        residues = list(ch.get_polymer())
        if not residues:
            continue
        for r in residues:
            key = (ch.name, r.seqid.num)
            if key in target_key:
                role = "target"
            elif key in guide_key:
                role = "guide"
            elif ch.name in record.protein_chain_ids:
                role = "protein"
            else:
                continue
            for a in r:
                coords.append([a.pos.x, a.pos.y, a.pos.z])
                radii.append(VDW_RADII_A.get(a.element.name.upper(), DEFAULT_VDW_RADIUS_A))
                roles.append(role)
                seqids.append(r.seqid.num)
    return AtomSet(
        coords=np.asarray(coords, dtype=np.float64),
        radii=np.asarray(radii, dtype=np.float64),
        role=np.asarray(roles, dtype="<U8"),
        seqid=np.asarray(seqids, dtype=int),
    )


# --- Shadow profiles --------------------------------------------------------


@dataclass(frozen=True)
class ShadowProfile:
    """Per-nucleotide steric shadow over one target strand."""

    pdb_id: str
    state: str
    probe_radius_A: float
    seqid: np.ndarray  # (M,) author residue numbers, 5'->3'
    register: np.ndarray  # (M,) duplex-register coordinate
    hybridized: np.ndarray  # (M,) bool, inside the guide:target footprint
    sasa_free_A2: np.ndarray  # (M,) target alone
    sasa_bound_A2: np.ndarray  # (M,) target within the occluding assembly
    shadow: np.ndarray  # (M,) s(i) in [0, 1]
    n_points: int
    notes: tuple[str, ...] = ()

    def __len__(self) -> int:
        return len(self.seqid)

    def as_dict(self) -> dict[int, float]:
        """``s(i)`` keyed by target residue number."""
        return {int(s): float(v) for s, v in zip(self.seqid, self.shadow)}

    def by_register(self) -> dict[int, float]:
        """``s(i)`` keyed by duplex register, the cross-structure coordinate."""
        return {int(r): float(v) for r, v in zip(self.register, self.shadow)}


# SASA_free below this is treated as unmeasurable rather than divided into.
# 1 A^2 is far under one atom's worth of surface at any probe radius.
_MIN_FREE_SASA_A2: float = 1.0


def _residue_sasa(
    atoms: AtomSet,
    occluder_mask: np.ndarray,
    probe_radius_A: float,
    n_points: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-target-residue SASA against the atoms selected by ``occluder_mask``."""
    tgt = atoms.where("target")
    keep = np.flatnonzero(occluder_mask)
    # Re-index the target atoms into the reduced coordinate array.
    pos_in_keep = {int(a): k for k, a in enumerate(keep)}
    subset = np.array([pos_in_keep[int(a)] for a in tgt], dtype=int)
    per_atom = sasa_for_subset(
        atoms.coords[keep], atoms.radii[keep], subset, probe_radius_A, n_points
    )
    seqids = atoms.seqid[tgt]
    order = sorted(set(int(s) for s in seqids))
    totals = np.array([per_atom[seqids == s].sum() for s in order], dtype=np.float64)
    return np.asarray(order, dtype=int), totals


def shadow_profile(
    record: ComplexRecord,
    probe_radius_A: float,
    state: State = "rnp",
    n_points: int = DEFAULT_SPHERE_POINTS,
    atoms: AtomSet | None = None,
) -> ShadowProfile:
    """Compute ``s(i)`` over the target strand for one probe radius.

    ``state`` selects the occluding body:

    * ``'rnp'``    - protein + guide + target: the full ribonucleoprotein.
    * ``'duplex'`` - guide + target only, protein deleted.  This is the
      antisense-oligonucleotide comparison case, and it is *extracted from the
      same cryo-EM coordinates* rather than modelled as an idealised A-form
      helix, so that the only difference between ``'rnp'`` and ``'duplex'``
      is the presence of the protein.  Same target atoms, same conformation,
      same denominator - a genuine like-for-like contrast.
    * ``'free'``   - target alone; ``s(i)`` is identically 0 by construction
      and it exists only as a self-consistency check.
    """
    atoms = atoms if atoms is not None else build_atom_set(record)
    notes: list[str] = []

    is_target = atoms.role == "target"
    if state == "rnp":
        mask = np.ones(len(atoms.role), dtype=bool)
    elif state == "duplex":
        mask = is_target | (atoms.role == "guide")
    elif state == "free":
        mask = is_target.copy()
    else:
        raise ValueError(f"unknown state {state!r}")

    seqids, bound = _residue_sasa(atoms, mask, probe_radius_A, n_points)
    _, free = _residue_sasa(atoms, is_target, probe_radius_A, n_points)

    with np.errstate(invalid="ignore", divide="ignore"):
        s = 1.0 - bound / free
    unmeasurable = free < _MIN_FREE_SASA_A2
    if unmeasurable.any():
        s[unmeasurable] = np.nan
        notes.append(
            f"{int(unmeasurable.sum())} nucleotide(s) had free SASA < "
            f"{_MIN_FREE_SASA_A2} A^2 at probe {probe_radius_A} A; s(i) set to NaN "
            "rather than divided into noise"
        )
    s = np.clip(s, 0.0, 1.0)

    reg = duplex_register(record)
    span = set(record.hybridized_target_seqids)
    return ShadowProfile(
        pdb_id=record.pdb_id,
        state=state,
        probe_radius_A=float(probe_radius_A),
        seqid=seqids,
        register=np.array([reg.get(int(x), 0) for x in seqids], dtype=int),
        hybridized=np.array([int(x) in span for x in seqids], dtype=bool),
        sasa_free_A2=free,
        sasa_bound_A2=bound,
        shadow=s,
        n_points=n_points,
        notes=tuple(notes) + record.flags,
    )


def shadow_series(
    record: ComplexRecord,
    probe_radii_A: Sequence[float] = DEFAULT_PROBE_RADII_A,
    state: State = "rnp",
    n_points: int = DEFAULT_SPHERE_POINTS,
) -> dict[float, ShadowProfile]:
    """``s(i)`` at several probe radii, sharing one atom set."""
    atoms = build_atom_set(record)
    return {
        float(r): shadow_profile(record, r, state=state, n_points=n_points, atoms=atoms)
        for r in probe_radii_A
    }


# --- ASO comparison ---------------------------------------------------------


def aso_atom_set(record: ComplexRecord, length_nt: int, atoms: AtomSet | None = None) -> AtomSet:
    """Atom set for an ``length_nt``-mer oligonucleotide on the same target.

    Built by deleting the protein and truncating the guide down to the
    ``length_nt`` duplex-forming nucleotides centred on the guide:target
    ladder.  Those coordinates are real, experimentally determined A-form
    guide:target duplex - not a modelled helix - so an ASO of length L and the
    RNP are compared on identical target atoms in an identical conformation.
    """
    atoms = atoms if atoms is not None else build_atom_set(record)
    if record.guide is None:
        raise ValueError(f"{record.pdb_id}: no guide segment")
    ladder_guide = [g for _, g in record.duplex_pairs]
    if not ladder_guide:
        raise ValueError(f"{record.pdb_id}: no guide:target ladder to truncate")
    order = {s: i for i, s in enumerate(record.guide.seqids)}
    lo, hi = min(order[g] for g in ladder_guide), max(order[g] for g in ladder_guide)
    duplex_guide = [s for s in record.guide.seqids if lo <= order[s] <= hi]
    if length_nt > len(duplex_guide):
        raise ValueError(
            f"{record.pdb_id}: requested ASO length {length_nt} exceeds the "
            f"{len(duplex_guide)} nt of resolved guide:target duplex available"
        )
    start = (len(duplex_guide) - length_nt) // 2
    keep_guide = set(duplex_guide[start : start + length_nt])

    keep = (atoms.role == "target") | (
        (atoms.role == "guide") & np.isin(atoms.seqid, list(keep_guide))
    )
    idx = np.flatnonzero(keep)
    return AtomSet(atoms.coords[idx], atoms.radii[idx], atoms.role[idx], atoms.seqid[idx])


def aso_shadow(
    record: ComplexRecord,
    length_nt: int,
    probe_radius_A: float,
    n_points: int = DEFAULT_SPHERE_POINTS,
    atoms: AtomSet | None = None,
) -> ShadowProfile:
    """``s(i)`` for a bare ``length_nt``-mer antisense oligonucleotide."""
    sub = aso_atom_set(record, length_nt, atoms=atoms)
    prof = shadow_profile(record, probe_radius_A, state="duplex", n_points=n_points, atoms=sub)
    return ShadowProfile(**{**prof.__dict__, "state": f"aso{length_nt}"})


@lru_cache(maxsize=8)
def cached_record(path: str) -> ComplexRecord:
    """Classification cache; structure parsing dominates short scripts."""
    return classify_structure(path)
