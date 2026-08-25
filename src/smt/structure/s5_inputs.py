"""S5 inputs: the U1 snRNA 5' end and the risdiplam ligand, both pinned.

WHAT S5 IS FOR, and what it is NOT. R1's readout assumes that
``P_unpaired(32052)`` -- a 2D ensemble probability -- tracks whether the -1
adenosine is actually PRESENTED as the stacking bulge risdiplam binds. That
2D -> 3D link is assumed, not tested, and it is the weakest joint in R1. S5
tests it by co-folding the U1 : 5'ss duplex with the ligand.

S5 does NOT produce a synergy number. Epistasis stays a dG_open quantity on
the pre-mRNA ensemble and docking cannot compute it.

RISDIPLAM IS A READOUT, NEVER A FOOTPRINT. ``hc_add_up`` means "this base
cannot pair" -- an ASO or a bound protein. Risdiplam does the OPPOSITE: it
stabilises the U1:5'ss helix by stacking on the bulged -1A. Modelling it as an
occupied window would model a drug that prevents the very duplex the drug
exists to promote. Same for branaplam / SMN-C5. Only occluders are footprints.

PRIORS, STATED BEFORE ANY PREDICTION IS RUN. This project has three 3D arms
whose gates PASSED (A4 4.63-5.71 A; S1 F1 0.971; S3 F1 0.857/0.923) and whose
results were all null, inconclusive, or unresolved. S5 enters with a WORSE
prior than any of them:

  1. Training-set contamination is worse here than in S1/S3. Campagne et al.
     2019 deposited SMN-C5 bound to the U1:SMN2 5'ss duplex, and Boltz-2
     trained on the PDB, so this complex is plausibly in training. Gating on
     that entry and then predicting it would be a memorisation test, not a
     gate -- the S3 lesson inverted. The held-out gate therefore EXCLUDES all
     risdiplam / SMN-C5 / branaplam-class entries.
  2. RNA-small-molecule structure prediction is substantially less reliable
     than protein-ligand, and n = 3 systems. A pass is a consistency check,
     not a test. The honest prior is MARGINAL or FAIL.

Nothing here is typed from memory. The ligand is pinned from PubChem by CID
and validated through rdkit (canonical round-trip, formula, MW, InChIKey); the
U1 sequence is read from a pinned RefSeq GenBank record with its length
asserted, exactly as ``smt.seq.reference`` does for NG_008728.
"""

from __future__ import annotations

from pathlib import Path

from smt.seq.reference import DATA

# --------------------------------------------------------------------------
# Risdiplam, pinned from PubChem by CID. Verified against the record 2026-08-24.
# --------------------------------------------------------------------------

RISDIPLAM_CID = 118513932
RISDIPLAM_SMILES = "CC1=CC(=NN2C1=NC(=C2)C)C3=CC(=O)N4C=C(C=CC4=N3)N5CCNC6(C5)CC6"
RISDIPLAM_FORMULA = "C22H23N7O"
RISDIPLAM_MW = 401.47
RISDIPLAM_INCHIKEY = "ASKZRYGFUPSJPN-UHFFFAOYSA-N"

# --------------------------------------------------------------------------
# U1 snRNA, pinned RefSeq record.
# --------------------------------------------------------------------------

U1_ACCESSION = "NR_004430.4"
U1_GB = DATA / "RNU1-1_NR_004430.gb"
U1_LENGTH = 164

#: U1 nucleotides 3-11. These are the bases that pair with the 9-nt 5' splice
#: site in the canonical U1:pre-mRNA duplex, antiparallel. NOT hardcoded as a
#: string -- ``u1_duplex_arm()`` slices them out of the pinned record, and
#: ``pairs_with`` re-derives the pairing rather than asserting it.
U1_ARM = (3, 11)


class LigandError(RuntimeError):
    """The pinned ligand no longer round-trips or no longer matches the record."""


def load_u1() -> str:
    """Full U1 snRNA as RNA (uppercase, U not T), from the pinned record."""
    if not U1_GB.exists():
        raise FileNotFoundError(
            f"{U1_GB} missing. Fetch with NCBI efetch db=nuccore "
            f"id={U1_ACCESSION.split('.')[0]} rettype=gbwithparts."
        )
    from Bio import SeqIO

    rec = SeqIO.read(U1_GB, "genbank")
    seq = str(rec.seq).upper().replace("T", "U")
    if len(seq) != U1_LENGTH:
        raise ValueError(
            f"{U1_ACCESSION} length {len(seq)} != pinned {U1_LENGTH}; "
            "record changed, re-verify the duplex arm before proceeding."
        )
    return seq


def u1_duplex_arm() -> str:
    """U1 nt 3-11, 5'->3'. The arm that reads the 5' splice site."""
    seq = load_u1()
    lo, hi = U1_ARM
    return seq[lo - 1 : hi]


_COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G"}


def pairs_with(u1_arm: str, five_ss: str, allow_wobble: bool = True) -> list[bool]:
    """Per-position Watson-Crick (optionally GU wobble) pairing of the duplex.

    The U1 arm runs 5'->3' and the splice site runs 5'->3', so they anneal
    ANTIPARALLEL: ``u1_arm`` reversed aligns with ``five_ss``. Returns one bool
    per splice-site position so a caller can see WHICH bases pair, rather than
    a single yes/no that hides a mismatch.
    """
    if len(u1_arm) != len(five_ss):
        raise ValueError(f"length mismatch: {len(u1_arm)} vs {len(five_ss)}")
    out = []
    for u, m in zip(reversed(u1_arm), five_ss):
        wc = _COMPLEMENT.get(u) == m
        wob = allow_wobble and {u, m} == {"G", "U"}
        out.append(wc or wob)
    return out


def validate_ligand() -> dict:
    """Round-trip the SMILES through rdkit and check it against the pinned record.

    Raises ``LigandError`` rather than returning a mismatch, because a ligand
    that is not what it claims to be poisons every downstream pose metric.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, rdMolDescriptors
    except ImportError as exc:  # pragma: no cover - env-dependent
        raise LigandError("rdkit unavailable; S5 cannot validate its ligand") from exc

    mol = Chem.MolFromSmiles(RISDIPLAM_SMILES)
    if mol is None:
        raise LigandError("pinned risdiplam SMILES does not parse")

    canonical = Chem.MolToSmiles(mol)
    if Chem.MolToSmiles(Chem.MolFromSmiles(canonical)) != canonical:
        raise LigandError("canonical SMILES is not round-trip stable")

    formula = rdMolDescriptors.CalcMolFormula(mol)
    mw = Descriptors.MolWt(mol)
    key = Chem.MolToInchiKey(mol)

    if formula != RISDIPLAM_FORMULA:
        raise LigandError(f"formula {formula} != pinned {RISDIPLAM_FORMULA}")
    if abs(mw - RISDIPLAM_MW) > 0.05:
        raise LigandError(f"MW {mw:.2f} != pinned {RISDIPLAM_MW}")
    if key != RISDIPLAM_INCHIKEY:
        raise LigandError(f"InChIKey {key} != pinned {RISDIPLAM_INCHIKEY}")

    return {"canonical_smiles": canonical, "formula": formula,
            "mw": mw, "inchikey": key,
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "rings": rdMolDescriptors.CalcNumRings(mol)}


# --------------------------------------------------------------------------
# Held-out gate: deposited RNA : small-molecule complexes.
# --------------------------------------------------------------------------

S5_DIR = Path("/home/dashii/orcd/pool/smt-isef-data/boltz_s5")

#: Found by QUERYING RCSB (protein-free RNA, non-polymer entity, <2.5 A,
#: >20 residues; 510 hits), not recalled. The top 40 by resolution were
#: screened for SMN / risdiplam / branaplam / "spinal muscular" in the title:
#: NONE were present, so the held-out exclusion is verifiable rather than
#: delicate.
#:
#: WHY EXCLUSION MATTERS. Campagne et al. 2019 deposited SMN-C5 bound to the
#: U1:SMN2 5'ss duplex and Boltz-2 trained on the PDB. Gating on that entry and
#: then predicting it would be a memorisation test, not a gate -- the S3 lesson
#: inverted.
#:
#: ligand character is recorded because it is NOT uniform, and the report must
#: say so rather than averaging over it:
#:   W53  430 Da  synthetic diamino-quinoline -- the closest analogue to
#:                risdiplam (401 Da heteroaromatic) and the most informative
#:   V8A  805 Da  large synthetic polyamide
#:   BTN  244 Da  biotin; a genuine specific small-molecule:RNA recognition
#:   HYJ  364 Da  a GMP DERIVATIVE, i.e. a nucleotide, not a drug-like
#:                molecule. Kept as a placement test, but it is the weakest
#:                evidence here for the drug-like regime S5 cares about.
GATE_SYSTEMS = {
    "8QMH": ("W53", "c1cc(nc2c1ccc(n2)N)CC(Cc3ccc4ccc(nc4n3)N)C(=O)NCCCN"),
    "6XRQ": ("V8A", "c1c2c3c4c(cc(c3C(=O)N(C2=O)CCCN5CCOCC5)NCCCN6CCOCC6)"
                    "C(=O)N(C(=O)c4c1NCCCN7CCOCC7)CCCN8CCOCC8"),
    "1F27": ("BTN", "C1[C@H]2[C@@H]([C@@H](S1)CCCCC(=O)O)NC(=O)N2"),
    "7QSH": ("HYJ", "c1nc2c(n1[C@H]3[C@@H]([C@@H]([C@H](O3)COP(=O)(O)O)O)O)"
                    "N=C(NC2=O)O"),
}

#: Minimum MODELED nucleotides for a gate complex to test anything -- the 9VVA
#: lesson from S3, where deposition metadata described a 26-mer but only three
#: nucleotides were resolved. Checked on coordinates, never on metadata.
#:
#: APPLIED TO THE COMPLEX TOTAL, NOT PER CHAIN, and the unit differs from S3's
#: deliberately. S3 gated a protein:RNA INTERFACE, which lives on one RNA, so
#: a per-chain minimum was the right unit there. A ligand POSE needs enough
#: modelled RNA to define a binding SITE, which several short chains can do
#: jointly: 8QMH is four 6-nt G2C4 strands totalling 24 nt and hosting a
#: measurable ligand site. Rejecting it per-chain would discard the closest
#: analogue to risdiplam for a reason that does not apply to this question.
#: Per-chain composition is still reported so a 4x6 assembly is never
#: presented as a 24-mer.
MIN_GATE_NT = 8


def rna_chains_from_cif(cif: Path) -> dict[str, str]:
    """Chain id -> RNA sequence, read from MODELLED residues only."""
    import gemmi

    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    out: dict[str, str] = {}
    for chain in st[0]:
        poly = chain.get_polymer()
        if len(poly) == 0:
            continue
        seq = ""
        for res in poly:
            info = gemmi.find_tabulated_residue(res.name.strip())
            code = (info.one_letter_code.upper().strip() if info else "")
            seq += {"T": "U"}.get(code, code) or "N"
        out[chain.name] = seq
    return out


def write_ligand_yaml(path: Path, rna: list[tuple[str, str]], smiles: str,
                      ligand_id: str = "L") -> None:
    """Boltz input: one or more RNA chains plus one SMILES ligand.

    ``smt.structure.inputs.write_yaml`` cannot express a ligand, which is why
    this exists rather than being folded into it. No ``msa: empty`` -- A4
    established that single-sequence mode fails the gate outright.
    """
    lines = ["version: 1", "sequences:"]
    for cid, seq in rna:
        lines += [f"  - rna:", f"      id: {cid}", f"      sequence: {seq}"]
    lines += ["  - ligand:", f"      id: {ligand_id}", f"      smiles: '{smiles}'"]
    path.write_text("\n".join(lines) + "\n")


def gate_inputs(out_dir: Path, cif_dir: Path | None = None) -> int:
    """Emit one Boltz YAML per held-out gate complex. Returns the count."""
    from smt.structure.inputs import PDB_DIR

    cif_dir = PDB_DIR if cif_dir is None else cif_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for pid, (lig, smiles) in sorted(GATE_SYSTEMS.items()):
        cif = cif_dir / f"{pid}.cif"
        if not cif.exists():
            raise FileNotFoundError(f"{cif} missing; fetch from RCSB first")
        chains = rna_chains_from_cif(cif)
        rna = [(cid, seq) for cid, seq in sorted(chains.items())]
        total = sum(len(s) for _, s in rna)
        if total < MIN_GATE_NT:
            raise ValueError(
                f"{pid}: only {total} modelled nt across {len(rna)} chain(s); "
                f"below MIN_GATE_NT={MIN_GATE_NT}. Metadata may claim more -- "
                "the 9VVA lesson. Excluded rather than graded.")
        write_ligand_yaml(out_dir / f"s5gate_{pid}.yaml", rna, smiles)
        sizes = "+".join(str(len(s)) for _, s in rna)
        print(f"  gate {pid}: {len(rna)} chain(s) [{sizes}] = {total} nt, "
              f"ligand {lig}")
        n += 1
    return n
