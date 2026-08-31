"""S3: hnRNP A1 steric parity, made structural.

THE UNMET SPEC BULLET. The project spec asks for "synthetic chemistries and
biological guide RNAs as identical steric blockades against hnRNP A1/A2". A4
answered only half of it: it measured the GEOMETRIC SHADOW of deposited dCas13
RNPs and concluded no 8-20 nt ASO reaches parity with one. It never modeled the
protein actually being competed with. hnRNP A1 is that protein -- it binds
ISS-N1 and represses SMN2 exon 7, and nusinersen works by displacing it.

WHAT IS AND IS NOT MEASURABLE HERE -- read before interpreting any output.
Boltz-2 has NO notion of binding equilibrium. Handing it protein + RNA + ASO
asks it to place all three; whatever it returns is not evidence about which
one wins, and reading "displacement" off such a prediction would be inventing
a competition assay out of a structure predictor. What IS measurable, and what
this module measures, is GEOMETRIC: which RNA nucleotides the UP1 interface
covers, and whether an ASO footprint overlaps them. That is the same steric
parity argument A4 made, now against the real repressor instead of a dCas13
stand-in. The ternary input is still built and predicted, but it is labeled a
CO-OCCUPANCY PROBE, not a competition experiment.

THE GATE, and why it is a real one. A4's precedent is binding: Boltz must
reproduce known complexes before any novel prediction is trusted, and in
single-sequence mode it FAILED, passing only with real MSAs and only at
fold level (4.6-5.7 A). S1's RNA-only gate was separate and passed at F1 0.971.
Neither covers PROTEIN-RNA interface placement, which is what S3 depends on.

So S3 gates on deposited hnRNP A1 : RNA complexes -- the same protein, the same
kind of interface, solved by X-ray:

    6DCL  2.50 A   UP1 + 12-nt RNA
    9VVA  1.75 A   hnRNP A1 + 26-mer RNA
    8X0N  2.80 A   hnRNP A1 + 12-nt RNA

These were found by querying RCSB for hnRNP A1 entries containing a nucleic
acid, not recalled -- and the ssDNA telomeric-repeat complexes that dominate
that list (2UP1, 1U1K, 1PO6, ...) were deliberately NOT used, because ssDNA is
not the polymer S3 reasons about.

THE PROTEIN SEQUENCE IS TAKEN FROM THE DEPOSITED FILE, never typed from
memory. A wrong UP1 sequence would produce a confident, wrong interface and
nothing downstream would catch it.

IF THE GATE FAILS, S3 IS REPORTED AS UNGATED AND NO STRUCTURAL CLAIM IS MADE.
That is pre-committed here, before any prediction exists.
"""

from __future__ import annotations

import argparse
import urllib.request
from pathlib import Path

import gemmi

from smt.seq.construct import compact_window, ISS_N1, SS5_EXON7
from smt.seq.reference import Interval
from smt.sieve.scan import NUSINERSEN
from smt.structure.inputs import PDB_DIR, RCSB_URL, write_yaml

S_DIR = Path("/home/dashii/orcd/pool/smt-isef-data/boltz_s")

#: Deposited hnRNP A1 : RNA complexes used as the S3 gate. RNA, not ssDNA.
#:
#: 9VVA is listed by RCSB as carrying a 26-mer RNA and was in this tuple until
#: its coordinates were read: only THREE nucleotides are actually modeled, the
#: rest being disordered in the crystal. Deposition metadata describes the
#: construct, not what was resolved. A 3-nt chain cannot carry an interface
#: test, so it is excluded by MIN_GATE_NT below rather than by editing it out
#: of the list -- the guard has to catch the next one too.
GATE_IDS = ("6DCL", "9VVA", "8X0N")

#: Minimum MODELED nucleotides for a gate complex to test anything. S1 rejected
#: 14- and 25-nt Cas13 chains as too short to carry a base-pairing test; 3 nt is
#: far below any interface footprint.
MIN_GATE_NT = 8

#: Contact distance for "this protein residue touches this nucleotide", any
#: heavy atom to any heavy atom. 4.0 A is the conventional non-bonded contact
#: cutoff and is used identically on deposited and predicted structures.
CONTACT_A = 4.0

#: Construct-space padding around the SMN2 region handed to Boltz. The window
#: must contain the exon 7 5' splice site AND ISS-N1, since the whole claim is
#: that occupancy at the latter changes usage of the former.
WINDOW_PAD = 45

INTRON6_TAIL = 200


def _fetch(pid: str, cif_dir: Path) -> Path:
    cif_dir.mkdir(parents=True, exist_ok=True)
    cif = cif_dir / f"{pid}.cif"
    if not cif.exists():
        urllib.request.urlretrieve(RCSB_URL.format(pid), cif)
    return cif


def chains_of(cif: Path) -> tuple[dict[str, str], dict[str, str]]:
    """(protein chains, rna chains) as one-letter sequences.

    Classified by residue composition, NOT by depositor annotation -- the same
    rule A4 adopted after 8FTI turned out to label its protein "Integrase".
    """
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    prot: dict[str, str] = {}
    rna: dict[str, str] = {}
    for chain in st[0]:
        poly = chain.get_polymer()
        if len(poly) == 0:
            continue
        kind = poly.check_polymer_type()
        seq = gemmi.one_letter_code(poly.extract_sequence())
        if not seq:
            continue
        if kind in (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD):
            prot[chain.name] = seq.upper()
        elif kind in (gemmi.PolymerType.Rna, gemmi.PolymerType.DnaRnaHybrid):
            rna[chain.name] = seq.upper().replace("T", "U")
    return prot, rna


def gate_inputs(out_dir: Path, cif_dir: Path = PDB_DIR) -> int:
    out_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for pid in GATE_IDS:
        cif = _fetch(pid, cif_dir)
        prot, rna = chains_of(cif)
        if not prot or not rna:
            raise ValueError(f"{pid}: expected one protein and one RNA chain, "
                             f"got protein={list(prot)} rna={list(rna)}")
        pseq = sorted(prot.values(), key=len)[-1]
        rseq = sorted(rna.values(), key=len)[-1]
        if len(rseq) < MIN_GATE_NT:
            print(f"  gate {pid}: SKIPPED -- only {len(rseq)} nt modeled "
                  f"({rseq}); deposition metadata overstates what was resolved")
            continue
        write_yaml(out_dir / f"s3gate_{pid}.yaml",
                   [("protein", "A", pseq), ("rna", "B", rseq)])
        print(f"  gate {pid}: protein {len(pseq)} aa + RNA {len(rseq)} nt ({rseq})")
        n += 1
    return n


def up1_sequence(cif_dir: Path = PDB_DIR) -> str:
    """UP1 taken from the deposited 6DCL coordinates, not from memory."""
    prot, _ = chains_of(_fetch("6DCL", cif_dir))
    return sorted(prot.values(), key=len)[-1]


def smn2_window() -> tuple[str, Interval, Interval, Interval]:
    """(RNA window, window iv, ISS-N1 iv, 5'ss iv) -- all CONSTRUCT space."""
    c = compact_window(intron6_tail=INTRON6_TAIL)
    iss = c.map_interval(ISS_N1)
    ss5 = c.map_interval(SS5_EXON7)
    lo = max(1, min(iss.start, ss5.start) - WINDOW_PAD)
    hi = min(len(c.sequence), max(iss.end, ss5.end) + WINDOW_PAD)
    win = Interval(lo, hi)
    return c.rna_sub(win), win, iss, ss5


def s3_inputs(out_dir: Path, cif_dir: Path = PDB_DIR) -> int:
    """Four systems. The ternary one is a CO-OCCUPANCY PROBE, not a competition
    assay -- see the module docstring."""
    out_dir.mkdir(parents=True, exist_ok=True)
    rna, win, iss, ss5 = smn2_window()
    up1 = up1_sequence(cif_dir)
    aso = NUSINERSEN.replace("T", "U")

    systems = {
        "s3_rna": [("rna", "A", rna)],
        "s3_rna_up1": [("rna", "A", rna), ("protein", "B", up1)],
        "s3_rna_aso": [("rna", "A", rna), ("rna", "C", aso)],
        "s3_rna_up1_aso": [("rna", "A", rna), ("protein", "B", up1), ("rna", "C", aso)],
    }
    for name, chains in systems.items():
        write_yaml(out_dir / f"{name}.yaml", chains)
    print(f"  SMN2 window {win} = {len(rna)} nt; ISS-N1 {iss}; 5'ss {ss5}")
    print(f"  UP1 {len(up1)} aa from 6DCL; ASO {aso} ({len(aso)} nt)")
    print(f"  {len(systems)} systems: {', '.join(systems)}")
    return len(systems)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gate-out", type=Path, default=S_DIR / "s3gate_in")
    ap.add_argument("--s3-out", type=Path, default=S_DIR / "s3_in")
    ap.add_argument("--cif", type=Path, default=PDB_DIR)
    a = ap.parse_args()
    print("S3 gate inputs (hnRNP A1 : RNA, deposited):")
    gate_inputs(a.gate_out, a.cif)
    print("S3 systems:")
    s3_inputs(a.s3_out, a.cif)
    print("\nGATE MUST PASS BEFORE ANY S3 RESULT IS BELIEVED.")


if __name__ == "__main__":
    main()
