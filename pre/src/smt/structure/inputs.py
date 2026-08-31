"""Generate Boltz-2 inputs for S (orthogonal 3D validation of A_5ss).

WHY THIS EXISTS. `A_5ss` comes from ViennaRNA nearest-neighbor thermodynamics
with a `max_bp_span` parameter, and 87% of A7 footprints flip the SIGN of their
predicted effect with that parameter. Boltz-2 is a learned potential with no
concept of max_bp_span at all, so it is an independent instrument on the same
question -- available months before the wet lab.

THE GATE IS NOT TSL2/ISTL1. The plan originally proposed gating on those, but
`smt.seq.construct` carries an explicit TODO(SHAPE) recording that this
project's TSL2 and ISTL1 coordinates are DERIVED FROM FOLDING, not from
literature. Gating one model on another model's output is circular and would
have quietly invalidated the entire exercise.

THE CAS13 RNA CHAINS ARE ALSO A BAD GATE, and that was checked rather than
assumed. Those guide RNAs are PROTEIN-BOUND: their deposited conformation is
induced by Cas13, so predicting them standalone tests whether Boltz can
reproduce a protein-imposed fold without the protein -- which it should not be
able to do, and failing that says nothing about S1. Two of the four chains are
also only 14 and 25 nt, too short to carry a base-pairing test at all.

The gate therefore uses DEPOSITED RNAs THAT FOLD ON THEIR OWN: tRNA-Phe
(1EHZ, 1.93 A X-ray), the UUCG tetraloop hairpin (2KOC, NMR), and the adenine
riboswitch aptamer (1Y26). Each has an experimentally determined secondary
structure that the molecule adopts unaided, which is exactly the regime S1
runs in.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import gemmi

from smt.seq.reference import Interval
from smt.validate.cryptic_sweep import TARGETS

PDB_DIR = Path("/home/dashii/orcd/pool/smt-isef-data/pdb")

#: Deposited RNA-only structures that fold WITHOUT a protein partner.
#: tRNA-Phe (cloverleaf + tertiary), UUCG tetraloop (the canonical hairpin
#: benchmark), adenine riboswitch aptamer (three-way junction + loop-loop kiss).
GATE_PDB_IDS = ("1EHZ", "2KOC", "1Y26")

#: Modified nucleotides gemmi's table leaves uncoded. Both are tRNA-Phe
#: hypermodifications; each is written back to the parent base it derives from.
_EXOTIC = {"YYG": "G", "5MC": "C", "YG": "G"}
RCSB_URL = "https://files.rcsb.org/download/{}.cif"

#: Seeds per prediction. Seed variance is what makes a within-pair comparison
#: statistical rather than anecdotal -- it is where the GPU budget goes.
N_SEEDS = 5

#: Window (nt) of construct sequence handed to Boltz around the region of
#: interest. Boltz cannot fold the full 784-827 nt construct; this keeps the
#: splice site and the footprint together in one chain with real context.
WINDOW_PAD = 60


def rna_chains(cif: Path) -> dict[str, str]:
    """One-letter RNA sequences per polymer chain in a deposited structure."""
    st = gemmi.read_structure(str(cif))
    st.setup_entities()
    out: dict[str, str] = {}
    for chain in st[0]:
        # get_polymer() drops waters, ions and ligands. Iterating the raw chain
        # instead hits HOH / MG / MN and aborts the sequence mid-way.
        poly = chain.get_polymer()
        seq = ""
        ok = True
        for res in poly:
            # gemmi's chemical-component table resolves MODIFIED nucleotides
            # to their parent base. tRNA-Phe (1EHZ) is ~10% modified residues
            # (PSU, H2U, 2MG, 7MG, 5MC, ...), so a naive A/U/G/C membership
            # test rejects the single best RNA gate target in the PDB.
            name = res.name.strip()
            info = gemmi.find_tabulated_residue(name)
            code = (info.one_letter_code.upper().strip() if info else "")
            if not code:
                code = _EXOTIC.get(name, "")
            if code in ("A", "U", "G", "C", "T"):
                seq += "U" if code == "T" else code
            else:
                ok = False
                break
        if ok and len(seq) >= 12:
            out[chain.name] = seq
    return out


def write_yaml(path: Path, chains: list[tuple[str, str, str]]) -> None:
    """chains: list of (kind, chain_id, sequence). No `msa: empty` -- A4 showed
    single-sequence mode fails the gate outright."""
    lines = ["version: 1", "sequences:"]
    for kind, cid, seq in chains:
        lines += [f"  - {kind}:", f"      id: {cid}", f"      sequence: {seq}"]
    path.write_text("\n".join(lines) + "\n")


def gate_inputs(out_dir: Path, cif_dir: Path) -> int:
    """Download (once) and emit the RNA-only gate targets."""
    import urllib.request
    cif_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    for pid in GATE_PDB_IDS:
        cif = cif_dir / f"{pid}.cif"
        if not cif.exists():
            urllib.request.urlretrieve(RCSB_URL.format(pid), cif)
        chains = rna_chains(cif)
        if not chains:
            raise ValueError(f"{pid}: no RNA chain >= 12 nt found")
        cid, seq = sorted(chains.items())[0]
        write_yaml(out_dir / f"gate_{pid}_{cid}.yaml", [("rna", "A", seq)])
        print(f"  gate {pid} chain {cid}: {len(seq)} nt")
        n += 1
    return n


def panel_inputs(panel_csv: Path, out_dir: Path) -> int:
    """S1: every prereg oligo, as ASO + its RNA window; plus unbound controls."""
    rows = list(csv.DictReader(panel_csv.open()))
    constructs = {k: t.construct_fn() for k, t in TARGETS.items()}

    # ONE window per target, spanning every footprint in the panel for it.
    #
    # A per-oligo window would give the two members of a matched pair DIFFERENT
    # RNA chains, so any measured difference could be the window rather than the
    # ASO -- which would silently destroy the whole within-pair design. Sharing
    # one chain per target makes the only difference inside a pair the oligo.
    windows: dict[str, tuple[Interval, str]] = {}
    for tk in {r["target"] for r in rows}:
        c, t = constructs[tk], TARGETS[tk]
        site = c.map_interval(t.site)
        fps = [c.map_interval(Interval(int(r["start"]), int(r["end"])))
               for r in rows if r["target"] == tk]
        lo = max(1, min([site.start] + [f.start for f in fps]) - WINDOW_PAD)
        hi = min(len(c.sequence), max([site.end] + [f.end for f in fps]) + WINDOW_PAD)
        iv = Interval(lo, hi)
        windows[tk] = (iv, c.rna_sub(iv))
        print(f"  {tk}: window {lo}-{hi} ({len(iv)} nt), site at {site.start}-{site.end}")

    n = 0
    for tk, (iv, window) in windows.items():
        write_yaml(out_dir / f"unbound_{tk}.yaml", [("rna", "A", window)])
        n += 1
    for r in rows:
        _, window = windows[r["target"]]
        name = (f"s1_{r['target']}_p{r['pair']}_{r['role']}"
                f"_{r['start']}_{r['length']}")
        write_yaml(out_dir / f"{name}.yaml",
                   [("rna", "A", window), ("rna", "B", r["aso_5to3"])])
        n += 1
    return n


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--panel", type=Path, default=Path("docs/prereg-panel.csv"))
    ap.add_argument("--stage", choices=["gate", "s1"], required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    n = (gate_inputs(args.out, args.out.parent / "gate_cif") if args.stage == "gate"
         else panel_inputs(args.panel, args.out))
    print(f"{args.stage}: wrote {n} Boltz inputs to {args.out}")


if __name__ == "__main__":
    main()
