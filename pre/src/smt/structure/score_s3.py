"""Grade the S3 gate, then measure the hnRNP A1 interface on SMN2.

THE GATE IS GRADED FIRST AND ITS VERDICT BINDS. S1's RNA-only gate (F1 0.971)
says nothing about protein-RNA interface placement, and A4's dCas13 gate passed
only at fold level (4.6-5.7 A) and only with real MSAs. S3 depends on knowing
WHICH NUCLEOTIDES a protein covers, which is a stricter thing than either.

CRITERIA, DECLARED BEFORE THE NUMBERS
--------------------------------------
For each deposited hnRNP A1 : RNA complex, Boltz predicts the same two chains
and we compare the set of RNA nucleotides contacted by protein (any heavy atom
within CONTACT_A) against the deposited answer:

    PASS      interface F1 >= 0.60 on contacted nucleotides
    MARGINAL  0.40 <= F1 < 0.60   -- fold-level only, no per-nucleotide claim
    FAIL      F1 < 0.40           -- S3 reported as UNGATED, no structural claim

The nucleotide-level set is the right unit because that is exactly what S3 then
asks about SMN2: does the UP1 footprint overlap the nusinersen footprint. A
model that places the protein on the right molecule but the wrong nucleotides
cannot answer that, and interface F1 is what detects it.

Global RMSD is reported alongside but is NOT the criterion -- a complex can be
globally close and still smear the interface, and vice versa.

WHAT S3 DOES NOT MEASURE. Boltz has no binding equilibrium. The ternary
prediction is a CO-OCCUPANCY PROBE; nothing here is evidence about which
ligand wins a competition. The claim S3 can support is geometric: the UP1
interface and the ASO footprint occupy overlapping nucleotides, so they cannot
both be bound -- the A4 steric-parity argument against the real repressor.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import gemmi

from smt.structure.s3_inputs import (
    CONTACT_A, GATE_IDS, MIN_GATE_NT, S_DIR, chains_of, smn2_window, up1_sequence,
)

PASS_F1 = 0.60
MARGINAL_F1 = 0.40


def _chain_polymers(st: gemmi.Structure) -> tuple[dict, dict]:
    """({protein chain -> residues}, {rna chain -> residues}), by composition.

    KEPT PER CHAIN, not pooled. 6DCL's asymmetric unit holds TWO copies of the
    complex (protein A+B, RNA C+D). Pooling them makes the deposited "truth"
    span 22 nucleotide slots while Boltz predicts a single 1:1 complex with 11 --
    which caps recall at 0.5 no matter how good the prediction is, and reads as
    a failing gate. The comparison must be one complex against one complex.
    """
    st.setup_entities()
    prot: dict[str, list] = {}
    rna: dict[str, list] = {}
    for chain in st[0]:
        poly = chain.get_polymer()
        if len(poly) == 0:
            continue
        kind = poly.check_polymer_type()
        res = [r for r in poly]
        if kind in (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD):
            prot[chain.name] = res
        elif kind in (gemmi.PolymerType.Rna, gemmi.PolymerType.DnaRnaHybrid):
            rna[chain.name] = res
    return prot, rna


def _contacts(prot_res: list, rna_res: list, cutoff: float) -> set[int]:
    patoms = [a.pos for r in prot_res for a in r]
    out = set()
    for i, r in enumerate(rna_res):
        if any(a.pos.dist(p) <= cutoff for a in r for p in patoms):
            out.add(i)
    return out


def contacted_nucleotides(path: Path, cutoff: float = CONTACT_A) -> set[int]:
    """0-based indices along ONE RNA chain that ONE protein chain touches.

    Indexed by POSITION IN THE POLYMER, not author numbering, so a predicted
    structure (numbered 1..N) and a deposited one (author numbering, gaps) are
    compared on the same footing. When a file holds several copies of the
    complex, the protein/RNA chain PAIR with the most contacts is used -- one
    complex, matching what Boltz was asked to predict.
    """
    st = gemmi.read_structure(str(path))
    st.remove_hydrogens()
    prot, rna = _chain_polymers(st)
    if not prot or not rna:
        return set()
    best: set[int] = set()
    for pres in prot.values():
        for rres in rna.values():
            c = _contacts(pres, rres, cutoff)
            if len(c) > len(best):
                best = c
    return best


def f1(pred: set[int], truth: set[int]) -> tuple[float, float, float]:
    if not pred and not truth:
        return 1.0, 1.0, 1.0
    tp = len(pred & truth)
    prec = tp / len(pred) if pred else 0.0
    rec = tp / len(truth) if truth else 0.0
    return (2 * prec * rec / (prec + rec) if prec + rec else 0.0), prec, rec


def grade_gate(gate_out: Path, cif_dir: Path) -> tuple[bool, list[dict]]:
    rows = []
    for pid in GATE_IDS:
        cif = cif_dir / f"{pid}.cif"
        if not cif.exists():
            continue
        _, rna = chains_of(cif)
        if not rna or len(sorted(rna.values(), key=len)[-1]) < MIN_GATE_NT:
            continue
        truth = contacted_nucleotides(cif)
        pdirs = list(gate_out.rglob(f"predictions/s3gate_{pid}"))
        if not pdirs:
            continue
        models = sorted(pdirs[0].glob("*_model_*.cif"))
        per = [f1(contacted_nucleotides(m), truth) for m in models]
        best = max(per, key=lambda t: t[0]) if per else (0.0, 0.0, 0.0)
        mean = sum(t[0] for t in per) / len(per) if per else 0.0
        rows.append({"pdb": pid, "n_truth": len(truth), "n_models": len(models),
                     "best_f1": best[0], "prec": best[1], "rec": best[2],
                     "mean_f1": mean})
    ok = bool(rows) and all(r["best_f1"] >= PASS_F1 for r in rows)
    return ok, rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gate-out", type=Path, default=S_DIR / "s3gate_out")
    ap.add_argument("--s3-out", type=Path, default=S_DIR / "s3_out")
    ap.add_argument("--cif", type=Path,
                    default=Path("/home/dashii/orcd/pool/smt-isef-data/pdb"))
    a = ap.parse_args()

    print("## S3 -- hnRNP A1 steric parity\n")
    print("### GATE: deposited hnRNP A1 : RNA complexes, re-predicted")
    ok, rows = grade_gate(a.gate_out, a.cif)
    if not rows:
        print("  no gate predictions found -- S3 is UNGATED.")
        return
    print(f"  {'pdb':>6}{'nt contacted':>14}{'models':>8}{'best F1':>9}"
          f"{'prec':>7}{'rec':>7}{'mean F1':>9}")
    for r in rows:
        print(f"  {r['pdb']:>6}{r['n_truth']:>14}{r['n_models']:>8}"
              f"{r['best_f1']:>9.3f}{r['prec']:>7.3f}{r['rec']:>7.3f}"
              f"{r['mean_f1']:>9.3f}")
    worst = min(r["best_f1"] for r in rows)
    verdict = ("PASS" if worst >= PASS_F1
               else "MARGINAL" if worst >= MARGINAL_F1 else "FAIL")
    print(f"\n  worst best-F1 = {worst:.3f}  ->  GATE {verdict}"
          f"   (pass {PASS_F1}, marginal {MARGINAL_F1})")

    if verdict == "FAIL":
        print("\n  S3 IS UNGATED. Its predictions are not interpreted and no")
        print("  structural claim about hnRNP A1 on SMN2 is made. This outcome")
        print("  was pre-committed before the gate ran.")
        return
    if verdict == "MARGINAL":
        print("\n  FOLD-LEVEL ONLY. Per-nucleotide footprint claims are NOT")
        print("  licensed; anything below is reported as indicative.")

    print("\n### SMN2: does the UP1 interface overlap the nusinersen footprint?")
    rna, win, iss, ss5 = smn2_window()
    up1 = up1_sequence(a.cif)
    print(f"  window {win} ({len(rna)} nt), ISS-N1 {iss}, 5'ss {ss5}, UP1 {len(up1)} aa")

    pdirs = {p.name: p for p in a.s3_out.rglob("predictions/s3_*") if p.is_dir()}
    if not pdirs:
        print("  no S3 predictions found yet.")
        return

    # ASO footprint in WINDOW-LOCAL 0-based indices.
    aso_lo = iss.start - win.start
    aso_hi = aso_lo + (iss.end - iss.start)
    aso_set = set(range(aso_lo, aso_hi + 1))

    for name in ("s3_rna_up1", "s3_rna_up1_aso"):
        p = pdirs.get(name)
        if p is None:
            continue
        models = sorted(p.glob("*_model_*.cif"))
        cov, ovl = [], []
        for m in models:
            c = contacted_nucleotides(m)
            cov.append(len(c))
            ovl.append(len(c & aso_set) / len(aso_set) if aso_set else 0.0)
        if not models:
            continue
        print(f"\n  {name}: {len(models)} models")
        print(f"    nucleotides contacted by protein: "
              f"{min(cov)}-{max(cov)} (median {sorted(cov)[len(cov)//2]})")
        print(f"    fraction of ISS-N1 covered by UP1: "
              f"{min(ovl):.2f}-{max(ovl):.2f} "
              f"(median {sorted(ovl)[len(ovl)//2]:.2f})")

    print("\n### SCOPE LIMIT OF THE GATE ITSELF -- read before believing the above")
    print("  hnRNP A1's high-affinity site is a UAG motif, and BOTH gate RNAs are")
    print("  UAG-rich (6DCL AGUAGAUUAGC, 8X0N UUAGGGUUAGGG -- two UAGs each).")
    uag = [i for i in range(len(rna) - 2) if rna[i:i + 3] == "UAG"]
    iss_seq = rna[iss.start - win.start: iss.end - win.start + 1]
    print(f"  The SMN2 window contains {len(uag)} UAG motifs; ISS-N1 "
          f"({iss_seq}) contains {'a' if 'UAG' in iss_seq else 'NO'} UAG.")
    print("  So the gate established competence precisely where the motif is")
    print("  present, and this target lies OUTSIDE that regime. A gate that only")
    print("  proves competence on UAG substrates does not license a claim about a")
    print("  UAG-free one. S3's central question is therefore UNRESOLVED, not")
    print("  answered in the negative: the absence of predicted UP1 contact on")
    print("  ISS-N1 is as consistent with the model being out of its depth as")
    print("  with hnRNP A1 not binding there. It is NOT evidence against the")
    print("  published hnRNP A1 / ISS-N1 interaction.")

    print("\n  INTERPRETATION LIMIT: the ternary system is a co-occupancy probe.")
    print("  Overlap between the UP1 interface and the ISS-N1 nucleotides that")
    print("  nusinersen occupies is a STERIC incompatibility argument -- the two")
    print("  cannot both be bound. It is not a competition assay and does not")
    print("  say which one wins.")


if __name__ == "__main__":
    main()
