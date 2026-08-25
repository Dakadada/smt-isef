"""Grade the S gate: can Boltz-2 reproduce RNA base-pairing at all?

S1 asks whether a splice site ends up PAIRED in a predicted structure. That
question is meaningless unless the predictor gets base-pairing right on RNAs
whose answer is already known, so this grades exactly that -- not global RMSD,
which can look respectable while every pair is wrong.

Pairs are detected geometrically from coordinates, identically in the deposited
and predicted structures, so the comparison cannot be biased by using a
different definition on each side. A Watson-Crick pair is called when the
purine N1 and pyrimidine N3 are within HBOND_MAX and the C1'-C1' distance is in
the canonical range -- the standard geometric signature.

Verdict thresholds are declared here, BEFORE looking at the numbers.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import gemmi

#: N1(purine)-N3(pyrimidine) hydrogen bond, angstroms.
HBOND_MAX = 3.5
#: C1'-C1' separation of a Watson-Crick pair, angstroms.
C1_MIN, C1_MAX = 8.5, 12.5

#: Declared before grading. F1 over base pairs, versus the deposited structure.
PASS_F1 = 0.70
MARGINAL_F1 = 0.50

_PURINE = {"A", "G"}
_WC = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}


def _code(res) -> str:
    info = gemmi.find_tabulated_residue(res.name.strip())
    c = (info.one_letter_code.upper().strip() if info else "")
    return {"T": "U", "YYG": "G", "5MC": "C"}.get(c or res.name.strip(), c)


def _atom(res, *names):
    for n in names:
        a = res.find_atom(n, "*")
        if a is not None:
            return a.pos
    return None


def base_pairs(path: Path, all_chains: bool = False) -> set[tuple[int, int]]:
    """Geometric Watson-Crick/wobble pairs as 0-based index pairs in the polymer.

    ``all_chains`` is OPT-IN so the S1/S3 single-chain behaviour is byte-for-byte
    unchanged. S5's gate complexes are multi-chain duplexes and assemblies
    (8QMH is four 6-nt strands), where stopping at the first chain would score
    an inter-chain duplex as having no pairs at all.
    """
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    residues = []
    for chain in st[0]:
        for res in chain.get_polymer():
            residues.append(res)
        if not all_chains:
            break                  # single-chain gate targets only
    pairs: set[tuple[int, int]] = set()
    for i in range(len(residues)):
        ci = _code(residues[i])
        for j in range(i + 2, len(residues)):
            cj = _code(residues[j])
            if (ci, cj) not in _WC:
                continue
            ni = _atom(residues[i], "N1" if ci in _PURINE else "N3")
            nj = _atom(residues[j], "N1" if cj in _PURINE else "N3")
            c1i, c1j = _atom(residues[i], "C1'"), _atom(residues[j], "C1'")
            if None in (ni, nj, c1i, c1j):
                continue
            if ni.dist(nj) <= HBOND_MAX and C1_MIN <= c1i.dist(c1j) <= C1_MAX:
                pairs.add((i, j))
    return pairs


def grade(deposited: Path, predicted: Path, all_chains: bool = False) -> dict:
    ref = base_pairs(deposited, all_chains)
    pred = base_pairs(predicted, all_chains)
    tp = len(ref & pred)
    prec = tp / len(pred) if pred else 0.0
    rec = tp / len(ref) if ref else float("nan")
    f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
    return {"ref": len(ref), "pred": len(pred), "tp": tp,
            "precision": prec, "recall": rec, "f1": f1}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cif-dir", type=Path,
                    default=Path("/home/dashii/orcd/pool/smt-isef-data/boltz_s/gate_cif"))
    ap.add_argument("--pred-dir", type=Path,
                    default=Path("/home/dashii/orcd/pool/smt-isef-data/boltz_s/gate_out/"
                                 "boltz_results_gate_in/predictions"))
    args = ap.parse_args()

    print("## S GATE -- base-pair recovery vs deposited structures")
    print(f"thresholds declared before grading: PASS F1 >= {PASS_F1}, "
          f"MARGINAL >= {MARGINAL_F1}\n")
    print(f"  {'target':<14}{'ref bp':>8}{'pred bp':>9}{'TP':>5}"
          f"{'prec':>8}{'recall':>8}{'F1':>8}")
    scores = []
    for pdir in sorted(args.pred_dir.iterdir()):
        model = pdir / f"{pdir.name}_model_0.cif"
        if not model.exists():
            continue
        pid = pdir.name.split("_")[1]
        dep = args.cif_dir / f"{pid}.cif"
        g = grade(dep, model)
        scores.append((pid, g["f1"]))
        print(f"  {pid:<14}{g['ref']:>8}{g['pred']:>9}{g['tp']:>5}"
              f"{g['precision']:>8.3f}{g['recall']:>8.3f}{g['f1']:>8.3f}")

    worst = min(f for _, f in scores)
    mean = sum(f for _, f in scores) / len(scores)
    print(f"\n  mean F1 {mean:.3f}   worst F1 {worst:.3f}")
    if worst >= PASS_F1:
        verdict = "PASS -- Boltz recovers RNA base-pairing; S1 may proceed."
    elif worst >= MARGINAL_F1:
        verdict = ("MARGINAL -- usable only for RANKING pairs against each other, "
                   "never for an absolute claim about any one structure.")
    else:
        verdict = ("FAIL -- S1-S4 are reported UNGATED and no structural claim "
                   "is made. This is a real result: the 2D observable could not "
                   "be corroborated by an orthogonal 3D model.")
    print(f"  VERDICT: {verdict}")


if __name__ == "__main__":
    main()
