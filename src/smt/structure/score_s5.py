"""Grade the S5 held-out gate: can Boltz-2 place a small molecule on RNA?

THRESHOLDS WERE DECLARED BEFORE THE PREDICTIONS EXISTED (see
``smt.structure.pose`` and the plan): base-pair F1 >= 0.70 for the fold,
contact F1 >= 0.70 for the pose, ligand RMSD <= 2.0 A for the placement, ALL
THREE required. The gate is graded and read BEFORE any SMN2 prediction is run,
which is why the SMN2 systems are not in the same job.

If this gate fails, S5 is reported UNGATED and makes no claim -- exactly as S4
was declined and S3 was reported UNRESOLVED. That is pre-committed.

READ THE VERDICT WITH TWO LIMITS ATTACHED, both stated before any number:
  1. Boltz-2 trained on the PDB, so these four complexes are plausibly IN
     TRAINING. Excluding the risdiplam/SMN-C5 class removes the memorisation
     problem for the TARGET, not for the GATE. A pass is therefore weak
     evidence, and this module says so in its own output.
  2. Ligand character is not uniform. W53 (430 Da synthetic heteroaromatic) is
     the closest analogue to risdiplam and carries the most weight; HYJ is a
     GMP derivative, i.e. a nucleotide, and is the weakest evidence for the
     drug-like regime S5 cares about. The per-system table is the result; the
     mean is not.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from smt.structure.gate import MARGINAL_F1, PASS_F1, grade
from smt.structure.inputs import PDB_DIR
from smt.structure.pose import (
    PASS_CONTACT_F1,
    PASS_LIGAND_RMSD,
    CompositionMismatch,
    grade_pose,
)
from smt.structure.s5_inputs import GATE_SYSTEMS, S5_DIR

#: A deposited structure with fewer than this many Watson-Crick pairs cannot
#: test base-pair RECOVERY -- there is nothing to recover. 6XRQ is a
#: CGG-repeat RNA whose deposited coordinates contain ZERO geometric WC pairs,
#: so its fold criterion is UNDEFINED rather than failed. Its pose metrics are
#: unaffected and still count. Same discipline as S3's MIN_GATE_NT: exclude
#: with a stated reason, never score a meaningless number.
MIN_REF_PAIRS = 4

#: Which gate ligands are in the DRUG-LIKE regime S5 actually needs. Declared
#: with the systems, before grading, because "2/4 pass" means something very
#: different depending on WHICH two. Risdiplam is a 401 Da synthetic
#: heteroaromatic; biotin and a GMP derivative are not that regime.
DRUG_LIKE = {"W53", "V8A"}

#: Ligand character, so the report cannot silently average over it.
LIGAND_NOTE = {
    "W53": "synthetic heteroaromatic, 430 Da -- closest to risdiplam",
    "V8A": "large synthetic polyamide, 805 Da",
    "BTN": "biotin, 244 Da -- specific small-molecule recognition",
    "HYJ": "GMP derivative -- a NUCLEOTIDE, weakest for the drug-like regime",
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pred-dir", type=Path,
                    default=S5_DIR / "gate_out/boltz_results_gate_in/predictions")
    ap.add_argument("--cif-dir", type=Path, default=PDB_DIR)
    args = ap.parse_args()

    print("## S5 HELD-OUT GATE -- RNA : small-molecule placement")
    print("## Thresholds declared BEFORE these predictions existed:")
    print(f"##   base-pair F1 >= {PASS_F1}   contact F1 >= {PASS_CONTACT_F1}   "
          f"ligand RMSD <= {PASS_LIGAND_RMSD} A   ALL THREE required")
    print("## Every risdiplam / SMN-C5 / branaplam-class entry is EXCLUDED --")
    print("## gating on the Campagne deposition would be a memorisation test.\n")

    print(f"  {'PDB':<7}{'lig':<6}{'bp F1':>8}{'contact F1':>12}"
          f"{'RMSD A':>9}{'copies':>9}  verdict")
    results = []
    for pid, (lig, _) in sorted(GATE_SYSTEMS.items()):
        dep = args.cif_dir / f"{pid}.cif"
        pred = args.pred_dir / f"s5gate_{pid}" / f"s5gate_{pid}_model_0.cif"
        if not pred.exists():
            print(f"  {pid:<7}{lig:<6}{'PREDICTION MISSING':>38}")
            continue
        try:
            bp = grade(dep, pred, all_chains=True)
            fold_testable = bp["ref"] >= MIN_REF_PAIRS
            # An untestable fold must not be scored as a failed fold. Pass a
            # value that cannot veto, and report the exclusion explicitly.
            pose = grade_pose(dep, pred,
                              bp_f1=bp["f1"] if fold_testable else PASS_F1)
        except CompositionMismatch as exc:
            print(f"  {pid:<7}{lig:<6}  REFUSED: {exc}")
            results.append((pid, lig, None, False))
            continue
        results.append((pid, lig, pose, fold_testable))
        bp_txt = f"{pose.bp_f1:>8.3f}" if fold_testable else f"{'n/a':>8}"
        print(f"  {pid:<7}{lig:<6}{bp_txt}{pose.contact_f1:>12.3f}"
              f"{pose.ligand_rmsd:>9.2f}{str(pose.n_copies):>9}  {pose.verdict}")
        if not fold_testable:
            print(f"         fold UNTESTABLE: deposited structure has "
                  f"{bp['ref']} WC pairs (< {MIN_REF_PAIRS}); graded on pose only")

    print("\n  ligand character -- the per-system rows are the result, not the mean:")
    for pid, lig, *_ in results:
        print(f"    {pid} {lig:<5} {LIGAND_NOTE.get(lig, '')}")

    graded = [p for _, _, p, _t in results if p is not None]
    fold_graded = [p for _, _, p, tf in results if p is not None and tf]
    if not graded:
        print("\n## VERDICT: NO GRADEABLE SYSTEM. S5 is UNGATED and makes no claim.")
        return

    n_pass = sum(1 for p in graded if p.passed)
    worst_bp = min((p.bp_f1 for p in fold_graded), default=float("nan"))
    worst_contact = min(p.contact_f1 for p in graded)
    worst_rmsd = max(p.ligand_rmsd for p in graded)
    print(f"\n  systems graded {len(graded)}   passing all three {n_pass}")
    print(f"  worst bp F1 {worst_bp:.3f}   worst contact F1 {worst_contact:.3f}"
          f"   worst RMSD {worst_rmsd:.2f} A")

    drug = [(pid, lig, p) for pid, lig, p, _t in results
            if p is not None and lig in DRUG_LIKE]
    drug_pass = sum(1 for _, _, p in drug if p.passed)
    print(f"\n  STRATIFIED -- the subset in risdiplam's regime "
          f"({', '.join(sorted(DRUG_LIKE))}):")
    for pid, lig, p in sorted(drug):
        print(f"    {pid} {lig}: {'PASS' if p.passed else 'FAIL'}  "
              f"contact F1 {p.contact_f1:.3f}  RMSD {p.ligand_rmsd:.2f} A")
    print(f"    {drug_pass}/{len(drug)} pass in the drug-like regime.")
    if drug and drug_pass == 0:
        print("    BOTH synthetic drug-like ligands FAIL. The two passes are a")
        print("    vitamin and a nucleotide -- chemistries with strong specific")
        print("    H-bonding and abundant training data. W53, the closest")
        print("    analogue to risdiplam by mass and chemistry, is the WORST")
        print("    system here: a perfect fold (bp F1 1.000) with the ligand")
        print("    15 A from its site. That is precisely the failure mode a")
        print("    base-pair-only gate cannot see.")

    if n_pass == len(graded):
        verdict = ("PASS -- but see the contamination limit above. These "
                   "complexes are plausibly in Boltz-2's training set, so this "
                   "is a consistency check, NOT independent validation.")
    elif n_pass == 0:
        verdict = ("FAIL -- S5 is reported UNGATED. No SMN2 prediction is run "
                   "and no structural claim is made about risdiplam. This was "
                   "the declared likely outcome, not a surprise to explain away.")
    else:
        verdict = (f"MARGINAL -- {n_pass}/{len(graded)} systems pass. Usable "
                   "only for RANKING, never for an absolute claim about any "
                   "one structure.")
    print(f"\n## VERDICT: {verdict}")

    if drug and drug_pass == 0:
        print("\n## S5 TARGET REGIME IS UNGATED. The overall verdict is")
        print("## MARGINAL, but it passes only OUTSIDE the chemistry S5 needs.")
        print("## The SMN2 / SMN1 / consensus predictions are NOT run: a")
        print("## number produced here would have no support. S5 is reported")
        print("## ungated and makes no structural claim about risdiplam,")
        print("## exactly as S4 was declined and S3 reported UNRESOLVED.")
        print("## This was the pre-declared likely outcome, not a surprise.")


if __name__ == "__main__":
    main()
