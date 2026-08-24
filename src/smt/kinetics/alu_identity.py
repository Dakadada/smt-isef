"""P: are C2B34's flanking inverted repeats actually ALU elements?

WHY IT MATTERS. A8 found inverted repeats in the introns flanking C2B34
(intron 1 <-> intron 4, z = -52.6, 209 revcomp vs 18 forward matches) and they
are ~280 nt -- Alu length. **Alu elements are primate-specific.** If these arms
are Alus, that supplies a concrete mechanism for why C2B34 is seen only in
primates: Alu-mediated backsplicing. That is the most genuinely novel claim
available to this project.

WHY IT IS HARD HERE. `smt.kinetics.alu` says plainly in its own docstring that
it CANNOT call Alu subfamilies -- there is no RepeatMasker or Repbase in this
environment. This module closes that gap WITHOUT a repeat database, using only
features internal to the sequence:

  * length in the Alu band (~250-320 nt);
  * the Pol III internal promoter A-box and B-box, at their canonical spacing,
    scored as degenerate motif matches rather than exact strings (real Alus
    diverge 5-20% from consensus, so exact matching finds nothing);
  * a 3' A-rich tail;
  * DIMERIC structure -- an Alu is a left arm and a right arm joined by an
    A-rich linker, and the two arms are homologous. This is the most
    Alu-specific feature and the hardest to hit by chance.

EVERY SCORE IS NULL-CONTROLLED. A motif this degenerate will match random
sequence at some rate, so the same scoring is run against composition-matched
shuffles and against random same-length windows drawn from the SAME introns.
The z-score is reported next to A8's -52.6 so the two are comparable.

HONEST CEILING, stated before the numbers: motif + dimer evidence is weaker
than RepeatMasker. If the hallmarks do not clear the null, the conclusion is
"Alu-scale inverted repeats, subfamily unresolved" -- which is still A8's
result, just not upgraded.
"""

from __future__ import annotations

import argparse
import random
import re
import statistics
from dataclasses import dataclass

from smt.kinetics.alu import (
    DOWNSTREAM_FLANK, IR_NULL_SEED, IR_NULL_TRIALS, UPSTREAM_FLANK,
    scan_flanks, top_candidates,
)
from smt.seq.reference import Interval, load_smn2

#: Alu is ~300 nt; free-standing monomers and truncations run shorter.
ALU_LEN_MIN, ALU_LEN_MAX = 240, 340

#: Pol III type-2 internal promoter. IUPAC, matched with mismatch tolerance.
A_BOX = "TGGCTCACGCC"
B_BOX = "GTTCGANNC"
#: Canonical A-box->B-box spacing in Alu, generous bounds.
BOX_GAP_MIN, BOX_GAP_MAX = 30, 90

_IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "N": ".",
          "W": "[AT]", "R": "[AG]", "Y": "[CT]", "S": "[GC]"}


def _best_motif_hit(seq: str, motif: str, max_mm: int) -> tuple[int, int] | None:
    """Best (position, mismatches) for a degenerate motif, or None."""
    pat = [_IUPAC[c] for c in motif]
    best = None
    for i in range(len(seq) - len(motif) + 1):
        mm = sum(0 if re.match(p, seq[i + j]) else 1 for j, p in enumerate(pat))
        if mm <= max_mm and (best is None or mm < best[1]):
            best = (i, mm)
    return best


def dimer_score(seq: str) -> float:
    """Best ungapped identity between the left and right halves of the arm.

    An Alu is a dimer of two homologous arms; random sequence is not. Compares
    the first and second halves over a sliding offset and returns the best
    fractional identity.
    """
    h = len(seq) // 2
    left, right = seq[:h], seq[h:]
    n = min(len(left), len(right))
    best = 0.0
    for off in range(-30, 31):
        a0, b0 = max(0, -off), max(0, off)
        span = min(n - a0, n - b0)
        if span < 60:
            continue
        ident = sum(1 for k in range(span) if left[a0 + k] == right[b0 + k]) / span
        best = max(best, ident)
    return best


def a_tail(seq: str, window: int = 40) -> float:
    """A-fraction of the 3' terminal window."""
    tail = seq[-window:]
    return tail.count("A") / len(tail)


@dataclass(frozen=True)
class Hallmarks:
    length_ok: bool
    a_box_mm: int | None
    b_box_mm: int | None
    spacing_ok: bool
    dimer: float
    a_tail: float

    @property
    def score(self) -> float:
        """Composite, deliberately simple so the null comparison is meaningful."""
        s = 0.0
        s += 1.0 if self.length_ok else 0.0
        s += 1.0 if self.a_box_mm is not None else 0.0
        s += 1.0 if self.b_box_mm is not None else 0.0
        s += 1.0 if self.spacing_ok else 0.0
        s += 2.0 * self.dimer
        s += 1.0 * self.a_tail
        return s


def hallmarks(seq: str) -> Hallmarks:
    seq = seq.upper().replace("U", "T")
    a = _best_motif_hit(seq, A_BOX, max_mm=3)
    b = _best_motif_hit(seq, B_BOX, max_mm=2)
    spacing = False
    if a and b:
        gap = b[0] - (a[0] + len(A_BOX))
        spacing = BOX_GAP_MIN <= gap <= BOX_GAP_MAX
    return Hallmarks(
        length_ok=ALU_LEN_MIN <= len(seq) <= ALU_LEN_MAX,
        a_box_mm=a[1] if a else None,
        b_box_mm=b[1] if b else None,
        spacing_ok=spacing,
        dimer=dimer_score(seq),
        a_tail=a_tail(seq),
    )


def null_scores(seq: str, full: str, trials: int, rng: random.Random) -> tuple[list, list]:
    """Composition-matched shuffles, and random same-length intron windows."""
    shuffled = []
    chars = list(seq)
    for _ in range(trials):
        rng.shuffle(chars)
        shuffled.append(hallmarks("".join(chars)).score)
    windows = []
    for iv in (UPSTREAM_FLANK, DOWNSTREAM_FLANK):
        span = iv.end - iv.start + 1 - len(seq)
        if span <= 0:
            continue
        for _ in range(trials):
            st = rng.randint(iv.start, iv.start + span)
            windows.append(hallmarks(full[st - 1: st - 1 + len(seq)]).score)
    return shuffled, windows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=3)
    ap.add_argument("--trials", type=int, default=IR_NULL_TRIALS)
    args = ap.parse_args()

    full = load_smn2()
    print("## P -- Alu identity of the C2B34 flanking inverted repeats")
    print("   (scanning flanks; ~3 min, per smt.kinetics.alu.scan_flanks)\n")
    hits = top_candidates(scan_flanks(full), args.top)
    rng = random.Random(IR_NULL_SEED)

    for i, h in enumerate(hits):
        for side, iv in (("upstream", h.matched_upstream),
                         ("downstream", h.matched_downstream)):
            seq = full[iv.start - 1: iv.end]
            hm = hallmarks(seq)
            sh, wi = null_scores(seq, full, args.trials, rng)
            pool = sh + wi
            mu, sd = statistics.mean(pool), statistics.pstdev(pool)
            z = (hm.score - mu) / sd if sd else float("nan")
            print(f"  hit {i} {side:<10} {iv.start}-{iv.end} ({len(seq)} nt)")
            print(f"    length_ok={hm.length_ok}  A-box mm={hm.a_box_mm}  "
                  f"B-box mm={hm.b_box_mm}  spacing_ok={hm.spacing_ok}")
            print(f"    dimer identity {hm.dimer:.3f}   3' A-tail {hm.a_tail:.3f}")
            print(f"    score {hm.score:.3f}  vs null {mu:.3f}+/-{sd:.3f}  "
                  f"z = {z:+.2f}  (A8 energy z was -52.6)")
    print("\n  A z near 0 means these are Alu-SCALE inverted repeats with the")
    print("  subfamily UNRESOLVED -- still A8's finding, just not upgraded.")


if __name__ == "__main__":
    main()
