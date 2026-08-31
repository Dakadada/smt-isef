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


def revcomp(seq: str) -> str:
    return seq.upper().replace("U", "T").translate(
        str.maketrans("ACGTN", "TGCAN"))[::-1]


def best_orientation(seq: str) -> tuple[Hallmarks, str]:
    """Hallmarks in whichever orientation scores higher, and which one that was.

    THIS IS A CORRECTION TO THE FIRST VERSION OF P, and it matters.

    An INVERTED repeat pair is, by definition, one arm plus the reverse
    complement of the other -- that is what makes the two arms able to base
    pair, and it is what `alu.scan_flanks` selected on (209 revcomp vs 18
    forward matches). So exactly ONE of the two arms is antisense relative to
    the element it derives from.

    Three of the six hallmarks are strand-specific: the A-box and B-box are
    motifs, their SPACING depends on A preceding B, and a 3' poly-A tail reads
    as a 5' poly-T head when reversed. Only `length_ok` and `dimer` survive a
    flip. Scoring both arms on the same strand therefore guarantees the
    inverted one scores near null NO MATTER WHAT IT IS -- which is exactly the
    pattern the first run reported (hit 1 upstream z = -0.54 against hit 1
    downstream z = +2.78), and that verdict cannot be trusted.

    Both orientations are scored and the better is taken. The null MUST be
    computed the same way -- see `null_scores` -- because a maximum over two
    draws beats a single draw on average, and comparing a max-of-2 observation
    to a max-of-1 null would manufacture significance out of the fix itself.
    """
    fwd, rev = hallmarks(seq), hallmarks(revcomp(seq))
    return (fwd, "forward") if fwd.score >= rev.score else (rev, "revcomp")


def _null_one(seq: str) -> float:
    """Null draw scored under the SAME max-over-orientations rule as the data."""
    return best_orientation(seq)[0].score


def null_scores(seq: str, full: str, trials: int, rng: random.Random) -> tuple[list, list]:
    """Composition-matched shuffles, and random same-length intron windows.

    Every draw is scored with `_null_one`, i.e. the same max-over-both-
    orientations rule applied to the observation. Without that the fix to
    `best_orientation` would inflate every z-score by construction.
    """
    shuffled = []
    chars = list(seq)
    for _ in range(trials):
        rng.shuffle(chars)
        shuffled.append(_null_one("".join(chars)))
    windows = []
    for iv in (UPSTREAM_FLANK, DOWNSTREAM_FLANK):
        span = iv.end - iv.start + 1 - len(seq)
        if span <= 0:
            continue
        for _ in range(trials):
            st = rng.randint(iv.start, iv.start + span)
            windows.append(_null_one(full[st - 1: st - 1 + len(seq)]))
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

    print("  ORIENTATION-AWARE. An inverted-repeat pair has one arm antisense by")
    print("  construction, and three of six hallmarks are strand-specific, so each")
    print("  arm is scored in BOTH orientations and so is every null draw.\n")

    #: z at which an arm is called Alu. Two-sided 0.01 under a normal null;
    #: declared here rather than chosen after seeing the numbers.
    Z_CALL = 2.58
    per_hit: dict[int, dict[str, tuple[float, str]]] = {}

    for i, h in enumerate(hits):
        for side, iv in (("upstream", h.matched_upstream),
                         ("downstream", h.matched_downstream)):
            seq = full[iv.start - 1: iv.end]
            hm, orient = best_orientation(seq)
            sh, wi = null_scores(seq, full, args.trials, rng)
            pool = sh + wi
            mu, sd = statistics.mean(pool), statistics.pstdev(pool)
            z = (hm.score - mu) / sd if sd else float("nan")
            per_hit.setdefault(i, {})[side] = (z, orient)
            print(f"  hit {i} {side:<10} {iv.start}-{iv.end} ({len(seq)} nt) "
                  f"[{orient}]")
            print(f"    length_ok={hm.length_ok}  A-box mm={hm.a_box_mm}  "
                  f"B-box mm={hm.b_box_mm}  spacing_ok={hm.spacing_ok}")
            print(f"    dimer identity {hm.dimer:.3f}   3' A-tail {hm.a_tail:.3f}")
            print(f"    score {hm.score:.3f}  vs null {mu:.3f}+/-{sd:.3f}  "
                  f"z = {z:+.2f}  (A8 energy z was -52.6)")

    print("\n### Verdict -- BOTH arms must be Alu, not one")
    print("  Alu-mediated backsplicing requires an inverted Alu PAIR: a single")
    print(f"  Alu arm juxtaposes nothing. Calling an arm Alu at z >= {Z_CALL}.")
    established = False
    for i, arms in sorted(per_hit.items()):
        zs = {k: v[0] for k, v in arms.items()}
        ok = all(z >= Z_CALL for z in zs.values()) and len(zs) == 2
        established = established or ok
        detail = "  ".join(f"{k} z={v:+.2f} [{arms[k][1]}]" for k, v in zs.items())
        print(f"    hit {i}: {detail}  ->  {'ALU PAIR' if ok else 'not a pair'}")
    if established:
        print("\n  ESTABLISHED: at least one hit is an inverted Alu PAIR. The")
        print("  primate-specificity mechanism for C2B34 backsplicing is supported.")
    else:
        print("\n  NOT ESTABLISHED: no hit has both arms clearing the null.")
        print("  Verdict stays the pre-declared fallback -- Alu-SCALE inverted")
        print("  repeats, subfamily UNRESOLVED. Still A8's finding, not upgraded.")


if __name__ == "__main__":
    main()
