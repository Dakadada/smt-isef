"""G1: the off-target axis of the Goldilocks length study.

`smt.validate.goldilocks` produced the ON-TARGET axis -- the best span-robust
footprint at each length 6-20, from the A3 sweep already on disk. Its own
printed conclusion says what is missing:

    "Off-target scoring (G1, targeted panel) attaches to these winners; until
    it runs, this is the on-target axis ALONE and implies nothing about a
    trade-off."

This module supplies the missing axis, and answers the spec's
"without inducing off-target affinity" / "zero friendly fire on POLR2H".

WHY NOT JUST COUNT STAGE-1 HITS
--------------------------------
The obvious implementation -- run `smt.sieve.scan` at each length and count
survivors -- is WRONG, and quietly so. `scan.DEFAULT_MAX_MISMATCHES` is 4,
an ABSOLUTE budget chosen for an 18-mer (~20% of its length). Applied across
lengths it is not one criterion but fifteen: 4 mismatches in a 6-mer is 33%
identity, in a 20-mer it is 80%. A count taken that way falls with length
mostly because the implied identity threshold RISES with length, and the
resulting "curve" would be an artifact of a constant that was never meant to
travel. Nothing about the molecule would be measured.

So G1 uses two axes, each length-comparable by construction:

AXIS 1 -- UNIQUENESS (combinatorial, exhaustive, transcriptome-wide).
    For an oligo of length L, the full mismatch histogram over EVERY window of
    EVERY scan unit in the transcriptome index: how many sites match at 0
    mismatches, at 1, at 2, ... No cutoff is applied, so no cutoff can be
    smuggled in. The comparison across lengths is then made at matched
    FRACTIONAL identity, not matched absolute mismatches. This is exhaustive,
    not seeded -- see `scan._scan_unit` on why seeding loses real hits.

AXIS 2 -- AFFINITY (thermodynamic, targeted panel).
    Mismatch counts are combinatorics, not chemistry: a 6-mer with zero
    mismatches is still a weak duplex. So each length's ASO is also scored
    with the project's own thermodynamics (IntaRNA dG_hybrid + RNAplfold
    dG_unfold, via `scan.score_candidate`) against a panel of the four
    validated rMATS hits -- POLR2H, MAP1S, ANKRD28, CCSER2 -- plus
    length-matched controls. For each transcript the BEST windows are taken
    by mismatch count, however many mismatches those happen to carry, so
    again no threshold enters the cross-length comparison.

THE ACTUAL TRADE-OFF QUANTITY
------------------------------
Neither axis alone is a specificity measure, because a longer ASO binds
EVERYTHING more tightly, its true target included. What matters is the
SEPARATION between the intended duplex and the best unintended one:

    margin = ddG_net(best off-target) - ddG_net(on-target)

in kcal/mol, positive = the off-target duplex is weaker = safer. This is the
quantity that is comparable across lengths, and the one the Goldilocks curve
is drawn on. Nusinersen's 18-mer is scored by the identical code path as a
reference point, though note it has NO span-robust footprint in the A3 scan
(see `goldilocks`), so for 18 the on-target site is taken as nusinersen's
published footprint rather than a scan winner -- flagged in the output.

HONEST EXPECTATION, STATED BEFORE THE NUMBERS
----------------------------------------------
Short ASOs should look bad on Axis 1 for a purely combinatorial reason: a
6-mer has 4^6 = 4096 possible sequences against ~10^8 nt of transcriptome, so
tens of thousands of exact sites are EXPECTED with no biology involved. If the
margin curve is monotonic in length with no interior optimum, THERE IS NO
GOLDILOCKS POINT, and that is the result. It does not get bins tuned until one
appears.

Two of the on-target winners (lengths 6 and 7) target AAUUAA / AAAUUAA -- low
complexity by inspection. If those turn out to be the promiscuous ones, that is
a sequence-composition finding, not a length finding, and must be reported as
such rather than folded into a length curve.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
from collections import defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional

import numpy as np

from smt.seq.construct import compact_window
from smt.seq.reference import Interval
from smt.sieve.scan import (
    NUSINERSEN,
    IntaRNA_PAD,
    Candidate,
)
from smt.sieve.transcriptome import INDEX_PATH, Transcript, load_index, revcomp
from smt.thermo.accessibility import (
    DEFAULT_MAX_BP_SPAN,
    DEFAULT_WINDOW,
    Method,
    dG_unfold,
)
from smt.thermo.hybridization import IntaRNAError, run_intarna
from smt.validate.goldilocks import EDGE_MARGIN, SPANS, SWEEP_DIR, load, robust, worst_span

OUT_DIR = Path("/home/dashii/orcd/pool/smt-isef-data/goldilocks_offtarget")
ANNOTATED = Path("/home/dashii/orcd/pool/smt-isef-data/sieve/annotated.csv")

#: The four rMATS-validated nusinersen off-target genes. POLR2H is the field's
#: canonical one; MAP1S was this project's novel blind prediction (A5).
PANEL_HITS = ("POLR2H", "MAP1S", "ANKRD28", "CCSER2")

#: Control transcripts per hit, matched on TOTAL SCANNED LENGTH. Length is the
#: confound that matters here: a longer transcript offers more windows, so it
#: accumulates more and better off-target hits for reasons having nothing to do
#: with being a real target. Controls are drawn from the same annotated sieve
#: output as the hits, so they are already transcripts the sieve flagged.
CONTROLS_PER_HIT = 5

#: Best windows per transcript promoted to the thermodynamic stage. Taken by
#: rank (lowest mismatch first), NOT by threshold, so the cross-length
#: comparison never depends on a cutoff.
TOP_WINDOWS_PER_UNIT = 3

#: Nusinersen's published intron 7 footprint (+10..+27), reference coords.
#: Used only as the length-18 reference point; see module docstring.
NUSINERSEN_SITE = (32062, 18)

INTRON6_TAIL = 200


# ---------------------------------------------------------------------------
# on-target winners, from the A3 sweep already on disk
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Winner:
    length: int
    start: int
    worst_fom: float
    target_rna: str
    aso: str
    edge_excluded: bool


def _complement(rna: str) -> str:
    return rna.translate(str.maketrans("ACGUacgu", "UGCAugca"))


def candidates_by_length(sweep_dir: Path = SWEEP_DIR) -> dict[int, list[Winner]]:
    """All EDGE-EXCLUDED span-robust footprints per length, best worst-span first.

    Edge exclusion is not optional here. `goldilocks` measured that every
    length 6-10 had its nominal best at 32040, the first scanned position, and
    that dropping three boundary positions cut those scores by 2-3 orders of
    magnitude. Scoring an ASO built from a window artifact would put a fake
    point on the trade-off curve at exactly the lengths where the curve is
    most interesting.

    A RANKED LIST, not a single winner, because the top-ranked footprint is not
    always a physically realisable ASO -- see `feasible`. The list order is the
    A3 ranking and is never re-sorted by anything measured downstream.
    """
    data = load(sweep_dir)
    lo = min(s for s, _ in data)
    by_len: dict[int, list[tuple[float, int]]] = defaultdict(list)
    for (start, length), ladder in data.items():
        if not robust(ladder):
            continue
        if start <= lo + EDGE_MARGIN - 1:
            continue
        by_len[length].append((worst_span(ladder), start))

    c = compact_window(intron6_tail=INTRON6_TAIL)
    out: dict[int, list[Winner]] = {}
    for length in sorted(by_len):
        ws = []
        for fom, start in sorted(by_len[length], reverse=True):
            rna = c.rna_sub(c.map_interval(Interval(start, start + length - 1)))
            ws.append(Winner(length, start, fom, rna, _complement(rna)[::-1], True))
        out[length] = ws
    return out


def winners(sweep_dir: Path = SWEEP_DIR) -> list[Winner]:
    """Top-ranked span-robust footprint per length (feasibility NOT applied)."""
    return [v[0] for _, v in sorted(candidates_by_length(sweep_dir).items())]


def feasible(w: Winner) -> Optional[tuple[float, float]]:
    """(dG_hybrid, dG_unfold) on-target, or None if NO DUPLEX FORMS AT ALL.

    THE GATE THIS EXISTS FOR, measured not assumed. The A3 on-target winners at
    lengths 6 and 7 are UUAAUU and UUAAUUU -- and IntaRNA finds no duplex for
    either against its OWN perfect complement, at any seed setting. A 6-mer
    over the GC core (UGGCAG) binds at -9.09 kcal/mol, while an 8-mer over the
    AU-rich site (CCUUAAUU) still forms nothing, so this is COMPOSITION, not
    length.

    That is a real limit of the framework and it must not be papered over:
    `A_5ss` scores the CONSEQUENCE of occluding a window, never the FEASIBILITY
    of occupying it. An all-AU hexamer has no stable duplex to form, so the
    accessibility gain the scan credits it with is unreachable by any real
    oligo of that sequence. Feasibility is therefore applied as a HARD PRE-GATE
    on the A3 ranking -- declared here, before any off-target number is looked
    at -- rather than as a post-hoc filter on results that came out wrong.
    """
    c = compact_window(intron6_tail=INTRON6_TAIL)
    seq = c.sequence
    if not seq:
        raise ValueError("empty construct sequence")
    iv = c.map_interval(Interval(w.start, w.start + w.length - 1))
    if c.sub(iv) != w.target_rna.replace("U", "T"):
        raise ValueError(
            f"construct slice {c.sub(iv)} != recorded target {w.target_rna}; "
            "coordinate mapping is wrong, do not trust any margin computed here"
        )
    return _ddG_net(seq, iv, w.aso)


def first_feasible(cands: list[Winner]) -> tuple[Optional[Winner], list[Winner]]:
    """(highest-ranked ASO that actually forms a duplex, those skipped above it)."""
    skipped = []
    for w in cands:
        if feasible(w) is not None:
            return w, skipped
        skipped.append(w)
    return None, skipped


def nusinersen_reference() -> Winner:
    """Length-18 reference point. NOT a scan winner -- 18 has no span-robust
    footprint at all (see `goldilocks`), so this uses the published site and is
    flagged `edge_excluded=False` to keep it visually distinct in the output."""
    c = compact_window(intron6_tail=INTRON6_TAIL)
    start, length = NUSINERSEN_SITE
    rna = c.rna_sub(c.map_interval(Interval(start, start + length - 1)))
    aso = _complement(rna)[::-1]
    if aso.replace("U", "T") != NUSINERSEN:
        raise ValueError(
            f"reference site {NUSINERSEN_SITE} yields {aso.replace('U','T')}, "
            f"not the published nusinersen sequence {NUSINERSEN}. "
            "Coordinates or construct changed; fix before trusting anything."
        )
    return Winner(length, start, float("nan"), rna, aso, False)


# ---------------------------------------------------------------------------
# AXIS 1: exhaustive mismatch histogram, transcriptome-wide
# ---------------------------------------------------------------------------


def _units(t: Transcript) -> list[tuple[str, str]]:
    out = [("mrna", t.mrna_seq)]
    for f in t.flanks:
        out.append((f"intron:{f.intron_index}:{f.side}", f.seq))
    return out


def mismatch_histogram(oligo: str, transcripts: list[Transcript]) -> np.ndarray:
    """Counts of transcriptome windows at each mismatch distance from ``oligo``.

    Returns an array of length ``len(oligo)+1``; index m is the number of
    windows matching with exactly m mismatches. EXHAUSTIVE -- every window of
    every scan unit is measured, so this is a complete census and not a
    cutoff-dependent count.

    The searched string is the reverse complement of the oligo, since an ASO
    binds the revcomp of its own sequence (same convention as `scan`).
    """
    target = revcomp(oligo.upper().replace("U", "T"))
    n = len(target)
    q = np.frombuffer(target.encode(), dtype=np.uint8)
    hist = np.zeros(n + 1, dtype=np.int64)
    for t in transcripts:
        for _, seq in _units(t):
            if len(seq) < n:
                continue
            arr = np.frombuffer(seq.encode(), dtype=np.uint8)
            mism = (np.lib.stride_tricks.sliding_window_view(arr, n) != q).sum(1)
            hist += np.bincount(mism, minlength=n + 1)[: n + 1]
    return hist


# ---------------------------------------------------------------------------
# AXIS 2: thermodynamic affinity on a length-matched panel
# ---------------------------------------------------------------------------


def scanned_length(t: Transcript) -> int:
    return sum(len(s) for _, s in _units(t))


def build_panel(
    transcripts: list[Transcript],
    annotated: Path = ANNOTATED,
    controls_per_hit: int = CONTROLS_PER_HIT,
) -> tuple[list[Transcript], list[Transcript]]:
    """(hits, controls). Controls are nearest-in-scanned-length non-hit genes
    drawn from the annotated sieve output, chosen deterministically."""
    with annotated.open(newline="") as fh:
        names = {r["gene_name"] for r in csv.DictReader(fh)}
    by_name = {t.gene_name: t for t in transcripts}

    hits = []
    for g in PANEL_HITS:
        if g not in by_name:
            raise KeyError(f"panel hit {g} absent from the transcriptome index")
        hits.append(by_name[g])

    pool = sorted(
        (t for t in transcripts if t.gene_name in names and t.gene_name not in PANEL_HITS),
        key=lambda t: t.gene_name,
    )
    chosen: dict[str, Transcript] = {}
    for h in hits:
        L = scanned_length(h)
        ranked = sorted(pool, key=lambda t: (abs(scanned_length(t) - L), t.gene_name))
        taken = 0
        for t in ranked:
            if t.gene_name in chosen:
                continue
            chosen[t.gene_name] = t
            taken += 1
            if taken == controls_per_hit:
                break
    return hits, [chosen[k] for k in sorted(chosen)]


def best_windows(t: Transcript, oligo: str, top: int = TOP_WINDOWS_PER_UNIT) -> list[Candidate]:
    """The ``top`` lowest-mismatch windows per scan unit, WHATEVER their
    mismatch count. Rank-based, not threshold-based, so no cutoff leaks into
    the cross-length comparison."""
    target = revcomp(oligo.upper().replace("U", "T"))
    n = len(target)
    q = np.frombuffer(target.encode(), dtype=np.uint8)
    out: list[Candidate] = []
    for unit, seq in _units(t):
        if len(seq) < n:
            continue
        arr = np.frombuffer(seq.encode(), dtype=np.uint8)
        mism = (np.lib.stride_tricks.sliding_window_view(arr, n) != q).sum(1)
        k = min(top, mism.size)
        idx = np.argpartition(mism, k - 1)[:k]
        for i in sorted(idx, key=lambda j: (int(mism[j]), int(j))):
            ls, le = int(i) + 1, int(i) + n
            if unit == "mrna":
                ga, gb = t.mrna_pos_to_genomic(ls), t.mrna_pos_to_genomic(le)
            else:
                _, si, side = unit.split(":")
                f = next(x for x in t.flanks if x.intron_index == int(si) and x.side == side)
                ga, gb = t.flank_pos_to_genomic(f, ls), t.flank_pos_to_genomic(f, le)
            out.append(Candidate(
                t.transcript_id, t.gene_name, t.chrom, t.strand, unit,
                ls, le, int(mism[i]), min(ga, gb), max(ga, gb),
            ))
    return out


#: IntaRNA's default minimum seed is 7 contiguous base pairs. A 6-mer CANNOT
#: form one, so IntaRNA exits non-zero and scores nothing -- which, left alone,
#: would make length 6 appear to have ZERO off-target affinity and rank as the
#: safest length on the trade-off curve. That is the tool refusing to answer,
#: not the molecule being specific. Short oligos therefore get seedBP lowered
#: to their own length. seedBP is a search constraint, not an energy parameter:
#: it widens where a duplex may be found, it does not change what a found
#: duplex is worth. The value used is recorded per length in the output.
INTARNA_DEFAULT_SEED_BP = 7


def seed_bp_for(oligo: str) -> int:
    return min(INTARNA_DEFAULT_SEED_BP, len(oligo))


def _ddG_net(seq: str, iv: Interval, oligo: str) -> Optional[tuple[float, float]]:
    """(dG_hybrid, dG_unfold) for ``oligo`` against ``seq[iv]`` in context.

    Mirrors `scan.score_candidate` -- IntaRNA E_hybrid plus an independently
    computed RNAplfold opening cost -- but routes short oligos past the seed
    limit above. ON-TARGET AND OFF-TARGET BOTH GO THROUGH THIS FUNCTION, so
    the margin between them is a difference of like quantities rather than of
    two separately-parameterised pipelines.
    """
    n = len(seq)
    pad_start = max(0, iv.start - 1 - IntaRNA_PAD)
    pad_end = min(n, iv.end + IntaRNA_PAD)
    try:
        recs = run_intarna(
            seq[pad_start:pad_end], oligo, target_id="t", query_id="oligo",
            n_best=1, extra_args=(f"--seedBP={seed_bp_for(oligo)}",),
        )
    except IntaRNAError:
        return None
    if not recs:
        return None
    win_size = min(DEFAULT_WINDOW, n)
    span = min(DEFAULT_MAX_BP_SPAN, win_size)
    try:
        du = dG_unfold(seq, iv, method=Method.LOCAL, window_size=win_size, max_bp_span=span)
    except Exception:
        return None
    return recs[0].dG_hybrid, du


def best_offtarget(t: Transcript, oligo: str, top: int = TOP_WINDOWS_PER_UNIT) -> Optional[dict]:
    """Most favourable (most negative ddG_net) scored window on one transcript."""
    best = None
    for c in best_windows(t, oligo, top):
        got = _ddG_net(_unit_seq_of(t, c.unit), Interval(c.local_start, c.local_end), oligo)
        if got is None:
            continue
        dh, du = got
        row = {
            "transcript_id": c.transcript_id, "gene_name": c.gene_name,
            "chrom": c.chrom, "strand": c.strand, "unit": c.unit,
            "local_start": c.local_start, "local_end": c.local_end,
            "mismatches": c.mismatches, "g_start": c.g_start, "g_end": c.g_end,
            "dG_hybrid": round(dh, 4), "dG_unfold": round(du, 4),
            "ddG_net": round(dh + du, 4),
        }
        if best is None or row["ddG_net"] < best["ddG_net"]:
            best = row
    return best


def _unit_seq_of(t: Transcript, unit: str) -> str:
    for name, seq in _units(t):
        if name == unit:
            return seq
    raise KeyError(f"unit {unit} not found on {t.transcript_id}")


def on_target_ddG(w: Winner) -> Optional[float]:
    """ddG_net of the ASO against its OWN site in the SMN2 construct.

    Same two terms and the same combination as `scan.score_candidate`
    (IntaRNA E_hybrid + independent RNAplfold opening cost), so the on-target
    and off-target numbers are on one scale and their difference means
    something.
    """
    got = feasible(w)
    return None if got is None else got[0] + got[1]


# ---------------------------------------------------------------------------
# per-length unit of work (idempotent, atomic)
# ---------------------------------------------------------------------------


def run_length(length: int, cands: list[Winner], transcripts: list[Transcript],
               out_dir: Path) -> Path:
    """One length = one unit of work. Idempotent, atomic."""
    out_dir.mkdir(parents=True, exist_ok=True)
    final = out_dir / f"len_{length:02d}.json"
    if final.exists():
        print(f"length {length}: exists, skipping -> {final}")
        return final

    w, skipped = first_feasible(cands)
    payload = {
        "length": length,
        "n_span_robust": len(cands),
        "infeasible_skipped": [
            {"start": x.start, "target_rna": x.target_rna, "aso": x.aso}
            for x in skipped
        ],
    }
    if w is None:
        payload["feasible"] = False
        print(f"length {length}: NO feasible ASO among {len(cands)} span-robust "
              "footprints -- recorded as such", flush=True)
        tmp = final.with_suffix(f".json.tmp.{os.getpid()}")
        tmp.write_text(json.dumps(payload, indent=1))
        os.replace(tmp, final)
        return final

    dh, du = feasible(w)
    print(f"length {length} @ {w.start}: aso={w.aso} target={w.target_rna} "
          f"({len(skipped)} infeasible skipped above it)", flush=True)

    hist = mismatch_histogram(w.aso, transcripts)
    print(f"  axis 1: exact={int(hist[0])} 1mm={int(hist[1])}", flush=True)

    hits, controls = build_panel(transcripts)
    panel = {}
    for group, ts in (("hit", hits), ("control", controls)):
        for t in ts:
            row = best_offtarget(t, w.aso)
            if row is not None:
                row["group"] = group
                panel[t.gene_name] = row
    print(f"  axis 2: {len(panel)}/{len(hits)+len(controls)} panel transcripts scored",
          flush=True)

    payload.update({
        "feasible": True,
        "start": w.start,
        "target_rna": w.target_rna,
        "aso": w.aso,
        "gc_fraction": round(sum(c in "GC" for c in w.aso) / w.length, 4),
        "worst_fom": w.worst_fom,
        "is_scan_winner": w.edge_excluded,
        "seed_bp": seed_bp_for(w.aso),
        "mismatch_histogram": [int(x) for x in hist],
        "on_target_dG_hybrid": round(dh, 4),
        "on_target_dG_unfold": round(du, 4),
        "on_target_ddG_net": round(dh + du, 4),
        "panel": panel,
    })
    tmp = final.with_suffix(f".json.tmp.{os.getpid()}")
    tmp.write_text(json.dumps(payload, indent=1))
    os.replace(tmp, final)
    print(f"  -> {final}", flush=True)
    return final


def all_units() -> list[Winner]:
    """Top-ranked footprint per length plus the nusinersen reference point."""
    return winners() + [nusinersen_reference()]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--task", type=int, required=True,
                    help="0..N-1 over sorted lengths; last task is nusinersen (18)")
    ap.add_argument("--out", type=Path, default=OUT_DIR)
    ap.add_argument("--index", type=Path, default=INDEX_PATH)
    a = ap.parse_args()

    cb = candidates_by_length()
    lengths = sorted(cb)
    tasks: list[tuple[int, list[Winner]]] = [(L, cb[L]) for L in lengths]
    tasks.append((NUSINERSEN_SITE[1], [nusinersen_reference()]))  # length 18 ref

    if not 0 <= a.task < len(tasks):
        raise SystemExit(f"task {a.task} out of range 0..{len(tasks)-1}")
    length, cands = tasks[a.task]
    final = a.out / f"len_{length:02d}.json"
    if final.exists():
        print(f"length {length}: exists, skipping -> {final}")
        return
    print(f"loading transcriptome index {a.index}", flush=True)
    transcripts = load_index(a.index)
    print(f"  {len(transcripts)} transcripts", flush=True)
    run_length(length, cands, transcripts, a.out)


if __name__ == "__main__":
    main()
