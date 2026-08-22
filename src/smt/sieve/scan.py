"""A5 step 2: transcriptome-wide off-target hybridization site scan.

Given an oligo (default: nusinersen), find every place in the
:mod:`smt.sieve.transcriptome` index it could plausibly hybridize, and score
each surviving site with the SAME thermodynamic primitives the rest of the
project uses -- IntaRNA for duplex energy, RNAplfold-style local accessibility
for opening cost. No energies are reimplemented here.

Two-stage design (for tractability)
------------------------------------
Stage 1 -- seed (cheap, pure Python, this module's ``seed_scan_transcript``):
    An ASO/oligo hybridizes to its target by Watson-Crick pairing, antiparallel,
    so the transcript sequence a perfectly-matched oligo would bind is the
    REVERSE COMPLEMENT of the oligo -- that is the ``target`` string searched
    for below. The target is split into non-overlapping tiling seeds of length
    ``seed_k`` (7-10 nt, see :func:`default_seed_k`); every EXACT occurrence of
    any seed in a transcript unit's sequence anchors a candidate window, which
    is then extended to the oligo's full length and scored by plain Hamming
    distance (no indels -- an ASO duplex is not gapped). Windows with at most
    ``max_mismatches`` survive to stage 2. Non-overlapping tiling means a
    target with <= 1 mismatch is GUARANTEED to be found (the mismatch can only
    land in one of the two tiles), which is what makes "near-exact" meaningful
    rather than aspirational.

Stage 2 -- thermodynamics (expensive, only run on stage-1 survivors):
    ``dG_hybrid`` from :func:`smt.thermo.hybridization.dG_hybrid` (IntaRNA's
    duplex-only ``E_hybrid`` -- no accessibility penalty, see that module's
    docstring on why ``E_hybrid`` and not ``E``), and ``dG_unfold`` from
    :func:`smt.thermo.accessibility.dG_unfold` (RNAplfold-style local
    accessibility of the candidate window in ITS transcript unit's sequence,
    i.e. the opening cost of the site given whatever else surrounds it). The
    reported quantity is

        ddG_net = dG_hybrid + dG_unfold

    which is directly analogous to how :mod:`smt.thermo.constrained` combines
    a duplex term with an independent opening-cost term (see that module's
    docstring on why ``E_hybrid`` and not IntaRNA's own ``E`` must be used
    when adding an independently computed accessibility penalty, to avoid
    double-counting). More negative ddG_net = more thermodynamically
    favorable off-target site.

Scan units and coordinates
---------------------------
Each transcript contributes scan units: one ``"mrna"`` unit (the full spliced
mRNA, so any hit here is by construction exonic or UTR) and one unit per
intron flank (``"intron:{index}:{side}"``, see
:mod:`smt.sieve.transcriptome`). All positions reported here are 1-based
inclusive, both in unit-local (transcript-oriented) coordinates and in
GENOMIC coordinates (via :meth:`Transcript.mrna_pos_to_genomic` /
:meth:`Transcript.flank_pos_to_genomic`), so :mod:`smt.sieve.polarity` can
classify every hit without re-touching the genome.

Sharding is idempotent, copied from :mod:`smt.validate.length_sweep`: the
transcript list order is deterministic (index file order), ``shard_of`` is a
stride partition over it, each shard skips its output file if already
present, and writes atomically via a ``.tmp`` -> rename.
"""

from __future__ import annotations

import argparse
import csv
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from smt.seq.reference import Interval
from smt.sieve.transcriptome import INDEX_PATH, Transcript, load_index, revcomp
from smt.thermo.accessibility import (
    DEFAULT_MAX_BP_SPAN,
    DEFAULT_WINDOW,
    Method,
    dG_unfold,
)
from smt.thermo.hybridization import IntaRNAError, run_intarna

# Nusinersen (Spinraza), the SMN2 ISS-N1-targeting 18-mer 2'-MOE/PS ASO.
# Sequence per Hua et al. 2010 (Genes Dev 24:1634) / Rigo et al. 2014 (J
# Pharmacol Exp Ther 350:46) / the approved drug label. Given in DNA letters
# by project convention (matches how ISS_N1 etc. are written in
# smt.seq.construct); the real molecule is 2'-MOE, phosphorothioate-linked.
# Its reverse complement (CCAGCATTATGAAAGTGA) covers construct.ISS_N1
# (CCAGCATTATGAAAG, 15 nt) plus 3 more nt of intron 7 -- consistent with the
# published footprint (intron 7 +10..+27 for the full 18-mer vs. the 15-nt
# ISS_N1 annotation of +10..+24).
NUSINERSEN = "TCACTTTCATAATGCTGG"

#: Seed length bounds, nt. Task spec asks for "k around 7-10"; the actual
#: value used for a given oligo is picked by `default_seed_k` to exactly tile
#: the oligo length in two (occasionally unequal) non-overlapping pieces.
SEED_K_MIN = 7
SEED_K_MAX = 10

#: Mismatches (Hamming distance) tolerated in the full-length extension before
#: a stage-1 candidate is dropped without ever reaching IntaRNA. Chosen as a
#: generous upper bound (~20% of an 18-mer) so the THERMODYNAMIC stage -- not
#: an arbitrary alignment cutoff -- makes the final call on marginal sites.
#: This is a tractability knob, not a biological threshold.
DEFAULT_MAX_MISMATCHES = 3

#: Extra sequence context (nt) given to IntaRNA on each side of a stage-1
#: candidate window, so it can explore nearby duplex registers rather than
#: being force-fit to our exact seed-derived window. Small because the
#: duplex itself is only as long as the oligo; this just avoids an edge
#: artifact at the window boundary.
IntaRNA_PAD = 15


def default_seed_k(oligo_len: int) -> int:
    """Seed length that tiles ``oligo_len`` into exactly two non-overlapping
    pieces, clamped to [SEED_K_MIN, SEED_K_MAX].

    Two non-overlapping tiles covering the whole oligo guarantee that any
    true site with at most 1 mismatch is found by at least one exact seed
    (pigeonhole: one mismatch can only fall in one of the two halves).
    """
    k = oligo_len // 2
    return max(SEED_K_MIN, min(SEED_K_MAX, k))


def tiling_seeds(target: str, k: int) -> list[tuple[int, str]]:
    """Non-overlapping seeds of length ``k`` tiling ``target`` left to right.

    Returns ``(offset, seed_seq)`` with ``offset`` 0-based into ``target``.
    The last seed is shifted left (not shortened) to still be length ``k`` and
    end exactly at ``len(target)``, so it may overlap its predecessor by a few
    nt when ``len(target)`` is not a multiple of ``k`` -- full coverage is
    kept, at the cost of the pigeonhole guarantee applying loosely at that
    one boundary.
    """
    n = len(target)
    if k > n:
        raise ValueError(f"seed length {k} exceeds target length {n}")
    seeds = []
    offset = 0
    while offset + k <= n:
        seeds.append((offset, target[offset : offset + k]))
        offset += k
    if offset < n:
        last_off = n - k
        if not seeds or seeds[-1][0] != last_off:
            seeds.append((last_off, target[last_off : last_off + k]))
    return seeds


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        raise ValueError(f"length mismatch {len(a)} vs {len(b)}")
    return sum(1 for x, y in zip(a, b) if x != y)


# ---------------------------------------------------------------------------
# stage 1: seed + extend
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Candidate:
    """One stage-1 survivor: a Hamming-scored window against the oligo target."""

    transcript_id: str
    gene_name: str
    chrom: str
    strand: str
    unit: str  # "mrna" or "intron:{index}:{side}"
    local_start: int  # 1-based inclusive, unit-local, transcript orientation
    local_end: int
    mismatches: int
    g_start: int  # genomic, 1-based inclusive, ascending regardless of strand
    g_end: int


def _find_all(haystack: str, needle: str) -> list[int]:
    """All 0-based start positions of exact occurrences of ``needle``."""
    out = []
    start = 0
    while True:
        i = haystack.find(needle, start)
        if i == -1:
            return out
        out.append(i)
        start = i + 1


def _scan_unit(
    unit_seq: str,
    target: str,
    seeds: list[tuple[int, str]],
    max_mismatches: int,
) -> list[tuple[int, int]]:
    """Return ``(window_start_0based, mismatches)`` for surviving windows.

    Deduplicated across seeds that anchor the same window (a window can be
    found by more than one seed when both happen to match).
    """
    n = len(target)
    seen: dict[int, int] = {}
    for seed_off, seed_seq in seeds:
        for hit in _find_all(unit_seq, seed_seq):
            win_start = hit - seed_off
            if win_start < 0 or win_start + n > len(unit_seq):
                continue
            if win_start in seen:
                continue
            mm = hamming(unit_seq[win_start : win_start + n], target)
            if mm <= max_mismatches:
                seen[win_start] = mm
    return sorted(seen.items())


def seed_scan_transcript(
    t: Transcript,
    oligo: str,
    seed_k: Optional[int] = None,
    max_mismatches: int = DEFAULT_MAX_MISMATCHES,
) -> list[Candidate]:
    """Stage 1 for one transcript, across its mRNA unit and all intron flanks."""
    target = revcomp(oligo.upper().replace("U", "T"))
    k = seed_k if seed_k is not None else default_seed_k(len(target))
    seeds = tiling_seeds(target, k)
    n = len(target)
    out: list[Candidate] = []

    for win_start, mm in _scan_unit(t.mrna_seq, target, seeds, max_mismatches):
        local_start, local_end = win_start + 1, win_start + n
        g_a = t.mrna_pos_to_genomic(local_start)
        g_b = t.mrna_pos_to_genomic(local_end)
        out.append(
            Candidate(
                t.transcript_id, t.gene_name, t.chrom, t.strand, "mrna",
                local_start, local_end, mm, min(g_a, g_b), max(g_a, g_b),
            )
        )

    for f in t.flanks:
        unit = f"intron:{f.intron_index}:{f.side}"
        for win_start, mm in _scan_unit(f.seq, target, seeds, max_mismatches):
            local_start, local_end = win_start + 1, win_start + n
            g_a = t.flank_pos_to_genomic(f, local_start)
            g_b = t.flank_pos_to_genomic(f, local_end)
            out.append(
                Candidate(
                    t.transcript_id, t.gene_name, t.chrom, t.strand, unit,
                    local_start, local_end, mm, min(g_a, g_b), max(g_a, g_b),
                )
            )
    return out


# ---------------------------------------------------------------------------
# stage 2: thermodynamics on survivors only
# ---------------------------------------------------------------------------


def _unit_seq(t: Transcript, unit: str) -> str:
    if unit == "mrna":
        return t.mrna_seq
    _, idx_s, side = unit.split(":")
    idx = int(idx_s)
    for f in t.flanks:
        if f.intron_index == idx and f.side == side:
            return f.seq
    raise KeyError(f"unit {unit} not found on {t.transcript_id}")


def score_candidate(t: Transcript, c: Candidate, oligo: str) -> Optional[dict]:
    """Stage 2: dG_hybrid + dG_unfold for one candidate. None if IntaRNA finds
    no interaction in the padded window (site is not hybridization-competent)."""
    seq = _unit_seq(t, c.unit)
    n = len(seq)
    pad_start = max(0, c.local_start - 1 - IntaRNA_PAD)
    pad_end = min(n, c.local_end + IntaRNA_PAD)
    window = seq[pad_start:pad_end]

    try:
        recs = run_intarna(window, oligo, target_id=c.unit, query_id="oligo", n_best=1)
    except IntaRNAError:
        return None
    if not recs:
        return None
    dg_hybrid = recs[0].dG_hybrid

    win_size = min(DEFAULT_WINDOW, n)
    span = min(DEFAULT_MAX_BP_SPAN, win_size)
    iv = Interval(c.local_start, c.local_end)
    try:
        dg_unfold = dG_unfold(seq, iv, method=Method.LOCAL, window_size=win_size, max_bp_span=span)
    except Exception:
        return None

    return {
        "transcript_id": c.transcript_id,
        "gene_name": c.gene_name,
        "chrom": c.chrom,
        "strand": c.strand,
        "unit": c.unit,
        "local_start": c.local_start,
        "local_end": c.local_end,
        "mismatches": c.mismatches,
        "g_start": c.g_start,
        "g_end": c.g_end,
        "dG_hybrid": f"{dg_hybrid:.4f}",
        "dG_unfold": f"{dg_unfold:.4f}",
        "ddG_net": f"{dg_hybrid + dg_unfold:.4f}",
    }


FIELDS = [
    "transcript_id", "gene_name", "chrom", "strand", "unit",
    "local_start", "local_end", "mismatches", "g_start", "g_end",
    "dG_hybrid", "dG_unfold", "ddG_net",
]


# ---------------------------------------------------------------------------
# sharding (copied pattern from smt.validate.length_sweep)
# ---------------------------------------------------------------------------


def shard_of(items: list, shard: int, nshards: int) -> list:
    """Stride partition. Union over all shards == items, with no overlap."""
    if not 0 <= shard < nshards:
        raise ValueError(f"shard {shard} out of range for nshards {nshards}")
    return items[shard::nshards]


def run_shard(
    oligo: str,
    shard: int,
    nshards: int,
    out_dir: Path,
    index_path: Path = INDEX_PATH,
    seed_k: Optional[int] = None,
    max_mismatches: int = DEFAULT_MAX_MISMATCHES,
) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    final = out_dir / f"scan_{shard:04d}_of_{nshards:04d}.csv"
    if final.exists():
        print(f"shard {shard}: output exists, skipping -> {final}")
        return final

    transcripts = load_index(index_path)
    items = shard_of(transcripts, shard, nshards)
    print(f"shard {shard}/{nshards}: {len(items)} transcripts")

    tmp = final.with_suffix(".csv.tmp")
    n_cand = 0
    n_scored = 0
    with tmp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        for i, t in enumerate(items, 1):
            cands = seed_scan_transcript(t, oligo, seed_k, max_mismatches)
            n_cand += len(cands)
            for c in cands:
                row = score_candidate(t, c, oligo)
                if row is not None:
                    w.writerow(row)
                    n_scored += 1
            if i % 500 == 0:
                fh.flush()
                print(f"  {i}/{len(items)} transcripts, {n_cand} stage-1 candidates, "
                      f"{n_scored} scored", flush=True)
    os.replace(tmp, final)  # atomic
    print(f"shard {shard}: {n_cand} stage-1 candidates, {n_scored} scored -> {final}")
    return final


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--shard", type=int, required=True)
    p.add_argument("--nshards", type=int, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--oligo", type=str, default=NUSINERSEN)
    p.add_argument("--index", type=Path, default=INDEX_PATH)
    p.add_argument("--seed-k", type=int, default=None)
    p.add_argument("--max-mismatches", type=int, default=DEFAULT_MAX_MISMATCHES)
    a = p.parse_args()
    run_shard(a.oligo, a.shard, a.nshards, a.out, a.index, a.seed_k, a.max_mismatches)


if __name__ == "__main__":
    main()
