"""A5 step 3: positional polarity -- turning a binding site into a SIGNED
predicted splicing effect.

The biology this module encodes
--------------------------------
The same short RNA motif can act as an intronic SPLICING SILENCER or an
exonic SPLICING ENHANCER depending only on which side of a splice site it
sits on (this is exactly nusinersen's own mechanism: it silences ISS-N1, an
intronic silencer in SMN2 intron 7, which is therapeutic; the SAME oligo also
binds an exonic splicing enhancer in POLR2H exon 2, which is the documented
off-target). A hybridization-competent site (:mod:`smt.sieve.scan`'s
``ddG_net``) is therefore NOT itself a splicing prediction -- it is a binding
affinity. This module is the step that turns binding + position into a signed
direction of effect.

Sign convention -- MUST match the rest of the project
-------------------------------------------------------
Identical to :class:`smt.thermo.constrained.OcclusionDelta`: POSITIVE means
the containing exon's splicing machinery is FREED, i.e. INCLUSION is
promoted. NEGATIVE means SKIPPING is promoted.

    site in an INTRON, near a splice site (silencer-like, e.g. ISS-N1)
        -> occluding it REMOVES a silencer -> INCLUSION promoted -> POSITIVE.
    site in an EXON (enhancer-like, e.g. the POLR2H ESE)
        -> occluding it REMOVES an enhancer -> SKIPPING promoted -> NEGATIVE.
    site in a UTR
        -> this project has no enhancer/silencer model for UTR-bound sites;
           reported with sign 0 (an explicit non-prediction, not a guessed
           small effect -- see ``Category.FIVE_UTR`` / ``THREE_UTR`` below).

This is a MODELING CHOICE (a "when in an exon, assume ESE; when in an intron
near a splice site, assume ISS-like silencer" prior), not a measurement. Every
real regulatory element has to be experimentally typed as an enhancer or
silencer; this module assumes the typing from POSITION ALONE, which is the
same simplifying assumption the field calls "the splicing code's positional
bias" (canonical intronic splicing silencers and exonic splicing enhancers
are both enriched near splice sites, but the enrichment is statistical, not
absolute -- there are excpetions in both directions). This is flagged here
explicitly, not asserted as fact.

Distance-to-splice-site magnitude modulation -- A MODELING ASSUMPTION
------------------------------------------------------------------------
Both classes of element (ISS-N1-like intronic silencers, ESE-like exonic
enhancers) are enriched close to the splice site they regulate; potency is
observed, across the splicing literature broadly, to fall off with distance.
We model that fall-off as a single exponential:

    weight(d) = exp(-d / DECAY_NT)

with ``d`` the distance in nt from the CENTER of the occluded site to the
nearest splice site of its exon, and ``DECAY_NT`` (default 30 nt, see the
constant below) a length scale picked to be roughly the span over which this
project's OWN blind falsification sweep (:mod:`smt.validate.length_sweep`,
gated in ``docs/HANDOFF.md``) found the exon 7 5'ss accessibility effect of an
occluded footprint concentrated -- NOT fitted to any measured PSI value (this
project has none, and placeholder PSI calibration is explicitly forbidden).
Treat ``DECAY_NT`` as an assumed parameter, reported alongside every score,
not a derived constant. A different, equally defensible choice (e.g. a step
function, or a slower power-law decay) would change absolute magnitudes;
sign and rank order among sites of the same category are far more robust to
this choice than the magnitude is.

Final score
-----------
    signed_effect = sign(category) * weight(distance_nt) * max(0, -ddG_net)

``max(0, -ddG_net)`` means only thermodynamically FAVORABLE sites (negative
ddG_net) produce a nonzero prediction; a site IntaRNA/accessibility together
score as net-unfavorable is treated as non-binding, not as a weak effect in
some direction. All quantities: ddG_net in kcal/mol (from
:mod:`smt.sieve.scan`), distance in nt (integer, genomic), signed_effect in
kcal/mol-scaled but UNITLESS as a splicing-outcome predictor -- it is a
RANKING SCORE, not a PSI or fold-change prediction, because this project has
no PSI data to calibrate a physical outcome scale against.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional

from smt.sieve.transcriptome import INDEX_PATH, Transcript, load_index

#: Exponential decay length, nt, for splice-site-distance modulation of
#: magnitude. ASSUMPTION -- see module docstring. Not fitted to any PSI data.
DEFAULT_DECAY_NT = 30.0

#: Distance (nt) inside which a hit is additionally labeled "proximal" for
#: reporting purposes. This threshold does NOT change the continuous
#: exponential weighting above; it is an informational category boundary
#: only, chosen to match the span ISS-N1 itself occupies (intron 7
#: +10..+24 of SMN2, i.e. 10-24 nt from the 5'ss -- see smt.seq.construct).
PROXIMITY_NT = 30


class Category(str, Enum):
    EXONIC_CORE = "exonic_core"
    EXONIC_PROXIMAL = "exonic_proximal"
    INTRONIC_CORE = "intronic_core"
    INTRONIC_PROXIMAL = "intronic_proximal"
    FIVE_UTR = "five_utr"
    THREE_UTR = "three_utr"


#: sign(category): +1 promotes inclusion (silencer disrupted), -1 promotes
#: skipping (enhancer disrupted), 0 = no polarity model (UTR). See module
#: docstring for the biological reasoning and its caveats.
CATEGORY_SIGN: dict[Category, int] = {
    Category.EXONIC_CORE: -1,
    Category.EXONIC_PROXIMAL: -1,
    Category.INTRONIC_CORE: +1,
    Category.INTRONIC_PROXIMAL: +1,
    Category.FIVE_UTR: 0,
    Category.THREE_UTR: 0,
}


@dataclass(frozen=True)
class PolarityCall:
    category: Category
    distance_nt: int  # nt from hit center to nearest splice site of its exon/intron
    sign: int
    weight: float  # exp(-distance/decay), in (0, 1]
    signed_effect: float  # sign * weight * max(0, -ddG_net)


def _exon_containing(t: Transcript, t_pos: int):
    for ex in t.exons:
        if ex.t_start <= t_pos <= ex.t_end:
            return ex
    return None


def _mrna_distance_to_splice_site(t: Transcript, t_pos: int) -> int:
    """Distance (nt, transcript space) from ``t_pos`` to the nearer splice
    site of the exon containing it. The first exon's upstream edge (TSS) and
    the last exon's downstream edge (poly-A site) are not splice sites and
    are excluded from the comparison.
    """
    ex = _exon_containing(t, t_pos)
    if ex is None:
        raise ValueError(f"{t_pos} is not inside any exon of {t.transcript_id}")
    candidates = []
    if ex.order > 1:  # upstream edge is a real 3'ss (acceptor), not the TSS
        candidates.append(t_pos - ex.t_start)
    if ex.order < len(t.exons):  # downstream edge is a real 5'ss (donor), not poly-A
        candidates.append(ex.t_end - t_pos)
    if not candidates:
        # Single-exon transcript: no splice sites at all. Treat as maximally
        # far so the exponential weight -> 0 rather than raising.
        return 10 ** 9
    return min(candidates)


def classify_hit(
    t: Transcript,
    unit: str,
    local_start: int,
    local_end: int,
    proximity_nt: int = PROXIMITY_NT,
) -> tuple[Category, int]:
    """Category + distance-to-nearest-splice-site (nt) for one scan hit.

    ``unit``/``local_start``/``local_end`` are exactly the fields
    :mod:`smt.sieve.scan` writes: ``local_start``/``local_end`` are 1-based
    inclusive in the unit's own local (transcript-oriented) coordinates.
    """
    center = (local_start + local_end) // 2
    if unit == "mrna":
        kind = t.genomic_region_kind(t.mrna_pos_to_genomic(center))
        if kind == "five_utr":
            return Category.FIVE_UTR, 0
        if kind == "three_utr":
            return Category.THREE_UTR, 0
        d = _mrna_distance_to_splice_site(t, center)
        return (Category.EXONIC_PROXIMAL if d <= proximity_nt else Category.EXONIC_CORE), d

    # intron flank unit: "intron:{index}:{side}"
    _, idx_s, side = unit.split(":")
    idx = int(idx_s)
    flank = next(f for f in t.flanks if f.intron_index == idx and f.side == side)
    span = local_end - local_start + 1
    d = flank.dist_to_splice_site(local_start, span)
    return (Category.INTRONIC_PROXIMAL if d <= proximity_nt else Category.INTRONIC_CORE), d


def score_polarity(
    category: Category,
    distance_nt: int,
    ddG_net: float,
    decay_nt: float = DEFAULT_DECAY_NT,
) -> PolarityCall:
    sign = CATEGORY_SIGN[category]
    weight = math.exp(-distance_nt / decay_nt)
    magnitude = max(0.0, -ddG_net)  # only favorable (negative ddG_net) binding counts
    return PolarityCall(category, distance_nt, sign, weight, sign * weight * magnitude)


def annotate_row(row: dict, index_by_id: dict[str, Transcript], decay_nt: float = DEFAULT_DECAY_NT) -> dict:
    """Add polarity columns to one row produced by :mod:`smt.sieve.scan`."""
    t = index_by_id[row["transcript_id"]]
    local_start, local_end = int(row["local_start"]), int(row["local_end"])
    category, distance_nt = classify_hit(t, row["unit"], local_start, local_end)
    call = score_polarity(category, distance_nt, float(row["ddG_net"]), decay_nt)
    out = dict(row)
    out["category"] = call.category.value
    out["distance_to_splice_site_nt"] = call.distance_nt
    out["sign"] = call.sign
    out["weight"] = f"{call.weight:.6f}"
    out["signed_effect"] = f"{call.signed_effect:.6f}"
    return out


ANNOTATED_FIELDS = [
    "transcript_id", "gene_name", "chrom", "strand", "unit",
    "local_start", "local_end", "mismatches", "g_start", "g_end",
    "dG_hybrid", "dG_unfold", "ddG_net",
    "category", "distance_to_splice_site_nt", "sign", "weight", "signed_effect",
]


def annotate_scan_dir(
    scan_dir: Path,
    out_path: Path,
    index_path: Path = INDEX_PATH,
    decay_nt: float = DEFAULT_DECAY_NT,
) -> Path:
    """Merge every ``scan_*.csv`` shard in ``scan_dir``, add polarity columns.

    Idempotent (skip if ``out_path`` exists) and atomic (``.tmp`` -> rename),
    same discipline as :mod:`smt.validate.length_sweep` and
    :mod:`smt.sieve.scan`.
    """
    if out_path.exists():
        print(f"annotated output exists, skipping -> {out_path}")
        return out_path
    transcripts = load_index(index_path)
    index_by_id = {t.transcript_id: t for t in transcripts}

    shard_files = sorted(scan_dir.glob("scan_*_of_*.csv"))
    if not shard_files:
        raise FileNotFoundError(f"no scan_*_of_*.csv shards in {scan_dir}")

    tmp = out_path.with_suffix(".csv.tmp")
    n = 0
    with tmp.open("w", newline="") as out_fh:
        w = csv.DictWriter(out_fh, fieldnames=ANNOTATED_FIELDS)
        w.writeheader()
        for shard_path in shard_files:
            with shard_path.open() as in_fh:
                for row in csv.DictReader(in_fh):
                    w.writerow(annotate_row(row, index_by_id, decay_nt))
                    n += 1
    os.replace(tmp, out_path)
    print(f"annotated {n} rows from {len(shard_files)} shards -> {out_path}")
    return out_path


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--scan-dir", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--index", type=Path, default=INDEX_PATH)
    p.add_argument("--decay-nt", type=float, default=DEFAULT_DECAY_NT)
    a = p.parse_args()
    annotate_scan_dir(a.scan_dir, a.out, a.index, a.decay_nt)


if __name__ == "__main__":
    main()
