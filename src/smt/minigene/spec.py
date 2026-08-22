"""A6: machine- and human-readable construct spec for a compressed minigene.

Emits two things a wet lab actually needs:

* a FASTA record of the compressed sequence;
* a block table -- every retained block in 5'->3' order, with both its
  REFERENCE (natural-locus) and CONSTRUCT (compressed-sequence) coordinates,
  plus the junction coordinates between consecutive blocks (i.e. exactly
  where a restriction site would need to go to reproduce this construct by
  cloning).

Operates on ``smt.minigene.compress.GenePlan`` (for labels/provenance) and
the ``Construct`` built from it (for the actual sequence + coordinate
mapping) -- no new coordinate logic is introduced; block boundaries are
derived directly from ``Construct.window``/``Construct.deletions``, the same
representation ``smt.seq.construct.Construct`` already uses everywhere else
in this codebase.
"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from smt.minigene.compress import GenePlan
from smt.seq.construct import Construct
from smt.seq.reference import Interval


@dataclass(frozen=True)
class Block:
    """One contiguous retained segment, in 5'->3' order."""

    label: str
    ref: Interval  # reference space
    con: Interval  # construct space (contiguous by construction)


def _label_for(iv: Interval, plan: Optional[GenePlan]) -> str:
    """Best-effort human label for a retained segment, using the GenePlan's
    own exon/intron feature names where the segment matches one; falls back
    to a generic "retained_<start>_<end>" label otherwise (e.g. for a segment
    that straddles an exon and a donor/acceptor keep block, which is common
    since exons and their flanking minimal blocks are contiguous)."""
    if plan is None:
        return f"retained_{iv.start}_{iv.end}"
    for name, exon_iv in plan.exons.items():
        if iv == exon_iv:
            return f"exon{name}"
    for iname, ip in plan.introns.items():
        if iv == ip.donor.keep:
            return f"intron{iname}_donor"
        if iv == ip.acceptor.keep:
            return f"intron{iname}_acceptor"
        for extra in ip.extra:
            if iv == extra.keep:
                return extra.name
    if iv == plan.upstream:
        return "promoter_upstream"
    return f"retained_{iv.start}_{iv.end}"


def blocks(construct: Construct, plan: Optional[GenePlan] = None) -> list[Block]:
    """Retained blocks in 5'->3' order, reference + construct coordinates.

    Recomputes segments directly from ``window``/``deletions`` (the same
    interval arithmetic ``Construct._segments`` performs internally) rather
    than depending on that private method, since this module lives outside
    ``smt.seq.construct``.
    """
    out: list[Block] = []
    cursor = construct.window.start
    con_offset = 1
    dels = sorted(construct.deletions, key=lambda d: d.start)
    for d in dels:
        if d.start > cursor:
            ref = Interval(cursor, d.start - 1)
            con = Interval(con_offset, con_offset + len(ref) - 1)
            out.append(Block(_label_for(ref, plan), ref, con))
            con_offset += len(ref)
        cursor = d.end + 1
    if cursor <= construct.window.end:
        ref = Interval(cursor, construct.window.end)
        con = Interval(con_offset, con_offset + len(ref) - 1)
        out.append(Block(_label_for(ref, plan), ref, con))
    return out


@dataclass(frozen=True)
class Junction:
    """Boundary between two consecutive retained blocks -- where a
    restriction site would go to reproduce this construct by cloning."""

    upstream_block: str
    downstream_block: str
    con_pos: int  # construct-space position immediately 5' of the junction
    ref_upstream_end: int  # last retained reference base of the upstream block
    ref_downstream_start: int  # first retained reference base of the downstream block
    excised_ref_span: int  # nt of natural sequence removed at this junction


def junctions(bs: list[Block]) -> list[Junction]:
    out = []
    for a, b in zip(bs, bs[1:]):
        out.append(
            Junction(
                upstream_block=a.label,
                downstream_block=b.label,
                con_pos=a.con.end,
                ref_upstream_end=a.ref.end,
                ref_downstream_start=b.ref.start,
                excised_ref_span=b.ref.start - a.ref.end - 1,
            )
        )
    return out


def to_fasta(construct: Construct, line_width: int = 70) -> str:
    header = f">{construct.name} length={len(construct)} window={construct.window.start}-{construct.window.end} ndeletions={len(construct.deletions)}"
    seq = construct.sequence
    lines = [seq[i : i + line_width] for i in range(0, len(seq), line_width)]
    return header + "\n" + "\n".join(lines) + "\n"


BLOCK_FIELDS = ["order", "label", "ref_start", "ref_end", "ref_length", "con_start", "con_end"]
JUNCTION_FIELDS = [
    "upstream_block", "downstream_block", "con_pos",
    "ref_upstream_end", "ref_downstream_start", "excised_ref_span",
]


def block_table_csv(bs: list[Block]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=BLOCK_FIELDS)
    w.writeheader()
    for i, b in enumerate(bs, 1):
        w.writerow(
            {
                "order": i, "label": b.label,
                "ref_start": b.ref.start, "ref_end": b.ref.end, "ref_length": len(b.ref),
                "con_start": b.con.start, "con_end": b.con.end,
            }
        )
    return buf.getvalue()


def junction_table_csv(js: list[Junction]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=JUNCTION_FIELDS)
    w.writeheader()
    for j in js:
        w.writerow(
            {
                "upstream_block": j.upstream_block, "downstream_block": j.downstream_block,
                "con_pos": j.con_pos, "ref_upstream_end": j.ref_upstream_end,
                "ref_downstream_start": j.ref_downstream_start, "excised_ref_span": j.excised_ref_span,
            }
        )
    return buf.getvalue()


def write_spec(construct: Construct, out_dir: Path, plan: Optional[GenePlan] = None) -> dict[str, Path]:
    """Write <name>.fasta, <name>_blocks.csv, <name>_junctions.csv into
    ``out_dir``. Returns the paths written."""
    out_dir.mkdir(parents=True, exist_ok=True)
    bs = blocks(construct, plan)
    js = junctions(bs)
    paths = {
        "fasta": out_dir / f"{construct.name}.fasta",
        "blocks": out_dir / f"{construct.name}_blocks.csv",
        "junctions": out_dir / f"{construct.name}_junctions.csv",
    }
    paths["fasta"].write_text(to_fasta(construct))
    paths["blocks"].write_text(block_table_csv(bs))
    paths["junctions"].write_text(junction_table_csv(js))
    return paths
