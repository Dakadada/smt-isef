"""A2: effector-target hybridization energies via the IntaRNA binary.

We do not reimplement duplex thermodynamics. IntaRNA 3.x is the primitive; this
module is a typed subprocess wrapper plus a CSV parser.

Energy bookkeeping (IntaRNA 3 conventions, all kcal/mol, all negative-is-good)
-----------------------------------------------------------------------------
``E_hybrid``
    Duplex formation energy alone: stacking + loops + init + dangles + helix
    ends, for the intermolecular helix. Contains NO accessibility penalty.
``ED1`` / ``ED2``
    Accessibility (opening) penalties for making the interacting site of the
    target (1) and the query (2) single stranded. Positive numbers. These are
    the same quantity as :func:`smt.thermo.accessibility.dG_unfold`, computed
    by IntaRNA with its own RNAplfold-style windows.
``E``
    The total interaction energy IntaRNA optimises: ``E ~= E_hybrid + ED1 +
    ED2``. It ALREADY CONTAINS the opening penalties.

So: if you are combining an IntaRNA number with an independently computed
opening cost (as :mod:`smt.thermo.constrained` does), use ``E_hybrid`` and add
your own penalty. Using ``E`` there double-counts ED1.
"""

from __future__ import annotations

import csv
import io
import shutil
import subprocess
from dataclasses import dataclass
from typing import Optional, Sequence

# Columns requested from IntaRNA. Names come from `IntaRNA --fullhelp`
# (--outCsvCols); IntaRNA emits them ';'-separated with a header line.
CSV_COLS: tuple[str, ...] = (
    "id1",
    "start1",
    "end1",
    "id2",
    "start2",
    "end2",
    "subseqDP",
    "hybridDP",
    "E",
    "E_hybrid",
    "ED1",
    "ED2",
    "Pu1",
    "Pu2",
)

CSV_SEP = ";"  # IntaRNA's fixed CSV separator in --outMode=C
DEFAULT_TIMEOUT_S = 600


class IntaRNAError(RuntimeError):
    pass


@dataclass(frozen=True)
class Hybridization:
    """One IntaRNA interaction record.

    Coordinates ``start1/end1`` are 1-based inclusive on the target sequence as
    handed to IntaRNA (i.e. construct space if you passed a construct window),
    ``start2/end2`` likewise on the query (the effector oligo).
    """

    target_id: str
    query_id: str
    start1: int
    end1: int
    start2: int
    end2: int
    hybrid_dp: str
    E_total: float
    E_hybrid: float
    ED1: float
    ED2: float
    Pu1: float
    Pu2: float

    @property
    def dG_hybrid(self) -> float:
        """Duplex-only energy, no accessibility penalty. kcal/mol."""
        return self.E_hybrid

    @property
    def dG_total(self) -> float:
        """Full interaction energy including ED1 and ED2. kcal/mol."""
        return self.E_total

    @property
    def dG_opening(self) -> float:
        """ED1 + ED2, the accessibility penalty already inside E_total."""
        return self.ED1 + self.ED2


def _to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def _exe() -> str:
    exe = shutil.which("IntaRNA")
    if not exe:
        raise IntaRNAError("IntaRNA binary not on PATH")
    return exe


def parse_intarna_csv(text: str) -> list[Hybridization]:
    """Parse ``--outMode=C`` output. Empty (no interaction found) -> []."""
    text = text.strip()
    if not text:
        return []
    rows = list(csv.DictReader(io.StringIO(text), delimiter=CSV_SEP))
    out: list[Hybridization] = []
    for r in rows:
        missing = [c for c in CSV_COLS if c not in r]
        if missing:
            raise IntaRNAError(f"IntaRNA CSV missing columns {missing}; got {list(r)}")
        out.append(
            Hybridization(
                target_id=r["id1"],
                query_id=r["id2"],
                start1=int(r["start1"]),
                end1=int(r["end1"]),
                start2=int(r["start2"]),
                end2=int(r["end2"]),
                hybrid_dp=r["hybridDP"],
                E_total=float(r["E"]),
                E_hybrid=float(r["E_hybrid"]),
                ED1=float(r["ED1"]),
                ED2=float(r["ED2"]),
                Pu1=float(r["Pu1"]),
                Pu2=float(r["Pu2"]),
            )
        )
    return out


def run_intarna(
    target: str,
    query: str,
    target_id: str = "target",
    query_id: str = "query",
    n_best: int = 1,
    extra_args: Sequence[str] = (),
    timeout_s: int = DEFAULT_TIMEOUT_S,
) -> list[Hybridization]:
    """Call IntaRNA and return parsed interactions, best first.

    ``target`` and ``query`` may be DNA or RNA; both are uppercased and T->U'd
    because IntaRNA is an RNA-RNA predictor.
    """
    cmd = [
        _exe(),
        "-t",
        _to_rna(target),
        "-q",
        _to_rna(query),
        "--outMode=C",
        f"--outCsvCols={','.join(CSV_COLS)}",
        f"--outNumber={n_best}",
        *extra_args,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout_s)
    if proc.returncode != 0:
        raise IntaRNAError(f"IntaRNA failed ({proc.returncode}): {proc.stderr.strip()}")
    recs = parse_intarna_csv(proc.stdout)
    # IntaRNA labels sequences with generated ids; rename to caller's labels.
    return [
        Hybridization(
            target_id, query_id, r.start1, r.end1, r.start2, r.end2,
            r.hybrid_dp, r.E_total, r.E_hybrid, r.ED1, r.ED2, r.Pu1, r.Pu2,
        )
        for r in recs
    ]


def dG_hybrid(
    target: str,
    query: str,
    include_accessibility: bool = False,
    **kwargs,
) -> Optional[float]:
    """Best interaction energy for ``query`` against ``target``, kcal/mol.

    ``include_accessibility=False`` (default) returns the duplex-only
    ``E_hybrid``; True returns IntaRNA's total ``E`` (duplex + ED1 + ED2).
    Returns None when IntaRNA predicts no interaction.
    """
    recs = run_intarna(target, query, **kwargs)
    if not recs:
        return None
    best = recs[0]
    return best.dG_total if include_accessibility else best.dG_hybrid
