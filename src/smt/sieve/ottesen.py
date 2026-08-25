"""Ottesen 2021 (PMC8395096) splicing tables -- the independent replication set.

WHY THIS DATASET EXISTS IN THE PROJECT. Every sieve grade produced so far came
from GSE307431 alone, which is why ``docs/a5-enrichment-pooled.txt`` had to
write "two groupings of ONE dataset are not independent replication, so the
denominator is 36, not 12." This is a DIFFERENT experiment -- a dose series
rather than a chemistry series -- in the SAME cell line (GM03813, the SMN1-null
line characterised in R3). ``docs/geo-validation.md`` catalogued it in detail
and then explicitly did not download it; this module is that download.

THE SIGN CONVENTION IS INVERTED, AND IT IS VERIFIED ON EVERY LOAD.
The published ``∆PSI`` is **control - treated**. Everywhere else in this project
dPSI is **treated - control** (``smt.rnaseq`` uses F18MOE - ScrMOE; rMATS uses
IncLevelDifference the same way). Scoring signed predictions against the
published column unflipped inverts every conclusion, silently -- the same class
of failure as P's single-strand scan and R3's SMN1/SMN2 tie-break.

Checked over all 4,286 rows carrying the three columns needed to test it:

    dPSI == control - treated :  100.0% within 1e-3  (max residual 3e-4, rounding)
    dPSI == treated - control :    0.0% within 5e-3  (median residual 0.396)

and confirmed semantically as well as arithmetically: in the exon-INCLUSION
sheet a negative published ∆PSI means treated PSI exceeded control PSI, i.e.
inclusion rose, which is what that sheet is supposed to contain.

``load_events`` re-runs that check on every call and raises
``SignConventionError`` rather than returning silently mis-signed data. The
constant is not trusted; it is re-derived from the file each time.

THE CONTROL DESIGN FAVOURS THE DOSE TEST, and it is the authors' design rather
than an assumption of ours: ``10mer`` is a non-targeting control ASO with no
effect on SMN2 exon 7 splicing, and LoA (5 nM Anti-N1) was **topped up to 100 nM
total oligo with the 10mer**. Sequence-independent oligo load is therefore
constant across doses, so the LoA-vs-HiA contrast is purely sequence-dependent.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

DEFAULT_XLSX = Path(
    "/home/dashii/orcd/pool/smt-isef-data/ottesen2021/"
    "Supplementary Tables_8-3-2021.xlsx"
)

#: sheet -> (event type, the dose the sheet was called significant at).
#: Every sheet carries BOTH doses on every row; the split is which dose the
#: authors called the event on, not which dose the row describes.
EVENT_SHEETS = {
    "Supplementary Table 3": ("EIN", "HiA"),
    "Supplementary Table 4": ("EIN", "LoA"),
    "Supplementary Table 5": ("ESK", "HiA"),
    "Supplementary Table 6": ("ESK", "LoA"),
}

#: Rounding slack in the published tables. Residuals max out at 3e-4.
_SIGN_TOL = 1e-3

#: Fraction of checkable rows that must match ``control - treated``.
_SIGN_MIN_AGREement = 0.99


class SignConventionError(RuntimeError):
    """The published dPSI stopped being control - treated. Stop, do not guess."""


@dataclass(frozen=True)
class Event:
    """One splicing event, with dPSI already converted to PROJECT convention."""

    gene: str
    ensembl: str
    etype: str
    called_at: str
    fdr: float | None
    psi_control: float | None
    dpsi_hia: float | None
    dpsi_loa: float | None

    @property
    def retention(self) -> float | None:
        """Fraction of the high-dose effect still present at 5 nM.

        Sign-convention independent: flipping numerator and denominator
        together leaves the ratio unchanged.
        """
        if self.dpsi_hia in (None, 0.0) or self.dpsi_loa is None:
            return None
        return self.dpsi_loa / self.dpsi_hia


def _num(x) -> float | None:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def verify_sign_convention(xlsx: Path = DEFAULT_XLSX) -> tuple[int, float]:
    """Re-derive the convention from the file. Returns (n_checked, agreement)."""
    import openpyxl

    wb = openpyxl.load_workbook(xlsx, read_only=True)
    checked = agree = 0
    for sheet in EVENT_SHEETS:
        ws = wb[sheet]
        it = ws.iter_rows(values_only=True)
        idx = {str(h): i for i, h in enumerate(next(it))}
        for row in it:
            ctl = _num(row[idx["10mer PSI"]])
            hia = _num(row[idx["HiA PSI"]])
            dps = _num(row[idx["HiA ∆PSI"]])
            if None in (ctl, hia, dps):
                continue
            checked += 1
            if abs((ctl - hia) - dps) < _SIGN_TOL:
                agree += 1
    if not checked:
        raise SignConventionError("no checkable rows -- schema changed")
    frac = agree / checked
    if frac < _SIGN_MIN_AGREement:
        raise SignConventionError(
            f"published dPSI matches 'control - treated' for only {frac:.1%} of "
            f"{checked} rows; it was 100%. Re-derive before scoring anything."
        )
    return checked, frac


def load_events(xlsx: Path = DEFAULT_XLSX, verify: bool = True) -> list[Event]:
    """All EIN/ESK events, dPSI flipped into PROJECT convention.

    ``dpsi_hia`` / ``dpsi_loa`` are **treated - control**: positive means the
    ASO INCREASED inclusion. This is the negation of the published column.
    """
    import openpyxl

    if verify:
        verify_sign_convention(xlsx)

    wb = openpyxl.load_workbook(xlsx, read_only=True)
    out: list[Event] = []
    for sheet, (etype, called_at) in EVENT_SHEETS.items():
        ws = wb[sheet]
        it = ws.iter_rows(values_only=True)
        idx = {str(h): i for i, h in enumerate(next(it))}
        for row in it:
            if row[idx["Gene name"]] is None:
                continue
            hia = _num(row[idx["HiA ∆PSI"]])
            loa = _num(row[idx["LoA ∆PSI"]])
            out.append(Event(
                gene=str(row[idx["Gene name"]]).strip(),
                ensembl=str(row[idx["Ensembl ID"]]),
                etype=etype,
                called_at=called_at,
                fdr=_num(row[idx["FDR"]]),
                psi_control=_num(row[idx["10mer PSI"]]),
                # THE FLIP. Published is control - treated.
                dpsi_hia=None if hia is None else -hia,
                dpsi_loa=None if loa is None else -loa,
            ))
    return out


def best_per_gene(events: list[Event]) -> dict[str, Event]:
    """Strongest |dpsi_hia| event per gene, among those carrying both doses."""
    best: dict[str, Event] = {}
    for ev in events:
        if ev.dpsi_hia is None or ev.dpsi_loa is None:
            continue
        cur = best.get(ev.gene)
        if cur is None or abs(ev.dpsi_hia) > abs(cur.dpsi_hia):
            best[ev.gene] = ev
    return best
