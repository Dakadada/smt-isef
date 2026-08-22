"""A2: cost of forcing a window open.

``dG_unfold(iv)`` is the free energy that must be paid to hold every base of
an interval unpaired:

    dG_unfold = -RT * ln P_unpaired(iv)

with ``P_unpaired`` the equilibrium probability that the *entire* interval is
simultaneously single-stranded. Two independent estimators are provided; the
caller picks:

``Method.LOCAL``
    RNAplfold-style. Sliding windows of bounded size and bounded base-pair
    span, via ``RNA.pfl_fold_up``. O(N * W^2). The only option for multi-kb
    windows. Ignores pairing beyond ``max_bp_span``, so it *underestimates*
    the opening cost of anything held by a long-range interaction.

``Method.GLOBAL``
    Exact ratio of constrained to unconstrained partition functions over the
    whole sequence:  P = Z(iv forced unpaired) / Z. O(N^3); use on windows of
    a few kb at most. This is the same quantity RNAup's opening energy reports,
    and is the estimator :mod:`smt.thermo.constrained` builds on.

All energies are kcal/mol. All intervals are 1-based inclusive on the sequence
passed in (construct space, not reference space).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from typing import Optional

from smt.seq.reference import Interval

# ViennaRNA's gas constant, cal/(mol K), and 0 C in K (RNA.GASCONST, RNA.K0).
_GASCONST_CAL = 1.98717
_K0 = 273.15

#: Folding temperature, deg C. 37 C throughout the project (ViennaRNA default).
TEMPERATURE_C = 37.0

#: RT at 37 C, kcal/mol = 1.98717 cal/(mol K) * 310.15 K / 1000.
RT = _GASCONST_CAL * (TEMPERATURE_C + _K0) / 1000.0

# RNAplfold defaults used by the ViennaRNA distribution (-W 80 -L 40 in the
# RNAplfold man page; 80/40 is also what the RNAup/CopraRNA pipelines use).
DEFAULT_WINDOW = 80
DEFAULT_MAX_BP_SPAN = 40

# Probabilities below this are treated as underflow; -RT*ln of it is ~+18
# kcal/mol, already far outside anything an effector can pay for.
_P_FLOOR = 1e-13


class Method(str, Enum):
    LOCAL = "local"
    GLOBAL = "global"


@dataclass(frozen=True)
class Accessibility:
    """Result of one accessibility query.

    p_unpaired:
        Probability the whole interval is simultaneously unpaired (0..1].
    dG_unfold:
        -RT ln p_unpaired, kcal/mol, >= 0 by construction.
    """

    interval: Interval
    p_unpaired: float
    dG_unfold: float
    method: str


def _to_rna(seq: str) -> str:
    return seq.upper().replace("T", "U")


def _dg_from_p(p: float) -> float:
    return -RT * math.log(max(p, _P_FLOOR))


# ---------------------------------------------------------------------------
# local (RNAplfold-style)
# ---------------------------------------------------------------------------


def p_unpaired_local(
    seq: str,
    iv: Interval,
    window_size: int = DEFAULT_WINDOW,
    max_bp_span: int = DEFAULT_MAX_BP_SPAN,
) -> float:
    """P(entire ``iv`` unpaired) from sliding-window partition functions."""
    import RNA

    u = len(iv)
    if u > window_size:
        raise ValueError(
            f"interval length {u} exceeds plfold window {window_size}; "
            "raise window_size or use Method.GLOBAL"
        )
    if iv.end > len(seq) or iv.start < 1:
        raise ValueError(f"{iv} outside sequence of length {len(seq)}")
    # pfl_fold_up returns a 1-based [position][unpaired_length] matrix where
    # entry [i][u] is P(the u-mer ENDING at i is unpaired).
    up = RNA.pfl_fold_up(_to_rna(seq), u, window_size, max_bp_span)
    return float(up[iv.end][u])


# ---------------------------------------------------------------------------
# global (constrained partition function ratio)
# ---------------------------------------------------------------------------


def _pf_kcal(seq: str, unpaired: Optional[Interval] = None) -> float:
    """Ensemble free energy (kcal/mol), optionally with ``unpaired`` forced open."""
    import RNA

    fc = RNA.fold_compound(seq)
    # Rescale Boltzmann factors around the MFE or long sequences overflow.
    _, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    if unpaired is not None:
        for i in range(unpaired.start, unpaired.end + 1):
            # hc_add_up(i, ALL_LOOPS): position i must be unpaired but may sit
            # in any loop type. That is exactly "occupied, cannot base pair".
            fc.hc_add_up(i, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
    _, g = fc.pf()
    return float(g)


def p_unpaired_global(seq: str, iv: Interval) -> float:
    """P(entire ``iv`` unpaired) as Z(constrained)/Z(free) over the whole seq."""
    if iv.end > len(seq) or iv.start < 1:
        raise ValueError(f"{iv} outside sequence of length {len(seq)}")
    rna = _to_rna(seq)
    g_free = _pf_kcal(rna)
    g_open = _pf_kcal(rna, unpaired=iv)
    # G_open >= G_free always; the difference IS the opening energy.
    return math.exp(-(g_open - g_free) / RT)


# ---------------------------------------------------------------------------
# public entry point
# ---------------------------------------------------------------------------


def accessibility(
    seq: str,
    iv: Interval,
    method: Method | str = Method.LOCAL,
    window_size: int = DEFAULT_WINDOW,
    max_bp_span: int = DEFAULT_MAX_BP_SPAN,
) -> Accessibility:
    """P_unpaired and dG_unfold for ``iv`` on ``seq`` (construct coordinates)."""
    m = Method(method)
    if m is Method.LOCAL:
        p = p_unpaired_local(seq, iv, window_size, max_bp_span)
    else:
        p = p_unpaired_global(seq, iv)
    return Accessibility(iv, p, _dg_from_p(p), m.value)


def dG_unfold(
    seq: str,
    iv: Interval,
    method: Method | str = Method.LOCAL,
    window_size: int = DEFAULT_WINDOW,
    max_bp_span: int = DEFAULT_MAX_BP_SPAN,
) -> float:
    """Opening cost of ``iv``, kcal/mol, >= 0."""
    return accessibility(seq, iv, method, window_size, max_bp_span).dG_unfold
