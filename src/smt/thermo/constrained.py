"""A2 core: occlusion-constrained accessibility of the exon 7 5' splice site.

The claim this module implements
--------------------------------
An effector (ASO, Cas13 spacer, small protein) bound to an intronic site does
two thermodynamically distinct things. It (i) pays an opening cost to unfold
its own footprint, and (ii) *removes that footprint from the folding ensemble*
of the rest of the transcript. Effect (ii) is the mechanism of interest: with
ISS-N1 no longer available to pair, whatever it was sequestering is released.

We model (ii) exactly. A bound effector is a HARD CONSTRAINT that every base of
its footprint is unpaired and unavailable. Under that constraint we re-run the
partition function and ask for

    A_5ss = P(all 9 nt of the exon 7 5' splice site are simultaneously unpaired)

which is the quantity U1 snRNP has to compete for. The primary output is

    delta_A = A_5ss(effector bound) - A_5ss(baseline)

Estimator
---------
``A_5ss`` is a *joint* probability over a 9-nt block, so base-pair
probabilities alone cannot give it: ``1 - sum_j p_ij`` is a marginal, and the
product over positions assumes independence that does not hold inside a helix.
The exact value is a ratio of constrained partition functions,

    A = Z(footprint blocked AND 5'ss forced unpaired) / Z(footprint blocked)
      = exp( -(G_num - G_den) / RT )

which is what ``Estimator.ENSEMBLE`` (the default) computes. ``Estimator.BPPM``
is provided only as a cross-check: it multiplies per-position marginal unpaired
probabilities from the base-pair probability matrix and is an *upper*-biased
approximation whenever the 5'ss positions pair cooperatively. Both are reported
so the gap is visible rather than hidden.

Units: kcal/mol, probabilities dimensionless. Intervals passed in are REFERENCE
space (NG_008728, 1-based inclusive); they are mapped into construct space
internally.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from typing import Iterable, Optional, Sequence

from smt.seq.construct import GLOBAL_PF_SOFT_LIMIT, Construct, compact_window
from smt.seq.reference import Interval
from smt.thermo.accessibility import RT

#: Feature name of the readout window. Exon 7 last 3 nt + intron 7 first 6 nt.
SS5_FEATURE = "SS5_EXON7"


class Estimator(str, Enum):
    ENSEMBLE = "ensemble"  # exact Z-ratio
    BPPM = "bppm"  # per-position marginals, independence approximation


@dataclass(frozen=True)
class A5ssResult:
    """Accessibility of the 5'ss under one occlusion state."""

    construct: str
    footprint: Optional[Interval]  # reference space, None = baseline
    ss5: Interval  # construct space
    A_5ss: float
    dG_open_5ss: float  # -RT ln A_5ss, kcal/mol
    estimator: str
    max_bp_span: int  # 0 == unrestricted global folding


@dataclass(frozen=True)
class OcclusionDelta:
    """Primary model output: effect of one effector footprint on the 5'ss."""

    construct: str
    footprint: Interval  # reference space
    baseline: A5ssResult
    blocked: A5ssResult

    @property
    def delta_A(self) -> float:
        """A_5ss(effector) - A_5ss(baseline). Positive == 5'ss freed."""
        return self.blocked.A_5ss - self.baseline.A_5ss

    @property
    def delta_dG_open(self) -> float:
        """Change in 5'ss opening cost, kcal/mol. Negative == easier to open."""
        return self.blocked.dG_open_5ss - self.baseline.dG_open_5ss

    @property
    def fold_change(self) -> float:
        return self.blocked.A_5ss / self.baseline.A_5ss


# ---------------------------------------------------------------------------
# fold_compound plumbing
# ---------------------------------------------------------------------------


def _blocked_fc(seq: str, blocked: Sequence[Interval], max_bp_span: int = 0):
    """fold_compound over ``seq`` with every base of ``blocked`` forced unpaired.

    ``hc_add_up(i, CONSTRAINT_CONTEXT_ALL_LOOPS)`` means: position i is unpaired
    in every loop context. That is the correct representation of a site occupied
    by an effector -- it can neither pair nor be removed from the chain.

    ``max_bp_span`` caps the sequence separation of allowed base pairs (0 = no
    cap). Unrestricted global folding of an ~800 nt pre-mRNA over-predicts
    long-range helices relative to what a co-transcriptionally folding
    transcript can reach, so results should be reported over a span ladder
    rather than at a single value. See ``docs`` / the A2 report.
    """
    import RNA

    md = RNA.md()
    if max_bp_span:
        md.max_bp_span = max_bp_span
    fc = RNA.fold_compound(seq, md)
    _, mfe = fc.mfe()
    # Boltzmann rescaling around the MFE; without it an ~800 nt PF overflows.
    fc.exp_params_rescale(mfe)
    for iv in blocked:
        for i in range(iv.start, iv.end + 1):
            fc.hc_add_up(i, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
    return fc


def _ensemble_G(seq: str, blocked: Sequence[Interval], max_bp_span: int = 0) -> float:
    """Ensemble free energy (kcal/mol) of ``seq`` with ``blocked`` forced open."""
    fc = _blocked_fc(seq, blocked, max_bp_span)
    _, g = fc.pf()
    return float(g)


def _merge(intervals: Iterable[Interval]) -> list[Interval]:
    """Merge overlapping/adjacent construct-space intervals."""
    ivs = sorted(intervals, key=lambda x: x.start)
    out: list[Interval] = []
    for iv in ivs:
        if out and iv.start <= out[-1].end + 1:
            out[-1] = Interval(out[-1].start, max(out[-1].end, iv.end))
        else:
            out.append(iv)
    return out


# ---------------------------------------------------------------------------
# estimators
# ---------------------------------------------------------------------------


def _a5ss_ensemble(
    seq: str, ss5: Interval, footprint: Optional[Interval], max_bp_span: int = 0
) -> float:
    blocked = [footprint] if footprint is not None else []
    g_den = _ensemble_G(seq, blocked, max_bp_span)
    g_num = _ensemble_G(seq, _merge(blocked + [ss5]), max_bp_span)
    # G_num >= G_den by construction (constraints only shrink the ensemble).
    return math.exp(-(g_num - g_den) / RT)


def _a5ss_bppm(
    seq: str, ss5: Interval, footprint: Optional[Interval], max_bp_span: int = 0
) -> float:
    blocked = [footprint] if footprint is not None else []
    fc = _blocked_fc(seq, blocked, max_bp_span)
    fc.pf()
    bpp = fc.bpp()  # 1-based upper-triangular matrix, bpp[i][j] for i < j
    n = len(seq)
    p = 1.0
    for i in range(ss5.start, ss5.end + 1):
        paired = 0.0
        for j in range(1, n + 1):
            if j == i:
                continue
            a, b = (i, j) if i < j else (j, i)
            paired += bpp[a][b]
        p *= max(0.0, 1.0 - paired)
    return p


# ---------------------------------------------------------------------------
# public API
# ---------------------------------------------------------------------------


def a_5ss(
    construct: Optional[Construct] = None,
    footprint: Optional[Interval] = None,
    estimator: Estimator | str = Estimator.ENSEMBLE,
    ss5_feature: str = SS5_FEATURE,
    max_bp_span: int = 0,
) -> A5ssResult:
    """P(exon 7 5'ss entirely unpaired), optionally with a footprint occluded.

    Parameters
    ----------
    construct:
        Defaults to :func:`smt.seq.construct.compact_window` (798 nt).
    footprint:
        REFERENCE-space interval occupied by the effector. ``None`` gives the
        unconstrained baseline.
    """
    c = construct if construct is not None else compact_window()
    if len(c) > GLOBAL_PF_SOFT_LIMIT:
        raise ValueError(
            f"construct {c.name} is {len(c)} nt, above the {GLOBAL_PF_SOFT_LIMIT} nt "
            "global partition-function limit; use compact_window() or a local method"
        )
    ss5 = c.feature(ss5_feature)
    fp_local: Optional[Interval] = None
    if footprint is not None:
        fp_local = c.map_interval(footprint)
        if not (fp_local.end < ss5.start or fp_local.start > ss5.end):
            raise ValueError(
                f"footprint {footprint} overlaps the 5'ss readout window; "
                "A_5ss would be trivially 1 and the comparison meaningless"
            )
    seq = c.rna()
    est = Estimator(estimator)
    p = (
        _a5ss_ensemble(seq, ss5, fp_local, max_bp_span)
        if est is Estimator.ENSEMBLE
        else _a5ss_bppm(seq, ss5, fp_local, max_bp_span)
    )
    p = min(max(p, 0.0), 1.0)
    dg = -RT * math.log(p) if p > 0 else float("inf")
    return A5ssResult(c.name, footprint, ss5, p, dg, est.value, max_bp_span)


def occlusion_delta(
    footprint: Interval,
    construct: Optional[Construct] = None,
    estimator: Estimator | str = Estimator.ENSEMBLE,
    ss5_feature: str = SS5_FEATURE,
    max_bp_span: int = 0,
) -> OcclusionDelta:
    """Primary output: A_5ss with and without the effector, and their delta."""
    c = construct if construct is not None else compact_window()
    base = a_5ss(c, None, estimator, ss5_feature, max_bp_span)
    blocked = a_5ss(c, footprint, estimator, ss5_feature, max_bp_span)
    return OcclusionDelta(c.name, footprint, base, blocked)


def scan_footprints(
    footprints: dict[str, Interval],
    construct: Optional[Construct] = None,
    estimator: Estimator | str = Estimator.ENSEMBLE,
    max_bp_span: int = 0,
) -> dict[str, OcclusionDelta]:
    """Run :func:`occlusion_delta` for several named footprints on one construct.

    Shares nothing between runs except the construct; each baseline is
    recomputed, which is cheap relative to the constrained pass and keeps the
    results independently reproducible.
    """
    c = construct if construct is not None else compact_window()
    return {
        name: occlusion_delta(iv, c, estimator, max_bp_span=max_bp_span)
        for name, iv in footprints.items()
    }
