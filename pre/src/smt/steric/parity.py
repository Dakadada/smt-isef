"""Parity metric: put an ASO and an RNP on one measurable axis.

One-sentence definition (also the docstring of :func:`aso_length_for_parity`):
two effectors are AT PARITY when, at a given probe radius, they occlude the
same number of target nucleotides above a shadow threshold and reach equally
far past the hybridized guide:target footprint on each side.

This module never re-derives structure or SASA; it only summarizes
:class:`smt.steric.shadow.ShadowProfile` objects that :mod:`smt.steric.shadow`
already computed, and sweeps :func:`smt.steric.shadow.aso_shadow` over ASO
length to find where that summary matches a reference (typically the RNP).

Units and conventions
----------------------
* All nucleotide counts are counts of *resolved* target residues (see
  :func:`smt.steric.structures.duplex_register`); unresolved density is never
  counted as occluded or unoccluded, it is simply absent from the sum.
* ``register`` is the cross-structure duplex-register coordinate: 0 is the
  5'-most base-paired target nucleotide, so "5' overhang" and "3' overhang"
  below are directional in that common coordinate, not raw seqid.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .shadow import ComplexRecord, ShadowProfile, aso_shadow

#: A nucleotide is counted "occluded" once its shadow reaches this fraction.
#: 0.5 is the natural midpoint of s(i) in [0, 1] and is used only as a
#: reporting threshold -- it does not feed back into how s(i) itself is
#: computed. APPROXIMATION: an arbitrary but standard choice; the sweep
#: functions below accept a different value if a study needs one.
DEFAULT_SHADOW_THRESHOLD: float = 0.5

#: ASO lengths swept by :func:`aso_length_for_parity` by default, matching the
#: 8-20 nt window named in the project brief for a therapeutically plausible
#: antisense oligonucleotide.
DEFAULT_ASO_LENGTHS: tuple[int, ...] = tuple(range(8, 21))


@dataclass(frozen=True)
class ParitySummary:
    """Scalar occlusion summary of one :class:`ShadowProfile`.

    ``n_occluded`` / ``occluded_fraction`` capture *how much* of the target is
    blocked; ``overhang_5p_nt`` / ``overhang_3p_nt`` capture *how far past the
    physical hybridization footprint* the shadow reaches on each side, in
    resolved nucleotides in duplex-register space. Both halves matter: a
    small, deeply-buried footprint and a shallow-but-wide one can have the
    same ``n_occluded`` and very different overhangs.
    """

    pdb_id: str
    state: str
    probe_radius_A: float
    threshold: float
    n_target: int
    n_measurable: int
    n_occluded: int
    occluded_fraction: float
    footprint_len: int
    overhang_5p_nt: int
    overhang_3p_nt: int

    @property
    def extent_nt(self) -> int:
        """Total occluded nucleotides outside the hybridized footprint."""
        return self.overhang_5p_nt + self.overhang_3p_nt


def summarize(
    profile: ShadowProfile, threshold: float = DEFAULT_SHADOW_THRESHOLD
) -> ParitySummary:
    """Reduce a per-nucleotide :class:`ShadowProfile` to one comparable scalar.

    ``NaN`` entries (unmeasurable SASA_free, see
    :mod:`smt.steric.shadow`) are excluded from every count rather than
    treated as occluded or clear.
    """
    s = profile.shadow
    measurable = ~np.isnan(s)
    occluded = measurable & (s >= threshold)

    hyb_reg = profile.register[profile.hybridized]
    if hyb_reg.size:
        lo, hi = int(hyb_reg.min()), int(hyb_reg.max())
    else:
        # No hybridized span resolved at all (should not happen for these
        # structures, but handled rather than raising, per the "flag rather
        # than guess" rule).
        lo, hi = 0, -1

    outside = occluded & ~profile.hybridized
    outside_reg = profile.register[outside]
    overhang_5p = int((outside_reg < lo).sum())
    overhang_3p = int((outside_reg > hi).sum())

    n_measurable = int(measurable.sum())
    n_occluded = int(occluded.sum())
    return ParitySummary(
        pdb_id=profile.pdb_id,
        state=profile.state,
        probe_radius_A=profile.probe_radius_A,
        threshold=threshold,
        n_target=len(profile),
        n_measurable=n_measurable,
        n_occluded=n_occluded,
        occluded_fraction=(n_occluded / n_measurable) if n_measurable else float("nan"),
        footprint_len=int(profile.hybridized.sum()),
        overhang_5p_nt=overhang_5p,
        overhang_3p_nt=overhang_3p,
    )


@dataclass(frozen=True)
class ParityGap:
    """Difference between two :class:`ParitySummary` objects, ``b`` minus ``a``."""

    a: ParitySummary
    b: ParitySummary
    delta_n_occluded: int
    delta_extent_nt: int
    delta_overhang_5p_nt: int
    delta_overhang_3p_nt: int

    @property
    def b_wider(self) -> bool:
        """True when ``b``'s shadow is strictly wider than ``a``'s."""
        return self.delta_n_occluded > 0


def compare(a: ParitySummary, b: ParitySummary) -> ParityGap:
    """Compare two summaries computed at the SAME probe radius and threshold.

    Raises rather than silently comparing apples to oranges.
    """
    if a.probe_radius_A != b.probe_radius_A:
        raise ValueError(
            f"cannot compare summaries at different probe radii "
            f"({a.probe_radius_A} A vs {b.probe_radius_A} A)"
        )
    if a.threshold != b.threshold:
        raise ValueError(
            f"cannot compare summaries at different thresholds "
            f"({a.threshold} vs {b.threshold})"
        )
    return ParityGap(
        a=a,
        b=b,
        delta_n_occluded=b.n_occluded - a.n_occluded,
        delta_extent_nt=b.extent_nt - a.extent_nt,
        delta_overhang_5p_nt=b.overhang_5p_nt - a.overhang_5p_nt,
        delta_overhang_3p_nt=b.overhang_3p_nt - a.overhang_3p_nt,
    )


@dataclass(frozen=True)
class ParitySweepResult:
    """Result of sweeping ASO length against one reference shadow."""

    reference: ParitySummary
    by_length: dict[int, ParitySummary]
    best_length_nt: int | None
    #: |n_occluded(L) - n_occluded(reference)| for the chosen length, or None
    #: if no length in the sweep produced a measurable ASO shadow.
    best_gap: int | None


def aso_length_for_parity(
    record: ComplexRecord,
    reference: ShadowProfile,
    probe_radius_A: float,
    lengths_nt: Sequence[int] = DEFAULT_ASO_LENGTHS,
    threshold: float = DEFAULT_SHADOW_THRESHOLD,
    n_points: int = 200,
) -> ParitySweepResult:
    """What ASO length L produces a shadow comparable to ``reference``?

    One sentence: sweep ``lengths_nt``, compute each length's bare
    guide:target-duplex shadow with :func:`smt.steric.shadow.aso_shadow`
    (real cryo-EM duplex coordinates, protein deleted -- never an idealized
    A-form helix), summarize it with :func:`summarize`, and return the length
    whose occluded-nucleotide count is closest to the reference's.

    ``reference`` and every ASO shadow are evaluated at the same
    ``probe_radius_A`` and ``threshold`` so the comparison is meaningful; a
    ``ValueError`` from :func:`~smt.steric.shadow.aso_shadow` for a length
    that exceeds the resolved guide:target duplex is caught and that length is
    simply absent from ``by_length`` rather than aborting the sweep.
    """
    ref_summary = summarize(reference, threshold)
    by_length: dict[int, ParitySummary] = {}
    best_length: int | None = None
    best_gap: int | None = None
    for length in lengths_nt:
        try:
            prof = aso_shadow(record, length, probe_radius_A, n_points=n_points)
        except ValueError:
            continue
        summ = summarize(prof, threshold)
        by_length[length] = summ
        gap = abs(summ.n_occluded - ref_summary.n_occluded)
        if best_gap is None or gap < best_gap:
            best_gap = gap
            best_length = length
    return ParitySweepResult(
        reference=ref_summary,
        by_length=by_length,
        best_length_nt=best_length,
        best_gap=best_gap,
    )
