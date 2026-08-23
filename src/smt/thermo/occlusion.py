"""A7: generalized, signed occlusion readout for effector footprints.

Retargets ``smt.thermo.constrained.a_5ss`` -- built for SMN2 exon 7 -- at
arbitrary named splice sites, with an explicit DIRECTION, so the same
thermodynamic machinery scores both regimes this project now covers:

  Direction.FREE     SMN2 exon 7 5'ss: an effector binds a distal silencer
                      element (ISS-N1 / the GC core) and the therapeutic goal
                      is to INCREASE accessibility of the 5'ss so U1 snRNP can
                      bind and the exon is INCLUDED.
  Direction.OCCLUDE   STMN2 / UNC13A cryptic splice sites: the therapeutic
                      goal is the opposite -- DECREASE accessibility of the
                      cryptic site so it is NOT recognised and the cryptic
                      exon is EXCLUDED.

No new thermodynamics is added here. Every number comes from
``smt.thermo.constrained.a_5ss`` / ``occlusion_delta``, which already accepts
an arbitrary ``ss5_feature`` name -- this module only adds the sign
convention and documents why one guard in that function must NOT be routed
around.

The signed figure of merit
---------------------------
    FOM = direction.sign * delta_A

where ``delta_A = A_site(effector bound) - A_site(baseline)`` exactly as
``OcclusionDelta.delta_A`` defines it. ``direction.sign`` is +1 for FREE
(wanting delta_A > 0) and -1 for OCCLUDE (wanting delta_A < 0). FOM is
POSITIVE, in BOTH regimes, exactly when the effector does the therapeutically
desired thing at that site, and its magnitude is directly comparable to
SMN2's existing delta_A numbers (same units: a probability difference on the
same kind of joint-unpaired readout). This is the single number
``cryptic_sweep`` ranks footprints by.

THE SUBTLE PART -- why a footprint occluding the site directly is refused
here too, in both regimes, and what that means for real ASO designs
---------------------------------------------------------------------------
``a_5ss`` raises ``ValueError`` if the footprint overlaps the readout window
AT ALL. That guard was written for the FREE regime ("forcing the readout
window itself unpaired makes A trivially 1"), but it is not a FREE-regime
quirk -- it is a hard boundary of what the hard-constraint model can
represent, and it applies with equal force, for a different reason, in the
OCCLUDE regime.

``fc.hc_add_up(i, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)`` forces base ``i`` to be
UNPAIRED and says nothing about what (if anything) is bound to it. That is
the right representation of "this intronic silencer base cannot form an
intramolecular pair any more" (mechanism 1 below). It is the WRONG
representation of "this exact base is occupied by an intermolecular RNA:ASO
duplex" (mechanism 2): forcing the splice site itself unpaired via
``hc_add_up`` says the site is perfectly single-stranded, which by the letter
of ``A = P(unpaired)`` reads as perfectly ACCESSIBLE -- backwards for a site
that is supposed to be sterically blocked by a drug sitting on top of it.

Two physically distinct occlusion mechanisms, only one of which this
framework models:

  1. INDIRECT / STRUCTURAL (what this module scores). The effector binds a
     DIFFERENT sequence element -- e.g. it locks one arm of a competing
     hairpin in place, or removes a single-stranded stretch a distal
     structure needed to stay open -- and the change propagates through the
     folding ensemble to the splice site. This is exactly the SMN2 mechanism
     run with the opposite desired sign: ``a_5ss``'s hard-constraint-plus-
     partition-function-ratio computation is agnostic to which way the
     ensemble moves, so it is equally valid whether the propagated effect
     frees or sequesters the readout window, AS LONG AS the footprint and the
     readout window are physically disjoint (the constraint acts only on the
     footprint; the readout is measured on an untouched window, so whatever
     the ensemble does there is a real consequence, not a constraint
     artifact). This is the actual published STMN2 strategy: Baughn et al.
     2023 (Science 379:1140) mask the TDP-43 GU-repeat region immediately 3'
     of the cryptic acceptor -- not the AG dinucleotide itself -- to restore
     correct splicing, i.e. a flanking, not overlapping, footprint.

  2. DIRECT / COMPETITIVE OCCUPANCY (NOT modelled here; a documented
     limitation, not an oversight). A splice-switching ASO that anneals
     directly across the splice site does not make that site "unpaired" in
     any sense ``hc_add_up`` can express -- it makes it double-stranded with
     the drug. Modelling this correctly means asking whether the ASO:site
     duplex free energy (``smt.thermo.hybridization``) outcompetes the
     intramolecular fold at that window (``smt.thermo.accessibility``'s
     ``dG_unfold``), i.e. the ``ddG_net`` construction already used in
     ``smt.sieve.scan`` for a different purpose. That is a materially
     different calculation from a hard-constraint partition-function ratio,
     and building a second, competitive-hybridization occlusion estimator
     was out of scope for this module. Rather than fake a number by routing
     around ``a_5ss``'s guard (e.g. by constraining the footprint bases
     PAIRED instead of unpaired, which has no principled interpretation
     either -- paired to WHAT is exactly the question mechanism 2 needs to
     answer and mechanism 1's machinery cannot), this module lets the
     ``ValueError`` propagate, and ``cryptic_sweep`` records any such
     footprint as skipped with ``direct_overlap=True`` rather than scoring
     it. This is the same class of documented boundary the HANDOFF already
     records for ISS-N1 (a purely-structural model that would need explicit
     hnRNP-competitive terms to go further) -- restated here for the
     opposite regime rather than re-derived.

Practical consequence for real ASO design, stated plainly: this framework
predicts and ranks FLANKING occlusion footprints (mechanism 1) for STMN2 and
UNC13A. It does not, and does not claim to, predict the effect of an ASO
tiled directly across the exact splice-site dinucleotides (mechanism 2) --
those footprints are visible in a tiling sweep only as "excluded, would need
a competitive-hybridization model."
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Optional

from smt.seq.construct import Construct
from smt.seq.reference import Interval
from smt.thermo.constrained import (
    SS5_FEATURE,
    A5ssResult,
    Estimator,
    OcclusionDelta,
    a_5ss,
    occlusion_delta,
)


class Direction(str, Enum):
    """Which way the effector is supposed to move accessibility."""

    FREE = "free"  # increase accessibility -- SMN2 exon7 5'ss regime
    OCCLUDE = "occlude"  # decrease accessibility -- cryptic splice site regime

    @property
    def sign(self) -> int:
        return 1 if self is Direction.FREE else -1


@dataclass(frozen=True)
class SignedOcclusionResult:
    """One footprint scored against one named site, signed by ``direction``.

    ``fom`` (figure of merit) is positive exactly when the effector does the
    therapeutically desired thing for that direction. See module docstring.
    """

    direction: Direction
    delta: OcclusionDelta

    @property
    def fom(self) -> float:
        return self.direction.sign * self.delta.delta_A

    @property
    def desired(self) -> bool:
        """True iff the effector moved accessibility the therapeutically
        desired way (FOM > 0), i.e. the sign, not the magnitude."""
        return self.fom > 0

    # -- passthroughs, for convenience when reporting -----------------------

    @property
    def delta_A(self) -> float:
        return self.delta.delta_A

    @property
    def fold_change(self) -> float:
        return self.delta.fold_change

    @property
    def baseline(self) -> A5ssResult:
        return self.delta.baseline

    @property
    def blocked(self) -> A5ssResult:
        return self.delta.blocked


def score(
    footprint: Interval,
    direction: Direction | str,
    site_feature: str = SS5_FEATURE,
    construct: Optional[Construct] = None,
    estimator: Estimator | str = Estimator.ENSEMBLE,
    max_bp_span: int = 0,
    baseline: Optional[A5ssResult] = None,
) -> SignedOcclusionResult:
    """Score one footprint against one named splice site, signed by direction.

    Thin wrapper around ``smt.thermo.constrained.occlusion_delta`` -- see
    module docstring for why the identical computation is correct for both
    directions, and why a footprint overlapping ``site_feature`` is refused
    (the underlying ``ValueError`` from ``a_5ss`` propagates unchanged; this
    function does not catch or reinterpret it).

    ``baseline`` is an optional precomputed unoccluded result, passed straight
    through to ``occlusion_delta`` (which validates it). It exists so a sweep
    over many footprints on one construct does not recompute the identical
    baseline partition function per footprint -- that is 50% of the cost, not
    a rounding error. See ``smt.validate.cryptic_sweep.build_baselines``.
    """
    d = occlusion_delta(
        footprint, construct, estimator, site_feature, max_bp_span, baseline=baseline
    )
    return SignedOcclusionResult(Direction(direction), d)


def overlaps(footprint: Interval, site: Interval) -> bool:
    """True if ``footprint`` and ``site`` (both REFERENCE space) share any
    position. Mirrors the check inside ``a_5ss`` so callers (e.g.
    ``cryptic_sweep``) can filter/flag before calling ``score`` rather than
    catching its ``ValueError``."""
    return not (footprint.end < site.start or footprint.start > site.end)


def mutate(construct: Construct, ref_pos: int, alt_base: str) -> Construct:
    """Return a copy of ``construct`` with one base substituted.

    ``ref_pos`` is REFERENCE space (same convention as a footprint
    ``Interval`` -- e.g. a gene-local coordinate from ``smt.seq.targets``),
    mapped to construct space internally via ``Construct.to_construct``, the
    same mapping ``Construct.feature``/``map_interval`` use. Used for scoring
    point variants (e.g. UNC13A rs12608932 / rs12973192) through exactly the
    same readout used for footprint scoring -- same window, same features,
    only the sequence differs by one base. DNA alphabet in, DNA alphabet
    stored (``Construct.rna()`` T->U's at fold time, same as everywhere else
    in this codebase).
    """
    con_pos = construct.to_construct(ref_pos)
    if con_pos is None:
        raise ValueError(f"reference position {ref_pos} outside construct {construct.name} (window {construct.window})")
    seq = construct.sequence
    new_seq = seq[: con_pos - 1] + alt_base.upper() + seq[con_pos:]
    return Construct(
        f"{construct.name}_variant_{ref_pos}{alt_base.upper()}",
        construct.window,
        construct.deletions,
        new_seq,
        construct.features,
    )


@dataclass(frozen=True)
class VariantEffect:
    """Baseline (no effector) accessibility of one named site, ref vs. variant
    background. Answers the literature's actual question for a risk SNP: does
    the substitution alone shift the intrinsic splicing decision, independent
    of any drug."""

    site_feature: str
    max_bp_span: int
    ref: A5ssResult
    alt: A5ssResult

    @property
    def delta_A(self) -> float:
        """A(variant) - A(reference). Positive == variant makes the site MORE
        accessible (for a cryptic site, more accessible baseline is the
        pro-pathogenic direction: it is easier for the spliceosome to find
        without needing TDP-43 loss to help)."""
        return self.alt.A_5ss - self.ref.A_5ss

    @property
    def fold_change(self) -> float:
        return self.alt.A_5ss / self.ref.A_5ss if self.ref.A_5ss > 0 else float("inf")


def variant_effect(
    ref_construct: Construct,
    ref_pos: int,
    alt_base: str,
    site_feature: str,
    estimator: Estimator | str = Estimator.ENSEMBLE,
    max_bp_span: int = 0,
) -> VariantEffect:
    """Baseline site accessibility, reference vs. single-base variant.

    ``ref_pos`` is REFERENCE space (see ``mutate``). No footprint is involved
    -- this scores the intrinsic effect of the sequence change itself on
    ``site_feature``'s accessibility, using ``a_5ss`` with ``footprint=None``
    on each background in turn. Reuses the exact same estimator as effector
    scoring; only the input sequence differs.
    """
    alt_construct = mutate(ref_construct, ref_pos, alt_base)
    ref_result = a_5ss(ref_construct, None, estimator, site_feature, max_bp_span)
    alt_result = a_5ss(alt_construct, None, estimator, site_feature, max_bp_span)
    return VariantEffect(site_feature, max_bp_span, ref_result, alt_result)
