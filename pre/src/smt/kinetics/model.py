"""A8: continuous-time kinetic model of C2B34 backsplicing vs. forward splicing.

The competition this model represents
--------------------------------------
A single nascent SMN2 transcript, once RNA Pol II has transcribed through
exon 4, presents TWO mutually exclusive fates for its exon 4 donor / exon 2a
acceptor pair (see ``smt.kinetics.alu`` for the sequence evidence that a
physical inverted-repeat bridge across this region exists and is not merely
GC-content noise):

  * BACKSPLICE: exon 4's 5' donor joins exon 2a's 3' acceptor, excising and
    circularizing exons 2a-2b-3-4 as the C2B34 circle.
  * FORWARD: exon 4's 5' donor joins exon 5's 3' acceptor (and, symmetrically,
    exon 1's 5' donor joins exon 2a's 3' acceptor) the ordinary, linear way,
    committing that transcript to the canonical splicing pathway.

These two fates draw on the same pool of nascent transcript and the same
donor/acceptor pair, so they are modelled as competing pseudo-first-order
reactions out of one shared state, N.

THE STATE THAT ENCODES THE MECHANISTIC TIMING PROBLEM
-------------------------------------------------------
Exons 2a-4 are transcribed, and the backsplice-vs-forward decision at their
flanking splice sites is made, BEFORE RNA Pol II has even reached intron 7
(where ISS-N1 and the nusinersen binding site live -- see
``smt.thermo.constrained``). An ASO occupying ISS-N1 cannot retroactively
un-make a decision that already happened upstream. This model does not
paper over that: transcripts that have already resolved their upstream
donor/acceptor pairing (whichever way) are moved into an explicit
intermediate state, P ("forward-committed precursor"), BEFORE the exon 7
inclusion/skipping decision is applied. P has no route back to N and no
route to C -- by construction, nothing that happens downstream (at intron 7)
can influence whether a given transcript is backspliced or not, because
that question is already settled by the time a transcript is in state P.
This is what lets the three hypotheses in ``smt.kinetics.hypotheses`` differ
in a mechanistically honest way: H1 and H3 change quantities that act on N
(the pre-decision pool) or on C's own decay, never smuggling an intron-7
signal backwards into the N -> C branch; only H2 is allowed to touch the
N -> C branch directly, and it is required to justify that with a physically
separate upstream binding event, not the ISS-N1 occupancy itself.

States (order used throughout this module and by
``smt.kinetics.hypotheses``): ``STATE_NAMES = ("N", "P", "C", "L_incl",
"L_skip")``

  N       nascent transcript, upstream (exon 2a/exon 4) fate undecided
  P       forward-committed precursor: upstream fate resolved to FORWARD,
          exon 7 inclusion/skipping still undecided
  C       C2B34 circle (backsplice fate)
  L_incl  mature linear mRNA, exon 7 included
  L_skip  mature linear mRNA, exon 7 skipped (the predominant SMN2 product)

Reactions (all pseudo-first-order except transcription, which is zero-order):

  (none) --k_txn--> N                    transcription
  N      --k_back--> C                   backsplicing
  N      --k_fwd-->  P                   forward (canonical) splicing of the
                                          exon 4 donor / exon 2a acceptor pair
  P      --k_incl--> L_incl              exon 7 inclusion
  P      --k_skip--> L_skip              exon 7 skipping
  N      --gamma_N--> (decayed)          nascent transcript turnover
  P      --gamma_P--> (decayed)          precursor turnover
  C      --gamma_C--> (decayed)          circle turnover
  L_incl --gamma_incl--> (decayed)       mRNA turnover
  L_skip --gamma_skip--> (decayed)       mRNA turnover

On "decay" and the no-catalysis constraint
--------------------------------------------
``tests/test_no_catalysis.py`` AST-scans ``smt.steric`` and ``smt.thermo``
(not ``smt.kinetics``) for catalytic/degradation vocabulary, because this
project's steric-occlusion and thermodynamic-accessibility modules must never
model the effector (ASO / dCas13x) itself as catalytic -- dCas13x is
catalytically dead (k_cat = 0 is not even representable there; the concept is
absent) and an ASO recruits no RNase H in this project's scope. That
constraint is about the EFFECTOR. It says nothing about the RNA species
themselves, which do turn over biologically (mRNA decay, nuclear RNA
surveillance) regardless of any effector. The ``gamma_*`` terms in this
module are ordinary first-order RNA decay/turnover of a model species, not
effector-catalyzed cleavage; no effector identity (ASO, dCas13x, guide RNA)
appears anywhere in this module's rate constants or state variables, and
nothing here gives the effector a rate constant of any kind, catalytic or
otherwise. Kept in a separate module (``smt.kinetics``, outside the scanned
directories) precisely so that legitimate RNA-species decay is never
conflated with, or forced to justify itself against, the effector-catalysis
prohibition that governs ``smt.steric``/``smt.thermo``.

Units
-----
Time: arbitrary consistent time unit (e.g. minutes or hours); this module
never fixes one, since no experimental timescale has been measured for this
system (see docs/HANDOFF.md -- no PSI or kinetic data exists yet). All rate
constants (``k_*``, ``gamma_*``) are in units of inverse time; ``k_txn`` is in
units of concentration (or copy number) per unit time.

Concentration / copy number: this module does not distinguish concentration
from absolute copy number -- there is no nuclear volume specified anywhere in
this project, so "state value" is treated as a dimensionless
quantity that means concentration in :func:`simulate_ode` and integer copy
number in :func:`simulate_gillespie`. Consistency between the two only
requires using the same numeric rate constants in both; no unit conversion
is performed or needed for the tests in this module (mean-field vs.
stochastic agreement at high copy number), but converting to a true per-cell
copy number for a real experiment would require an assumed nuclear volume,
which is explicitly NOT supplied here.

ALL RATE CONSTANTS ARE PARAMETRIC
-----------------------------------
:class:`RateConstants` has no default values for any field. There is no
experimental PSI or kinetic data for this system (docs/HANDOFF.md), so
placeholder calibration is forbidden project-wide. Every rate constant must
be supplied explicitly by the caller; omitting one raises ``TypeError``
(the ordinary dataclass behaviour for a missing required field -- see
``tests/test_kinetics.py::test_rate_constants_cannot_be_defaulted``).
:func:`illustrative_rates` provides one example parameter set for
demonstration and is labelled, in its own docstring and return value, as
illustrative -- explicitly NOT a measurement and NOT a fitted value.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np
from scipy.integrate import solve_ivp

#: Canonical state ordering used by every function in this module and by
#: smt.kinetics.hypotheses.
STATE_NAMES: tuple[str, str, str, str, str] = ("N", "P", "C", "L_incl", "L_skip")
N_IDX, P_IDX, C_IDX, LINCL_IDX, LSKIP_IDX = range(5)


@dataclass(frozen=True)
class RateConstants:
    """One instantaneous parameterisation of the reaction network.

    No field has a default value -- every rate constant must be supplied
    explicitly. This is deliberate (see module docstring, "ALL RATE
    CONSTANTS ARE PARAMETRIC"): there is no measured or fitted value for any
    of these in this project, and a silently-defaulted rate constant would
    be exactly the placeholder calibration this project forbids.

    All rates are inverse-time except ``k_txn`` (concentration or copy
    number per unit time). See the module docstring for the reaction each
    symbol drives.
    """

    k_txn: float
    k_back: float
    k_fwd: float
    k_incl: float
    k_skip: float
    gamma_N: float
    gamma_P: float
    gamma_C: float
    gamma_incl: float
    gamma_skip: float

    def __post_init__(self) -> None:
        for name in self.__dataclass_fields__:
            v = getattr(self, name)
            if v < 0:
                raise ValueError(f"rate constant {name}={v!r} must be >= 0")


def illustrative_rates() -> RateConstants:
    """One ILLUSTRATIVE parameter set. NOT a measurement, NOT a fit.

    Used only to demonstrate that the model runs and to give
    :mod:`smt.kinetics.hypotheses` a common baseline to perturb. Two
    qualitative choices are literature-motivated even though the numeric
    values themselves are arbitrary:

    * ``k_skip >> k_incl`` at baseline: SMN2's hallmark defect is that the
      large majority of its transcripts skip exon 7 (Lorson et al. 1999,
      PNAS 96:6307; Lorson & Androphy 2000, Hum Mol Genet 9:259) -- the whole
      reason nusinersen exists is to shift this ratio.
    * ``gamma_C << gamma_incl, gamma_skip``: circRNAs lack free 5'/3' ends
      and are broadly reported to be resistant to exonucleolytic decay
      relative to their linear counterparts (Jeck et al. 2013, RNA
      19:141; Suzuki & Tsukahara 2014, Int J Mol Sci 15:9331) -- this is the
      textbook basis of H3 in ``smt.kinetics.hypotheses``.

    Every numeric value is otherwise an arbitrary round number chosen only so
    the model reaches a nonzero, finite steady state; none should be read as
    a rate this project has measured.
    """
    return RateConstants(
        k_txn=1.0,  # ILLUSTRATIVE: arbitrary unit rate, not measured
        k_back=0.01,  # ILLUSTRATIVE: backsplicing assumed rare at baseline
        k_fwd=1.0,  # ILLUSTRATIVE: forward splicing assumed to dominate at baseline
        k_incl=0.1,  # ILLUSTRATIVE: baseline low inclusion (SMN2 hallmark direction, not a measured value)
        k_skip=0.9,  # ILLUSTRATIVE: baseline high skipping (SMN2 hallmark direction, not a measured value)
        gamma_N=0.05,  # ILLUSTRATIVE
        gamma_P=0.05,  # ILLUSTRATIVE
        gamma_C=0.01,  # ILLUSTRATIVE: small (circRNA stability direction, not a measured half-life)
        gamma_incl=1.0,  # ILLUSTRATIVE
        gamma_skip=1.0,  # ILLUSTRATIVE
    )


# ---------------------------------------------------------------------------
# deterministic (ODE) model
# ---------------------------------------------------------------------------


def flux(state: np.ndarray, rates: RateConstants) -> np.ndarray:
    """Instantaneous d(state)/dt for one set of (possibly time-sliced) rates.

    Pure function of the current state and a single ``RateConstants``; has no
    knowledge of time itself. ``smt.kinetics.hypotheses`` reuses this
    function with a DIFFERENT ``RateConstants`` recomputed at each timestep
    to implement time- and dose-dependent rates, so the reaction network
    itself (this function) is defined in exactly one place.
    """
    N, P, C, L_incl, L_skip = state
    r = rates
    dN = r.k_txn - (r.k_back + r.k_fwd + r.gamma_N) * N
    dP = r.k_fwd * N - (r.k_incl + r.k_skip + r.gamma_P) * P
    dC = r.k_back * N - r.gamma_C * C
    dL_incl = r.k_incl * P - r.gamma_incl * L_incl
    dL_skip = r.k_skip * P - r.gamma_skip * L_skip
    return np.array([dN, dP, dC, dL_incl, dL_skip])


def steady_state(rates: RateConstants) -> np.ndarray:
    """Exact analytic steady state of the (linear, constant-rate) ODE system.

    Solving dN/dt=dP/dt=dC/dt=dL_incl/dt=dL_skip/dt=0 in sequence (the system
    is a linear cascade, so each equation only needs the previous state).
    Useful as an initial condition (the state "before treatment", i.e. the
    system's own equilibrium under baseline rates) and as an exact reference
    value for tests.
    """
    r = rates
    N_ss = r.k_txn / (r.k_back + r.k_fwd + r.gamma_N)
    P_ss = r.k_fwd * N_ss / (r.k_incl + r.k_skip + r.gamma_P)
    C_ss = r.k_back * N_ss / r.gamma_C
    L_incl_ss = r.k_incl * P_ss / r.gamma_incl
    L_skip_ss = r.k_skip * P_ss / r.gamma_skip
    return np.array([N_ss, P_ss, C_ss, L_incl_ss, L_skip_ss])


def simulate_ode(
    rates: RateConstants,
    y0: Optional[np.ndarray] = None,
    t_span: tuple[float, float] = (0.0, 100.0),
    t_eval: Optional[np.ndarray] = None,
    rhs: Optional[Callable[[float, np.ndarray], np.ndarray]] = None,
    **solve_ivp_kwargs,
):
    """Integrate the ODE system with ``scipy.integrate.solve_ivp``.

    Parameters
    ----------
    rates:
        Constant-rate parameterisation. Ignored if ``rhs`` is supplied
        (see below).
    y0:
        Initial state, ``[N, P, C, L_incl, L_skip]``. Defaults to the
        analytic steady state of ``rates`` (i.e. "the system at its own
        baseline equilibrium"), which is a principled default, not a
        placeholder value -- it falls out exactly from ``rates`` itself.
    rhs:
        Optional override of the right-hand side, ``f(t, y) -> dy/dt``, for
        callers (``smt.kinetics.hypotheses``) that need time/dose-dependent
        rates. If omitted, the constant-rate system built from ``rates`` is
        used.
    **solve_ivp_kwargs:
        Forwarded to ``scipy.integrate.solve_ivp`` (e.g. ``method``,
        ``max_step``, ``rtol``, ``atol``).

    Returns
    -------
    The ``scipy.integrate.OdeResult`` from ``solve_ivp`` (``.t``, ``.y``).
    """
    if y0 is None:
        y0 = steady_state(rates)
    f = rhs if rhs is not None else (lambda t, y: flux(y, rates))
    return solve_ivp(f, t_span, y0, t_eval=t_eval, **solve_ivp_kwargs)


# ---------------------------------------------------------------------------
# stochastic (Gillespie) model -- for the low-copy-number regime
# ---------------------------------------------------------------------------

#: (state_index_changes, propensity_fn) for each of the 10 reactions in the
#: module docstring's reaction list, in the same order they are listed there.
#: ``delta`` is added to the state vector when the reaction fires.
_ReactionSpec = tuple[np.ndarray, Callable[[np.ndarray, RateConstants], float]]


def _reactions() -> list[_ReactionSpec]:
    def d(**kw) -> np.ndarray:
        v = np.zeros(5)
        for name, val in kw.items():
            v[STATE_NAMES.index(name)] = val
        return v

    return [
        (d(N=+1), lambda s, r: r.k_txn),  # transcription (zero-order)
        (d(N=-1, C=+1), lambda s, r: r.k_back * s[N_IDX]),  # backsplicing
        (d(N=-1, P=+1), lambda s, r: r.k_fwd * s[N_IDX]),  # forward splicing
        (d(N=-1), lambda s, r: r.gamma_N * s[N_IDX]),  # N decay
        (d(P=-1, L_incl=+1), lambda s, r: r.k_incl * s[P_IDX]),  # inclusion
        (d(P=-1, L_skip=+1), lambda s, r: r.k_skip * s[P_IDX]),  # skipping
        (d(P=-1), lambda s, r: r.gamma_P * s[P_IDX]),  # P decay
        (d(C=-1), lambda s, r: r.gamma_C * s[C_IDX]),  # C decay
        (d(L_incl=-1), lambda s, r: r.gamma_incl * s[LINCL_IDX]),  # L_incl decay
        (d(L_skip=-1), lambda s, r: r.gamma_skip * s[LSKIP_IDX]),  # L_skip decay
    ]


@dataclass(frozen=True)
class GillespieResult:
    """One stochastic trajectory: jump times and the state after each jump.

    ``t[0] == t_start`` with ``y[0]`` the initial condition; every subsequent
    row is the state immediately after the reaction that fired at ``t[k]``.
    """

    t: np.ndarray  # shape (n_events + 1,)
    y: np.ndarray  # shape (n_events + 1, 5)

    def sample_at(self, times: np.ndarray) -> np.ndarray:
        """Right-continuous step-function evaluation at arbitrary ``times``.

        For each query time, returns the state after the most recent event
        at or before it (i.e. the state that trajectory was actually in).
        """
        idx = np.searchsorted(self.t, times, side="right") - 1
        idx = np.clip(idx, 0, len(self.t) - 1)
        return self.y[idx]


def simulate_gillespie(
    rates: RateConstants,
    y0: np.ndarray,
    t_max: float,
    rng: Optional[np.random.Generator] = None,
    max_events: int = 2_000_000,
) -> GillespieResult:
    """Exact stochastic simulation (Gillespie 1977, J. Phys. Chem. 81:2340).

    Direct method: draw the waiting time to the next reaction from an
    exponential with rate equal to the total propensity, then choose which
    reaction fires with probability proportional to its own propensity.
    Copy numbers in ``y0`` are treated as integers (rounded if not already);
    this is the low-copy-number regime the mean-field ODE approximates.

    Parameters
    ----------
    y0:
        Initial copy numbers, ``[N, P, C, L_incl, L_skip]``. Unlike
        :func:`simulate_ode`, there is no principled default here (a
        "steady state" copy number is generally not an integer), so it must
        be supplied.
    max_events:
        Safety cap so a mis-parameterised system (e.g. all decay rates zero
        with nonzero transcription) cannot run forever; raises
        ``RuntimeError`` if hit.
    """
    if rng is None:
        rng = np.random.default_rng()
    reactions = _reactions()
    state = np.array(y0, dtype=float).round()
    t = 0.0
    ts = [t]
    ys = [state.copy()]
    n_events = 0
    while t < t_max:
        props = np.array([p(state, rates) for _, p in reactions])
        a0 = props.sum()
        if a0 <= 0:
            break  # absorbing: nothing can happen from here (e.g. all-zero, all decay)
        tau = rng.exponential(1.0 / a0)
        t += tau
        if t > t_max:
            break
        u = rng.random() * a0
        cum = np.cumsum(props)
        idx = int(np.searchsorted(cum, u))
        delta, _ = reactions[idx]
        state = state + delta
        ts.append(t)
        ys.append(state.copy())
        n_events += 1
        if n_events > max_events:
            raise RuntimeError(
                f"simulate_gillespie exceeded max_events={max_events}; "
                "check for a mis-parameterised (e.g. non-decaying, growing) system"
            )
    return GillespieResult(np.array(ts), np.array(ys))
