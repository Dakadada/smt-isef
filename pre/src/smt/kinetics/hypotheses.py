"""A8: three competing mechanisms for nusinersen-induced C2B34 overexpression,
encoded as three parameterisations of the shared model in
:mod:`smt.kinetics.model`, plus a discriminator over their predicted
signatures.

Why circle SYNTHESIS is the thing in dispute, not exon 7 inclusion
--------------------------------------------------------------------
Nusinersen occupying ISS-N1 directly frees the exon 7 5' splice site for U1
(this is the established, already-modelled mechanism in
``smt.thermo.constrained``, not something this module re-litigates). All
three hypotheses below apply that SAME direct-occlusion boost to ``k_incl``,
identically, as a function of dose -- see :func:`_incl_rate`. What the three
hypotheses actually disagree about is why the C2B34 CIRCLE goes up by as
much as ~1000-fold, since (per the model's own construction in
``smt.kinetics.model``) nothing at ISS-N1 can directly touch the N -> C
branch: that branch is decided before intron 7 is even transcribed.

  H1  SMN FEEDBACK. Restored SMN protein (translated only from the
      exon-7-included isoform) feeds back on snRNP biogenesis and shifts
      total spliceosome availability, including the N -> C branching ratio.
      Because this must run through translation and snRNP re-import before
      it can touch splicing, it is implemented with an explicit intermediate
      state S (a proxy for available SMN-dependent spliceosome capacity) with
      its own production/decay timescale -- this is what MANUFACTURES A LAG
      that the other two hypotheses do not have. H1's effect on k_back
      depends on the protein level S(t), not on ASO occupancy directly.

  H2  DIRECT UPSTREAM EFFECT. At high concentration, free (unbound) ASO
      also occupies a site in the flanking introns near the inverted-repeat
      arms found by ``smt.kinetics.alu`` (see that module's docstring and
      report for the actual sequence evidence) and directly competes with
      or stabilises the Alu-Alu duplex, changing k_back immediately (no
      protein-synthesis delay) but only once occupancy at this SEPARATE,
      lower-affinity site becomes appreciable -- i.e. only at high dose.

  H3  STABILITY, NOT SYNTHESIS. Neither the backsplicing branch ratio nor
      any splicing rate changes at all. The circle's own decay rate,
      gamma_C, decreases with dose. (The prompt's phrasing, "circles
      accumulate because linear decay changes", is ambiguous on a literal
      reading -- decreasing LINEAR mRNA decay would, if anything, shrink the
      circle:linear ratio, not grow it. The reading implemented here --
      decay of the CIRCLE itself changes, differentially from the linear
      species -- is the only one consistent with "flat synthesis flux"
      actually producing circle accumulation, and is also the textbook
      circRNA-stability argument (circles resist exonucleolytic decay
      because they lack free ends). This is a disclosed INTERPRETIVE
      assumption, not a verbatim restatement of the prompt.)

All parameters below are symbolic/illustrative (see
``smt.kinetics.model.illustrative_rates`` docstring for what "illustrative"
means in this project) -- there is no dose-response or time-course data for
C2B34 yet (docs/HANDOFF.md), so none of EC50, gain, or lag-timescale
parameters are fitted to anything. They exist to be swept, and the
discriminator below is explicitly about how well-separated the three
hypotheses' PREDICTED SHAPES are, which does not require the absolute
parameter values to be correct -- only the qualitative dose/time dependence
each hypothesis is mechanistically committed to.

Units: dose is in arbitrary concentration units (never fixed to nM/uM
because no real potency has been measured for either the ISS-N1 site or the
hypothetical H2 upstream site). Time is in the same arbitrary unit as
``smt.kinetics.model``.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from smt.kinetics.model import LINCL_IDX, N_IDX, RateConstants, flux, steady_state

HYPOTHESES: tuple[str, str, str] = ("H1_feedback", "H2_direct", "H3_stability")


def occupancy(dose: np.ndarray | float, ec50: float) -> np.ndarray | float:
    """Simple Hill (n=1) receptor-occupancy curve, dose/(dose+EC50).

    ASSUMPTION: a minimal saturating binding curve, not a fitted
    pharmacological model -- there is no measured Hill coefficient or
    cooperativity for either binding site in this project.
    """
    return dose / (dose + ec50)


@dataclass(frozen=True)
class HypothesisParams:
    """Shared + hypothesis-specific parameters. No field has a default value.

    baseline:
        The untreated (dose=0) RateConstants, e.g.
        ``smt.kinetics.model.illustrative_rates()`` -- but callers must pass
        one explicitly (see module-level docstring: nothing here is silently
        defaulted).
    ec50_intron7:
        Occupancy EC50 (dose units) for the established, high-affinity
        ISS-N1 site. Drives ``k_incl`` identically in all three hypotheses.
    g_incl:
        Maximum fold-change gain applied to k_incl at saturating ISS-N1
        occupancy.
    ec50_upstream:
        H2-only: occupancy EC50 for the hypothesised separate, lower-affinity
        upstream (Alu-region) site. Mechanistically required to be >>
        ec50_intron7 ("only at high concentration" per the H2 description).
    g_back_direct:
        H2-only: maximum fold-change gain applied to k_back at saturating
        upstream occupancy.
    g_back_feedback:
        H1-only: maximum fold-change gain applied to k_back at saturating
        SMN-feedback signal S.
    alpha_translate:
        H1-only: rate at which L_incl produces the feedback signal S
        (a lumped proxy for translation + nuclear re-import of new snRNP
        capacity, not a mechanistic translation-rate measurement).
    delta_S:
        H1-only: decay/turnover rate of the feedback signal S. Together with
        ``alpha_translate`` this sets the LAG timescale ~1/delta_S (roughly)
        that is H1's defining, testable signature.
    s50:
        H1-only: half-saturation of the feedback signal's effect on k_back.
    g_decay_stability:
        H3-only: maximum fractional REDUCTION of gamma_C at saturating
        ISS-N1 occupancy (uses the same ec50_intron7 curve as k_incl, on the
        assumption that whatever stabilises the circle is downstream of the
        same ASO-binding event, not a separate site -- an explicit
        assumption, since H3 does not claim a distinct binding site the way
        H2 does).
    """

    baseline: RateConstants
    ec50_intron7: float
    g_incl: float
    ec50_upstream: float
    g_back_direct: float
    g_back_feedback: float
    alpha_translate: float
    delta_S: float
    s50: float
    g_decay_stability: float

    def __post_init__(self) -> None:
        if self.ec50_upstream <= self.ec50_intron7:
            raise ValueError(
                "H2 requires ec50_upstream > ec50_intron7 (the upstream site is "
                "hypothesised to be lower-affinity, engaged only at high dose); "
                f"got ec50_upstream={self.ec50_upstream} <= ec50_intron7={self.ec50_intron7}"
            )


def illustrative_hypothesis_params() -> HypothesisParams:
    """One ILLUSTRATIVE parameter set for all three hypotheses. Not a fit.

    ``baseline`` overrides ``smt.kinetics.model.illustrative_rates()``'s own
    gamma_C (0.01 -> 0.1): a faster baseline circle turnover than that
    module's own illustrative choice, picked purely so that H2's rate-change
    relaxation timescale, H3's decay-driven relaxation timescale, and H1's
    feedback-driven lag timescale (set by ``delta_S`` below) sit far enough
    apart to be visually and numerically distinguishable in the demonstration
    below. Both values are illustrative; neither is measured.
    """
    from dataclasses import replace as _replace

    from smt.kinetics.model import illustrative_rates

    return HypothesisParams(
        baseline=_replace(illustrative_rates(), gamma_C=0.1),  # ILLUSTRATIVE override, see docstring
        ec50_intron7=10.0,  # ILLUSTRATIVE dose unit
        g_incl=8.0,  # ILLUSTRATIVE: up to 9x k_incl at saturating ISS-N1 occupancy
        ec50_upstream=200.0,  # ILLUSTRATIVE: 20x higher than ec50_intron7 ("high concentration" per H2)
        g_back_direct=15.0,  # ILLUSTRATIVE: H2 can push k_back up ~16x at saturating upstream occupancy
        g_back_feedback=15.0,  # ILLUSTRATIVE: H1 can push k_back up ~16x at saturating feedback signal
        alpha_translate=0.5,  # ILLUSTRATIVE: sets H1's lag timescale together with delta_S
        delta_S=0.02,  # ILLUSTRATIVE: slow turnover -> the feedback signal, and hence k_back, moves slowly
        s50=40.0,  # ILLUSTRATIVE: set well above the untreated feedback-signal level (~2.2, see
        # _h1_self_consistent_baseline) so the pre-treatment (disease-state, low-SMN) baseline
        # sits on the shallow part of the feedback curve, as SMA biology motivates -- not fitted.
        g_decay_stability=0.98,  # ILLUSTRATIVE: gamma_C can shrink up to 98% at saturating occupancy
    )


# ---------------------------------------------------------------------------
# shared component: the direct, dose-tracking, no-lag exon-7 inclusion boost
# ---------------------------------------------------------------------------


def _incl_rate(params: HypothesisParams, dose: float) -> float:
    """k_incl(dose): identical across all three hypotheses (see module docstring)."""
    theta = occupancy(dose, params.ec50_intron7)
    return params.baseline.k_incl * (1.0 + params.g_incl * theta)


# ---------------------------------------------------------------------------
# H2: direct upstream effect -- constant rates after t=0, no lag
# ---------------------------------------------------------------------------


def _rates_h2(params: HypothesisParams, dose: float) -> RateConstants:
    theta_up = occupancy(dose, params.ec50_upstream)
    return replace(
        params.baseline,
        k_incl=_incl_rate(params, dose),
        k_back=params.baseline.k_back * (1.0 + params.g_back_direct * theta_up),
    )


# ---------------------------------------------------------------------------
# H3: stability, not synthesis -- k_back untouched, gamma_C decreases
# ---------------------------------------------------------------------------


def _rates_h3(params: HypothesisParams, dose: float) -> RateConstants:
    theta = occupancy(dose, params.ec50_intron7)
    return replace(
        params.baseline,
        k_incl=_incl_rate(params, dose),
        gamma_C=params.baseline.gamma_C * (1.0 - params.g_decay_stability * theta),
    )


# ---------------------------------------------------------------------------
# H1: SMN feedback -- extended state [N, P, C, L_incl, L_skip, S], has a lag
# ---------------------------------------------------------------------------


def _h1_rhs(t: float, y: np.ndarray, params: HypothesisParams, dose: float) -> np.ndarray:
    state, S = y[:5], y[5]
    k_incl = _incl_rate(params, dose)
    feedback = S / (S + params.s50)
    k_back = params.baseline.k_back * (1.0 + params.g_back_feedback * feedback)
    rates = replace(params.baseline, k_incl=k_incl, k_back=k_back)
    dstate = flux(state, rates)
    dS = params.alpha_translate * state[LINCL_IDX] - params.delta_S * S
    return np.concatenate([dstate, [dS]])


def _h1_fixed_point(params: HypothesisParams, dose: float) -> np.ndarray:
    """The self-consistent (N,P,C,L_incl,L_skip,S) fixed point of H1 at a
    given, CONSTANT dose (dose=0 gives the pre-treatment baseline).

    H1's feedback is a real, ALWAYS-ON physiological loop, not something that
    switches on only with drug -- so even the untreated state must be the
    self-consistent equilibrium of the whole feedback loop, not
    ``smt.kinetics.model.steady_state(params.baseline)`` (which holds k_back
    fixed at its unmodulated value and would silently ignore whatever
    feedback the pre-treatment SMN level already exerts). At any FIXED dose
    (hence fixed k_incl), solved with a 1-D root find on the feedback signal
    S: for a fixed S the rest of the state has the model's ordinary
    closed-form steady state (``smt.kinetics.model.steady_state``), and S
    itself is determined by that same state through
    ``alpha_translate``/``delta_S`` -- a fixed point, not a free parameter.
    This also gives an exact, cheap (root-find, not long-time ODE
    integration) new-steady-state value for :func:`dose_response`, which
    needs one evaluation per dose.

    Returns the 6-vector ``[N, P, C, L_incl, L_skip, S]``.
    """
    k_incl = _incl_rate(params, dose)

    def residual(S: float) -> float:
        feedback = S / (S + params.s50)
        k_back = params.baseline.k_back * (1.0 + params.g_back_feedback * feedback)
        ss = steady_state(replace(params.baseline, k_incl=k_incl, k_back=k_back))
        return params.alpha_translate * ss[LINCL_IDX] / params.delta_S - S

    # Bracket: S=0 has residual >= 0 (translation term is non-negative); a
    # generous upper bound comes from the maximum possible L_incl (achieved
    # at maximal k_back, i.e. maximal feedback) -- use the no-feedback
    # L_incl_ss scaled up by (1 + g_back_feedback) as a safe multiple of the
    # true bound, then widen further if that is not enough.
    naive_S = params.alpha_translate * steady_state(replace(params.baseline, k_incl=k_incl))[LINCL_IDX] / params.delta_S
    hi = max(10.0 * (1.0 + params.g_back_feedback) * max(naive_S, 1e-6), 1.0)
    while residual(hi) > 0:
        hi *= 10.0
    S0 = brentq(residual, 0.0, hi)
    feedback0 = S0 / (S0 + params.s50)
    k_back0 = params.baseline.k_back * (1.0 + params.g_back_feedback * feedback0)
    ss0 = steady_state(replace(params.baseline, k_incl=k_incl, k_back=k_back0))
    return np.concatenate([ss0, [S0]])


def simulate(
    hyp: str,
    params: HypothesisParams,
    dose: float,
    t_span: tuple[float, float],
    t_eval: np.ndarray,
) -> dict[str, np.ndarray]:
    """Simulate one hypothesis from the pre-treatment baseline steady state.

    ``dose`` is applied as a step function at t=0 (held constant for
    t_span); this models a single bolus/steady dosing regimen, not a
    pharmacokinetic decay of free drug -- an explicit simplification, since
    no PK data for this system is being modelled here.

    Returns a dict with keys 'N','P','C','L_incl','L_skip' (and 'S' for H1),
    each an array aligned with ``t_eval``, plus 'ratio' = C / (L_incl+L_skip)
    and 'synthesis_flux' = instantaneous k_back(t) * N(t) (a diagnostic, not
    something assumed measurable in a real experiment -- see
    docs/c2b34-discriminating-experiment.md for what IS measurable).
    """
    if hyp not in HYPOTHESES:
        raise ValueError(f"unknown hypothesis {hyp!r}; must be one of {HYPOTHESES}")

    if hyp == "H1_feedback":
        y0 = _h1_fixed_point(params, 0.0)  # pre-treatment (dose=0) baseline
        sol = solve_ivp(
            lambda t, y: _h1_rhs(t, y, params, dose), t_span, y0, t_eval=t_eval, rtol=1e-8, atol=1e-10
        )
        state = sol.y[:5]
        S = sol.y[5]
        k_back_t = params.baseline.k_back * (1.0 + params.g_back_feedback * (S / (S + params.s50)))
        out = {name: state[i] for i, name in enumerate(("N", "P", "C", "L_incl", "L_skip"))}
        out["S"] = S
        out["synthesis_flux"] = k_back_t * state[N_IDX]
    else:
        y0_base = steady_state(params.baseline)
        rates_fn = _rates_h2 if hyp == "H2_direct" else _rates_h3
        rates = rates_fn(params, dose)
        sol = solve_ivp(lambda t, y: flux(y, rates), t_span, y0_base, t_eval=t_eval, rtol=1e-8, atol=1e-10)
        state = sol.y
        out = {name: state[i] for i, name in enumerate(("N", "P", "C", "L_incl", "L_skip"))}
        out["synthesis_flux"] = rates.k_back * state[N_IDX]

    out["ratio"] = out["C"] / (out["L_incl"] + out["L_skip"])
    return out


# ---------------------------------------------------------------------------
# predicted signatures
# ---------------------------------------------------------------------------


def dose_response(hyp: str, params: HypothesisParams, doses: np.ndarray, t_late: float) -> np.ndarray:
    """Circle:linear ratio at the NEW STEADY STATE reached after treatment,
    for each dose.

    Computed directly from the (exact, closed-form for H2/H3; root-found for
    H1) new steady state at each dose rather than by numerically integrating
    the ODE out to ``t_late`` -- for a constant post-treatment dose the two
    are the same thing once the system has relaxed, but the direct steady
    state is exact and orders of magnitude cheaper, which matters because
    :func:`smt.kinetics.hypotheses.separability_power` calls this once per
    dose per synthetic replicate. ``t_late`` is kept as a parameter (used by
    callers to also validate against :func:`time_course` at a matching late
    time, e.g. in tests) but no longer drives an integration length here.
    """
    if hyp not in HYPOTHESES:
        raise ValueError(f"unknown hypothesis {hyp!r}; must be one of {HYPOTHESES}")
    del t_late  # kept for API/back-compat with callers that pass it; see docstring
    out = np.empty_like(doses, dtype=float)
    for i, dose in enumerate(doses):
        if hyp == "H1_feedback":
            ss = _h1_fixed_point(params, float(dose))[:5]
        else:
            rates_fn = _rates_h2 if hyp == "H2_direct" else _rates_h3
            ss = steady_state(rates_fn(params, float(dose)))
        out[i] = ss[2] / (ss[3] + ss[4])  # C / (L_incl + L_skip)
    return out


def time_course(hyp: str, params: HypothesisParams, dose: float, times: np.ndarray) -> np.ndarray:
    """Circle:linear ratio vs. time after a single fixed dose."""
    t_eval = np.concatenate([[0.0], times]) if times[0] != 0.0 else times
    r = simulate(hyp, params, dose, (0.0, float(times[-1])), t_eval)
    return r["ratio"][-len(times) :]


def recommended_t_late(params: HypothesisParams) -> float:
    """A 'late enough' time: several H1 feedback timescales (the slowest one)."""
    return 20.0 / params.delta_S


# ---------------------------------------------------------------------------
# discriminator
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Discrimination:
    """Result of scoring one observed dataset against all three hypotheses."""

    chi2_by_hypothesis: dict[str, float]
    best_hypothesis: str
    delta_chi2: float  # chi2 of best runner-up minus chi2 of the best (>=0; bigger = more decisive)


def _chi2(predicted: np.ndarray, observed: np.ndarray, sigma: np.ndarray) -> float:
    return float(np.sum(((predicted - observed) / sigma) ** 2))


def discriminate(
    params_by_hypothesis: dict[str, HypothesisParams],
    doses: np.ndarray,
    t_late: float,
    dose_fixed: float,
    times: np.ndarray,
    observed_dose_response: np.ndarray,
    observed_time_course: np.ndarray,
    sigma_dose_response: np.ndarray,
    sigma_time_course: np.ndarray,
) -> Discrimination:
    """Score each hypothesis's PREDICTED shape against one observed dataset.

    This does not fit any hypothesis-specific free parameter to the data --
    it compares the shapes each hypothesis is mechanistically committed to
    (given ``params_by_hypothesis``) directly, by chi-square. That is
    deliberate: the discriminator in this project is about whether the
    hypotheses' predicted SIGNATURES (dose-tracking-with-no-lag vs.
    lagged-and-protein-dependent vs. flat-flux-with-slow-accumulation) are
    distinguishable in principle, not a parameter-fitting exercise for which
    there is no data to fit against (docs/HANDOFF.md).
    """
    chi2: dict[str, float] = {}
    for hyp, params in params_by_hypothesis.items():
        pred_dr = dose_response(hyp, params, doses, t_late)
        pred_tc = time_course(hyp, params, dose_fixed, times)
        chi2[hyp] = _chi2(pred_dr, observed_dose_response, sigma_dose_response) + _chi2(
            pred_tc, observed_time_course, sigma_time_course
        )
    ordered = sorted(chi2.items(), key=lambda kv: kv[1])
    best_hyp, best_chi2 = ordered[0]
    runner_up_chi2 = ordered[1][1]
    return Discrimination(chi2, best_hyp, runner_up_chi2 - best_chi2)


def separability_power(
    true_hyp: str,
    params_by_hypothesis: dict[str, HypothesisParams],
    doses: np.ndarray,
    t_late: float,
    dose_fixed: float,
    times: np.ndarray,
    sigma_rel: float,
    n_reps: int,
    rng: Optional[np.random.Generator] = None,
) -> float:
    """Fraction of noisy synthetic replicates for which the discriminator
    correctly recovers ``true_hyp`` as the best fit.

    Generates ``n_reps`` synthetic datasets by simulating ``true_hyp``'s own
    prediction and adding Gaussian noise with standard deviation
    ``sigma_rel * |signal|`` (a relative-precision noise model, e.g. a
    fractional RT-qPCR/ddPCR technical CV -- ILLUSTRATIVE, no measured assay
    precision for this system exists), then runs :func:`discriminate`
    against all three hypotheses and checks whether it recovers the true
    one. This is a real (if small) Monte Carlo power analysis, not a
    plausibility argument -- see docs/c2b34-discriminating-experiment.md for
    the concrete recommendation it produces.
    """
    if rng is None:
        rng = np.random.default_rng()
    true_params = params_by_hypothesis[true_hyp]
    true_dr = dose_response(true_hyp, true_params, doses, t_late)
    true_tc = time_course(true_hyp, true_params, dose_fixed, times)
    sigma_dr = np.maximum(sigma_rel * np.abs(true_dr), 1e-12)
    sigma_tc = np.maximum(sigma_rel * np.abs(true_tc), 1e-12)

    correct = 0
    for _ in range(n_reps):
        obs_dr = true_dr + rng.normal(0.0, sigma_dr)
        obs_tc = true_tc + rng.normal(0.0, sigma_tc)
        result = discriminate(
            params_by_hypothesis, doses, t_late, dose_fixed, times, obs_dr, obs_tc, sigma_dr, sigma_tc
        )
        if result.best_hypothesis == true_hyp:
            correct += 1
    return correct / n_reps
