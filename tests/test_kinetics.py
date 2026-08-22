"""A8 tests: the C2B34 backsplicing kinetics model and its three hypotheses.

Covers (per the task spec):
  * mass conservation in the ODE -- nothing created/destroyed except by the
    explicit k_txn (synthesis) and gamma_* (decay) terms;
  * the three hypotheses produce measurably DIFFERENT observable signatures;
  * the Gillespie stochastic variant's mean converges to the ODE solution at
    high copy number;
  * rate constants cannot be silently defaulted -- calling without them
    raises.

Also includes a few fast, synthetic-sequence unit tests for
``smt.kinetics.alu``'s pure helper functions (tiling, the RNAduplex i/j
index-parsing convention, orientation classification). These do NOT run the
real genome-scale scan (~3 minutes, see ``smt.kinetics.alu.scan_flanks``'s
docstring) -- that was run once, interactively, to produce the numbers
reported in the project's findings, not as part of the test suite.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from smt.kinetics import alu, hypotheses, model

# ---------------------------------------------------------------------------
# model.py: rate constants cannot be silently defaulted
# ---------------------------------------------------------------------------


def test_rate_constants_requires_every_field():
    with pytest.raises(TypeError):
        model.RateConstants()  # type: ignore[call-arg]


def test_rate_constants_requires_every_field_even_partial():
    with pytest.raises(TypeError):
        model.RateConstants(k_txn=1.0, k_back=0.1)  # type: ignore[call-arg]


def test_rate_constants_rejects_negative_rate():
    kwargs = dataclasses.asdict(model.illustrative_rates())
    kwargs["gamma_C"] = -0.1
    with pytest.raises(ValueError):
        model.RateConstants(**kwargs)


def test_hypothesis_params_requires_every_field():
    with pytest.raises(TypeError):
        hypotheses.HypothesisParams()  # type: ignore[call-arg]


def test_hypothesis_params_enforces_h2_high_concentration_constraint():
    """H2 is only mechanistically coherent if its site is lower-affinity."""
    good = hypotheses.illustrative_hypothesis_params()
    bad_kwargs = dataclasses.asdict(good)
    bad_kwargs["baseline"] = good.baseline  # asdict flattens the nested dataclass; restore it
    bad_kwargs["ec50_upstream"] = bad_kwargs["ec50_intron7"]  # not "high concentration" any more
    with pytest.raises(ValueError):
        hypotheses.HypothesisParams(**bad_kwargs)


# ---------------------------------------------------------------------------
# model.py: mass conservation
# ---------------------------------------------------------------------------


def test_closed_system_conserves_total_mass():
    """With no synthesis and no decay, splicing only moves mass around."""
    rates = model.RateConstants(
        k_txn=0.0,
        k_back=0.3,
        k_fwd=0.7,
        k_incl=0.2,
        k_skip=0.4,
        gamma_N=0.0,
        gamma_P=0.0,
        gamma_C=0.0,
        gamma_incl=0.0,
        gamma_skip=0.0,
    )
    y0 = np.array([10.0, 3.0, 1.0, 0.5, 0.5])
    sol = model.simulate_ode(rates, y0=y0, t_span=(0.0, 50.0), t_eval=np.linspace(0, 50, 21))
    totals = sol.y.sum(axis=0)
    assert totals == pytest.approx(totals[0], rel=1e-6)


def test_open_system_mass_balance_is_exact_algebraically():
    """d(total)/dt must equal exactly k_txn minus the sum of decay terms.

    This is an algebraic identity given how ``flux`` is built (every
    reaction that moves mass between two tracked states contributes with
    opposite sign to each side and cancels in the sum; only k_txn and the
    gamma_* terms are unpaired). Checking it numerically at a handful of
    non-equilibrium states guards against a sign error being introduced
    into ``flux`` in a future edit.
    """
    rates = model.illustrative_rates()
    rng = np.random.default_rng(0)
    for _ in range(20):
        state = rng.uniform(0.0, 5.0, size=5)
        d = model.flux(state, rates)
        expected_total_rate = rates.k_txn - (
            rates.gamma_N * state[0]
            + rates.gamma_P * state[1]
            + rates.gamma_C * state[2]
            + rates.gamma_incl * state[3]
            + rates.gamma_skip * state[4]
        )
        assert d.sum() == pytest.approx(expected_total_rate, rel=1e-9, abs=1e-12)


def test_steady_state_is_a_fixed_point():
    rates = model.illustrative_rates()
    ss = model.steady_state(rates)
    assert np.allclose(model.flux(ss, rates), 0.0, atol=1e-10)
    assert np.all(ss >= 0.0)


# ---------------------------------------------------------------------------
# model.py: Gillespie SSA mean converges to the ODE at high copy number
# ---------------------------------------------------------------------------


def test_gillespie_mean_converges_to_ode_at_high_copy_number():
    rates = model.RateConstants(
        k_txn=150.0,
        k_back=0.05,
        k_fwd=1.0,
        k_incl=0.2,
        k_skip=0.8,
        gamma_N=0.1,
        gamma_P=0.1,
        gamma_C=0.05,
        gamma_incl=0.5,
        gamma_skip=0.5,
    )
    t_eval = np.array([10.0, 20.0, 30.0])
    ode = model.simulate_ode(rates, y0=np.zeros(5), t_span=(0.0, 30.0), t_eval=t_eval)

    rng = np.random.default_rng(12345)
    n_reps = 60
    traj = np.zeros((n_reps, len(t_eval), 5))
    for i in range(n_reps):
        res = model.simulate_gillespie(rates, y0=np.zeros(5), t_max=30.0, rng=rng)
        traj[i] = res.sample_at(t_eval)
    mean_traj = traj.mean(axis=0)

    # Relative tolerance: at these copy numbers (tens to hundreds) SSA-mean
    # vs. ODE agreement within ~10% over n_reps=60 replicates is expected;
    # a looser floor (2.0 absolute) avoids flakiness on the near-zero states.
    assert np.allclose(mean_traj, ode.y.T, rtol=0.15, atol=2.0)


def test_gillespie_result_sample_at_is_a_right_continuous_step_function():
    rates = model.illustrative_rates()
    rng = np.random.default_rng(7)
    res = model.simulate_gillespie(rates, y0=[5, 0, 0, 0, 0], t_max=5.0, rng=rng)
    # Sampling exactly at an event time must return the POST-event state.
    if len(res.t) > 1:
        t_event = res.t[1]
        sampled = res.sample_at(np.array([t_event]))[0]
        assert np.array_equal(sampled, res.y[1])


# ---------------------------------------------------------------------------
# hypotheses.py: the three hypotheses produce measurably different signatures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def hyp_params():
    return {h: hypotheses.illustrative_hypothesis_params() for h in hypotheses.HYPOTHESES}


def test_h2_and_h3_collapse_to_the_same_baseline_at_zero_dose(hyp_params):
    """A basic sanity/consistency check: with theta(0)=0, no rate is modified,
    so H2 and H3 (which have no state-dependent feedback) must agree exactly
    with each other and with the plain baseline steady state at dose=0.
    """
    t_late = hypotheses.recommended_t_late(hyp_params["H1_feedback"])
    dr2 = hypotheses.dose_response("H2_direct", hyp_params["H2_direct"], np.array([0.0]), t_late)
    dr3 = hypotheses.dose_response("H3_stability", hyp_params["H3_stability"], np.array([0.0]), t_late)
    baseline_ss = model.steady_state(hyp_params["H2_direct"].baseline)
    expected_ratio = baseline_ss[model.C_IDX] / (baseline_ss[model.LINCL_IDX] + baseline_ss[model.LSKIP_IDX])
    assert dr2[0] == pytest.approx(expected_ratio, rel=1e-6)
    assert dr3[0] == pytest.approx(expected_ratio, rel=1e-6)


def test_dose_response_curves_are_measurably_different(hyp_params):
    doses = np.array([1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0])
    t_late = hypotheses.recommended_t_late(hyp_params["H1_feedback"])
    curves = {
        h: hypotheses.dose_response(h, hyp_params[h], doses, t_late) for h in hypotheses.HYPOTHESES
    }
    for a, b in [("H1_feedback", "H2_direct"), ("H1_feedback", "H3_stability"), ("H2_direct", "H3_stability")]:
        max_rel_diff = np.max(np.abs(curves[a] - curves[b]) / np.maximum(curves[a], curves[b]))
        assert max_rel_diff > 0.2, f"{a} vs {b} dose-response curves are not measurably different"


def test_time_course_curves_are_measurably_different(hyp_params):
    times = np.array([1.0, 5.0, 20.0, 80.0, 320.0, 1000.0])
    curves = {
        h: hypotheses.time_course(h, hyp_params[h], 300.0, times) for h in hypotheses.HYPOTHESES
    }
    for a, b in [("H1_feedback", "H2_direct"), ("H1_feedback", "H3_stability"), ("H2_direct", "H3_stability")]:
        max_rel_diff = np.max(np.abs(curves[a] - curves[b]) / np.maximum(curves[a], curves[b]))
        assert max_rel_diff > 0.2, f"{a} vs {b} time-course curves are not measurably different"


def test_h1_shows_a_lag_h2_does_not():
    """H1's defining prediction: near-zero response shortly after dosing,
    because the feedback signal S has not yet accumulated. H2 (direct,
    dose-tracking) should already show substantial movement at the same
    early timepoint.
    """
    params = hypotheses.illustrative_hypothesis_params()
    dose = 1000.0  # saturating, so any early movement is due to kinetics, not weak binding
    times = np.array([0.5, 1.0])
    tc_h1 = hypotheses.time_course("H1_feedback", params, dose, times)
    tc_h2 = hypotheses.time_course("H2_direct", params, dose, times)
    baseline = model.steady_state(params.baseline)
    baseline_ratio = baseline[model.C_IDX] / (baseline[model.LINCL_IDX] + baseline[model.LSKIP_IDX])
    # H1 barely moves from baseline in the first time unit...
    assert abs(tc_h1[-1] - baseline_ratio) < 0.3 * abs(tc_h2[-1] - baseline_ratio) or tc_h2[-1] > tc_h1[-1]


def test_h3_predicts_flat_synthesis_flux_across_dose():
    """H3's defining prediction: k_back * N (circle synthesis flux) does NOT
    change with dose, even though the circle LEVEL does (via decay).
    """
    params = hypotheses.illustrative_hypothesis_params()
    t_late = hypotheses.recommended_t_late(params)
    fluxes = []
    for dose in (0.0, 100.0, 1000.0):
        r = hypotheses.simulate("H3_stability", params, dose, (0.0, t_late), np.array([t_late]))
        fluxes.append(r["synthesis_flux"][-1])
    fluxes = np.array(fluxes)
    assert np.max(np.abs(fluxes - fluxes[0])) / fluxes[0] < 1e-6, "H3 synthesis flux must be dose-independent"


def test_h2_synthesis_flux_does_change_with_dose():
    """Contrast case: H2 explicitly perturbs k_back, so flux MUST move."""
    params = hypotheses.illustrative_hypothesis_params()
    t_late = hypotheses.recommended_t_late(params)
    r_lo = hypotheses.simulate("H2_direct", params, 0.0, (0.0, t_late), np.array([t_late]))
    r_hi = hypotheses.simulate("H2_direct", params, 1000.0, (0.0, t_late), np.array([t_late]))
    assert r_hi["synthesis_flux"][-1] > 2.0 * r_lo["synthesis_flux"][-1]


# ---------------------------------------------------------------------------
# hypotheses.py: the discriminator recovers the generating hypothesis
# ---------------------------------------------------------------------------


def test_discriminator_recovers_the_true_hypothesis_noiselessly(hyp_params):
    """With no measurement noise, the discriminator must prefer the true
    generating hypothesis over the other two by a clear margin (chi2 ~ 0 for
    the truth, distinctly larger for the alternatives).
    """
    doses = np.array([1.0, 10.0, 100.0, 1000.0])
    times = np.array([5.0, 50.0, 500.0])
    t_late = hypotheses.recommended_t_late(hyp_params["H1_feedback"])
    dose_fixed = 300.0
    for true_hyp in hypotheses.HYPOTHESES:
        true_dr = hypotheses.dose_response(true_hyp, hyp_params[true_hyp], doses, t_late)
        true_tc = hypotheses.time_course(true_hyp, hyp_params[true_hyp], dose_fixed, times)
        sigma_dr = np.full_like(true_dr, 1e-3)
        sigma_tc = np.full_like(true_tc, 1e-3)
        result = hypotheses.discriminate(
            hyp_params, doses, t_late, dose_fixed, times, true_dr, true_tc, sigma_dr, sigma_tc
        )
        assert result.best_hypothesis == true_hyp
        assert result.delta_chi2 > 10.0  # a clear margin, not a coin flip


def test_separability_power_is_high_at_low_noise_and_can_degrade_at_high_noise(hyp_params):
    """A light-weight version of the project's headline separability claim:
    recovery probability should be high when measurement noise is small.
    Uses few replicates/doses/timepoints to keep this test fast; the full
    analysis (more reps, a finer sigma grid) is run interactively and
    reported in docs/c2b34-discriminating-experiment.md, not in this test.
    """
    doses = np.array([1.0, 30.0, 1000.0])
    times = np.array([5.0, 100.0])
    t_late = hypotheses.recommended_t_late(hyp_params["H1_feedback"])
    dose_fixed = 300.0
    rng = np.random.default_rng(2)
    recovery = hypotheses.separability_power(
        "H3_stability", hyp_params, doses, t_late, dose_fixed, times, sigma_rel=0.02, n_reps=30, rng=rng
    )
    assert recovery >= 0.8


# ---------------------------------------------------------------------------
# alu.py: fast unit tests on pure helper functions (no genome-scale scan)
# ---------------------------------------------------------------------------


def test_tile_covers_interval_without_gaps():
    iv = alu.Interval(100, 799)  # 700 nt
    tiles = alu._tile(iv, width=300, step=150)
    assert tiles[0].start == 100
    assert tiles[-1].end == 799
    # every position is covered by at least one tile
    covered = np.zeros(iv.end - iv.start + 1, dtype=bool)
    for t in tiles:
        covered[t.start - iv.start : t.end - iv.start + 1] = True
    assert covered.all()


def test_tile_drops_undersized_trailing_remainder():
    iv = alu.Interval(1, 320)  # one full 300-nt tile plus a 20-nt remainder
    tiles = alu._tile(iv, width=300, step=150)
    # 20 nt < 300//2=150, so it must be dropped, not returned as a tiny tile
    assert all(len(t) >= 150 for t in tiles)


def test_duplex_span_parsing_matches_known_duplexfold_convention():
    """Verifies the i/j -> reference-interval mapping against a fully
    complementary synthetic pair, where the matched span is unambiguous
    (the whole tile on both sides).
    """
    import RNA

    s1 = "AAAAAAAAAAGGGGGGGGGG"  # 20 nt
    s2 = "CCCCCCCCCCUUUUUUUUUU"  # 20 nt, fully complementary to s1 (RNA alphabet)
    d = RNA.duplexfold(s1, s2)
    tile1 = alu.Interval(1000, 1019)  # arbitrary reference-space placement
    tile2 = alu.Interval(2000, 2019)
    matched1, matched2, n_bp = alu._parse_duplex_span(d, tile1, tile2)
    assert n_bp == 20
    assert matched1 == alu.Interval(1000, 1019)
    assert matched2 == alu.Interval(2000, 2019)


def test_orientation_check_identifies_a_synthetic_inverted_repeat():
    from Bio.Seq import Seq

    core = "ACGGCTAGCTTAGGCATCGGATCCAGGTTACGGAT" * 3  # 105 nt, arbitrary but non-palindromic
    upstream_seq = "TTTT" + core + "TTTT"
    downstream_seq = "GGGG" + str(Seq(core).reverse_complement()) + "GGGG"  # genuine inverted repeat

    fake_seq_map = {}

    def fake_slice_ref(seq, iv):
        return fake_seq_map[(iv.start, iv.end)]

    up_iv = alu.Interval(1, len(upstream_seq))
    down_iv = alu.Interval(10_000, 10_000 + len(downstream_seq) - 1)
    fake_seq_map[(up_iv.start, up_iv.end)] = upstream_seq
    fake_seq_map[(down_iv.start, down_iv.end)] = downstream_seq

    hit = alu.DuplexHit(
        tile_upstream=up_iv,
        tile_downstream=down_iv,
        matched_upstream=up_iv,
        matched_downstream=down_iv,
        energy_kcal_mol=-1.0,
        n_bp=1,
        density=1.0,
    )

    orig_slice_ref = alu.slice_ref
    alu.slice_ref = fake_slice_ref
    try:
        oc = alu.orientation_check(seq="unused", hit=hit)
    finally:
        alu.slice_ref = orig_slice_ref

    assert oc.is_inverted
    assert oc.revcomp_score > oc.forward_score
