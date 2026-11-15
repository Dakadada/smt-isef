import pytest

from conference.aso_pipeline.item6 import (
    censored_no_duplex_lower_tail_p,
    summarize_null,
)


def test_null_summary_uses_population_sd_linear_quantiles_and_lower_tail():
    summary = summarize_null(-1.0, [-5.0, -4.0, -3.0, -2.0, -1.0], lower_is_stronger=True)
    assert summary["null_n_numeric"] == 5
    assert summary["null_mean"] == -3.0
    assert summary["null_sd_population"] == pytest.approx(2**0.5)
    assert summary["null_q25"] == -4.0
    assert summary["null_median"] == -3.0
    assert summary["null_q75"] == -2.0
    assert summary["real_minus_null_mean"] == 2.0
    assert summary["one_sided_enrichment_p"] == pytest.approx(1.0)


def test_null_summary_higher_is_stronger_reverses_enrichment_tail():
    summary = summarize_null(6.0, [1.0, 2.0, 3.0, 4.0, 5.0], lower_is_stronger=False)
    assert summary["real_percentile_midrank"] == 100.0
    assert summary["one_sided_enrichment_p"] == pytest.approx(1 / 6)


def test_null_summary_refuses_empty_numeric_values():
    with pytest.raises(ValueError, match="no numeric"):
        summarize_null(-1.0, [], lower_is_stronger=True)


def test_no_duplex_draws_enter_tail_denominator_without_fake_energy():
    assert censored_no_duplex_lower_tail_p(-10.0, [-12.0, -8.0], 8) == pytest.approx(
        2 / 11
    )
