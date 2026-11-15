import pytest

from conference.aso_pipeline.psi_stub import PsiDataUnavailable, fit_dg_psi


def test_future_psi_interface_validates_required_schema_before_refusing_fit():
    with pytest.raises(ValueError, match="replicate"):
        fit_dg_psi(
            [{"aso_id": "F18", "length": 18, "PSI_percent": 50.0}],
            {"F18": -23.1},
        )


def test_future_psi_interface_never_computes_without_later_explicit_enablement():
    rows = [
        {"aso_id": "F18", "length": 18, "replicate": "r1", "PSI_percent": 50.0}
    ]
    with pytest.raises(PsiDataUnavailable, match="not available"):
        fit_dg_psi(rows, {"F18": -23.1})

