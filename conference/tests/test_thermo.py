from pathlib import Path

import pytest

import math

from conference.aso_pipeline.thermo import (
    DuplexError,
    effective_monovalent_from_mixed_salt,
    intarna_many_targets,
    parse_rnaduplex,
    rnaduplex_score,
    score_pair,
)
from feng.titration.evaluate import intarna_on_target


INTARNA = Path("/home/dashii/orcd/pool/feng-env/bin/IntaRNA")
RNADUPLEX = Path("/home/dashii/orcd/pool/feng-env/bin/RNAduplex")


def test_von_ahsen_effective_monovalent_folds_in_magnesium():
    # 50 mM mono + 1.5 mM Mg2+ -> 0.050 + 120*sqrt(1.5) mM = 0.19697 M
    eq = effective_monovalent_from_mixed_salt(0.050, 0.0015)
    assert math.isclose(eq, 0.050 + 120.0 * math.sqrt(1.5) / 1000.0, rel_tol=1e-9)
    assert math.isclose(eq, 0.196969, abs_tol=1e-5)
    # zero divalent collapses to the plain monovalent concentration
    assert math.isclose(effective_monovalent_from_mixed_salt(0.050, 0.0), 0.050)
    # monotonic in Mg2+
    assert effective_monovalent_from_mixed_salt(0.050, 0.003) > eq
    with pytest.raises(ValueError):
        effective_monovalent_from_mixed_salt(0.0, 0.0)
    with pytest.raises(ValueError):
        effective_monovalent_from_mixed_salt(-0.01, 0.0015)


def test_rnaduplex_parser_keeps_structure_and_literal_energy_separate():
    parsed = parse_rnaduplex("((((.&.))))  1,5  :  2,6  (-7.30)\n")
    assert parsed.structure == "((((.&.))))"
    assert parsed.energy_kcal_mol == -7.30
    assert parsed.paired_bases == 4


def test_rnaduplex_parser_rejects_missing_energy_instead_of_fabricating_zero():
    with pytest.raises(DuplexError, match="parse"):
        parse_rnaduplex("no favorable structure")


@pytest.mark.skipif(not RNADUPLEX.is_file(), reason="RNAduplex unavailable")
def test_rnaduplex_real_binary_uses_explicit_temperature_and_salt():
    result = rnaduplex_score(
        "AUAAUGCUGG", "CCAGCAUUAU", rnaduplex=RNADUPLEX, salt_molar=0.089
    )
    assert result.energy_kcal_mol < 0
    assert result.raw_reference_dg_kcal_mol == -12.8
    assert result.salt_correction_kcal_mol == pytest.approx(2.50, abs=0.01)
    assert result.energy_kcal_mol == pytest.approx(-10.30, abs=0.01)
    assert result.paired_bases == 10
    assert result.temperature_c == 37.0
    assert result.salt_molar == 0.089


@pytest.mark.skipif(not RNADUPLEX.is_file(), reason="RNAduplex unavailable")
def test_lower_monovalent_salt_makes_duplex_less_favorable():
    low = rnaduplex_score(
        "AUAAUGCUGG", "CCAGCAUUAU", rnaduplex=RNADUPLEX, salt_molar=0.025
    )
    high = rnaduplex_score(
        "AUAAUGCUGG", "CCAGCAUUAU", rnaduplex=RNADUPLEX, salt_molar=0.150
    )
    assert low.raw_reference_dg_kcal_mol == high.raw_reference_dg_kcal_mol
    assert low.energy_kcal_mol > high.energy_kcal_mol


def test_dna_rna_rerun_requires_an_explicit_validated_hybrid_parameter_file():
    with pytest.raises(ValueError, match="DNA:RNA"):
        rnaduplex_score(
            "ATAATGCTGG",
            "CCAGCAUUAU",
            rnaduplex=RNADUPLEX,
            salt_molar=0.089,
            chemistry="DNA:RNA",
        )


@pytest.mark.skipif(
    not INTARNA.is_file() or not RNADUPLEX.is_file(), reason="thermo binaries unavailable"
)
def test_score_pair_reports_both_raw_models_without_combining_them():
    result = score_pair(
        "AUAAUGCUGG",
        "CCAGCAUUAU",
        intarna=INTARNA,
        rnaduplex=RNADUPLEX,
        salt_molar=0.089,
    )
    assert result.intarna_dg_kcal_mol < 0
    assert result.intarna_no_duplex is False
    assert result.nn_dg_kcal_mol < 0
    assert result.nn_salt_molar == 0.089


@pytest.mark.skipif(not INTARNA.is_file(), reason="IntaRNA unavailable")
def test_batched_intarna_matches_existing_wrapper_with_same_orientation():
    targets = {
        "exact": "CCAGCAUUAU",
        "one_mismatch": "CCAGCAUUAA",
    }

    batched = intarna_many_targets(
        "AUAAUGCUGG", targets, intarna=INTARNA, threads=2
    )

    assert batched == {
        target_id: intarna_on_target("AUAAUGCUGG", target, intarna=INTARNA)
        for target_id, target in targets.items()
    }
