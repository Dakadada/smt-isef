"""Tests for the Arm E per-mer CPU scorer layer (feng.titration.evaluate)."""

from __future__ import annotations

from pathlib import Path

import pytest

from feng.bench.scores import reverse_complement
from feng.honesty.states import ResultState
from feng.titration import evaluate as ev

# A known perfectly-complementary ASO/target pair (nusinersen 18-mer, DNA space).
ASO_18 = "TCACTTTCATAATGCTGG"
TARGET_18 = reverse_complement(ASO_18)  # CCAGCATTATGAAAGTGA


def _rung(length: int, aso: str, target: str) -> dict:
    return {
        "length": length,
        "chrom": "chr5",
        "start": 0,
        "end": length,
        "strand": "+",
        "target_5to3": target,
        "aso_5to3": aso,
        "state": ResultState.CANDIDATE.value,
    }


# --- binary resolution -------------------------------------------------------


def _have(tool: str) -> Path | None:
    return ev.find_vienna(tool)


def test_find_vienna_resolves_env_binaries():
    # In feng-env at least one Vienna tool must resolve via the interpreter sibling.
    assert any(_have(t) is not None for t in ev.VIENNA_TOOLS)


# --- A. IntaRNA 0.0 trap -----------------------------------------------------


def test_intarna_no_duplex_flag_on_noncomplementary():
    intarna = ev.resolve_intarna(None)
    if intarna is None:
        pytest.skip("IntaRNA not available")
    # A valid-length but non-complementary pair -> IntaRNA emits no data line and
    # scores.py returns its "no favourable duplex" 0.0. This is the trap the flag
    # exists to catch (distinct from a real weak-but-negative energy).
    poly = "AAAAAAAAAAAA"  # 12-mer, above IntaRNA's seed length
    dg, no_duplex = ev.intarna_on_target(poly, poly, intarna=intarna)
    assert dg == 0.0
    assert no_duplex is True


def test_intarna_errors_on_subseed_mer_is_not_fabricated():
    intarna = ev.resolve_intarna(None)
    if intarna is None:
        pytest.skip("IntaRNA not available")
    # A 2-mer is below IntaRNA's seed length: IntaRNA hard-errors (exit 1). It must
    # surface as IntaRNAError, never a fabricated 0.0.
    from feng.bench.scores import IntaRNAError

    with pytest.raises(IntaRNAError):
        ev.intarna_on_target("AT", "AT", intarna=intarna)


def test_intarna_real_duplex_is_negative_and_not_flagged():
    intarna = ev.resolve_intarna(None)
    if intarna is None:
        pytest.skip("IntaRNA not available")
    dg, no_duplex = ev.intarna_on_target(ASO_18, TARGET_18, intarna=intarna)
    assert dg < 0.0
    assert no_duplex is False


# --- B/C/D. Vienna wrappers parse real output --------------------------------


def test_rnaduplex_parses_strong_complementary_energy():
    exe = _have("RNAduplex")
    if exe is None:
        pytest.skip("RNAduplex not available")
    dg = ev.rnaduplex_dg(ASO_18, TARGET_18, rnaduplex=exe)
    assert isinstance(dg, float)
    assert dg < -10.0  # a perfect 18-bp duplex is strongly negative


def test_rnacofold_parses_energy():
    exe = _have("RNAcofold")
    if exe is None:
        pytest.skip("RNAcofold not available")
    dg = ev.rnacofold_dg(ASO_18, TARGET_18, rnacofold=exe)
    assert isinstance(dg, float)
    assert dg < -10.0


def test_rnaup_parses_total_energy():
    exe = _have("RNAup")
    if exe is None:
        pytest.skip("RNAup not available")
    dg = ev.rnaup_dg(ASO_18, TARGET_18, rnaup=exe)
    assert isinstance(dg, float)
    assert dg < 0.0


def _rc_rna(seq: str) -> str:
    comp = {"A": "U", "U": "A", "G": "C", "C": "G"}
    return "".join(comp[b] for b in reversed(seq))


def test_rnaplfold_accessibility_context_lowers_with_pairing_flank():
    # MED-1 regression: accessibility must reflect the target's genomic context,
    # not an isolated fragment. Embed the SAME 12-nt target at the SAME offset in
    # two 40-nt flanks: a pairing flank (the target's reverse complement placed
    # adjacent) vs a non-pairing poly-A flank. The pairing context lets the target
    # base-pair into its flank, so its mean unpaired probability must drop below
    # the neutral context. If accessibility ignored context the two would match.
    exe = _have("RNAplfold")
    if exe is None:
        pytest.skip("RNAplfold not available")
    target = "GGAGGAGGAGGA"
    rc = _rc_rna(target)  # UCCUCCUCCUCC -- complementary to the target
    pad = 40
    tlen = len(target)
    neutral = "A" * pad + target + "A" * pad
    # rc immediately abutting the target on each side (rest poly-A) so it pairs.
    left = "A" * (pad - tlen) + rc
    right = rc + "A" * (pad - tlen)
    pairing = left + target + right
    assert len(neutral) == len(pairing) == tlen + 2 * pad
    acc_neutral = ev.mean_unpaired_over_window(neutral, pad, tlen, rnaplfold=exe)
    acc_pairing = ev.mean_unpaired_over_window(pairing, pad, tlen, rnaplfold=exe)
    assert 0.0 <= acc_pairing <= 1.0
    assert 0.0 <= acc_neutral <= 1.0
    assert acc_pairing < acc_neutral


def test_rnaplfold_accessibility_context_offset_selects_target_only():
    # The averaged rows must be exactly the target bases, not the pad: a fully
    # poly-A window is ~1.0 unpaired everywhere, so the mean is ~1.0 and finite.
    exe = _have("RNAplfold")
    if exe is None:
        pytest.skip("RNAplfold not available")
    pad, tlen = 20, 10
    seq = "A" * (2 * pad + tlen)
    acc = ev.mean_unpaired_over_window(seq, pad, tlen, rnaplfold=exe)
    assert 0.9 <= acc <= 1.0


def test_evaluate_ladder_accessibility_none_without_fasta():
    # No genome_fasta -> accessibility stays None rather than falling back to the
    # isolated-window bug.
    rungs = [_rung(18, ASO_18, TARGET_18)]
    out = ev.evaluate_ladder(rungs)
    assert out[0]["rnaplfold_accessibility"] is None


_V50_FASTA = Path(
    "/home/dashii/orcd/scratch/feng/gencode/v50/GRCh38.primary_assembly.genome.fa"
)


@pytest.mark.skipif(not _V50_FASTA.is_file(), reason="v50 genome FASTA not staged")
def test_rnaplfold_accessibility_real_genomic_context():
    exe = _have("RNAplfold")
    if exe is None:
        pytest.skip("RNAplfold not available")
    from feng.titration.ladder import build_ladder

    rungs = build_ladder(genome_fasta=_V50_FASTA, min_len=18, max_len=18)
    acc = ev.rnaplfold_accessibility(rungs[0], rnaplfold=exe, genome_fasta=_V50_FASTA)
    assert isinstance(acc, float)
    assert 0.0 <= acc <= 1.0


def test_rnafold_parses_mfe():
    exe = _have("RNAfold")
    if exe is None:
        pytest.skip("RNAfold not available")
    mfe = ev.rnafold_mfe(ASO_18, rnafold=exe)
    assert isinstance(mfe, float)
    assert mfe <= 0.0  # MFE is never positive


# --- E. off-target scan on a synthetic reference -----------------------------

# Filler with a composition that cannot contain the target window within budget.
_FILLER = "AAAACCCCAAAACCCCAAAACCCC"


def test_off_target_counts_planted_exact_hit():
    reference = _FILLER + TARGET_18 + _FILLER
    hits = ev.off_target_hits(ASO_18, reference, max_mismatch=0)
    assert len(hits) == 1


def test_off_target_no_hit_reference_returns_zero():
    reference = _FILLER + _FILLER
    hits = ev.off_target_hits(ASO_18, reference, max_mismatch=0)
    assert len(hits) == 0


def test_off_target_near_miss_excluded_at_zero_included_at_one():
    # One substitution in the planted target: not a hit at mm=0, a hit at mm=1.
    near_miss = "C" + TARGET_18[1:] if TARGET_18[0] != "C" else "A" + TARGET_18[1:]
    assert near_miss != TARGET_18
    reference = _FILLER + near_miss + _FILLER
    assert len(ev.off_target_hits(ASO_18, reference, max_mismatch=0)) == 0
    assert len(ev.off_target_hits(ASO_18, reference, max_mismatch=1)) == 1


# --- assembler ---------------------------------------------------------------


def test_evaluate_ladder_schema_and_state():
    rungs = [_rung(2, "AT", "AT"), _rung(18, ASO_18, TARGET_18)]
    out = ev.evaluate_ladder(rungs)
    assert len(out) == 2
    for row in out:
        for key in ev._SCORE_KEYS:
            assert key in row
        assert row["state"] == ResultState.CANDIDATE.value
        # Off-target column is None unless a reference is supplied.
        assert row["off_target_hits"] is None
    # Original ladder fields survive.
    assert out[1]["aso_5to3"] == ASO_18


def test_evaluate_ladder_intarna_three_states():
    intarna = ev.resolve_intarna(None)
    if intarna is None:
        pytest.skip("IntaRNA not available")
    poly = "AAAAAAAAAAAA"
    rungs = [
        _rung(2, "AT", "AT"),               # sub-seed -> IntaRNA errors -> None
        _rung(12, poly, poly),              # valid length, no duplex -> 0.0 flagged
        _rung(18, ASO_18, TARGET_18),       # real complementary -> negative
    ]
    subseed, no_duplex, full = ev.evaluate_ladder(rungs)
    # Sub-seed: could not compute -> None, never a fabricated 0.0.
    assert subseed["intarna_dg"] is None and subseed["intarna_no_duplex"] is None
    # No favourable duplex: 0.0 correctly flagged, not read as weak binding.
    assert no_duplex["intarna_dg"] == 0.0 and no_duplex["intarna_no_duplex"] is True
    # Real duplex: negative and unflagged.
    assert full["intarna_dg"] < 0.0 and full["intarna_no_duplex"] is False


def test_evaluate_ladder_off_target_column_populates_with_reference():
    reference = _FILLER + TARGET_18 + _FILLER
    rungs = [_rung(18, ASO_18, TARGET_18)]
    out = ev.evaluate_ladder(rungs, off_target_reference=reference, off_target_max_mismatch=0, off_target_max_edit=None)
    assert out[0]["off_target_hits"] == 1


# --- never-fabricate discipline ---------------------------------------------


def test_module_has_no_hardcoded_energy_literals():
    # Guard against baking a fixture's printed energy into the module.
    src = Path(ev.__file__).read_text()
    for forbidden in ("-27.90", "-26.99", "-27.91"):
        assert forbidden not in src
