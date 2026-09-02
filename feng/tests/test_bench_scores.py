"""Frozen sequence-score tests: Hamming baseline + IntaRNA wrapper."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from feng.bench.scores import (
    DEFAULT_WINDOW,
    NUSINERSEN_5P_TO_3P,
    NUSINERSEN_RC,
    IntaRNAError,
    baseline_hamming_scores,
    exon_window_specs,
    find_intarna,
    intarna_min_energy,
    min_hamming_to_target,
    reverse_complement,
)

ENV_PREFIX = Path("/home/dashii/orcd/pool/feng-env")


def _samtools() -> Path | None:
    candidate = ENV_PREFIX / "bin" / "samtools"
    if candidate.is_file():
        return candidate
    found = shutil.which("samtools")
    return Path(found) if found else None


def _intarna() -> Path | None:
    return find_intarna(ENV_PREFIX / "bin" / "IntaRNA") or find_intarna()


def _write_indexed_fasta(path: Path, contigs: dict[str, str], samtools: Path) -> None:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in contigs.items()))
    subprocess.run([str(samtools), "faidx", str(path)], check=True, capture_output=True)


# --- pure sequence helpers ---------------------------------------------------


def test_reverse_complement_roundtrip() -> None:
    assert reverse_complement("ACGT") == "ACGT"
    assert reverse_complement("AAAA") == "TTTT"
    assert reverse_complement(reverse_complement(NUSINERSEN_5P_TO_3P)) == NUSINERSEN_5P_TO_3P


def test_nusinersen_rc_is_expected() -> None:
    assert NUSINERSEN_5P_TO_3P == "TCACTTTCATAATGCTGG"
    assert NUSINERSEN_RC == "CCAGCATTATGAAAGTGA"
    assert len(NUSINERSEN_RC) == 18


def test_min_hamming_exact_and_mismatch() -> None:
    target = NUSINERSEN_RC
    assert min_hamming_to_target("AAAA" + target + "GGGG", target) == 0
    one_off = "T" + target[1:] if target[0] != "T" else "A" + target[1:]
    assert min_hamming_to_target(one_off, target) == 1


def test_min_hamming_skips_windows_with_N_and_short_seq() -> None:
    target = NUSINERSEN_RC
    # a window with N must not be scored as a match
    assert min_hamming_to_target("N" * len(target), target) == len(target)
    # too-short sequence returns the worst distance
    assert min_hamming_to_target("ACGT", target) == len(target)


def test_exon_window_specs_clamps_and_pads() -> None:
    events = pd.DataFrame(
        [{"exon_key": "chr1:5-20:+", "chrom": "chr1", "exon_start_0base": 5, "exon_end": 20, "strand": "+"}]
    )
    spec = exon_window_specs(events, window=10)[0]
    assert spec.start_0base == 0  # 5 - 10 clamped to 0
    assert spec.end == 30
    assert spec.strand == "+"


# --- baseline over a real (synthetic) FASTA ----------------------------------


def test_baseline_hamming_plus_and_minus_strand(tmp_path: Path) -> None:
    samtools = _samtools()
    if samtools is None:
        pytest.skip("samtools required for FASTA extraction")
    rc = NUSINERSEN_RC
    # + strand: transcript sequence == genomic; embed RC directly
    plus = "A" * 20 + rc + "C" * 20
    # - strand: transcript is RC of genomic, so embed RC-of-RC on the genome
    minus = "G" * 15 + reverse_complement(rc) + "T" * 15
    fasta = tmp_path / "mini.fa"
    _write_indexed_fasta(fasta, {"chr1": plus, "chr2": minus}, samtools)

    events = pd.DataFrame(
        [
            {"exon_key": "chr1:20-38:+", "chrom": "chr1", "exon_start_0base": 20, "exon_end": 38, "strand": "+"},
            {"exon_key": "chr2:15-33:-", "chrom": "chr2", "exon_start_0base": 15, "exon_end": 33, "strand": "-"},
        ]
    )
    scored = baseline_hamming_scores(events, fasta, window=0, samtools=samtools)
    assert scored["baseline_min_hamming"].tolist() == [0, 0]


# --- IntaRNA wrapper ---------------------------------------------------------


def test_intarna_wrapper_missing_binary_raises() -> None:
    with pytest.raises(IntaRNAError):
        intarna_min_energy("ACGUACGU", intarna=Path("/no/such/IntaRNA"))


def test_intarna_energy_sign_on_tiny_fixture() -> None:
    exe = _intarna()
    if exe is None:
        pytest.skip("IntaRNA not available in this environment")
    # a target complementary to nusinersen must yield a favourable (negative) energy
    complementary = reverse_complement(NUSINERSEN_5P_TO_3P)
    e_bind = intarna_min_energy(complementary, intarna=exe)
    assert e_bind < 0.0
    # a poly-A target has no favourable duplex -> IntaRNA reports none (0.0)
    e_none = intarna_min_energy("A" * 40, intarna=exe)
    assert e_none == 0.0


# --- match -> scores seam (end to end) ---------------------------------------


def _seam_universe():
    """Small universe with covariate spread and modest coordinates for the seam."""
    import pandas as pd

    rows = []
    # positive at chr1 1-based 200..249 -> exon_key chr1:199-249:+
    def row(start0, end, psi, jr, cpm, gc):
        return {
            "exon_key": f"chr1:{start0}-{end}:+", "chrom": "chr1",
            "exon_start_0base": start0, "exon_end": end, "strand": "+",
            "detectable": True, "baseline_psi": psi,
            "mean_total_junction_reads": jr, "control_mean_cpm": cpm,
            "exon_length": end - start0, "gc_fraction": gc,
        }
    rows.append(row(199, 249, 0.5, 50.0, 10.0, 0.5))            # positive (len 50)
    for i in range(6):                                          # backgrounds, real spread
        s0 = 400 + i * 100
        rows.append(row(s0, s0 + 48 + i * 3, 0.5 + i * 0.01, 50 + i, 10 + i, 0.5 + i * 0.005))
    universe = pd.DataFrame(rows)
    positives = pd.DataFrame(
        [{"event_id": "POSseam", "chrom": "chr1", "start": 200, "end": 249, "strand": "+"}]
    )
    return universe, positives


def test_build_strata_output_feeds_scores_without_missing_columns(tmp_path: Path) -> None:
    samtools = _samtools()
    if samtools is None:
        pytest.skip("samtools required for FASTA extraction")
    from feng.bench.match import MatchConfig, build_strata

    universe, positives = _seam_universe()
    strata, summary = build_strata(
        universe, positives, MatchConfig(n_background=5, caliper_sd=5.0, min_background=1)
    )
    # positives appear in the strata, flagged, with coordinates
    assert int(strata["is_positive"].sum()) == 1
    anchor = strata[strata["is_positive"]].iloc[0]
    assert anchor["exon_key"] == "chr1:199-249:+"
    assert {"chrom", "exon_start_0base", "exon_end", "strand"} <= set(strata.columns)

    # a FASTA long enough to cover every exon +/- window
    fasta = tmp_path / "seam.fa"
    _write_indexed_fasta(fasta, {"chr1": "ACGT" * 300}, samtools)  # 1200 nt

    # baseline path: no missing-column error, score column attached, positive kept
    scored = baseline_hamming_scores(strata, fasta, window=50, samtools=samtools)
    assert "baseline_min_hamming" in scored.columns
    assert int(scored["is_positive"].sum()) == 1

    # IntaRNA-spec path: window specs build for every strata row (incl. positive)
    specs = exon_window_specs(strata, window=DEFAULT_WINDOW)
    assert len(specs) == len(strata)
    assert any(s.exon_key == "chr1:199-249:+" for s in specs)
