"""Frozen sequence-score tests: Hamming baseline + IntaRNA wrapper."""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from feng.bench.scores import (
    DEFAULT_WINDOW,
    NUSINERSEN_5P_TO_3P,
    NUSINERSEN_RC,
    IntaRNAError,
    WindowSpec,
    baseline_hamming_scores,
    exon_window_specs,
    find_intarna,
    intarna_min_energy,
    join_score_tables,
    min_hamming_to_target,
    reverse_complement,
    score_cache_key,
    _fetch_transcript_sequences,
    _fetch_with_pysam,
    _fetch_with_samtools,
    _pad_window_to_length,
)
from feng.bench.strata import GenomeFastaError
from feng.bench.strata import FreezeGuardError, UNIVERSE_PROVENANCE

ENV_PREFIX = Path("/home/dashii/orcd/pool/feng-env")
SCORE_SBATCH = Path(__file__).resolve().parents[2] / "feng" / "scripts" / "bench_intarna_scores.sbatch"


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
            {"exon_key": "chr1:20-38:+", "chrom": "chr1", "exon_start_0base": 20, "exon_end": 38, "strand": "+", **UNIVERSE_PROVENANCE},
            {"exon_key": "chr2:15-33:-", "chrom": "chr2", "exon_start_0base": 15, "exon_end": 33, "strand": "-", **UNIVERSE_PROVENANCE},
        ]
    )
    scored = baseline_hamming_scores(events, fasta, window=0, samtools=samtools)
    assert scored["baseline_min_hamming"].tolist() == [0, 0]


# --- contig-edge N-padding (Arm E hoff off-target windows) -------------------
#
# Off-target context windows near the edge of a SHORT scaffold run off the
# contig (start < 1 or end > contig_length). The fetch must clamp to the contig
# and N-pad the missing flank(s) back to the exact requested length -- the
# Pangolin/SpliceAI convention for chromosome ends -- rather than raising and
# dropping the site (which would change the per-mer burden denominator).

_EDGE_CONTIG = "ACGTACGTAC" * 3  # 30 nt, len == 30


def _pysam_available() -> bool:
    try:
        import pysam  # noqa: F401
    except ImportError:
        return False
    return True


def test_pad_window_to_length_left_right_and_both_edges() -> None:
    seq = _EDGE_CONTIG
    n = len(seq)  # 30
    # left overflow only: window [-5, 15) -> "NNNNN" + seq[0:15]
    left = _pad_window_to_length(seq[0:15], -5, 15, n)
    assert left == "N" * 5 + seq[0:15]
    assert len(left) == 20
    # right overflow only: window [20, 40) -> seq[20:30] + "N"*10
    right = _pad_window_to_length(seq[20:30], 20, 40, n)
    assert right == seq[20:30] + "N" * 10
    assert len(right) == 20
    # both edges: window [-4, 34) -> "NNNN" + seq + "NNNN"
    both = _pad_window_to_length(seq, -4, 34, n)
    assert both == "N" * 4 + seq + "N" * 4
    assert len(both) == 38
    # no overflow: identity
    assert _pad_window_to_length(seq[5:25], 5, 25, n) == seq[5:25]


def test_pad_window_to_length_rejects_non_overlapping_window() -> None:
    # A window entirely off the contig is a genuinely bad spec, not an edge
    # condition -- fail closed instead of returning all-N.
    with pytest.raises(GenomeFastaError, match="does not overlap"):
        _pad_window_to_length("", 40, 60, 30)


def _edge_specs(n: int) -> list[WindowSpec]:
    """A start<1 case and an end>contig_length case on a length-``n`` contig."""
    return [
        WindowSpec(exon_key="left", chrom="scaf1", start_0base=-4, end=12, strand="+"),
        WindowSpec(exon_key="right", chrom="scaf1", start_0base=n - 8, end=n + 6, strand="+"),
        WindowSpec(exon_key="minus", chrom="scaf1", start_0base=-3, end=n + 3, strand="-"),
    ]


def _assert_edge_fetch(out: dict[str, str], seq: str) -> None:
    n = len(seq)
    # left: [-4, 12) -> requested 16, "NNNN" + seq[0:12]
    assert out["left"] == "N" * 4 + seq[0:12]
    assert len(out["left"]) == 16
    assert out["left"][4] == seq[0]  # true center offset preserved
    # right: [n-8, n+6) -> requested 14, seq[n-8:n] + "N"*6
    assert out["right"] == seq[n - 8 : n] + "N" * 6
    assert len(out["right"]) == 14
    # minus: genomic [-3, n+3) padded then reverse-complemented
    genomic = "N" * 3 + seq + "N" * 3
    assert out["minus"] == reverse_complement(genomic)
    assert len(out["minus"]) == n + 6


def test_fetch_with_samtools_n_pads_contig_edges(tmp_path: Path) -> None:
    samtools = _samtools()
    if samtools is None:
        pytest.skip("samtools required for FASTA extraction")
    seq = _EDGE_CONTIG
    fasta = tmp_path / "scaf.fa"
    _write_indexed_fasta(fasta, {"scaf1": seq}, samtools)
    out = _fetch_with_samtools(_edge_specs(len(seq)), fasta, samtools)
    _assert_edge_fetch(out, seq)


def test_fetch_with_pysam_n_pads_contig_edges(tmp_path: Path) -> None:
    if not _pysam_available():
        pytest.skip("pysam required for this backend")
    seq = _EDGE_CONTIG
    fasta = tmp_path / "scaf.fa"
    fasta.write_text(f">scaf1\n{seq}\n")  # pysam auto-builds the .fai on open
    out = _fetch_with_pysam(_edge_specs(len(seq)), fasta)
    _assert_edge_fetch(out, seq)


def test_fetch_transcript_sequences_edge_end_to_end(tmp_path: Path) -> None:
    # The real dispatcher: pysam-preferred, samtools fallback. Whichever backend
    # is present in this env must produce the N-padded edge windows.
    samtools = _samtools()
    if samtools is None and not _pysam_available():
        pytest.skip("neither pysam nor samtools available")
    seq = _EDGE_CONTIG
    fasta = tmp_path / "scaf.fa"
    if samtools is not None:
        _write_indexed_fasta(fasta, {"scaf1": seq}, samtools)
    else:
        fasta.write_text(f">scaf1\n{seq}\n")
    out = _fetch_transcript_sequences(_edge_specs(len(seq)), fasta, samtools=samtools)
    _assert_edge_fetch(out, seq)


def test_offtarget_edge_window_preserves_genomic_padding_and_footprint_offset(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    from feng.titration import gpu_scores as gs

    samtools = _samtools()
    if samtools is None and not _pysam_available():
        pytest.skip("neither pysam nor samtools available")
    contig = "CAAAAGTC"
    fasta = tmp_path / "short-scaffold.fa"
    if samtools is not None:
        _write_indexed_fasta(fasta, {"scaf1": contig}, samtools)
    else:
        fasta.write_text(f">scaf1\n{contig}\n")
    fake_gpu_python = tmp_path / "python"
    fake_gpu_python.write_text("")
    captured: list[dict] = []

    def worker(jobs: list[dict], d: int, gpu_python: Path) -> list[dict]:
        captured.extend(jobs)
        return [{"id": job["id"], "max_gain": 0.0, "max_loss": 0.0} for job in jobs]

    monkeypatch.setattr(gs, "_run_worker", worker)
    result = gs.pangolin_off_target_burden(
        "AAAA",
        [{"chrom": "scaf1", "start": 2, "end": 5, "strand": "-", "target_5to3": "TTTT"}],
        genome_fasta=fasta,
        gpu_python=fake_gpu_python,
        d=1,
    )

    assert result["pangolin_off_target_incomplete"] is False
    assert len(captured) == 1
    job = captured[0]
    window_start_1b, window_end_1b = gs.window_bounds_1b(3, d=1)
    left_pad = 1 - window_start_1b
    right_pad = window_end_1b - len(contig)
    expected_genomic_window = "N" * left_pad + contig + "N" * right_pad
    assert len(job["ref_seq"]) == window_end_1b - window_start_1b + 1
    assert job["ref_seq"] == expected_genomic_window
    footprint_offset = 2 - window_start_1b
    assert job["ref_seq"][footprint_offset : footprint_offset + 4] == "AAAA"
    assert job["alt_seq"][footprint_offset : footprint_offset + 4] == "NNNN"
    assert job["alt_seq"][:footprint_offset] == job["ref_seq"][:footprint_offset]
    assert job["alt_seq"][footprint_offset + 4 :] == job["ref_seq"][footprint_offset + 4 :]
    assert job["hybridization_query_5to3_provenance"] == "TTTT"


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


def test_score_cache_key_changes_with_input_content(tmp_path: Path) -> None:
    events = tmp_path / "events.parquet"
    fasta = tmp_path / "genome.fa"
    events.write_bytes(b"events-v1")
    fasta.write_bytes(b">chr1\nAAAA\n")
    params = {"window": 100, "query": NUSINERSEN_5P_TO_3P}
    first = score_cache_key({"events": events, "fasta": fasta}, params)
    events.write_bytes(b"events-v2")
    assert score_cache_key({"events": events, "fasta": fasta}, params) != first


def test_score_cache_key_is_location_independent_and_role_framed(tmp_path: Path) -> None:
    left = tmp_path / "left"
    right = tmp_path / "right"
    left.mkdir()
    right.mkdir()
    (left / "a").write_bytes(b"same")
    (right / "renamed").write_bytes(b"same")
    assert score_cache_key({"events": left / "a"}, {}) == score_cache_key(
        {"events": right / "renamed"}, {}
    )
    assert score_cache_key({"events": left / "a"}, {}) != score_cache_key(
        {"reference": left / "a"}, {}
    )


def _score_join_frames():
    keys = {
        "positive_event_id": ["P1", "P1"],
        "positive_exon_key": ["p", "p"],
        "is_positive": [True, False],
        "exon_key": ["p", "b"],
    }
    provenance = {key: [value, value] for key, value in UNIVERSE_PROVENANCE.items()}
    baseline = pd.DataFrame({**keys, **provenance, "baseline_min_hamming": [1, 2]})
    test = pd.DataFrame({**keys, **provenance, "intarna_min_energy": [-4.0, -2.0]})
    return baseline, test


def test_score_join_is_exactly_one_to_one() -> None:
    baseline, test = _score_join_frames()
    joined = join_score_tables(baseline, test.iloc[::-1])
    assert len(joined) == 2
    assert joined.loc[joined["exon_key"] == "p", "intarna_min_energy"].item() == -4.0


def test_score_join_rejects_duplicate_or_missing_keys() -> None:
    baseline, test = _score_join_frames()
    with pytest.raises(ValueError, match="duplicate join keys"):
        join_score_tables(pd.concat([baseline, baseline.iloc[[0]]]), test)
    with pytest.raises(ValueError, match="non-identical keys"):
        join_score_tables(baseline, test.iloc[[0]])


def test_score_publication_switches_one_atomic_manifest_symlink() -> None:
    script = SCORE_SBATCH.read_text()
    assert "parquet_sha256" in script
    assert 'os.replace(tmp_name, manifest_path)' in script
    assert 'mktemp "${SHARD_DIR}/.intarna_${shard}.tmp.XXXXXXXX.parquet"' in script
    assert 'flock -x 9' in script
    assert 'validate_shard "${shard_out}" "${shard}"' in script
    assert '9>"${CACHE_DIR}/publish.lock"' in script
    assert 'os.replace(candidate, parquet)' in script
    assert 'os.replace(tmp, current)' in script
    assert '.current.${SLURM_JOB_ID:-manual}' not in script
    assert 'cp "${CACHE_DIR}/scored_strata.parquet"' not in script


def test_score_array_assigns_one_logical_shard_per_task_and_finalizes_once_complete() -> None:
    script = SCORE_SBATCH.read_text()
    assert "#SBATCH --array=0-31%16" in script
    assert 'SHARD_FIRST="${SLURM_ARRAY_TASK_ID}"' in script
    assert 'SHARD_STEP="${SLURM_ARRAY_TASK_COUNT}"' in script
    assert 'for (( shard=SHARD_FIRST; shard<N_SHARDS; shard+=SHARD_STEP )); do' in script
    assert 'echo "checkpoint complete; waiting for remaining shards"' in script
    assert script.index('waiting for remaining shards') < script.index('scored_tmp="$(mktemp')


def test_score_cache_binds_python_scorer_and_samtools_implementation() -> None:
    script = SCORE_SBATCH.read_text()
    assert '"samtools_executable": Path(sys.argv[4])' in script
    assert '"scorer_source": Path(scores.__file__)' in script
    assert '"cache_schema": 2' in script


def test_score_cache_bootstrap_heredoc_executes(tmp_path: Path) -> None:
    """Execute the script's cache-key Python, not merely inspect its text."""
    script = SCORE_SBATCH.read_text()
    match = re.search(
        r"CACHE_KEY=.*?<<'PYEOF'\n(?P<body>.*?)\nPYEOF\n\)\"",
        script,
        flags=re.DOTALL,
    )
    assert match is not None
    inputs = []
    for name, content in (
        ("strata.parquet", b"strata"),
        ("genome.fa", b">chr1\nA\n"),
        ("IntaRNA", b"intarna"),
        ("samtools", b"samtools"),
    ):
        path = tmp_path / name
        path.write_bytes(content)
        inputs.append(str(path))
    result = subprocess.run(
        [sys.executable, "-c", match.group("body"), *inputs, "100", "32"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert re.fullmatch(r"[0-9a-f]{64}", result.stdout.strip())


# --- match -> scores seam (end to end) ---------------------------------------


def _seam_universe():
    """Small universe with covariate spread and modest coordinates for the seam."""
    import pandas as pd

    input_provenance = {
        "gencode_gtf_sha256": "b003af56d9e8d210a11095891e9291779dbe9647bc9ca55df743110abe5b84d3",
        "scrmoe_sj_sha256": '{"SRR35298699":"f9e66850810b97749709e891fb175a11b19b350f57f833e462a18d3e51213d30","SRR35298704":"5e8fff9e575e947724a0814468dbc6ea9433c0e11979c06d5794bb85c3cc0cef","SRR35298715":"77f9e8b34348b52f09f05c5edf139084c7b0278255440bbb04514ee86873f321"}',
        "genome_fasta_sha256": "8937614b4e206071d1668802ba883ea7d509ef88998c5f5316755fe9fe897632",
        "scrmoe_reads_per_gene_sha256": '{"SRR35298699":"2bd726e264d4a2472058b046e0ba093834729baaea9a831ae8d9351b110ca14c","SRR35298704":"a0a27928874f30a0bbb935c08fa6a955ec1d5b5557318a97308644d5768203b8","SRR35298715":"06e7cb98a595664bf5cb8118b952e1514e2375027423c0bf80a9612e0dae15c4"}',
        "junction_read_length": 100,
        "junction_read_length_sha256": "eea8254c7500ba3de996aa8ad6af399183f04e17d4a8102fde539dbc93a90012",
    }
    rows = []
    # positive at chr1 1-based 200..249 -> exon_key chr1:199-249:+
    def row(start0, end, psi, jr, cpm, gc):
        return {
            "exon_key": f"chr1:{start0}-{end}:+", "chrom": "chr1",
            "exon_start_0base": start0, "exon_end": end, "strand": "+",
            "detectable": True, "baseline_psi": psi,
            "mean_total_junction_reads": jr, "control_mean_cpm": cpm,
            "exon_length": end - start0, "gc_fraction": gc,
            **UNIVERSE_PROVENANCE,
            **input_provenance,
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


def test_scoring_rejects_events_without_scrmoe_provenance(tmp_path: Path) -> None:
    events = pd.DataFrame([{
        "exon_key": "chr1:0-4:+", "chrom": "chr1",
        "exon_start_0base": 0, "exon_end": 4, "strand": "+",
    }])
    with pytest.raises(FreezeGuardError, match="ScrMOE provenance"):
        baseline_hamming_scores(events, tmp_path / "unused.fa", window=0)
