"""Tests for Arm E axis H (Pangolin deep-splice) -- feng.titration.gpu_scores.

The module must import cleanly under the CPU feng-env (torch/pangolin are only
touched inside the GPU worker). Pure-logic tests (mask/coords/readout) run
everywhere; scorer-gating tests prove no fabrication when the GPU env is absent;
the real Pangolin smoke is opt-in (FENG_RUN_PANGOLIN=1 + importable GPU env).
"""

from __future__ import annotations

import os
import hashlib
import json
from pathlib import Path
import shutil
import sys

import pytest

from feng.titration import gpu_scores as gs


def _tree_sha(entries):
    digest = hashlib.sha256()
    for entry in sorted(entries, key=lambda item: item["path"]):
        digest.update(entry["path"].encode("utf-8"))
        digest.update(b"\0")
        digest.update(entry["git_blob_sha256"].encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


# ---------------------------------------------------------------- pure helpers
def test_window_bounds_length_and_center():
    d, ctx = 50, gs.PANGOLIN_CONTEXT
    start, end = gs.window_bounds_1b(gs.EXON7_LAST_BASE_1B, d=d, context=ctx)
    assert end - start + 1 == 2 * (ctx + d) + 1
    # center sits exactly in the middle
    assert gs.EXON7_LAST_BASE_1B - start == ctx + d


def test_footprint_mask_places_Ns_and_preserves_length():
    ref = "ACGT" * 10  # len 40, window_start 1-based = 101 -> covers 101..140
    win_start = 101
    # mask genomic 111..115 (5 bases) -> offset 10..14
    out = gs.footprint_mask(ref, 111, 115, win_start)
    assert len(out) == len(ref)
    assert out[10:15] == "NNNNN"
    assert out[:10] == ref[:10] and out[15:] == ref[15:]


def test_footprint_mask_is_deterministic():
    ref = "GGGGCCCCAAAATTTT"
    a = gs.footprint_mask(ref, 5, 8, 1)
    b = gs.footprint_mask(ref, 5, 8, 1)
    assert a == b == "GGGGNNNNAAAATTTT"


def test_footprint_mask_rejects_out_of_window():
    ref = "ACGTACGT"
    with pytest.raises(ValueError):
        gs.footprint_mask(ref, 100, 105, 1)  # past the end
    with pytest.raises(ValueError):
        gs.footprint_mask(ref, 0, 2, 1)      # before the start (offset < 0)
    with pytest.raises(ValueError):
        gs.footprint_mask(ref, 5, 4, 1)      # end precedes start


def test_readout_index_center_and_bounds():
    assert gs.readout_index(1000, 1000, d=50) == 50   # center -> d
    assert gs.readout_index(1000, 990, d=50) == 40
    assert gs.readout_index(1000, 1010, d=50) == 60
    with pytest.raises(ValueError):
        gs.readout_index(1000, 1100, d=50)            # outside +/- d


def test_anchor_constants_are_the_verified_frame():
    # exon-7 donor is one base past the last exon base; ISS-N1 target is the 18-mer
    assert gs.EXON7_DONOR_1B == gs.EXON7_LAST_BASE_1B + 1
    assert gs.ISS_N1_END_1B - gs.ISS_N1_START_1B + 1 == 18


def test_exon7_default_readout_uses_last_exonic_base(monkeypatch, tmp_path):
    fake_py = tmp_path / "python"
    fake_py.write_text("")
    observed = {}

    def fetch(center_1b, d, genome_fasta):
        observed["center_1b"] = center_1b
        return "A" * 10101, center_1b - gs.PANGOLIN_CONTEXT - d

    def worker(jobs, d, gpu_python):
        observed["readout_idx"] = jobs[0]["readout_idx"]
        return [{"id": 18, "delta_at_readout": 0.1,
                 "max_gain": 0.1, "max_gain_offset": 0}]

    monkeypatch.setattr(gs, "_fetch_window", fetch)
    monkeypatch.setattr(gs, "_run_worker", worker)
    gs.pangolin_exon7_delta(
        [{"length": 18, "start": 70076584, "end": 70076601}],
        gpu_python=fake_py,
    )
    assert observed == {
        "center_1b": gs.EXON7_LAST_BASE_1B,
        "readout_idx": gs.DEFAULT_DISTANCE,
    }


def test_pangolin_outputs_name_masking_as_n_ablation_proxy():
    on_target = gs.pangolin_exon7_delta(
        [{"length": 18, "start": 70076584, "end": 70076601}],
        gpu_python=Path("/nonexistent/bin/python"),
    )[0]
    off_target = gs.pangolin_off_target_burden("AAAA", [])
    assert on_target["masking_method"] == "N-ablation sensitivity proxy"
    assert off_target["masking_method"] == "N-ablation sensitivity proxy"
    assert on_target["state"] == off_target["state"] == "candidate"


# --------------------------------------------------- gating: never fabricate
def test_gpu_ready_false_when_env_absent(monkeypatch):
    def fake_check(prefix=None, **kw):
        return {"env_exists": False, "packages": {}}
    monkeypatch.setattr("feng.gpu.env.check_gpu_env", fake_check)
    ready, why = gs._gpu_ready(Path("/nope/bin/python"))
    assert ready is False and "absent" in why.lower()


def test_gpu_ready_false_when_pangolin_missing(monkeypatch):
    def fake_check(prefix=None, **kw):
        return {"env_exists": True, "packages": {"torch": True, "pangolin": False}}
    monkeypatch.setattr("feng.gpu.env.check_gpu_env", fake_check)
    ready, why = gs._gpu_ready(Path("/x/bin/python"))
    assert ready is False and "pangolin" in why.lower()


def test_exon7_delta_incomplete_when_env_python_absent():
    # Cheap gate: no fabricated score when the GPU env python doesn't exist.
    rows = gs.pangolin_exon7_delta(
        [{"length": 18, "start": 70076584, "end": 70076601},
         {"length": 16, "start": 70076585, "end": 70076600}],
        gpu_python=Path("/nonexistent/bin/python"),
    )
    assert len(rows) == 2
    for r in rows:
        assert r["pangolin_exon7_delta"] is None
        assert r["pangolin_incomplete"] is True
        assert r["state"] == "candidate"
        assert "absent" in r["reason"].lower()


def test_exon7_delta_incomplete_on_worker_import_error(monkeypatch, tmp_path):
    # AUTHORITATIVE path: env python exists but the worker's real import fails ->
    # incomplete with the real reason, never a fabricated number.
    fake_py = tmp_path / "python"
    fake_py.write_text("")  # exists() is True
    monkeypatch.setattr(gs, "_fetch_window", lambda *a, **k: ("A" * 10101, 70071525))

    def boom(jobs, d, gpu_python):
        raise RuntimeError("worker error: gpu import failed: ModuleNotFoundError: No module named 'torch'")

    monkeypatch.setattr(gs, "_run_worker", boom)
    rows = gs.pangolin_exon7_delta([{"length": 18, "start": 70076584, "end": 70076601}], gpu_python=fake_py)
    assert rows[0]["pangolin_exon7_delta"] is None
    assert rows[0]["pangolin_incomplete"] is True
    assert "import failed" in rows[0]["reason"]


def test_gpu_worker_reauthenticates_source_tree_before_import(monkeypatch, tmp_path):
    source_root = tmp_path / "snapshot"
    package_root = source_root / "feng/src/feng"
    (package_root / "titration").mkdir(parents=True)
    repo_package = Path(gs.__file__).resolve().parents[1]
    copied = [
        (repo_package / "__init__.py", package_root / "__init__.py"),
        (repo_package / "source_auth.py", package_root / "source_auth.py"),
        (Path(gs.__file__), package_root / "titration/gpu_scores.py"),
        (Path(gs.__file__).parent / "__init__.py", package_root / "titration/__init__.py"),
    ]
    entries = []
    for source, destination in copied:
        shutil.copyfile(source, destination)
        raw = destination.read_bytes()
        entries.append({
            "path": destination.relative_to(source_root).as_posix(),
            "git_blob_sha256": hashlib.sha256(raw).hexdigest(),
        })
    commit = "b" * 40
    tree_sha = _tree_sha(entries)
    entries.sort(key=lambda item: item["path"])
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({
        "schema_version": 1,
        "executed_git_commit": commit,
        "source_tree_sha256": tree_sha,
        "entries": entries,
    }))
    monkeypatch.setenv("FENG_ARME_SOURCE_ROOT", str(source_root))
    monkeypatch.setenv("FENG_ARME_SOURCE_MANIFEST", str(manifest))
    monkeypatch.setenv("FENG_ARME_SOURCE_TREE_SHA256", tree_sha)
    monkeypatch.setenv("FENG_ARME_SEALED_COMMIT", commit)
    (package_root / "titration/gpu_scores.py").write_bytes(b"tampered\n")

    with pytest.raises(RuntimeError, match="source authentication failed.*gpu_scores.py"):
        gs._run_worker([], gs.DEFAULT_DISTANCE, Path(sys.executable))


def test_offtarget_burden_none_without_sites():
    out = gs.pangolin_off_target_burden("TCACTTTCATAATGCTGG", [])
    assert out["pangolin_off_target_burden"] is None
    assert out["pangolin_off_target_incomplete"] is True   # missing scan != zero burden


def test_offtarget_burden_incomplete_when_env_python_absent():
    out = gs.pangolin_off_target_burden(
        "TCACTTTCATAATGCTGG", [{"chrom": "chr5", "start": 1000, "end": 1017}],
        gpu_python=Path("/nonexistent/bin/python"),
    )
    assert out["pangolin_off_target_burden"] is None and out["pangolin_off_target_incomplete"] is True


@pytest.mark.parametrize(
    ("strand", "genomic_target"), [("+", "TTTT"), ("-", "AAAA")]
)
def test_offtarget_honours_aso_orientation_and_target_strand(
    monkeypatch, tmp_path, strand, genomic_target
):
    fake_py = tmp_path / "python"
    fake_py.write_text("")
    ref = list("C" * 101)
    ref[50:54] = genomic_target
    monkeypatch.setattr(gs, "_fetch_window_generic", lambda *a, **k: ("".join(ref), 50))
    captured = []

    def worker(jobs, d, gpu_python):
        captured.extend(jobs)
        return [{"id": 0, "max_gain": 0.25, "max_loss": -0.1}]

    monkeypatch.setattr(gs, "_run_worker", worker)
    out = gs.pangolin_off_target_burden(
        "AAAA",
        [{"chrom": "chr1", "start": 100, "end": 103, "strand": strand,
          "target_5to3": "TTTT"}],
        gpu_python=fake_py,
    )
    assert out["pangolin_off_target_burden"] == pytest.approx(0.35)
    assert captured[0]["strand"] == strand
    assert captured[0]["hybridization_query_5to3_provenance"] == "TTTT"


def test_offtarget_missing_worker_value_is_incomplete_not_zero(monkeypatch, tmp_path):
    fake_py = tmp_path / "python"
    fake_py.write_text("")
    ref = "C" * 50 + "TTTT" + "C" * 47
    monkeypatch.setattr(gs, "_fetch_window_generic", lambda *a, **k: (ref, 50))
    monkeypatch.setattr(gs, "_run_worker", lambda *a, **k: [{"id": 0, "max_gain": None, "max_loss": -0.1}])
    out = gs.pangolin_off_target_burden(
        "AAAA",
        [{"chrom": "chr1", "start": 100, "end": 103, "strand": "+",
          "target_5to3": "TTTT"}],
        gpu_python=fake_py,
    )
    assert out["pangolin_off_target_burden"] is None
    assert out["pangolin_off_target_incomplete"] is True


def test_offtarget_burden_counts_only_directional_excursions(monkeypatch, tmp_path):
    fake_py = tmp_path / "python"
    fake_py.write_text("")
    ref = "C" * 50 + "TTTT" + "C" * 47
    monkeypatch.setattr(gs, "_fetch_window_generic", lambda *a, **k: (ref, 50))
    monkeypatch.setattr(
        gs, "_run_worker",
        lambda *a, **k: [{"id": 0, "max_gain": 0.2, "max_loss": 0.1}],
    )
    out = gs.pangolin_off_target_burden(
        "AAAA",
        [{"chrom": "chr1", "start": 100, "end": 103, "strand": "+",
          "target_5to3": "TTTT"}],
        gpu_python=fake_py,
    )
    assert out["pangolin_off_target_burden"] == pytest.approx(0.2)


def test_offtarget_indel_orientation_is_checked(monkeypatch, tmp_path):
    fake_py = tmp_path / "python"
    fake_py.write_text("")
    # ASO AAAA hybridizes in TTTT orientation; TTT is a one-base deletion hit.
    ref = "C" * 50 + "TTT" + "C" * 48
    monkeypatch.setattr(gs, "_fetch_window_generic", lambda *a, **k: (ref, 50))
    monkeypatch.setattr(
        gs, "_run_worker",
        lambda *a, **k: [{"id": 0, "max_gain": 0.2, "max_loss": -0.1}],
    )
    out = gs.pangolin_off_target_burden(
        "AAAA",
        [{"chrom": "chr1", "start": 100, "end": 102, "strand": "+",
          "target_5to3": "TTT"}],
        gpu_python=fake_py,
    )
    assert out["pangolin_off_target_burden"] == pytest.approx(0.3)


def test_offtarget_palindromic_orientation_tie_is_scored(monkeypatch, tmp_path):
    fake_py = tmp_path / "python"
    fake_py.write_text("")
    target = "CACTTTATAAAGTG"  # chr12:115999899-115999912:- (GRCh38)
    ref = "C" * 50 + target + "C" * 37
    monkeypatch.setattr(
        gs, "_fetch_window_generic", lambda *a, **k: (ref, 115999849)
    )
    monkeypatch.setattr(
        gs,
        "_run_worker",
        lambda *a, **k: [{"id": 0, "max_gain": 0.2, "max_loss": -0.1}],
    )
    out = gs.pangolin_off_target_burden(
        "CACTTTCATAATGCTG",
        [{"chrom": "chr12", "start": 115999899, "end": 115999912, "strand": "-",
          "target_5to3": target}],
        gpu_python=fake_py,
    )
    assert out["pangolin_off_target_burden"] == pytest.approx(0.3)
    assert out["pangolin_off_target_incomplete"] is False


def test_offtarget_nonpalindromic_orientation_tie_is_scored(monkeypatch, tmp_path):
    from feng.scan.seed_extend import _alignment_cost

    fake_py = tmp_path / "python"
    fake_py.write_text("")
    target = "TAAT"
    ref = "C" * 50 + target + "C" * 47
    captured = []
    monkeypatch.setattr(gs, "_fetch_window_generic", lambda *a, **k: (ref, 50))

    def worker(jobs, d, gpu_python):
        captured.extend(jobs)
        return [{"id": 0, "max_gain": 0.2, "max_loss": -0.1}]

    monkeypatch.setattr(gs, "_run_worker", worker)
    assert _alignment_cost("TTTT", target)[0] == 2
    assert _alignment_cost("AAAA", target)[0] == 2
    out = gs.pangolin_off_target_burden(
        "AAAA",
        [{"chrom": "chr1", "start": 100, "end": 103, "strand": "+",
          "target_5to3": target}],
        gpu_python=fake_py,
    )
    assert len(captured) == 1
    assert out["pangolin_off_target_incomplete"] is False


def test_offtarget_direct_aso_orientation_remains_incomplete(monkeypatch, tmp_path):
    fake_py = tmp_path / "python"
    fake_py.write_text("")
    ref = "C" * 50 + "AAAA" + "C" * 47
    monkeypatch.setattr(gs, "_fetch_window_generic", lambda *a, **k: (ref, 50))
    monkeypatch.setattr(gs, "_run_worker", lambda *a, **k: pytest.fail("worker must not run"))
    out = gs.pangolin_off_target_burden(
        "AAAA",
        [{"chrom": "chr1", "start": 100, "end": 103, "strand": "+",
          "target_5to3": "AAAA"}],
        gpu_python=fake_py,
    )
    assert out["pangolin_off_target_burden"] is None
    assert out["pangolin_off_target_incomplete"] is True
    assert "direct-ASO orientation" in out["reason"]


# ------------------------------------------------- real Pangolin smoke (opt-in)
# NOTE: the gate is a CHEAP env-flag check only -- the actual GPU readiness probe
# (which subprocess-imports torch+pangolin, ~25s) is done INSIDE the test so it
# never runs at collection time.
@pytest.mark.skipif(
    os.environ.get("FENG_RUN_PANGOLIN") != "1",
    reason="set FENG_RUN_PANGOLIN=1 (with an importable GPU env) to run the real Pangolin smoke",
)
def test_real_pangolin_exon7_delta_nusinersen_sign():
    from feng.titration.ladder import build_ladder

    ready, why = gs._gpu_ready(gs.DEFAULT_GPU_PYTHON)
    if not ready:
        pytest.skip(f"GPU env not ready: {why}")
    rungs = [r for r in build_ladder(genome_fasta=gs.DEFAULT_GENOME_FASTA) if r["length"] in (12, 18, 24)]
    rows = gs.pangolin_exon7_delta(rungs)
    by_len = {r["length"]: r for r in rows}
    assert by_len[18]["pangolin_exon7_delta"] is not None
    # masking ISS-N1 should not LOSE exon-7 inclusion (expected sign >= ~0)
    assert by_len[18]["pangolin_exon7_delta"] >= -0.05
