"""Tests for Arm E axis H (Pangolin deep-splice) -- feng.titration.gpu_scores.

The module must import cleanly under the CPU feng-env (torch/pangolin are only
touched inside the GPU worker). Pure-logic tests (mask/coords/readout) run
everywhere; scorer-gating tests prove no fabrication when the GPU env is absent;
the real Pangolin smoke is opt-in (FENG_RUN_PANGOLIN=1 + importable GPU env).
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from feng.titration import gpu_scores as gs


# ---------------------------------------------------------------- pure helpers
def test_window_bounds_length_and_center():
    d, ctx = 50, gs.PANGOLIN_CONTEXT
    start, end = gs.window_bounds_1b(70076575, d=d, context=ctx)
    assert end - start + 1 == 2 * (ctx + d) + 1
    # center sits exactly in the middle
    assert 70076575 - start == ctx + d


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


def test_offtarget_burden_none_without_sites():
    out = gs.pangolin_offtarget_burden("TCACTTTCATAATGCTGG", [])
    assert out["pangolin_offtarget_burden"] is None
    assert out["pangolin_offtarget_incomplete"] is True   # missing scan != zero burden


def test_offtarget_burden_incomplete_when_env_python_absent():
    out = gs.pangolin_offtarget_burden(
        "TCACTTTCATAATGCTGG", [{"chrom": "chr5", "start": 1000, "end": 1017}],
        gpu_python=Path("/nonexistent/bin/python"),
    )
    assert out["pangolin_offtarget_burden"] is None and out["pangolin_offtarget_incomplete"] is True


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
