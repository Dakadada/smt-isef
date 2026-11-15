"""CPU infrastructure tests for the J4 in-process scale driver and completion gate.

StubScorer outputs are NEVER model evidence; these tests check runner mechanics
(byte-equivalence with the per-config CLI path, idempotence, failure isolation,
one model load per task, preemption flush, gate accounting) and never a
biological claim. Real-model equivalence is a post-submit check, see
ops/briefs/J4_DRIVER_RESULT.md.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pyarrow.parquet as pq
import pytest

from conference.model_v2 import j4_build_configs as jb
from conference.model_v2.tests.test_gpu_perturbation_atlas import StubScorer

JOBS = Path(__file__).resolve().parents[1] / "jobs"
REPO = Path(__file__).resolve().parents[3]
CYCLE = "AACCGGTTAACCGGTTAACCGGTT"


def load_script(name: str):
    """jobs/ holds scripts, not a package; load by path like the sbatch does."""
    path = JOBS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"j4_jobs_{name}", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def driver():
    return load_script("j4_scale_driver")


@pytest.fixture
def gate():
    return load_script("j4_scale_gate")


@pytest.fixture
def atlas():
    from conference.model_v2 import gpu_perturbation_atlas
    return gpu_perturbation_atlas


@pytest.fixture
def indexed_reference(tmp_path):
    # Index-only fixture for assembly validation, never model evidence.
    lengths = [248956422, 242193529, 198295559, 190214555, 181538259,
               170805979, 159345973, 145138636, 138394717, 133797422,
               135086622, 133275309, 114364328, 107043718, 101991189,
               90338345, 83257441, 80373285, 58617616, 64444167,
               46709983, 50818468, 156040895, 57227415, 16569]
    names = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chr5\nACGT\n")
    Path(str(fasta) + ".fai").write_text("".join(
        f"{name}\t{length}\t6\t4\t5\n" for name, length in zip(names, lengths)))
    return fasta


@pytest.fixture
def stubbed(atlas, monkeypatch, indexed_reference):
    """Deterministic position-dependent reference slices + a stub scorer."""
    monkeypatch.setattr(atlas.ti.FaiFastaReader, "fetch",
                        lambda self, chrom, start, end: "".join(
                            CYCLE[(p - 1) % len(CYCLE)] for p in range(start, end + 1)))
    monkeypatch.setattr(atlas, "PangolinScorer", StubScorer)
    return indexed_reference


def toy_config(config_id: str, center: int) -> dict:
    """A builder-canonical config small enough for StubScorer.context == 2."""
    return jb.make_config(config_id, f"toy_{center}", "chr5", "+", center, 2,
                          center, "donor", center - 2, center + 2, (2, 3), 2)


TOY = [("toyA__donor", 12), ("toyB__donor", 20), ("toyC__donor", 30)]
TOY_ROWS = 16  # 4 footprints x 1 site x 4 tissue channels


def write_shard(directory: Path, task_id: int = 0, configs=TOY) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"shard_{task_id:04d}.jsonl"
    path.write_text("".join(jb.canonical(toy_config(cid, center)) + "\n"
                            for cid, center in configs))
    return path


def run_driver(driver, shard_dir: Path, atlas_dir: Path, fasta, *, task_id=0,
               batch_size=2, extra=()):
    argv = ["--shard-dir", str(shard_dir), "--atlas-dir", str(atlas_dir),
            "--task-id", str(task_id), "--fasta", str(fasta),
            "--batch-size", str(batch_size), *extra]
    return driver.main(argv)


# --------------------------------------------------------------------------
# byte-equivalence with the per-config CLI path
# --------------------------------------------------------------------------

def normalize(value, root: Path) -> str:
    return json.dumps(value, sort_keys=True).replace(str(root.resolve()), "<ROOT>")


def test_driver_matches_per_config_cli_output(driver, atlas, stubbed, tmp_path):
    shard_dir, by_driver, by_cli = tmp_path / "shards", tmp_path / "drv", tmp_path / "cli"
    write_shard(shard_dir)
    assert run_driver(driver, shard_dir, by_driver, stubbed) == 0

    for config_id, center in TOY:
        config_path = tmp_path / f"{config_id}.json"
        config_path.write_text(jb.canonical(toy_config(config_id, center)))
        assert atlas.main(["--config", str(config_path), "--batch-size", "2",
                           "--fasta", str(stubbed),
                           "--output-dir", str(by_cli / config_id)]) == 0

    for config_id, _ in TOY:
        drv = atlas.load_manifest(by_driver / config_id / "shard_0.done")
        cli = atlas.load_manifest(by_cli / config_id / "shard_0.done")
        assert drv["run_id"] == cli["run_id"]
        assert drv["expected_rows"] == cli["expected_rows"] == TOY_ROWS
        # Same bytes in the published map, not merely the same row values.
        assert drv["ablation"]["sha256"] == cli["ablation"]["sha256"]
        assert normalize(drv, by_driver) == normalize(cli, by_cli)
        assert (pq.read_table(drv["ablation"]["path"]).to_pylist()
                == pq.read_table(cli["ablation"]["path"]).to_pylist())
        meta = pq.read_metadata(drv["ablation"]["path"]).metadata
        assert json.loads(meta[b"provenance"])["model_evidence"] is False
        rows = pq.read_table(drv["ablation"]["path"]).to_pylist()
        assert {r["evidence_class"] for r in rows} == {"sequence_ablation_proxy"}


def test_output_tree_and_config_json_match_shell_plan_step(driver, stubbed, tmp_path):
    shard_dir, atlas_dir = tmp_path / "shards", tmp_path / "atlas"
    write_shard(shard_dir)
    assert run_driver(driver, shard_dir, atlas_dir, stubbed) == 0
    for config_id, center in TOY:
        directory = atlas_dir / config_id
        text = jb.canonical(toy_config(config_id, center))
        assert (directory / "config.json").read_text() == text + "\n"
        assert (directory / ".j4_config_sha256").read_text() == jb.sha256_text(text)
        assert (directory / "ablation_maps" / "shard_0.parquet").stat().st_size > 0
        assert (directory / "reference_scores").is_dir()
        assert sorted(p.name for p in directory.glob("*.tmp.*")) == []
    assert not list(atlas_dir.glob("failures_task_*"))


def test_published_sha_matches_builder_manifest_column(driver, stubbed, tmp_path):
    """The marker the gate reads is the builder's config_sha256, not a new hash."""
    shard_dir, atlas_dir = tmp_path / "shards", tmp_path / "atlas"
    shard = write_shard(shard_dir)
    records = driver.read_plan(shard)
    line = shard.read_text().splitlines()[0]
    assert records[0]["sha256"] == hashlib.sha256(line.encode()).hexdigest()
    assert run_driver(driver, shard_dir, atlas_dir, stubbed) == 0
    assert (atlas_dir / records[0]["config_id"] / ".j4_config_sha256").read_text() \
        == jb.sha256_text(jb.canonical(records[0]["config"]))


# --------------------------------------------------------------------------
# one model load per task
# --------------------------------------------------------------------------

def test_scorer_constructed_once_for_whole_shard(driver, atlas, stubbed, tmp_path, monkeypatch):
    builds = []

    class Counted(StubScorer):
        def __init__(self, *a, **k):
            builds.append(1)
            super().__init__(*a, **k)

    monkeypatch.setattr(atlas, "PangolinScorer", Counted)
    shard_dir = tmp_path / "shards"
    write_shard(shard_dir)
    assert run_driver(driver, shard_dir, tmp_path / "atlas", stubbed) == 0
    assert builds == [1], "one Pangolin load per array task, not per config"


def test_default_batch_size_is_32(driver, atlas, stubbed, tmp_path, monkeypatch):
    seen = []
    original = atlas.run_atlas

    def recording(*args, **kwargs):
        seen.append(kwargs["batch_size"])
        return original(*args, **kwargs)

    monkeypatch.setattr(atlas, "run_atlas", recording)
    shard_dir = tmp_path / "shards"
    write_shard(shard_dir)
    assert driver.DEFAULT_BATCH_SIZE == 32
    assert driver.main(["--shard-dir", str(shard_dir), "--atlas-dir", str(tmp_path / "atlas"),
                        "--task-id", "0", "--fasta", str(stubbed)]) == 0
    assert seen == [32, 32, 32]


def test_task_id_defaults_to_slurm_array_task_id(driver, stubbed, tmp_path, monkeypatch):
    shard_dir = tmp_path / "shards"
    write_shard(shard_dir, task_id=7)
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "7")
    assert driver.main(["--shard-dir", str(shard_dir), "--atlas-dir", str(tmp_path / "atlas"),
                        "--fasta", str(stubbed), "--batch-size", "2"]) == 0
    assert (tmp_path / "atlas" / "toyA__donor" / "shard_0.done").is_file()
    monkeypatch.delenv("SLURM_ARRAY_TASK_ID")
    with pytest.raises(SystemExit):
        driver.main(["--shard-dir", str(shard_dir), "--atlas-dir", str(tmp_path / "atlas"),
                     "--fasta", str(stubbed)])


def test_ready_marker_published_once_before_work(driver, stubbed, tmp_path, monkeypatch):
    ready = tmp_path / "ready"
    monkeypatch.setenv("J4_READY_FILE", str(ready))
    shard_dir = tmp_path / "shards"
    write_shard(shard_dir)
    assert run_driver(driver, shard_dir, tmp_path / "atlas", stubbed) == 0
    assert ready.is_file()
    # A second invocation must get a FRESH marker dir; a stale one is an error.
    with pytest.raises(FileExistsError):
        run_driver(driver, shard_dir, tmp_path / "atlas2", stubbed)


# --------------------------------------------------------------------------
# idempotence
# --------------------------------------------------------------------------

def test_second_run_skips_without_scoring(driver, atlas, stubbed, tmp_path, monkeypatch):
    shard_dir, atlas_dir = tmp_path / "shards", tmp_path / "atlas"
    write_shard(shard_dir)
    assert run_driver(driver, shard_dir, atlas_dir, stubbed) == 0
    before = {cid: (atlas_dir / cid / "ablation_maps" / "shard_0.parquet").read_bytes()
              for cid, _ in TOY}

    class Refuse(StubScorer):
        def predict(self, sequences):
            raise AssertionError("completed config re-scored on resume")

    monkeypatch.setattr(atlas, "PangolinScorer", Refuse)
    assert run_driver(driver, shard_dir, atlas_dir, stubbed) == 0
    checkpoint = json.loads((atlas_dir / ".j4_driver" / "task_0.checkpoint.json").read_text())
    assert (checkpoint["skipped"], checkpoint["done"], checkpoint["failed"]) == (3, 0, 0)
    assert {cid: (atlas_dir / cid / "ablation_maps" / "shard_0.parquet").read_bytes()
            for cid, _ in TOY} == before


@pytest.mark.parametrize("damage", ["done", "parquet", "sha", "sha_mismatch", "empty_done"])
def test_skip_gate_requires_marker_parquet_and_matching_sha(driver, stubbed, tmp_path, damage):
    shard_dir, atlas_dir = tmp_path / "shards", tmp_path / "atlas"
    write_shard(shard_dir)
    assert run_driver(driver, shard_dir, atlas_dir, stubbed) == 0
    directory = atlas_dir / "toyA__donor"
    sha = driver.read_plan(shard_dir / "shard_0000.jsonl")[0]["sha256"]
    assert driver.already_done(directory, sha)
    if damage == "done":
        (directory / "shard_0.done").unlink()
    elif damage == "empty_done":
        (directory / "shard_0.done").write_text("")
    elif damage == "parquet":
        (directory / "ablation_maps" / "shard_0.parquet").unlink()
    elif damage == "sha":
        (directory / ".j4_config_sha256").unlink()
    else:
        (directory / ".j4_config_sha256").write_text("0" * 64)
    assert not driver.already_done(directory, sha)


@pytest.mark.parametrize("damage", ["done", "sha", "sha_mismatch"])
def test_rerun_republishes_a_missing_sha_marker(driver, stubbed, tmp_path, damage):
    shard_dir, atlas_dir = tmp_path / "shards", tmp_path / "atlas"
    write_shard(shard_dir)
    assert run_driver(driver, shard_dir, atlas_dir, stubbed) == 0
    directory = atlas_dir / "toyA__donor"
    if damage == "done":
        (directory / "shard_0.done").unlink()  # engine checkpoint still carries the work
    elif damage == "sha":
        (directory / ".j4_config_sha256").unlink()
    else:
        (directory / ".j4_config_sha256").write_text("0" * 64)
    assert run_driver(driver, shard_dir, atlas_dir, stubbed) == 0
    assert (directory / "shard_0.done").is_file()
    assert (directory / ".j4_config_sha256").read_text() == jb.sha256_text(
        jb.canonical(toy_config("toyA__donor", 12)))


def test_corrupt_published_parquet_is_a_failure_not_a_silent_rerun(driver, stubbed, tmp_path):
    """A .done citing a vanished artifact must not be papered over, as in v1."""
    shard_dir, atlas_dir = tmp_path / "shards", tmp_path / "atlas"
    write_shard(shard_dir)
    assert run_driver(driver, shard_dir, atlas_dir, stubbed) == 0
    (atlas_dir / "toyA__donor" / "ablation_maps" / "shard_0.parquet").unlink()
    assert run_driver(driver, shard_dir, atlas_dir, stubbed) == 1
    failures = (atlas_dir / "failures_task_0.tsv").read_text().splitlines()
    assert failures == ["toyA__donor\t1"]
    reason = json.loads((atlas_dir / "failures_task_0.jsonl").read_text().splitlines()[0])
    assert "artifact integrity failed" in reason["reason"]


# --------------------------------------------------------------------------
# failure isolation / fail-closed parsing
# --------------------------------------------------------------------------

def test_failed_config_is_explicit_and_does_not_strand_the_shard(driver, atlas, stubbed,
                                                                 tmp_path, monkeypatch):
    original = atlas.run_atlas

    def failing(windows, grid, genome, scorer, output_dir, **kwargs):
        if Path(output_dir).name == "toyB__donor":
            raise RuntimeError("engine raised for this locus")
        return original(windows, grid, genome, scorer, output_dir, **kwargs)

    monkeypatch.setattr(atlas, "run_atlas", failing)
    shard_dir, atlas_dir = tmp_path / "shards", tmp_path / "atlas"
    write_shard(shard_dir)
    assert run_driver(driver, shard_dir, atlas_dir, stubbed) == 1
    failures = (atlas_dir / "failures_task_0.tsv").read_text().splitlines()
    assert failures == ["toyB__donor\t1"], "two columns, exactly as the v1 shell wrote"
    detail = json.loads((atlas_dir / "failures_task_0.jsonl").read_text().splitlines()[0])
    assert detail["config_id"] == "toyB__donor" and detail["task_id"] == 0
    assert "engine raised for this locus" in detail["reason"]
    assert not (atlas_dir / "toyB__donor" / ".j4_config_sha256").exists()
    for config_id in ("toyA__donor", "toyC__donor"):
        assert (atlas_dir / config_id / "shard_0.done").is_file()
        assert (atlas_dir / config_id / ".j4_config_sha256").is_file()


def test_failed_rows_block_the_sha_marker(driver, atlas, stubbed, tmp_path, monkeypatch):
    """complete_with_failures is a failure, matching the CLI's nonzero exit."""
    class AlwaysFail(StubScorer):
        def predict(self, sequences):
            raise RuntimeError("test scoring failure")

    monkeypatch.setattr(atlas, "PangolinScorer", AlwaysFail)
    shard_dir, atlas_dir = tmp_path / "shards", tmp_path / "atlas"
    write_shard(shard_dir)
    assert run_driver(driver, shard_dir, atlas_dir, stubbed) == 1
    assert len((atlas_dir / "failures_task_0.tsv").read_text().splitlines()) == 3
    for config_id, _ in TOY:
        done = json.loads((atlas_dir / config_id / "shard_0.done").read_text())
        assert done["status"] == "complete_with_failures" and done["failures"] > 0
        assert not (atlas_dir / config_id / ".j4_config_sha256").exists()
        rows = pq.read_table(done["ablation"]["path"]).to_pylist()
        assert all(r["status"] == "failed" and r["model_evidence"] is False for r in rows)


def test_noncanonical_or_fenced_config_fails_only_itself(driver, stubbed, tmp_path):
    shard_dir, atlas_dir = tmp_path / "shards", tmp_path / "atlas"
    shard_dir.mkdir(parents=True)
    bad_fence = dict(toy_config("toyBad__donor", 40), assembly="GRCh37")
    noncanonical = json.dumps(toy_config("toyLoose__donor", 50))  # spaces, key order
    (shard_dir / "shard_0000.jsonl").write_text(
        jb.canonical(toy_config("toyA__donor", 12)) + "\n"
        + jb.canonical(bad_fence) + "\n" + noncanonical + "\n")
    assert run_driver(driver, shard_dir, atlas_dir, stubbed) == 1
    failed = dict(line.split("\t") for line in
                  (atlas_dir / "failures_task_0.tsv").read_text().splitlines())
    assert set(failed) == {"toyBad__donor", "toyLoose__donor"}
    reasons = " ".join((atlas_dir / "failures_task_0.jsonl").read_text().splitlines())
    assert "GRCh38 / GENCODE v50 assembly fence required" in reasons
    assert "not canonical" in reasons
    assert (atlas_dir / "toyA__donor" / ".j4_config_sha256").is_file()


@pytest.mark.parametrize("config_id", ["../escape", "a/b", ".hidden", "", None])
def test_unsafe_config_id_rejected(driver, tmp_path, config_id):
    shard = tmp_path / "shard_0000.jsonl"
    config = dict(toy_config("ok__donor", 12), config_id=config_id)
    shard.write_text(jb.canonical(config) + "\n")
    with pytest.raises(ValueError, match="unsafe config_id"):
        driver.read_plan(shard)


def test_plan_rejects_duplicates_and_empty_shards(driver, tmp_path):
    duplicate = tmp_path / "shard_0000.jsonl"
    line = jb.canonical(toy_config("toyA__donor", 12))
    duplicate.write_text(f"{line}\n{line}\n")
    with pytest.raises(ValueError, match="duplicate config_id"):
        driver.read_plan(duplicate)
    empty = tmp_path / "shard_0001.jsonl"
    empty.write_text("\n\n")
    with pytest.raises(ValueError, match="no configs"):
        driver.read_plan(empty)


def test_shard_file_name_is_zero_padded(driver, tmp_path):
    assert driver.shard_file_for(tmp_path, 318).name == "shard_0318.jsonl"
    with pytest.raises(ValueError, match="non-negative"):
        driver.shard_file_for(tmp_path, -1)


# --------------------------------------------------------------------------
# preemption: payload flushes, shell owns the requeue
# --------------------------------------------------------------------------

def test_driver_flush_hook_writes_checkpoint_without_requeue(driver, stubbed, tmp_path,
                                                             monkeypatch):
    captured = {}
    original = driver.ck.checkpoint_on_signal

    def recording(on_signal, **kwargs):
        captured.setdefault("calls", []).append(kwargs)
        captured.setdefault("hooks", []).append(on_signal)
        return original(on_signal, **kwargs)

    monkeypatch.setattr(driver.ck, "checkpoint_on_signal", recording)
    monkeypatch.setattr(driver.ck, "requeue_self",
                        lambda **kw: pytest.fail("payload must not requeue"))
    shard_dir, atlas_dir = tmp_path / "shards", tmp_path / "atlas"
    write_shard(shard_dir)
    assert run_driver(driver, shard_dir, atlas_dir, stubbed) == 0
    # Every layer (driver readiness, task-level flush, run_atlas per config)
    # must be flush-only; the shell is the sole requeue owner.
    assert all(kwargs["requeue"] is False for kwargs in captured["calls"])
    # hooks[0] = readiness no-op, hooks[1] = the task-level flush, hooks[2:] =
    # run_atlas's own per-config flush nested inside it (one per config).
    assert len(captured["hooks"]) == 2 + len(TOY)
    checkpoint = atlas_dir / ".j4_driver" / "task_0.checkpoint.json"
    checkpoint.unlink()
    captured["hooks"][1](15)  # the task-level handler, outside any config
    state = json.loads(checkpoint.read_text())
    assert state["configs"] == 3 and state["next"] == 3 and state["done"] == 3


def test_sigterm_mid_config_flushes_and_resumes(driver, stubbed, tmp_path):
    shard_dir, atlas_dir = tmp_path / "shards", tmp_path / "atlas"
    write_shard(shard_dir)
    fake_bin = tmp_path / "bin"
    fake_bin.mkdir()
    (fake_bin / "scontrol").write_text(
        f'#!/bin/bash\n: > "{tmp_path}/requeued"\n')
    (fake_bin / "scontrol").chmod(0o755)
    child = f'''
import os, signal, sys
sys.path.insert(0, {str(REPO)!r})
import importlib.util
from conference.model_v2 import gpu_perturbation_atlas as a
from conference.model_v2 import trans_intervention as ti
from conference.model_v2.tests.test_gpu_perturbation_atlas import StubScorer
CYCLE = {CYCLE!r}
ATLAS = {str(atlas_dir)!r}
ti.FaiFastaReader.fetch = lambda self, chrom, start, end: "".join(
    CYCLE[(p - 1) % len(CYCLE)] for p in range(start, end + 1))
class Terminated(StubScorer):
    def predict(self, sequences):
        # Preempt mid-SECOND-config, after config A is durable and config B has
        # already published one part, so the resume is a real partial resume.
        done_a = os.path.join(ATLAS, "toyA__donor", "shard_0.done")
        parts_b = os.path.join(ATLAS, "toyB__donor", "parts")
        if os.path.exists(done_a) and len(os.listdir(parts_b) if os.path.isdir(parts_b) else []) == 1:
            os.kill(os.getpid(), signal.SIGTERM)
        return super().predict(sequences)
a.PangolinScorer = Terminated
spec = importlib.util.spec_from_file_location("drv", {str(JOBS / "j4_scale_driver.py")!r})
drv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(drv)
raise SystemExit(drv.main(["--shard-dir", {str(shard_dir)!r}, "--atlas-dir", {str(atlas_dir)!r},
                           "--task-id", "0", "--fasta", {str(stubbed)!r}, "--batch-size", "1"]))
'''
    env = dict(os.environ, PATH=f"{fake_bin}:{os.environ['PATH']}",
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
    result = subprocess.run([sys.executable, "-c", child], env=env, cwd=str(REPO),
                            capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stderr[-2000:]
    assert not (tmp_path / "requeued").exists(), "shell owns the single requeue"
    first = atlas_dir / "toyA__donor"
    assert (first / "shard_0.done").is_file() and (first / ".j4_config_sha256").is_file()
    interrupted = atlas_dir / "toyB__donor"
    checkpoint = json.loads((interrupted / "shard_0.checkpoint.json").read_text())
    assert checkpoint["next"] == len(checkpoint["parts"]) == 1, "partial work survived"
    assert not (interrupted / "shard_0.done").exists()
    assert not (atlas_dir / "toyC__donor").exists(), "preempted before the last config"
    # The requeued task resumes: finished configs are skipped, the rest complete.
    assert run_driver(driver, shard_dir, atlas_dir, stubbed, batch_size=1) == 0
    for config_id, _ in TOY:
        done = json.loads((atlas_dir / config_id / "shard_0.done").read_text())
        assert done["status"] == "complete" and done["expected_rows"] == TOY_ROWS


# --------------------------------------------------------------------------
# completion gate
# --------------------------------------------------------------------------

def manifest_text(rows) -> str:
    header = ["# test manifest", "#" + "\t".join(jb.MANIFEST_COLUMNS)]
    body = ["\t".join(str(row.get(c, "NA")) for c in jb.MANIFEST_COLUMNS) for row in rows]
    return "\n".join(header + body) + "\n"


def toy_manifest(tmp_path: Path, extra=()) -> Path:
    rows = []
    for config_id, center in TOY:
        rows.append({"config_id": config_id, "status": "emitted", "rows": TOY_ROWS,
                     "shard": 0, "shard_file": "shards/shard_0000.jsonl",
                     "config_sha256": jb.sha256_text(jb.canonical(toy_config(config_id, center)))})
    rows.extend(extra)
    path = tmp_path / "j4_configs_manifest.tsv"
    path.write_text(manifest_text(rows))
    return path


def test_gate_passes_on_complete_atlas_and_ignores_skipped_rows(gate, driver, stubbed, tmp_path):
    shard_dir, atlas_dir = tmp_path / "shards", tmp_path / "atlas"
    write_shard(shard_dir)
    assert run_driver(driver, shard_dir, atlas_dir, stubbed) == 0
    manifest = toy_manifest(tmp_path, extra=[
        {"config_id": "toyZ__donor", "status": "skipped", "reason": "context_out_of_chrom"}])
    report = gate.gate(manifest, atlas_dir, deep=True)
    assert (report["emitted"], report["complete"], report["gaps"]) == (3, 3, 0)
    assert gate.main(["--manifest", str(manifest), "--atlas-dir", str(atlas_dir), "--deep"]) == 0


@pytest.mark.parametrize("damage,reason", [
    ("dir", "missing_config_dir"),
    ("done", "missing_done"),
    ("parquet", "missing_parquet"),
    ("sha", "missing_config_sha"),
    ("sha_mismatch", "config_sha_mismatch"),
])
def test_gate_reports_each_gap_and_exits_1(gate, driver, stubbed, tmp_path, damage, reason):
    shard_dir, atlas_dir = tmp_path / "shards", tmp_path / "atlas"
    write_shard(shard_dir)
    assert run_driver(driver, shard_dir, atlas_dir, stubbed) == 0
    manifest = toy_manifest(tmp_path)
    directory = atlas_dir / "toyC__donor"
    if damage == "dir":
        for path in sorted(directory.rglob("*"), reverse=True):
            path.rmdir() if path.is_dir() else path.unlink()
        directory.rmdir()
    elif damage == "done":
        (directory / "shard_0.done").unlink()
    elif damage == "parquet":
        (directory / "ablation_maps" / "shard_0.parquet").unlink()
    elif damage == "sha":
        (directory / ".j4_config_sha256").unlink()
    else:
        (directory / ".j4_config_sha256").write_text("0" * 64)
    report = gate.gate(manifest, atlas_dir)
    assert report["gaps"] == 1 and report["complete"] == 2
    assert report["gap_rows"][0] == {"config_id": "toyC__donor", "shard": "0", "reason": reason}
    assert gate.main(["--manifest", str(manifest), "--atlas-dir", str(atlas_dir)]) == 1


def test_gate_deep_mode_catches_row_count_and_integrity_drift(gate, driver, stubbed, tmp_path):
    shard_dir, atlas_dir = tmp_path / "shards", tmp_path / "atlas"
    write_shard(shard_dir)
    assert run_driver(driver, shard_dir, atlas_dir, stubbed) == 0
    rows = [{"config_id": cid, "status": "emitted", "rows": TOY_ROWS + 1, "shard": 0,
             "config_sha256": jb.sha256_text(jb.canonical(toy_config(cid, center)))}
            for cid, center in TOY]
    manifest = tmp_path / "wrong_rows.tsv"
    manifest.write_text(manifest_text(rows))
    assert gate.gate(manifest, atlas_dir, deep=False)["gaps"] == 0
    deep = gate.gate(manifest, atlas_dir, deep=True)
    assert deep["gaps"] == 3 and deep["gap_reasons"] == {"row_count_mismatch": 3}
    marker = atlas_dir / "toyA__donor" / "shard_0.done"
    forged = json.loads(marker.read_text())
    forged["expected_rows"] = TOY_ROWS + 1
    marker.write_text(json.dumps(forged))
    assert gate.check_config(atlas_dir, rows[0], deep=True).startswith("done_manifest_invalid")


def test_gate_rejects_manifest_shape_drift(gate, tmp_path):
    good = toy_manifest(tmp_path)
    lines = good.read_text().splitlines()
    short = tmp_path / "short.tsv"
    short.write_text("\n".join([lines[0], lines[1], "only\ttwo"]) + "\n")
    with pytest.raises(ValueError, match="expected 24 fields"):
        gate.read_manifest(short)
    renamed = tmp_path / "renamed.tsv"
    renamed.write_text("\n".join(
        [lines[0], "#config_id\tstatus\tbogus"] + lines[2:]) + "\n")
    with pytest.raises(ValueError, match="unexpected manifest columns"):
        gate.read_manifest(renamed)
    headless = tmp_path / "headless.tsv"
    headless.write_text("\n".join(lines[2:]) + "\n")
    with pytest.raises(ValueError, match="no column header"):
        gate.read_manifest(headless)


def test_gate_matches_real_manifest_columns(gate):
    assert len(jb.MANIFEST_COLUMNS) == 24
    assert jb.MANIFEST_COLUMNS[0] == "config_id" and jb.MANIFEST_COLUMNS[-1] == "config_sha256"


# --------------------------------------------------------------------------
# runner wiring
# --------------------------------------------------------------------------

def test_array2_sbatch_is_single_process_and_requeue_safe():
    text = (JOBS / "j4_scale_array2.sbatch").read_text()
    assert "#SBATCH --time=00:45:00" in text
    assert "#SBATCH --mem=32G" in text
    assert "#SBATCH --gres=gpu:1" in text
    assert "#SBATCH --partition=mit_preemptable" in text
    assert "#SBATCH --requeue" in text and "#SBATCH --signal=B:SIGTERM@60" in text
    assert "#SBATCH --array" not in text, "array size is set at submit time"
    assert 'J4_BATCH_SIZE:-32' in text
    # Exactly ONE driver invocation per array task (comments may cite it freely).
    body = "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("#"))
    assert body.count("j4_scale_driver.py") == 1
    assert "scontrol requeue" in body and body.count("scontrol requeue") == 1
    assert "gpu_perturbation_atlas" not in body, "no per-config CLI spawn in v2"
