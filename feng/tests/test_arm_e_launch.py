"""Arm E transcriptome scan and H-off launch contracts."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
from contextlib import contextmanager
import fcntl
import shutil
import subprocess
import sys
from types import SimpleNamespace

import pytest

from feng.titration import arm_e_scan as scan
from feng.titration import arm_e_hoff as hoff
from feng.titration import gpu_scores as gs


_REAL_PANGOLIN_PROVENANCE = hoff._pangolin_provenance


@pytest.fixture(autouse=True)
def _stable_pangolin_provenance(monkeypatch):
    provenance = {
        "pangolin_package_version": "test-version",
        "pangolin_package_sha256": "b" * 64,
        "pangolin_package_pinned_sha256": "b" * 64,
        "pangolin_package_manifest": {
            "file_count": 1,
            "files": [{"path": "model.py", "sha256": "e" * 64}],
        },
        "pangolin_weight_sha256": {
            f"final.{j}.{i}.3.v2": f"{i}{j}" * 32
            for i in (0, 2, 4, 6) for j in (1, 2, 3)
        },
    }
    monkeypatch.setattr(
        hoff, "_pangolin_provenance", lambda _gpu_python: provenance,
        raising=False,
    )
    monkeypatch.setattr(
        hoff, "_source_digest_bindings", lambda: {
            "pangolin_distance": 50,
            "pangolin_context": 5000,
            "source_tree_sha256": "c" * 64,
            "executed_git_commit": "a" * 40,
        },
        raising=False,
    )
    monkeypatch.setenv("FENG_ARME_RUN_ID", "20260910T120000Z")
    monkeypatch.setenv("SLURM_ARRAY_JOB_ID", "12345")
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", "1")
    monkeypatch.setenv("FENG_ARME_LAUNCHER_BLOB_SHA256", "d" * 64)
    monkeypatch.setenv("FENG_ARME_SEALED_COMMIT", "a" * 40)


def _write_scan_burden(manifest: Path) -> Path:
    rows = list(scan.read_site_manifest(manifest))
    assert rows
    burden = manifest.with_suffix(".burden.json")
    burden.write_text(json.dumps({
        "schema_version": scan.SCAN_ARTIFACT_SCHEMA_VERSION,
        "algorithm_version": scan.SCAN_ALGORITHM_VERSION,
        "terminal": True,
        "complete": True,
        "input_digest": "0" * 64,
        "reference_sha256": "1" * 64,
        "manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
        "length": rows[0]["length"],
        "aso_5to3": rows[0]["aso_5to3"],
        "off_target_hits": len(rows),
    }))
    return burden


def test_gtf_exons_are_deduplicated_sorted_and_stranded(tmp_path: Path) -> None:
    gtf = tmp_path / "tiny.gtf"
    gtf.write_text(
        "chr2\ts\texon\t5\t8\t.\t-\t.\tgene_id \"g2\";\n"
        "chr1\ts\tgene\t1\t9\t.\t+\t.\tgene_id \"g\";\n"
        "chr1\ts\texon\t2\t4\t.\t+\t.\tgene_id \"g1\";\n"
        "chr1\ts\texon\t2\t4\t.\t+\t.\tgene_id \"g1b\";\n"
    )
    assert scan.load_unique_exons(gtf) == [
        scan.Exon("chr1", 2, 4, "+"),
        scan.Exon("chr2", 5, 8, "-"),
    ]


def test_gtf_union_merges_same_strand_despite_interleaved_opposite_strand(tmp_path: Path) -> None:
    gtf = tmp_path / "overlap.gtf"
    gtf.write_text(
        "chr1\ts\texon\t1\t10\t.\t+\t.\tx \"1\";\n"
        "chr1\ts\texon\t3\t4\t.\t-\t.\tx \"2\";\n"
        "chr1\ts\texon\t5\t12\t.\t+\t.\tx \"3\";\n"
    )
    assert scan.load_unique_exons(gtf) == [
        scan.Exon("chr1", 1, 12, "+"), scan.Exon("chr1", 3, 4, "-")
    ]


def test_build_reference_pins_sources_and_output(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nAACCGGTT\n")
    scan.write_fai(fasta)
    gtf = tmp_path / "a.gtf"
    gtf.write_text(
        "chr1\ts\texon\t2\t4\t.\t+\t.\tx \"1\";\n"
        "chr1\ts\texon\t5\t8\t.\t-\t.\tx \"2\";\n"
    )
    out = tmp_path / "ref.fa"
    pin = tmp_path / "ref.pin.json"
    scan.build_reference(gtf=gtf, genome_fasta=fasta, out=out, pin=pin)
    assert out.read_text() == ">chr1:2-4:+\nACC\n>chr1:5-8:-\nAACC\n"
    payload = json.loads(pin.read_text())
    assert payload["reference_path"] == str(out.resolve())
    assert payload["reference_sha256"] == hashlib.sha256(out.read_bytes()).hexdigest()
    assert payload["record_count"] == 2
    scan.validate_pinned_reference(out, pin)


def test_validate_pinned_reference_rejects_tamper(tmp_path: Path) -> None:
    ref = tmp_path / "r.fa"
    ref.write_text(">chr1:1-2:+\nAA\n")
    pin = tmp_path / "p.json"
    pin.write_text(json.dumps({"reference_path": str(ref.resolve()), "reference_sha256": hashlib.sha256(ref.read_bytes()).hexdigest()}))
    ref.write_text(ref.read_text() + "A")
    with pytest.raises(ValueError, match="SHA-256"):
        scan.validate_pinned_reference(ref, pin)


def test_scan_rung_emits_genomic_coords_and_excludes_intended_target(tmp_path: Path) -> None:
    ref = tmp_path / "r.fa"
    # ASO A*8 searches for T*8. Minus record is transcript-oriented sequence.
    ref.write_text(">chr1:101-110:+\nCTTTTTTTTC\n>chr2:201-210:-\nCTTTTTTTTC\n>chr5:300-307:+\nTTTTTTTT\n")
    pin = tmp_path / "p.json"
    pin.write_text(json.dumps({"reference_path": str(ref.resolve()), "reference_sha256": hashlib.sha256(ref.read_bytes()).hexdigest()}))
    manifest = tmp_path / "sites.tsv"
    burden = tmp_path / "burden.json"
    scan.scan_rung(
        {"length": 8, "aso_5to3": "AAAAAAAA", "chrom": "chr5", "start": 300, "end": 307},
        reference=ref, pin=pin, manifest=manifest, burden_out=burden,
        checkpoint=tmp_path / "c.json", max_edit=0,
    )
    rows = list(scan.read_site_manifest(manifest))
    assert [(r["chrom"], r["start"], r["end"], r["strand"], r["target_5to3"]) for r in rows] == [
        ("chr1", 102, 109, "+", "TTTTTTTT"),
        ("chr2", 202, 209, "-", "TTTTTTTT"),
    ]
    assert json.loads(burden.read_text())["off_target_hits"] == 2


def test_length_aware_edit_policy_marks_short_mers_uninformative() -> None:
    assert scan.effective_edit_policy(2, 5)["informative"] is False
    assert scan.effective_edit_policy(4, 5)["informative"] is False
    assert scan.effective_edit_policy(6, 5)["informative"] is False
    assert scan.effective_edit_policy(8, 5)["effective_max_edit"] == 1
    assert scan.effective_edit_policy(18, 5)["effective_max_edit"] == 3
    assert scan.effective_edit_policy(30, 5)["effective_max_edit"] == 5


def test_uninformative_short_rung_is_terminal_incomplete_not_zero(tmp_path: Path) -> None:
    ref = tmp_path / "r.fa"
    ref.write_text(">chr1:1-8:+\nTTTTTTTT\n")
    pin = tmp_path / "p.json"
    pin.write_text(json.dumps({
        "reference_path": str(ref.resolve()),
        "reference_sha256": hashlib.sha256(ref.read_bytes()).hexdigest(),
    }))
    manifest, burden = tmp_path / "sites.tsv", tmp_path / "burden.json"
    result = scan.scan_rung(
        {"length": 4, "aso_5to3": "AAAA", "chrom": "chr5", "start": 1, "end": 4},
        reference=ref, pin=pin, manifest=manifest, burden_out=burden,
        checkpoint=tmp_path / "c.json", max_edit=5,
    )
    assert list(scan.read_site_manifest(manifest)) == []
    assert result["off_target_hits"] is None
    assert result["complete"] is False
    assert result["terminal"] is True
    assert result["reason"] == "length_below_minimum_informative_length"


def test_scan_counts_each_genomic_site_once_within_one_exonic_locus(
    monkeypatch, tmp_path: Path,
) -> None:
    ref = tmp_path / "r.fa"
    ref.write_text(">chr1:101-109:+\nTTTTTTTTT\n")
    pin = tmp_path / "p.json"
    pin.write_text(json.dumps({
        "reference_path": str(ref.resolve()),
        "reference_sha256": hashlib.sha256(ref.read_bytes()).hexdigest(),
    }))
    manifest, burden = tmp_path / "sites.tsv", tmp_path / "burden.json"
    def hits(*_args, **_kwargs):
        yield scan.HybridizationHit("worse", 0, 8, 1, 1)
        yield scan.HybridizationHit("best", 0, 8, 0, 0)
        yield scan.HybridizationHit("second-site", 1, 9, 0, 0)

    monkeypatch.setattr(scan, "seed_and_extend", hits)
    scan.scan_rung(
        {"length": 8, "aso_5to3": "AAAAAAAA", "chrom": "chr5", "start": 1, "end": 8},
        reference=ref, pin=pin, manifest=manifest, burden_out=burden,
        checkpoint=tmp_path / "c.json", max_edit=5,
    )
    rows = list(scan.read_site_manifest(manifest))
    assert [(row["start"], row["end"], row["edit_distance"]) for row in rows] == [
        (101, 108, 0), (102, 109, 0),
    ]
    assert json.loads(burden.read_text())["off_target_hits"] == 2


def test_scan_rejects_stale_checkpoint_digest(tmp_path: Path) -> None:
    ref = tmp_path / "r.fa"
    ref.write_text(">chr1:1-8:+\nTTTTTTTT\n")
    pin = tmp_path / "p.json"
    pin.write_text(json.dumps({
        "reference_path": str(ref.resolve()),
        "reference_sha256": hashlib.sha256(ref.read_bytes()).hexdigest(),
    }))
    checkpoint = tmp_path / "c.json"
    checkpoint.write_text(json.dumps({"record_index": 0, "input_digest": "stale"}))
    with pytest.raises(ValueError, match="checkpoint input/config digest"):
        scan.scan_rung(
            {"length": 8, "aso_5to3": "AAAAAAAA", "chrom": "chr5", "start": 1, "end": 8},
            reference=ref, pin=pin, manifest=tmp_path / "sites.tsv",
            burden_out=tmp_path / "burden.json", checkpoint=checkpoint, max_edit=5,
        )


def test_scan_signal_checkpoint_uses_last_committed_snapshot(monkeypatch, tmp_path: Path) -> None:
    ref = tmp_path / "r.fa"
    ref.write_text(">chr1:1-8:+\nTTTTTTTT\n")
    pin = tmp_path / "p.json"
    pin.write_text(json.dumps({
        "reference_path": str(ref.resolve()),
        "reference_sha256": hashlib.sha256(ref.read_bytes()).hexdigest(),
    }))
    callback = {}

    @contextmanager
    def capture(on_signal, **_kwargs):
        callback["fn"] = on_signal
        yield

    def interrupted_hits(*_args, **_kwargs):
        yield scan.HybridizationHit("x", 0, 8, 0, 0)
        callback["fn"](15)
        raise RuntimeError("simulated eviction")

    monkeypatch.setattr(scan, "checkpoint_on_signal", capture)
    monkeypatch.setattr(scan, "seed_and_extend", interrupted_hits)
    checkpoint = tmp_path / "c.json"
    with pytest.raises(RuntimeError, match="simulated eviction"):
        scan.scan_rung(
            {"length": 8, "aso_5to3": "AAAAAAAA", "chrom": "chr5", "start": 20, "end": 27},
            reference=ref, pin=pin, manifest=tmp_path / "sites.tsv",
            burden_out=tmp_path / "burden.json", checkpoint=checkpoint, max_edit=0,
        )
    saved = json.loads(checkpoint.read_text())
    assert saved["record_index"] == 0
    assert saved["off_target_hits"] == 0


def test_reference_builder_signal_does_not_commit_written_unindexed_record(monkeypatch, tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nAACCGGTT\n")
    scan.write_fai(fasta)
    gtf = tmp_path / "a.gtf"
    gtf.write_text(
        "chr1\ts\texon\t1\t2\t.\t+\t.\tx \"1\";\n"
        "chr1\ts\texon\t5\t6\t.\t+\t.\tx \"2\";\n"
    )
    callback, calls = {}, 0

    @contextmanager
    def capture(on_signal, **_kwargs):
        callback["fn"] = on_signal
        yield

    real_fetch = scan._fetch

    def interrupted_fetch(*args):
        nonlocal calls
        calls += 1
        if calls == 2:
            callback["fn"](15)
            raise RuntimeError("simulated eviction")
        return real_fetch(*args)

    monkeypatch.setattr(scan, "checkpoint_on_signal", capture)
    monkeypatch.setattr(scan, "_fetch", interrupted_fetch)
    checkpoint = tmp_path / "ref.ckpt.json"
    with pytest.raises(RuntimeError, match="simulated eviction"):
        scan.build_reference(
            gtf=gtf, genome_fasta=fasta, out=tmp_path / "ref.fa",
            pin=tmp_path / "pin.json", checkpoint=checkpoint,
        )
    saved = json.loads(checkpoint.read_text())
    assert saved["record_index"] == 0
    assert saved["output_offset"] == 0


def test_launch_scripts_are_requeue_safe_and_arm_e_specific() -> None:
    root = Path(__file__).resolve().parents[1]
    cpu = (root / "scripts" / "arm_e_scan.sbatch").read_text()
    gpu = (root / "scripts" / "arm_e_hoff.sbatch").read_text()
    template = (root / "scripts" / "templates" / "gpu_preemptable.sbatch").read_text()
    ref = (root / "scripts" / "arm_e_reference.sbatch").read_text()
    for text in (cpu, gpu, ref):
        assert "--requeue" in text
        assert ".done" in text
        assert "mit_preemptable" in text
    assert "feng.titration.arm_e_scan scan" in cpu
    assert "feng.titration.arm_e_hoff" in gpu
    assert "__PAYLOAD_CMD__" not in gpu
    assert "--gres=gpu:1" in gpu
    assert "__PAYLOAD_CMD__" in template
    assert "feng.titration.arm_e_hoff" not in template


def _sealed_launcher_repo(root: Path, tmp_path: Path) -> tuple[Path, str, str]:
    repo = tmp_path / "sealed-repo"
    shutil.copytree(
        root / "src", repo / "feng" / "src",
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    scripts = repo / "feng" / "scripts"
    scripts.mkdir(parents=True)
    launcher = root / "scripts" / "arm_e_hoff.sbatch"
    shutil.copyfile(launcher, scripts / launcher.name)
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    subprocess.run(["git", "-C", str(repo), "add", "feng"], check=True)
    subprocess.run(
        [
            "git", "-C", str(repo), "-c", "user.name=Test",
            "-c", "user.email=test@example.invalid", "commit", "-qm", "sealed fixture",
        ],
        check=True,
    )
    commit = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    launcher_sha = hashlib.sha256(launcher.read_bytes()).hexdigest()
    return repo, commit, launcher_sha


def _launcher_env(
    root: Path, tmp_path: Path, *, repo: Path, commit: str, launcher_sha: str,
) -> dict[str, str]:
    env = dict(os.environ)
    for key in tuple(env):
        if key.startswith("BASH_FUNC_conda"):
            env.pop(key)
    env.update({
        "FENG_ARME_SEALED_COMMIT": commit,
        "FENG_ARME_LAUNCHER_BLOB_SHA256": launcher_sha,
        "FENG_ARME_RUN_ID": "20260910T120000Z",
        "FENG_GPU_ENV_PREFIX": str(Path(sys.executable).resolve().parents[1]),
        "SCRATCH": str(tmp_path / "scratch"),
        "SLURM_ARRAY_JOB_ID": "12345",
        "SLURM_ARRAY_TASK_ID": "1",
        "SLURM_JOB_ID": "12345_1",
        "SLURM_SUBMIT_DIR": str(repo),
    })
    return env


def test_hoff_launcher_rejects_running_launcher_blob_mismatch(
    tmp_path: Path,
) -> None:
    root = Path(__file__).resolve().parents[1]
    repo, commit, _launcher_sha = _sealed_launcher_repo(root, tmp_path)
    env = _launcher_env(
        root, tmp_path, repo=repo, commit=commit, launcher_sha="0" * 64,
    )
    completed = subprocess.run(
        ["bash", str(root / "scripts" / "arm_e_hoff.sbatch")],
        env=env, check=False, capture_output=True, text=True,
    )
    assert completed.returncode != 0
    assert "launcher SHA-256 mismatch" in completed.stderr
    assert "starting payload" not in completed.stdout


def test_hoff_launcher_rejects_expected_hash_not_bound_to_commit(
    tmp_path: Path,
) -> None:
    root = Path(__file__).resolve().parents[1]
    repo, _commit, launcher_sha = _sealed_launcher_repo(root, tmp_path)
    committed_launcher = repo / "feng/scripts/arm_e_hoff.sbatch"
    committed_launcher.write_bytes(committed_launcher.read_bytes() + b"# different blob\n")
    subprocess.run(["git", "-C", str(repo), "add", str(committed_launcher)], check=True)
    subprocess.run(
        [
            "git", "-C", str(repo), "-c", "user.name=Test",
            "-c", "user.email=test@example.invalid", "commit", "-qm", "different launcher",
        ],
        check=True,
    )
    commit = subprocess.run(
        ["git", "-C", str(repo), "rev-parse", "HEAD"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    env = _launcher_env(
        root, tmp_path, repo=repo, commit=commit, launcher_sha=launcher_sha,
    )

    completed = subprocess.run(
        ["bash", str(root / "scripts" / "arm_e_hoff.sbatch")],
        env=env, check=False, capture_output=True, text=True,
    )

    assert completed.returncode != 0
    assert "launcher hash is not the sealed commit blob" in completed.stderr
    assert "starting payload" not in completed.stdout


def test_hoff_launcher_materializes_and_authenticates_whole_committed_tree(
    tmp_path: Path,
) -> None:
    root = Path(__file__).resolve().parents[1]
    repo, commit, launcher_sha = _sealed_launcher_repo(root, tmp_path)
    env = _launcher_env(
        root, tmp_path, repo=repo, commit=commit, launcher_sha=launcher_sha,
    )
    fake_bin = _fake_hpc_commands(tmp_path)
    env["PATH"] = str(fake_bin) + os.pathsep + env["PATH"]
    completed = subprocess.run(
        ["bash", str(root / "scripts" / "arm_e_hoff.sbatch")],
        env=env, check=False, capture_output=True, text=True,
    )
    assert completed.returncode != 0  # controlled fixture has no scan inputs
    source_parent = (
        Path(env["SCRATCH"]) / "feng" / "arm_e" / "source"
        / env["FENG_ARME_SEALED_COMMIT"] / env["SLURM_ARRAY_JOB_ID"]
        / env["SLURM_ARRAY_TASK_ID"]
    )
    snapshots = sorted(source_parent.glob("snapshot.*"))
    assert len(snapshots) == 1
    manifest = snapshots[0] / "source-manifest.json"
    payload = json.loads(manifest.read_text())
    committed_paths = subprocess.run(
        ["git", "-C", str(repo), "ls-tree", "-r", "--name-only", commit, "--", "feng/src"],
        check=True, capture_output=True, text=True,
    ).stdout.splitlines()
    assert [entry["path"] for entry in payload["entries"]] == committed_paths
    assert len(payload["source_tree_sha256"]) == 64


def _fake_hpc_commands(tmp_path: Path) -> Path:
    fake_bin = tmp_path / "fake-bin"
    fake_bin.mkdir(exist_ok=True)
    module = fake_bin / "module"
    module.write_text("#!/usr/bin/env bash\nexit 0\n")
    module.chmod(0o755)
    conda_base = tmp_path / "conda-base"
    conda_init = conda_base / "etc" / "profile.d" / "conda.sh"
    conda_init.parent.mkdir(parents=True, exist_ok=True)
    conda_init.write_text("conda() { return 0; }\n")
    conda = fake_bin / "conda"
    conda.write_text(
        "#!/usr/bin/env bash\n"
        f"if [[ $1 == info ]]; then printf '%s\\n' '{conda_base}'; fi\n"
    )
    conda.chmod(0o755)
    return fake_bin


def test_hoff_supersession_recovers_mid_move_then_atomically_prepares(
    tmp_path: Path,
) -> None:
    root = Path(__file__).resolve().parents[1]
    repo, commit, launcher_sha = _sealed_launcher_repo(root, tmp_path)
    env = _launcher_env(
        root, tmp_path, repo=repo, commit=commit, launcher_sha=launcher_sha,
    )
    fake_bin = _fake_hpc_commands(tmp_path)
    env["PATH"] = str(fake_bin) + os.pathsep + env["PATH"]
    out_dir = tmp_path / "outputs"
    out_dir.mkdir()
    env["FENG_GPU_OUT_DIR"] = str(out_dir)
    stale = {
        "pangolin_hoff_L02.json": b"old-output\n",
        "pangolin_hoff_L02.ckpt.json": b"old-checkpoint\n",
    }
    for name, raw in stale.items():
        (out_dir / name).write_bytes(raw)
    superseded = out_dir / "superseded" / env["FENG_ARME_RUN_ID"] / "L02"
    superseded.mkdir(parents=True)
    (superseded / "pangolin_hoff_L02.done").write_bytes(b"")

    completed = subprocess.run(
        ["bash", str(root / "scripts" / "arm_e_hoff.sbatch")],
        env=env, check=False, capture_output=True, text=True,
    )
    assert completed.returncode != 0  # controlled fixture has no scan inputs
    for name, raw in stale.items():
        assert not (out_dir / name).exists()
        assert (superseded / name).read_bytes() == raw
    assert not (out_dir / "pangolin_hoff_L02.done").exists()
    assert (out_dir / ".prepared-20260910T120000Z-L02").is_file()


def test_hoff_supersession_lock_rejects_duplicate_length(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parents[1]
    repo, commit, launcher_sha = _sealed_launcher_repo(root, tmp_path)
    env = _launcher_env(
        root, tmp_path, repo=repo, commit=commit, launcher_sha=launcher_sha,
    )
    env["PATH"] = str(_fake_hpc_commands(tmp_path)) + os.pathsep + env["PATH"]
    out_dir = tmp_path / "outputs"
    lock_path = out_dir / ".locks" / "pangolin_hoff_L02.lock"
    lock_path.parent.mkdir(parents=True)
    env["FENG_GPU_OUT_DIR"] = str(out_dir)
    with lock_path.open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        completed = subprocess.run(
            ["bash", str(root / "scripts" / "arm_e_hoff.sbatch")],
            env=env, check=False, capture_output=True, text=True,
        )
    assert completed.returncode != 0
    assert "another H-off task owns length L02" in completed.stderr
    assert not (out_dir / ".prepared-20260910T120000Z-L02").exists()


def test_hoff_submit_wrapper_streams_committed_launcher_and_exports_its_hash(
    tmp_path: Path,
) -> None:
    root = Path(__file__).resolve().parents[1]
    repo, commit, launcher_sha = _sealed_launcher_repo(root, tmp_path)
    fake_bin = tmp_path / "submit-bin"
    fake_bin.mkdir()
    captured_script = tmp_path / "submitted.sbatch"
    captured_args = tmp_path / "sbatch.args"
    sbatch = fake_bin / "sbatch"
    sbatch.write_text(
        "#!/usr/bin/env bash\n"
        "printf '%s\\n' \"$@\" >\"$CAPTURED_ARGS\"\n"
        "cp /dev/stdin \"$CAPTURED_SCRIPT\"\n"
    )
    sbatch.chmod(0o755)
    env = dict(os.environ)
    env.update({
        "PATH": str(fake_bin) + os.pathsep + env["PATH"],
        "CAPTURED_ARGS": str(captured_args),
        "CAPTURED_SCRIPT": str(captured_script),
        "SCRATCH": "/orcd/scratch/orcd/011/dashii",
    })

    completed = subprocess.run(
        [
            str(root / "scripts" / "submit_arm_e_hoff.sh"), commit,
            "20260910T120000Z",
        ],
        cwd=repo, env=env, check=False, capture_output=True, text=True,
    )

    assert completed.returncode == 0
    committed_launcher = subprocess.run(
        ["git", "-C", str(repo), "show", f"{commit}:feng/scripts/arm_e_hoff.sbatch"],
        check=True, capture_output=True,
    ).stdout
    assert captured_script.read_bytes() == committed_launcher
    exports = captured_args.read_text()
    assert f"FENG_ARME_SEALED_COMMIT={commit}" in exports
    assert f"FENG_ARME_LAUNCHER_BLOB_SHA256={launcher_sha}" in exports
    assert "FENG_ARME_RUN_ID=20260910T120000Z" in exports


def test_hoff_accumulates_chunks_and_checkpoints(monkeypatch, tmp_path: Path) -> None:
    manifest = tmp_path / "sites.tsv"
    manifest.write_text(
        "\t".join(scan.SITE_COLUMNS) + "\n" +
        "4\tAAAA\tchr1\t10\t13\t+\tTTTT\tr1\t0\t0\n" +
        "4\tAAAA\tchr1\t20\t23\t+\tTTTT\tr2\t0\t0\n" +
        "4\tAAAA\tchr1\t30\t33\t+\tTTTT\tr3\t0\t0\n"
    )
    calls = []
    def fake(aso, sites, **kwargs):
        calls.append(len(sites))
        return {"pangolin_off_target_burden": float(len(sites)), "reason": None}
    monkeypatch.setattr(hoff, "pangolin_off_target_burden", fake)
    out = hoff.run_sweep(manifest=manifest, out=tmp_path / "out.json",
                         checkpoint=tmp_path / "ckpt.json", chunk_size=2,
                         scan_burden=_write_scan_burden(manifest))
    assert calls == [2, 1]
    assert out["pangolin_off_target_burden"] == 3.0
    assert json.loads((tmp_path / "ckpt.json").read_text())["site_index"] == 3


def test_hoff_rejects_unauthenticated_nonempty_manifest_before_pangolin(
    monkeypatch, tmp_path: Path,
) -> None:
    manifest = tmp_path / "sites.tsv"
    manifest.write_text(
        "\t".join(scan.SITE_COLUMNS) + "\n"
        "8\tAAAAAAAA\tchr1\t1\t8\t+\tTTTTTTTT\tr1\t0\t0\n"
    )
    burden = _write_scan_burden(manifest)
    payload = json.loads(burden.read_text())
    payload["manifest_sha256"] = "f" * 64
    burden.write_text(json.dumps(payload))
    called = False

    def forbidden(*_args, **_kwargs):
        nonlocal called
        called = True
        raise AssertionError("Pangolin must not run")

    monkeypatch.setattr(hoff, "pangolin_off_target_burden", forbidden)
    out = tmp_path / "out.json"
    done = tmp_path / "out.done"
    with pytest.raises(ValueError, match="scan manifest SHA-256 mismatch"):
        hoff.run_sweep(
            manifest=manifest, scan_burden=burden, out=out,
            checkpoint=tmp_path / "ckpt.json",
        )
    assert called is False
    assert not out.exists()
    assert not done.exists()
    completed = subprocess.run(
        [
            sys.executable, "-m", "feng.titration.arm_e_hoff",
            "--manifest", str(manifest), "--scan-burden", str(burden),
            "--out", str(out), "--checkpoint", str(tmp_path / "ckpt.json"),
        ],
        check=False, capture_output=True, text=True,
    )
    assert completed.returncode != 0
    assert "scan manifest SHA-256 mismatch" in completed.stderr
    assert not done.exists()


def test_hoff_distinguishes_complete_zero_sites_from_uninformative_short_rung(tmp_path: Path) -> None:
    ref = tmp_path / "r.fa"
    ref.write_text(">chr1:1-8:+\nCCCCCCCC\n")
    pin = tmp_path / "p.json"
    pin.write_text(json.dumps({
        "reference_path": str(ref.resolve()),
        "reference_sha256": hashlib.sha256(ref.read_bytes()).hexdigest(),
    }))
    for length, expected_incomplete in ((8, False), (4, True)):
        manifest = tmp_path / f"sites_L{length:02d}.tsv"
        burden = tmp_path / f"burden_L{length:02d}.json"
        scan.scan_rung(
            {"length": length, "aso_5to3": "A" * length,
             "chrom": "chr5", "start": 20, "end": 19 + length},
            reference=ref, pin=pin, manifest=manifest, burden_out=burden,
            checkpoint=tmp_path / f"scan_L{length:02d}.json", max_edit=0,
        )
        result = hoff.run_sweep(
            manifest=manifest, scan_burden=burden,
            out=tmp_path / f"hoff_L{length:02d}.json",
            checkpoint=tmp_path / f"hoff_L{length:02d}.ckpt.json",
            genome_fasta=ref,
        )
        assert result["pangolin_off_target_incomplete"] is expected_incomplete
        assert result["pangolin_off_target_burden"] == (None if expected_incomplete else 0.0)
        assert result["reason"] == (
            "length_below_minimum_informative_length"
            if expected_incomplete else "completed scan found zero off-target sites"
        )


def test_hoff_authenticates_live_genome_before_checkpoint_reuse(tmp_path: Path) -> None:
    manifest = tmp_path / "sites.tsv"
    manifest.write_text(
        "\t".join(scan.SITE_COLUMNS) + "\n"
        "8\tAAAAAAAA\tchr1\t1\t8\t+\tTTTTTTTT\tr1\t0\t0\n"
    )
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nTTTTTTTT\n")
    pin = tmp_path / "genome.pin.json"
    pin.write_text(json.dumps({
        "source_genome_fasta_sha256": hashlib.sha256(genome.read_bytes()).hexdigest(),
    }))
    checkpoint = tmp_path / "hoff.ckpt.json"
    checkpoint.write_text(json.dumps({"site_index": 0, "input_digest": "stale"}))
    genome.write_text(">chr1\nAAAAAAAA\n")
    with pytest.raises(ValueError, match="genome FASTA SHA-256 mismatch"):
        hoff.run_sweep(
            manifest=manifest, out=tmp_path / "out.json", checkpoint=checkpoint,
            genome_fasta=genome, genome_pin=pin,
            scan_burden=_write_scan_burden(manifest),
        )


def test_hoff_validate_output_rejects_stale_done_after_fasta_and_pin_change(
    monkeypatch, tmp_path: Path,
) -> None:
    manifest = tmp_path / "sites.tsv"
    manifest.write_text(
        "\t".join(scan.SITE_COLUMNS) + "\n"
        "8\tAAAAAAAA\tchr1\t1\t8\t+\tTTTTTTTT\tr1\t0\t0\n"
    )
    burden = _write_scan_burden(manifest)
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nTTTTTTTT\n")
    pin = tmp_path / "genome.pin.json"
    pin.write_text(json.dumps({
        "source_genome_fasta_sha256": hashlib.sha256(genome.read_bytes()).hexdigest(),
    }))
    out = tmp_path / "out.json"
    checkpoint = tmp_path / "ckpt.json"
    monkeypatch.setattr(
        hoff, "pangolin_off_target_burden",
        lambda *_args, **_kwargs: {"pangolin_off_target_burden": 1.0, "reason": None},
    )
    hoff.run_sweep(
        manifest=manifest, scan_burden=burden, out=out, checkpoint=checkpoint,
        genome_fasta=genome, genome_pin=pin,
    )
    done = tmp_path / "out.done"
    done.touch()

    genome.write_text(">chr1\nAAAAAAAA\n")
    argv = [
        "arm_e_hoff", "--manifest", str(manifest), "--scan-burden", str(burden),
        "--out", str(out), "--checkpoint", str(checkpoint),
        "--genome-fasta", str(genome), "--genome-pin", str(pin),
        "--validate-output",
    ]
    monkeypatch.setattr(sys, "argv", argv)
    with pytest.raises(ValueError, match="genome FASTA SHA-256 mismatch"):
        hoff._main()

    pin.write_text(json.dumps({
        "source_genome_fasta_sha256": hashlib.sha256(genome.read_bytes()).hexdigest(),
    }))
    with pytest.raises(ValueError, match="completed H-off input/config digest mismatch"):
        hoff._main()
    assert done.exists()  # Rejected by nonzero validation; wrapper cannot accept/touch it.


def test_hoff_propagates_none_without_fabrication(monkeypatch, tmp_path: Path) -> None:
    manifest = tmp_path / "sites.tsv"
    manifest.write_text("\t".join(scan.SITE_COLUMNS) + "\n4\tAAAA\tchr1\t10\t13\t+\tTTTT\tr1\t0\t0\n")
    monkeypatch.setattr(hoff, "pangolin_off_target_burden",
                        lambda *a, **k: {"pangolin_off_target_burden": None, "reason": "worker failed"})
    out_path = tmp_path / "out.json"
    with pytest.raises(hoff.HoffIncompleteError, match="worker failed"):
        hoff.run_sweep(manifest=manifest, out=out_path,
                       checkpoint=tmp_path / "ckpt.json",
                       scan_burden=_write_scan_burden(manifest))
    out = json.loads(out_path.read_text())
    assert out["pangolin_off_target_burden"] is None
    assert out["pangolin_off_target_incomplete"] is True
    assert out["terminal"] is False


def test_hoff_failed_diagnostic_does_not_block_checkpointed_retry(monkeypatch, tmp_path: Path) -> None:
    manifest = tmp_path / "sites.tsv"
    manifest.write_text("\t".join(scan.SITE_COLUMNS) + "\n8\tAAAAAAAA\tchr1\t10\t17\t+\tTTTTTTTT\tr1\t0\t0\n")
    out, checkpoint = tmp_path / "out.json", tmp_path / "ckpt.json"
    monkeypatch.setattr(hoff, "pangolin_off_target_burden",
                        lambda *a, **k: {"pangolin_off_target_burden": None, "reason": "failed once"})
    scan_burden = _write_scan_burden(manifest)
    with pytest.raises(hoff.HoffIncompleteError):
        hoff.run_sweep(manifest=manifest, out=out, checkpoint=checkpoint,
                       scan_burden=scan_burden)
    monkeypatch.setattr(hoff, "pangolin_off_target_burden",
                        lambda *a, **k: {"pangolin_off_target_burden": 1.25, "reason": None})
    result = hoff.run_sweep(manifest=manifest, out=out, checkpoint=checkpoint,
                            scan_burden=scan_burden)
    assert result["terminal"] is True
    assert result["pangolin_off_target_burden"] == 1.25


def test_hoff_checkpoint_is_bound_to_manifest_and_config(monkeypatch, tmp_path: Path) -> None:
    manifest = tmp_path / "sites.tsv"
    manifest.write_text("\t".join(scan.SITE_COLUMNS) + "\n4\tAAAA\tchr1\t10\t13\t+\tTTTT\tr1\t0\t0\n")
    checkpoint = tmp_path / "ckpt.json"
    checkpoint.write_text(json.dumps({"site_index": 0, "input_digest": "stale"}))
    with pytest.raises(ValueError, match="checkpoint input/config digest"):
        hoff.run_sweep(manifest=manifest, out=tmp_path / "out.json", checkpoint=checkpoint,
                       scan_burden=_write_scan_burden(manifest))


@pytest.mark.parametrize(
    "checkpoint_state",
    [
        {},
        {"site_index": 1, "input_digest": "a" * 64},
        {"site_index": 1, "input_digest": "a" * 64,
         "pangolin_off_target_burden": 0},
        {"site_index": 1, "input_digest": "a" * 64,
         "pangolin_off_target_burden": float("nan")},
        {"site_index": 1, "input_digest": "a" * 64,
         "pangolin_off_target_burden": -0.1},
    ],
)
def test_hoff_refuses_checkpoint_without_valid_prior_burden(
    monkeypatch, tmp_path: Path, checkpoint_state: dict,
) -> None:
    manifest = tmp_path / "sites.tsv"
    manifest.write_text(
        "\t".join(scan.SITE_COLUMNS) + "\n"
        "4\tAAAA\tchr1\t10\t13\t+\tTTTT\tr1\t0\t0\n"
        "4\tAAAA\tchr1\t20\t23\t+\tTTTT\tr2\t0\t0\n"
    )
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "T" * 40 + "\n")
    checkpoint = tmp_path / "ckpt.json"
    checkpoint.write_text(json.dumps(checkpoint_state))
    monkeypatch.setattr(hoff, "_input_digest", lambda _payload: "a" * 64)
    monkeypatch.setattr(hoff, "_digest", lambda _payload: "a" * 64)
    monkeypatch.setattr(
        hoff, "pangolin_off_target_burden",
        lambda *_a, **_k: pytest.fail("Pangolin must not run for invalid checkpoint"),
    )
    with pytest.raises(ValueError, match="checkpoint"):
        hoff.run_sweep(
            manifest=manifest, out=tmp_path / "out.json", checkpoint=checkpoint,
            scan_burden=_write_scan_burden(manifest), genome_fasta=genome,
        )


@pytest.mark.parametrize("site_index", [pytest.param(None, id="missing"), "1", 1.0, True])
def test_hoff_refuses_checkpoint_without_exact_integer_site_index(
    monkeypatch, tmp_path: Path, site_index: object,
) -> None:
    manifest = tmp_path / "sites.tsv"
    manifest.write_text(
        "\t".join(scan.SITE_COLUMNS) + "\n"
        "4\tAAAA\tchr1\t10\t13\t+\tTTTT\tr1\t0\t0\n"
        "4\tAAAA\tchr1\t20\t23\t+\tTTTT\tr2\t0\t0\n"
    )
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "T" * 40 + "\n")
    checkpoint_state = {
        "input_digest": "a" * 64,
        "pangolin_off_target_burden": 0.0,
    }
    if site_index is not None:
        checkpoint_state["site_index"] = site_index
    checkpoint = tmp_path / "ckpt.json"
    checkpoint.write_text(json.dumps(checkpoint_state))
    monkeypatch.setattr(hoff, "_input_digest", lambda _payload: "a" * 64, raising=False)
    monkeypatch.setattr(hoff, "_digest", lambda _payload: "a" * 64)
    monkeypatch.setattr(
        hoff, "pangolin_off_target_burden",
        lambda *_a, **_k: pytest.fail("Pangolin must not run for invalid checkpoint"),
    )
    with pytest.raises(ValueError, match="checkpoint site_index missing or invalid"):
        hoff.run_sweep(
            manifest=manifest, out=tmp_path / "out.json", checkpoint=checkpoint,
            scan_burden=_write_scan_burden(manifest), genome_fasta=genome,
        )


def test_hoff_input_identity_changes_with_pangolin_weight_digest(
    monkeypatch, tmp_path: Path,
) -> None:
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nTTTT\n")
    scan_burden = tmp_path / "scan.json"
    scan_burden.write_text("{}")
    kwargs = {
        "manifest_sha": "1" * 64,
        "scan_burden": scan_burden,
        "scan_payload": {"length": 4, "aso_5to3": "AAAA"},
        "genome_fasta": genome,
        "genome_pin": None,
        "gpu_python": Path("/gpu/bin/python"),
        "chunk_size": 2,
    }
    first_hashes = {name: "a" * 64 for name in hoff.PANGOLIN_WEIGHT_FILENAMES}
    second_hashes = dict(first_hashes)
    second_hashes[hoff.PANGOLIN_WEIGHT_FILENAMES[0]] = "b" * 64
    first = {
        "pangolin_package_version": "1.0",
        "pangolin_package_sha256": "c" * 64,
        "pangolin_package_pinned_sha256": "c" * 64,
        "pangolin_package_manifest": {
            "file_count": 1,
            "files": [{"path": "model.py", "sha256": "d" * 64}],
        },
        "pangolin_weight_sha256": first_hashes,
    }
    second = {
        "pangolin_package_version": "1.0",
        "pangolin_package_sha256": "c" * 64,
        "pangolin_package_pinned_sha256": "c" * 64,
        "pangolin_package_manifest": {
            "file_count": 1,
            "files": [{"path": "model.py", "sha256": "d" * 64}],
        },
        "pangolin_weight_sha256": second_hashes,
    }
    monkeypatch.setattr(hoff, "_pangolin_provenance", lambda _path: first)
    inputs_a = hoff._hoff_inputs(**kwargs)
    monkeypatch.setattr(hoff, "_pangolin_provenance", lambda _path: second)
    inputs_b = hoff._hoff_inputs(**kwargs)
    assert inputs_a["pangolin_package_version"] == "1.0"
    assert set(inputs_a["pangolin_weight_sha256"]) == set(hoff.PANGOLIN_WEIGHT_FILENAMES)
    assert len(inputs_a["pangolin_weight_sha256"]) == 12
    assert hoff._digest(inputs_a) != hoff._digest(inputs_b)


def test_pangolin_provenance_includes_python_source_identity(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    gpu_python = tmp_path / "python"
    gpu_python.write_text("")
    package_identity = {
        "pangolin_package_sha256": "b" * 64,
        "pangolin_package_manifest": {
            "file_count": 1,
            "files": [{"path": "model.py", "sha256": "e" * 64}],
        },
    }
    weight_hashes = {
        name: "a" * 64 for name in hoff.PANGOLIN_WEIGHT_FILENAMES
    }
    monkeypatch.setattr(
        hoff.source_auth,
        "PINNED_PANGOLIN_PACKAGE_SHA256",
        "b" * 64,
        raising=False,
    )
    monkeypatch.setattr(
        hoff.source_auth,
        "probe_pangolin_package_identity",
        lambda _python: package_identity,
    )
    monkeypatch.setattr(
        hoff.subprocess,
        "run",
        lambda *_args, **_kwargs: SimpleNamespace(
            returncode=0,
            stdout=json.dumps({
                "pangolin_package_version": "1.0",
                "pangolin_weight_sha256": weight_hashes,
            }) + "\n",
            stderr="",
        ),
    )

    provenance = _REAL_PANGOLIN_PROVENANCE(gpu_python)

    assert provenance == {
        "pangolin_package_version": "1.0",
        "pangolin_weight_sha256": weight_hashes,
        "pangolin_package_pinned_sha256": "b" * 64,
        **package_identity,
    }


def test_pangolin_provenance_rejects_live_identity_different_from_committed_pin(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    gpu_python = tmp_path / "python"
    gpu_python.write_text("")
    monkeypatch.setattr(
        hoff.source_auth,
        "PINNED_PANGOLIN_PACKAGE_SHA256",
        "a" * 64,
        raising=False,
    )
    monkeypatch.setattr(
        hoff.source_auth,
        "probe_pangolin_package_identity",
        lambda _python: {
            "pangolin_package_sha256": "b" * 64,
            "pangolin_package_manifest": {
                "file_count": 1,
                "files": [{"path": "model.py", "sha256": "e" * 64}],
            },
        },
    )

    with pytest.raises(
        hoff.source_auth.SourceAuthenticationError,
        match="committed pin",
    ):
        _REAL_PANGOLIN_PROVENANCE(gpu_python)


def test_hoff_input_digest_changes_with_pangolin_package_code_hash(
    monkeypatch, tmp_path: Path,
) -> None:
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nTTTT\n")
    scan_burden = tmp_path / "scan.json"
    scan_burden.write_text("{}")
    kwargs = {
        "manifest_sha": "1" * 64,
        "scan_burden": scan_burden,
        "scan_payload": {"length": 4, "aso_5to3": "AAAA"},
        "genome_fasta": genome,
        "genome_pin": None,
        "gpu_python": Path("/gpu/bin/python"),
        "chunk_size": 2,
    }
    first = hoff._pangolin_provenance(Path("/gpu/bin/python"))
    second = {
        **first,
        "pangolin_package_sha256": "f" * 64,
    }
    pinned_drift = {
        **first,
        "pangolin_package_pinned_sha256": "0" * 64,
    }
    monkeypatch.setattr(hoff, "_pangolin_provenance", lambda _path: first)
    inputs_a = hoff._hoff_inputs(**kwargs)
    monkeypatch.setattr(hoff, "_pangolin_provenance", lambda _path: second)
    inputs_b = hoff._hoff_inputs(**kwargs)
    monkeypatch.setattr(hoff, "_pangolin_provenance", lambda _path: pinned_drift)
    inputs_c = hoff._hoff_inputs(**kwargs)

    assert inputs_a["pangolin_package_sha256"] == "b" * 64
    assert inputs_a["pangolin_package_pinned_sha256"] == "b" * 64
    assert hoff._input_digest(inputs_a) != hoff._input_digest(inputs_b)
    assert hoff._input_digest(inputs_a) != hoff._input_digest(inputs_c)


def test_worker_rejects_missing_package_identity_before_runtime_imports(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    real_import = __import__

    def guarded_import(name, *args, **kwargs):
        if name.split(".", 1)[0] in {"numpy", "torch", "pangolin"}:
            pytest.fail(f"runtime dependency imported before package authentication: {name}")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(
        gs.source_auth,
        "PINNED_PANGOLIN_PACKAGE_SHA256",
        "0" * 64,
        raising=False,
    )
    monkeypatch.setattr(importlib.util, "find_spec", lambda _name: None)
    monkeypatch.setattr("builtins.__import__", guarded_import)

    with pytest.raises(
        gs.source_auth.SourceAuthenticationError,
        match="Pangolin package location unavailable",
    ):
        gs._worker_main()


def test_worker_rejects_package_identity_mismatch_before_runtime_imports(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    package = tmp_path / "pangolin"
    package.mkdir()
    (package / "model.py").write_text("MODEL = 1\n")
    monkeypatch.setattr(
        importlib.util,
        "find_spec",
        lambda _name: SimpleNamespace(submodule_search_locations=[str(package)]),
    )
    monkeypatch.setattr(
        gs.source_auth,
        "PINNED_PANGOLIN_PACKAGE_SHA256",
        "0" * 64,
        raising=False,
    )
    real_import = __import__

    def guarded_import(name, *args, **kwargs):
        if name.split(".", 1)[0] in {"numpy", "torch", "pangolin"}:
            pytest.fail(f"runtime dependency imported before package authentication: {name}")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr("builtins.__import__", guarded_import)

    with pytest.raises(
        gs.source_auth.SourceAuthenticationError,
        match="committed pin",
    ):
        gs._worker_main()


def test_hoff_input_digest_binds_distance_context_and_whole_source_tree(
    monkeypatch,
) -> None:
    bindings = hoff._source_digest_bindings()
    assert bindings["pangolin_distance"] == 50
    assert bindings["pangolin_context"] == 5000
    assert bindings["source_tree_sha256"] == "c" * 64
    assert bindings["executed_git_commit"] == "a" * 40


def test_hoff_inputs_expose_whole_source_identity_and_bind_it_into_digest(
    monkeypatch, tmp_path: Path,
) -> None:
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nTTTT\n")
    scan_burden = tmp_path / "scan.json"
    scan_burden.write_text("{}")
    kwargs = {
        "manifest_sha": "1" * 64,
        "scan_burden": scan_burden,
        "scan_payload": {"length": 4, "aso_5to3": "AAAA"},
        "genome_fasta": genome,
        "genome_pin": None,
        "gpu_python": Path("/gpu/bin/python"),
        "chunk_size": 2,
    }
    inputs_a = hoff._hoff_inputs(**kwargs)
    assert inputs_a["source_tree_sha256"] == "c" * 64
    assert inputs_a["executed_git_commit"] == "a" * 40
    monkeypatch.setattr(
        hoff, "_source_digest_bindings", lambda: {
            "pangolin_distance": 50,
            "pangolin_context": 5000,
            "source_tree_sha256": "e" * 64,
            "executed_git_commit": "a" * 40,
        },
    )
    inputs_b = hoff._hoff_inputs(**kwargs)
    assert hoff._input_digest(inputs_a) != hoff._input_digest(inputs_b)


def test_execution_provenance_is_recorded_but_not_bound_into_input_digest(
    monkeypatch, tmp_path: Path,
) -> None:
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nTTTT\n")
    scan_burden = tmp_path / "scan.json"
    scan_burden.write_text("{}")
    kwargs = {
        "manifest_sha": "1" * 64,
        "scan_burden": scan_burden,
        "scan_payload": {"length": 4, "aso_5to3": "AAAA"},
        "genome_fasta": genome,
        "genome_pin": None,
        "gpu_python": Path("/gpu/bin/python"),
        "chunk_size": 2,
    }
    inputs = hoff._hoff_inputs(**kwargs)
    digest = hoff._input_digest(inputs)
    first = hoff._execution_provenance()
    monkeypatch.setenv("FENG_ARME_RUN_ID", "20260911T120000Z")
    monkeypatch.setenv("SLURM_ARRAY_JOB_ID", "54321")
    second = hoff._execution_provenance()
    assert first != second
    assert hoff._input_digest(inputs) == digest
    assert set(first) == {
        "run_id", "slurm_array_job_id", "slurm_array_task_id",
        "launcher_blob_sha256",
    }


def _empty_hoff_fixture(tmp_path: Path, *, scan_complete: bool) -> dict:
    length = 8 if scan_complete else 4
    manifest = tmp_path / f"sites_L{length:02d}.tsv"
    manifest.write_text("\t".join(scan.SITE_COLUMNS) + "\n")
    scan_burden = tmp_path / f"burden_L{length:02d}.json"
    scan_burden.write_text(json.dumps({
        "schema_version": scan.SCAN_ARTIFACT_SCHEMA_VERSION,
        "algorithm_version": scan.SCAN_ALGORITHM_VERSION,
        "terminal": True,
        "complete": scan_complete,
        "input_digest": "0" * 64,
        "reference_sha256": "1" * 64,
        "manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
        "length": length,
        "aso_5to3": "A" * length,
        "off_target_hits": 0 if scan_complete else None,
        "reason": None if scan_complete else "length_below_minimum_informative_length",
    }))
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nCCCCCCCC\n")
    out = tmp_path / f"hoff_L{length:02d}.json"
    checkpoint = tmp_path / f"hoff_L{length:02d}.ckpt.json"
    hoff.run_sweep(
        manifest=manifest, scan_burden=scan_burden, out=out,
        checkpoint=checkpoint, genome_fasta=genome,
    )
    return {
        "path": out,
        "manifest": manifest,
        "scan_burden": scan_burden,
        "genome_fasta": genome,
    }


def _validate_fixture(fixture: dict) -> dict:
    return hoff.validate_output(
        fixture["path"], fixture["manifest"],
        scan_burden=fixture["scan_burden"],
        genome_fasta=fixture["genome_fasta"],
    )


def test_hoff_payload_records_execution_provenance_outside_digest_inputs(
    tmp_path: Path,
) -> None:
    fixture = _empty_hoff_fixture(tmp_path, scan_complete=True)
    payload = json.loads(fixture["path"].read_text())
    assert payload["execution_provenance"] == {
        "run_id": "20260910T120000Z",
        "slurm_array_job_id": "12345",
        "slurm_array_task_id": "1",
        "launcher_blob_sha256": "d" * 64,
    }
    assert not set(payload["execution_provenance"]) & set(payload["inputs"])
    _validate_fixture(fixture)


@pytest.mark.parametrize("bad_flag", [None, 0, 1, "false", [], {}])
def test_hoff_validate_output_requires_exact_boolean_incomplete_flag(
    tmp_path: Path, bad_flag: object,
) -> None:
    fixture = _empty_hoff_fixture(tmp_path, scan_complete=True)
    payload = json.loads(fixture["path"].read_text())
    payload["pangolin_off_target_incomplete"] = bad_flag
    fixture["path"].write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="incomplete flag must be a boolean"):
        _validate_fixture(fixture)


@pytest.mark.parametrize("bad_count", [None, True, 0.0, "0"])
def test_hoff_validate_output_requires_exact_integer_site_count_unconditionally(
    tmp_path: Path, bad_count: object,
) -> None:
    fixture = _empty_hoff_fixture(tmp_path, scan_complete=False)
    payload = json.loads(fixture["path"].read_text())
    if bad_count is None:
        payload.pop("n_sites")
    else:
        payload["n_sites"] = bad_count
    fixture["path"].write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="site count must be an exact integer"):
        _validate_fixture(fixture)


def test_hoff_validate_output_rejects_boolean_site_count_equal_to_one(
    monkeypatch, tmp_path: Path,
) -> None:
    manifest = tmp_path / "sitess.tsv"
    manifest.write_text(
        "\t".join(scan.SITE_COLUMNS) + "\n"
        "8\tAAAAAAAA\tchr1\t1\t8\t+\tTTTTTTTT\tr1\t0\t0\n"
    )
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nTTTTTTTT\n")
    scan_burden = _write_scan_burden(manifest)
    out = tmp_path / "hoff.json"
    monkeypatch.setattr(
        hoff, "pangolin_off_target_burden",
        lambda *_args, **_kwargs: {
            "pangolin_off_target_burden": 1.0,
            "reason": None,
        },
    )
    hoff.run_sweep(
        manifest=manifest, scan_burden=scan_burden, out=out,
        checkpoint=tmp_path / "hoff.ckpt.json", genome_fasta=genome,
    )
    payload = json.loads(out.read_text())
    payload["n_sites"] = True
    out.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="site count must be an exact integer"):
        hoff.validate_output(
            out, manifest, scan_burden=scan_burden, genome_fasta=genome,
        )


@pytest.mark.parametrize("scan_complete", [True, False])
def test_hoff_validate_output_cross_checks_empty_state_against_scan_completeness(
    tmp_path: Path, scan_complete: bool,
) -> None:
    fixture = _empty_hoff_fixture(tmp_path, scan_complete=scan_complete)
    payload = json.loads(fixture["path"].read_text())
    if scan_complete:
        payload["pangolin_off_target_incomplete"] = True
        payload["pangolin_off_target_burden"] = None
        payload["reason"] = "misrepresented as uninformative"
    else:
        payload["pangolin_off_target_incomplete"] = False
        payload["pangolin_off_target_burden"] = 0.0
        payload["reason"] = "misrepresented as complete zero"
    fixture["path"].write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="empty H-off state disagrees with authenticated scan"):
        _validate_fixture(fixture)


def test_hoff_validate_output_requires_candidate_state(tmp_path: Path) -> None:
    manifest = tmp_path / "sites.tsv"
    manifest.write_text("\t".join(scan.SITE_COLUMNS) + "\n")
    out = tmp_path / "hoff.json"
    out.write_text(json.dumps({
        "schema_version": hoff.HOFF_ARTIFACT_SCHEMA_VERSION,
        "algorithm_version": hoff.HOFF_ALGORITHM_VERSION,
        "terminal": True,
        "input_digest": "0" * 64,
        "manifest_sha256": hashlib.sha256(manifest.read_bytes()).hexdigest(),
        "state": "sealed",
    }))
    with pytest.raises(ValueError, match="candidate state"):
        hoff.validate_output(out, manifest)


def test_hoff_signal_checkpoint_uses_last_committed_batch(monkeypatch, tmp_path: Path) -> None:
    manifest = tmp_path / "sites.tsv"
    manifest.write_text(
        "\t".join(scan.SITE_COLUMNS) + "\n"
        "8\tAAAAAAAA\tchr1\t10\t17\t+\tTTTTTTTT\tr1\t0\t0\n"
        "8\tAAAAAAAA\tchr1\t20\t27\t+\tTTTTTTTT\tr2\t0\t0\n"
    )
    callback, snapshots = {}, []

    @contextmanager
    def capture(on_signal, **_kwargs):
        callback["fn"] = on_signal
        yield

    real_save = hoff.save_checkpoint

    def record_save(state, path):
        if "site_index" in state:
            snapshots.append(dict(state))
        return real_save(state, path)

    def fake(*_args, **_kwargs):
        callback["fn"](15)
        return {"pangolin_off_target_burden": 2.0, "reason": None}

    monkeypatch.setattr(hoff, "checkpoint_on_signal", capture)
    monkeypatch.setattr(hoff, "save_checkpoint", record_save)
    monkeypatch.setattr(hoff, "pangolin_off_target_burden", fake)
    hoff.run_sweep(
        manifest=manifest, out=tmp_path / "out.json",
        checkpoint=tmp_path / "ckpt.json", chunk_size=2,
        scan_burden=_write_scan_burden(manifest),
    )
    assert snapshots[0]["site_index"] == 0
    assert snapshots[0]["pangolin_off_target_burden"] == 0.0
    assert snapshots[1]["site_index"] == 2
    assert snapshots[1]["pangolin_off_target_burden"] == 2.0


def test_wrappers_only_touch_done_after_valid_complete_payload() -> None:
    root = Path(__file__).resolve().parents[1]
    cpu = (root / "scripts" / "arm_e_scan.sbatch").read_text()
    gpu = (root / "scripts" / "arm_e_hoff.sbatch").read_text()
    assert "validate-burden" in cpu
    assert "validate-output" in gpu
    assert cpu.index("validate-burden") < cpu.index('touch "${DONE}"')
    assert gpu.index("validate-output") < gpu.index('touch "${DONE}"')
