"""Phase-9 permutation-test seal + pure ranking-primitive tests."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from feng.bench.permute import (
    PHASE9_UNLOCK_ENV,
    assert_phase9_unlocked,
    percentile_within_stratum,
    phase9_is_unlocked,
    run_phase9_test,
    verify_freeze_manifest,
    verify_scored_strata,
)
from feng.bench.strata import FreezeGuardError


# --- pure ranking primitive (safe to test) -----------------------------------


def test_percentile_lower_is_better() -> None:
    bg = np.array([5.0, 6.0, 7.0, 8.0])
    # positive better (lower) than all background -> high percentile
    assert percentile_within_stratum(1.0, bg, lower_is_better=True) == pytest.approx(0.9)
    # positive worse (higher) than all -> low percentile
    assert percentile_within_stratum(9.0, bg, lower_is_better=True) == pytest.approx(0.1)


def test_percentile_higher_is_better_symmetry() -> None:
    bg = np.array([5.0, 6.0, 7.0, 8.0])
    assert percentile_within_stratum(9.0, bg, lower_is_better=False) == pytest.approx(0.9)


def test_percentile_ties_use_midrank() -> None:
    bg = np.array([1.0, 1.0, 3.0])
    # positive ties two backgrounds; mid-rank keeps it centred
    val = percentile_within_stratum(1.0, bg, lower_is_better=True)
    assert 0.0 < val < 1.0


def test_percentile_empty_stratum_raises() -> None:
    with pytest.raises(ValueError):
        percentile_within_stratum(1.0, np.array([]), lower_is_better=True)


# --- freeze seal on the outcome computation ----------------------------------


def test_assert_phase9_unlocked_raises_prefreeze(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(PHASE9_UNLOCK_ENV, raising=False)
    assert not phase9_is_unlocked()
    with pytest.raises(FreezeGuardError):
        assert_phase9_unlocked()


def test_phase9_unlock_via_env(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(PHASE9_UNLOCK_ENV, "1")
    monkeypatch.setattr("feng.bench.permute.verify_freeze_manifest", lambda path: True)
    manifest = Path("freeze-manifest.json")
    assert phase9_is_unlocked(freeze_manifest=manifest)
    assert_phase9_unlocked(freeze_manifest=manifest)


def test_phase9_unlock_via_sentinel(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(PHASE9_UNLOCK_ENV, raising=False)
    sentinel = tmp_path / "unlock.flag"
    monkeypatch.setattr("feng.bench.permute.verify_freeze_manifest", lambda path: True)
    manifest = tmp_path / "freeze-manifest.json"
    assert not phase9_is_unlocked(sentinel, freeze_manifest=manifest)  # absent unlock
    sentinel.write_text("ok")
    assert phase9_is_unlocked(sentinel, freeze_manifest=manifest)


def test_phase9_unlock_refused_without_valid_freeze_manifest(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv(PHASE9_UNLOCK_ENV, "1")
    monkeypatch.setattr("feng.bench.permute.verify_freeze_manifest", lambda path: False)
    assert not phase9_is_unlocked()
    with pytest.raises(FreezeGuardError, match="freeze"):
        assert_phase9_unlocked()


def _frozen_manifest_fixture(tmp_path, monkeypatch):
    from feng.bench import permute

    prereg = tmp_path / permute.EXPECTED_PREREGISTRATION_PATH
    prereg.parent.mkdir(parents=True, exist_ok=True)
    prereg.write_bytes(b"preregistered")
    payloads = {
        relative: (b"scores" if relative == permute.EXPECTED_SCORED_STRATA_PATH else relative.encode())
        for relative in permute.EXPECTED_ARTIFACT_PATHS
    }
    pointer_artifacts = {
        "feng/artifacts/bench/detectable_universe_v50.parquet",
        "feng/artifacts/bench/detectable_universe_v50_README.md",
    }
    tagged = {}
    for relative, content in payloads.items():
        artifact = tmp_path / relative
        artifact.parent.mkdir(parents=True, exist_ok=True)
        if relative in pointer_artifacts:
            bundle_relative = Path(
                "feng/artifacts/bench/detectable_universe_v50.bundles/frozen"
            )
            bundle_artifact = tmp_path / bundle_relative / artifact.name
            bundle_artifact.parent.mkdir(parents=True, exist_ok=True)
            bundle_artifact.write_bytes(content)
            artifact.unlink(missing_ok=True)
            artifact.symlink_to(
                Path("detectable_universe_v50.current") / artifact.name
            )
            tagged[relative] = os.readlink(artifact).encode()
            tagged[(bundle_relative / artifact.name).as_posix()] = content
        else:
            artifact.write_bytes(content)
            tagged[relative] = content
    current_relative = Path("feng/artifacts/bench/detectable_universe_v50.current")
    current = tmp_path / current_relative
    current.unlink(missing_ok=True)
    current.symlink_to(
        Path("detectable_universe_v50.bundles/frozen"), target_is_directory=True
    )
    tagged[current_relative.as_posix()] = os.readlink(current).encode()
    commit = "a" * 40
    tag_object = b"object tagged commit\ntag feng-freeze\ntagger approved\n\nfreeze\n"
    record = {
        "schema": permute.FREEZE_MANIFEST_SCHEMA,
        "tag": permute.APPROVED_FREEZE_TAG,
        "preregistration": {
            "path": permute.EXPECTED_PREREGISTRATION_PATH,
            "sha256": hashlib.sha256(b"preregistered").hexdigest(),
            "locked": True,
        },
        "artifacts": [
            {"path": relative, "sha256": hashlib.sha256(content).hexdigest()}
            for relative, content in sorted(payloads.items())
        ],
    }
    manifest = tmp_path / "feng/artifacts/bench/freeze_manifest.json"
    manifest.write_text(json.dumps(record))
    tagged.update({
        permute.EXPECTED_PREREGISTRATION_PATH: b"preregistered",
        "feng/artifacts/bench/freeze_manifest.json": manifest.read_bytes(),
    })

    def fake_git(args, *, text=True):
        if args == ["rev-parse", "--show-toplevel"]:
            return SimpleNamespace(stdout=str(tmp_path))
        if args == ["rev-parse", "feng-freeze^{commit}"]:
            return SimpleNamespace(stdout=commit)
        if args == ["cat-file", "tag", "refs/tags/feng-freeze"]:
            return SimpleNamespace(stdout=tag_object)
        if args[0] == "show":
            rel = args[1].split(":", 1)[1]
            return SimpleNamespace(stdout=tagged[rel])
        raise AssertionError(args)

    monkeypatch.setattr(permute, "FENG_REPO_ROOT", tmp_path)
    monkeypatch.setattr(permute, "FREEZE_MANIFEST_PATH", manifest)
    monkeypatch.setattr(
        permute,
        "APPROVED_FREEZE_TAG_OBJECT_SHA256",
        hashlib.sha256(tag_object).hexdigest(),
    )
    monkeypatch.setattr(permute, "_git", fake_git)
    return permute, manifest, record


def test_freeze_manifest_binds_exact_repo_schema_paths_and_file_hashes(tmp_path, monkeypatch):
    permute, manifest, _ = _frozen_manifest_fixture(tmp_path, monkeypatch)
    assert verify_freeze_manifest(manifest)
    artifact = tmp_path / permute.EXPECTED_SCORED_STRATA_PATH
    artifact.write_bytes(b"changed after freeze")
    assert not verify_freeze_manifest(manifest)


def test_freeze_manifest_resolves_versioned_bundle_symlinks_and_detects_tampering(
    tmp_path, monkeypatch
):
    permute, manifest, _ = _frozen_manifest_fixture(tmp_path, monkeypatch)
    assert verify_freeze_manifest(manifest)

    current = tmp_path / "feng/artifacts/bench/detectable_universe_v50.current"
    replacement = current.parent / "detectable_universe_v50.bundles/replacement"
    replacement.mkdir()
    for name in (
        "detectable_universe_v50.parquet",
        "detectable_universe_v50_README.md",
    ):
        (replacement / name).write_bytes((current / name).read_bytes())
    current.unlink()
    current.symlink_to(
        Path("detectable_universe_v50.bundles/replacement"), target_is_directory=True
    )

    assert not verify_freeze_manifest(manifest)


def test_freeze_manifest_rejects_arbitrary_prereg_and_incomplete_artifact_set(tmp_path, monkeypatch):
    _, manifest, record = _frozen_manifest_fixture(tmp_path, monkeypatch)
    record["preregistration"]["path"] = "feng/anything.md"
    manifest.write_text(json.dumps(record))
    assert not verify_freeze_manifest(manifest)

    _, manifest, record = _frozen_manifest_fixture(tmp_path, monkeypatch)
    record["artifacts"].pop()
    manifest.write_text(json.dumps(record))
    assert not verify_freeze_manifest(manifest)


def test_freeze_manifest_rejects_manifest_selected_tag_or_unpinned_tag_object(
    tmp_path, monkeypatch
):
    permute, manifest, record = _frozen_manifest_fixture(tmp_path, monkeypatch)
    record["tag"] = "feng-freeze-attacker"
    manifest.write_text(json.dumps(record))
    assert not verify_freeze_manifest(manifest)

    permute, manifest, _ = _frozen_manifest_fixture(tmp_path, monkeypatch)
    monkeypatch.setattr(permute, "APPROVED_FREEZE_TAG_OBJECT_SHA256", "f" * 64)
    assert not verify_freeze_manifest(manifest)

def test_scored_strata_path_and_bytes_must_match_manifest(tmp_path, monkeypatch):
    permute, manifest, _ = _frozen_manifest_fixture(tmp_path, monkeypatch)
    scored = tmp_path / permute.EXPECTED_SCORED_STRATA_PATH
    assert verify_scored_strata(scored, manifest)
    alias = tmp_path / "alias.parquet"
    alias.write_bytes(scored.read_bytes())
    assert not verify_scored_strata(alias, manifest)


def test_phase9_guard_binds_scored_bytes_before_unimplemented_body(tmp_path, monkeypatch):
    permute, manifest, _ = _frozen_manifest_fixture(tmp_path, monkeypatch)
    monkeypatch.setenv(PHASE9_UNLOCK_ENV, "1")
    scored = tmp_path / permute.EXPECTED_SCORED_STRATA_PATH
    with pytest.raises(NotImplementedError):
        run_phase9_test(scored, freeze_manifest=manifest)


def test_run_phase9_test_locked_raises_freezeguard(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv(PHASE9_UNLOCK_ENV, raising=False)
    with pytest.raises(FreezeGuardError):
        run_phase9_test(Path("scored_strata.parquet"))
