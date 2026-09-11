"""Fail-closed authentication tests for the committed ``feng/src`` snapshot."""

from __future__ import annotations

import hashlib
import importlib
import importlib.machinery
import importlib.util
import json
import os
from pathlib import Path
import py_compile
from types import SimpleNamespace

import pytest


def _tree_sha(entries: list[dict[str, str]]) -> str:
    digest = hashlib.sha256()
    for entry in sorted(entries, key=lambda item: item["path"]):
        digest.update(entry["path"].encode("utf-8"))
        digest.update(b"\0")
        digest.update(entry["git_blob_sha256"].encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def _snapshot(tmp_path: Path) -> tuple[Path, Path, str, str]:
    source_root = tmp_path / "snapshot"
    files = {
        "feng/src/feng/__init__.py": b"\n",
        "feng/src/feng/example.py": b"VALUE = 7\n",
    }
    entries = []
    for relative, raw in files.items():
        path = source_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(raw)
        entries.append({
            "path": relative,
            "git_blob_sha256": hashlib.sha256(raw).hexdigest(),
        })
    commit = "a" * 40
    tree_sha = _tree_sha(entries)
    manifest = tmp_path / "source-manifest.json"
    manifest.write_text(json.dumps({
        "schema_version": 1,
        "executed_git_commit": commit,
        "source_tree_sha256": tree_sha,
        "entries": entries,
    }))
    return source_root, manifest, commit, tree_sha


def test_whole_source_tree_verification_passes_for_exact_snapshot(tmp_path: Path) -> None:
    source_auth = importlib.import_module("feng.source_auth")
    source_root, manifest, commit, tree_sha = _snapshot(tmp_path)

    identity = source_auth.verify_source_tree(
        source_root=source_root,
        manifest_path=manifest,
        expected_commit=commit,
        expected_tree_sha256=tree_sha,
    )

    assert identity == {
        "source_tree_sha256": tree_sha,
        "executed_git_commit": commit,
    }


@pytest.mark.parametrize("mutation", ["tampered", "missing", "extra"])
def test_whole_source_tree_verification_fails_closed_on_file_set_or_bytes(
    tmp_path: Path, mutation: str,
) -> None:
    source_auth = importlib.import_module("feng.source_auth")
    source_root, manifest, commit, tree_sha = _snapshot(tmp_path)
    target = source_root / "feng/src/feng/example.py"
    if mutation == "tampered":
        target.write_bytes(b"VALUE = 8\n")
    elif mutation == "missing":
        target.unlink()
    else:
        (source_root / "feng/src/feng/extra.py").write_bytes(b"unexpected\n")

    with pytest.raises(source_auth.SourceAuthenticationError):
        source_auth.verify_source_tree(
            source_root=source_root,
            manifest_path=manifest,
            expected_commit=commit,
            expected_tree_sha256=tree_sha,
        )


def test_source_auth_bootstrap_verifies_before_importing_feng(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_auth = importlib.import_module("feng.source_auth")
    source_root, manifest, commit, tree_sha = _snapshot(tmp_path)
    (source_root / "feng/src/feng/example.py").write_bytes(b"raise AssertionError('imported')\n")
    monkeypatch.setenv("FENG_ARME_SOURCE_ROOT", str(source_root))
    monkeypatch.setenv("FENG_ARME_SOURCE_MANIFEST", str(manifest))
    monkeypatch.setenv("FENG_ARME_SOURCE_TREE_SHA256", tree_sha)
    monkeypatch.setenv("FENG_ARME_SEALED_COMMIT", commit)

    with pytest.raises(source_auth.SourceAuthenticationError):
        source_auth.run_authenticated_module("feng.example", [])


def test_pangolin_package_identity_hashes_sorted_python_sources(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_auth = importlib.import_module("feng.source_auth")
    package = tmp_path / "pangolin"
    (package / "nested").mkdir(parents=True)
    files = {
        "__init__.py": b"# package\n",
        "model.py": b"MODEL = 1\n",
        "nested/helper.py": b"HELPER = 2\n",
    }
    for relative, raw in files.items():
        (package / relative).write_bytes(raw)
    (package / "models").mkdir()
    (package / "models/weight.bin").write_bytes(b"not Python source")
    monkeypatch.setattr(
        importlib.util,
        "find_spec",
        lambda name: SimpleNamespace(submodule_search_locations=[str(package)])
        if name == "pangolin" else None,
    )

    entries = [
        {"path": relative, "sha256": hashlib.sha256(raw).hexdigest()}
        for relative, raw in sorted(files.items())
    ]
    digest = hashlib.sha256()
    for entry in entries:
        digest.update(entry["path"].encode("utf-8"))
        digest.update(b"\0")
        digest.update(entry["sha256"].encode("ascii"))
        digest.update(b"\n")

    identity = source_auth.pangolin_package_identity()

    assert identity == {
        "pangolin_package_sha256": digest.hexdigest(),
        "pangolin_package_manifest": {
            "file_count": 3,
            "files": entries,
        },
    }


def test_pangolin_package_identity_fails_when_package_is_absent(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_auth = importlib.import_module("feng.source_auth")
    monkeypatch.setattr(importlib.util, "find_spec", lambda _name: None)

    with pytest.raises(
        source_auth.SourceAuthenticationError,
        match="Pangolin package location unavailable",
    ):
        source_auth.pangolin_package_identity()


def test_pangolin_package_verification_rejects_source_mismatch(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_auth = importlib.import_module("feng.source_auth")
    package = tmp_path / "pangolin"
    package.mkdir()
    (package / "model.py").write_text("MODEL = 1\n")
    monkeypatch.setattr(
        importlib.util,
        "find_spec",
        lambda _name: SimpleNamespace(submodule_search_locations=[str(package)]),
    )
    monkeypatch.setattr(source_auth, "PINNED_PANGOLIN_PACKAGE_SHA256", "0" * 64)

    with pytest.raises(
        source_auth.SourceAuthenticationError,
        match="Pangolin package source identity mismatch",
    ):
        source_auth.verify_pangolin_package()


def test_pangolin_package_verification_rejects_valid_pyc_shadow(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_auth = importlib.import_module("feng.source_auth")
    package = tmp_path / "pangolin"
    package.mkdir()
    (package / "__init__.py").write_text("\n")
    (package / "model.py").write_text("MODEL = 'reviewed'\n")
    spec = importlib.util.spec_from_file_location(
        "pangolin", package / "__init__.py",
        submodule_search_locations=[str(package)],
    )
    monkeypatch.setattr(importlib.util, "find_spec", lambda _name: spec)
    pyc = package / "__pycache__" / "model.cpython-311.pyc"
    pyc.parent.mkdir()
    py_compile.compile(
        str(package / "model.py"), cfile=str(pyc),
        dfile=str(package / "model.py"),
        invalidation_mode=py_compile.PycInvalidationMode.UNCHECKED_HASH,
    )
    reviewed_body = pyc.read_bytes()[16:]
    reviewed = source_auth.pangolin_package_identity()

    hostile_source = tmp_path / "hostile_model.py"
    hostile_source.write_text("MODEL = 'substituted'\n")
    py_compile.compile(
        str(hostile_source), cfile=str(pyc), dfile=str(package / "model.py"),
        invalidation_mode=py_compile.PycInvalidationMode.UNCHECKED_HASH,
    )

    changed = source_auth.pangolin_package_identity()
    assert pyc.read_bytes()[16:] != reviewed_body
    assert changed["pangolin_package_sha256"] != reviewed["pangolin_package_sha256"]
    assert any(
        entry["path"] == "__pycache__/model.cpython-311.pyc"
        for entry in changed["pangolin_package_manifest"]["files"]
    )
    with pytest.raises(
        source_auth.SourceAuthenticationError,
        match="Pangolin package source identity mismatch",
    ):
        monkeypatch.setattr(
            source_auth, "PINNED_PANGOLIN_PACKAGE_SHA256",
            reviewed["pangolin_package_sha256"],
        )
        source_auth.verify_pangolin_package()


def test_pangolin_package_identity_normalizes_benign_pyc_recompile_header(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_auth = importlib.import_module("feng.source_auth")
    package = tmp_path / "pangolin"
    package.mkdir()
    (package / "__init__.py").write_text("\n")
    model = package / "model.py"
    model.write_text("MODEL = 'reviewed'\n")
    spec = importlib.util.spec_from_file_location(
        "pangolin", package / "__init__.py",
        submodule_search_locations=[str(package)],
    )
    monkeypatch.setattr(importlib.util, "find_spec", lambda _name: spec)
    pyc = package / "__pycache__" / "model.cpython-311.pyc"
    pyc.parent.mkdir()
    os.utime(model, (1_700_000_000, 1_700_000_000))
    py_compile.compile(
        str(model), cfile=str(pyc), dfile=str(model),
        invalidation_mode=py_compile.PycInvalidationMode.TIMESTAMP,
    )
    first_bytes = pyc.read_bytes()
    first = source_auth.pangolin_package_identity()

    os.utime(model, (1_700_000_060, 1_700_000_060))
    py_compile.compile(
        str(model), cfile=str(pyc), dfile=str(model),
        invalidation_mode=py_compile.PycInvalidationMode.TIMESTAMP,
    )
    second_bytes = pyc.read_bytes()
    second = source_auth.pangolin_package_identity()

    assert first_bytes[:16] != second_bytes[:16]
    assert first_bytes[16:] == second_bytes[16:]
    assert first["pangolin_package_sha256"] == second["pangolin_package_sha256"]


def test_pangolin_package_identity_rejects_truncated_bytecode(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_auth = importlib.import_module("feng.source_auth")
    package = tmp_path / "pangolin"
    package.mkdir()
    (package / "__init__.py").write_text("\n")
    (package / "model.pyc").write_bytes(b"x" * 15)
    spec = importlib.util.spec_from_file_location(
        "pangolin", package / "__init__.py",
        submodule_search_locations=[str(package)],
    )
    monkeypatch.setattr(importlib.util, "find_spec", lambda _name: spec)

    with pytest.raises(
        source_auth.SourceAuthenticationError,
        match="truncated Pangolin bytecode",
    ):
        source_auth.pangolin_package_identity()


def test_pangolin_package_verification_rejects_extension_shadow(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    source_auth = importlib.import_module("feng.source_auth")
    package = tmp_path / "pangolin"
    package.mkdir()
    (package / "__init__.py").write_text("\n")
    (package / "model.py").write_text("MODEL = 'reviewed'\n")
    spec = importlib.util.spec_from_file_location(
        "pangolin", package / "__init__.py",
        submodule_search_locations=[str(package)],
    )
    monkeypatch.setattr(importlib.util, "find_spec", lambda _name: spec)
    reviewed = source_auth.pangolin_package_identity()

    extension = package / f"model{importlib.machinery.EXTENSION_SUFFIXES[0]}"
    extension.write_bytes(b"substituted extension module")
    changed = source_auth.pangolin_package_identity()

    assert changed["pangolin_package_sha256"] != reviewed["pangolin_package_sha256"]
    assert any(
        entry["path"] == extension.name
        for entry in changed["pangolin_package_manifest"]["files"]
    )
    with pytest.raises(
        source_auth.SourceAuthenticationError,
        match="non-source loader",
    ):
        monkeypatch.setattr(
            source_auth, "PINNED_PANGOLIN_PACKAGE_SHA256",
            changed["pangolin_package_sha256"],
        )
        source_auth.verify_pangolin_package()
