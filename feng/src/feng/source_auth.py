"""Authenticate a committed ``feng/src`` snapshot before importing it.

This file is also a stdlib-only executable bootstrap.  The launch script runs it
by path, so no ``feng`` package is imported until :func:`verify_source_tree`
has re-hashed every materialized source entry.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.machinery
import importlib.util
import json
import os
from pathlib import Path
import re
import runpy
import stat
import subprocess
import sys
import tempfile
from typing import Iterable


SOURCE_MANIFEST_SCHEMA_VERSION = 1
_SHA256_RE = re.compile(r"[0-9a-f]{64}")
_COMMIT_RE = re.compile(r"(?:[0-9a-f]{40}|[0-9a-f]{64})")
# Seeded 2026-09-11 by calling pangolin_package_identity() under feng-gpu-env;
# this was a stdlib-only file scan and did not import Pangolin or Torch. Full
# seed manifest and reproduction command: feng/src/feng/pangolin_package_seed_provenance.json.
PINNED_PANGOLIN_PACKAGE_SHA256 = (
    "0af9857403226c5e72719e42d52fa5babfcda8fde9957ce7824df277e0051234"
)
_PANGOLIN_EXECUTABLE_SUFFIXES = tuple(sorted({
    ".py", ".pyc", ".pyo", ".pyd", *importlib.machinery.EXTENSION_SUFFIXES,
}, key=lambda suffix: (-len(suffix), suffix)))
_PANGOLIN_BYTECODE_SUFFIXES = (".pyc", ".pyo")


class SourceAuthenticationError(RuntimeError):
    """Runtime source does not match its authenticated expected identity."""


def _pangolin_candidate_sha256(path: Path) -> str:
    raw = path.read_bytes()
    if path.name.endswith(_PANGOLIN_BYTECODE_SUFFIXES):
        # SourceFileLoader is mandatory, so bytecode is only a cache. Normalize
        # CPython 3.7+'s volatile 16-byte header (magic, flags, mtime/hash+size)
        # while still authenticating the executable marshal/code body.
        if len(raw) < 16:
            raise SourceAuthenticationError(f"truncated Pangolin bytecode: {path}")
        raw = raw[16:]
    return hashlib.sha256(raw).hexdigest()


def _pangolin_executable_entries(package_dir: Path) -> list[dict[str, str]]:
    entries = []
    for directory, dirnames, filenames in os.walk(package_dir, followlinks=False):
        directory_path = Path(directory)
        symlink_dirs = sorted(
            name for name in dirnames if (directory_path / name).is_symlink()
        )
        if symlink_dirs:
            raise SourceAuthenticationError(
                f"unsupported Pangolin symlink directory: {directory_path / symlink_dirs[0]}"
            )
        dirnames[:] = sorted(dirnames)
        for name in sorted(filenames):
            if not name.endswith(_PANGOLIN_EXECUTABLE_SUFFIXES):
                continue
            path = directory_path / name
            if not path.is_file():
                raise SourceAuthenticationError(
                    f"unsupported Pangolin executable import candidate: {path}"
                )
            entries.append({
                "path": path.relative_to(package_dir).as_posix(),
                "sha256": _pangolin_candidate_sha256(path),
            })
    entries.sort(key=lambda item: item["path"])
    if not entries:
        raise SourceAuthenticationError("Pangolin package contains no executable import candidates")
    return entries


def _pangolin_entries_sha256(entries: Iterable[dict[str, str]]) -> str:
    digest = hashlib.sha256()
    for entry in sorted(entries, key=lambda item: item["path"]):
        digest.update(entry["path"].encode("utf-8"))
        digest.update(b"\0")
        digest.update(entry["sha256"].encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def pangolin_package_identity() -> dict[str, object]:
    """Hash every installed Pangolin executable import candidate without importing it."""
    spec = importlib.util.find_spec("pangolin")
    locations = list(spec.submodule_search_locations or ()) if spec is not None else []
    if len(locations) != 1:
        raise SourceAuthenticationError("Pangolin package location unavailable or ambiguous")
    package_dir = Path(locations[0]).resolve()
    if not package_dir.is_dir():
        raise SourceAuthenticationError("Pangolin package location unavailable or invalid")
    entries = _pangolin_executable_entries(package_dir)
    return {
        "pangolin_package_sha256": _pangolin_entries_sha256(entries),
        "pangolin_package_manifest": {
            "file_count": len(entries),
            "files": entries,
        },
    }


def _validated_pangolin_identity(identity: object) -> dict[str, object]:
    if type(identity) is not dict:
        raise SourceAuthenticationError("Pangolin package identity is malformed")
    package_sha = identity.get("pangolin_package_sha256")
    manifest = identity.get("pangolin_package_manifest")
    if type(package_sha) is not str or not _SHA256_RE.fullmatch(package_sha):
        raise SourceAuthenticationError("Pangolin package SHA-256 is malformed")
    if type(manifest) is not dict or set(manifest) != {"file_count", "files"}:
        raise SourceAuthenticationError("Pangolin package manifest is malformed")
    files = manifest["files"]
    if (type(manifest["file_count"]) is not int or type(files) is not list
            or not files or manifest["file_count"] != len(files)):
        raise SourceAuthenticationError("Pangolin package manifest file count is invalid")
    previous = None
    for entry in files:
        if type(entry) is not dict or set(entry) != {"path", "sha256"}:
            raise SourceAuthenticationError("Pangolin package manifest entry is malformed")
        path = entry["path"]
        file_sha = entry["sha256"]
        if (type(path) is not str or not path or path.startswith("/")
                or previous is not None and path <= previous):
            raise SourceAuthenticationError("Pangolin package manifest paths are invalid")
        if type(file_sha) is not str or not _SHA256_RE.fullmatch(file_sha):
            raise SourceAuthenticationError("Pangolin package file SHA-256 is malformed")
        previous = path
    if _pangolin_entries_sha256(files) != package_sha:
        raise SourceAuthenticationError("Pangolin package manifest digest mismatch")
    return identity


def _candidate_module(relative: str, package_dir: Path) -> tuple[str, Path] | None:
    parts = Path(relative).parts
    if "__pycache__" in parts:
        return None
    filename = parts[-1]
    suffix = next(
        (candidate for candidate in _PANGOLIN_EXECUTABLE_SUFFIXES
         if filename.endswith(candidate)),
        None,
    )
    if suffix is None:
        return None
    stem = filename[:-len(suffix)]
    module_parts = list(parts[:-1]) if stem == "__init__" else [*parts[:-1], stem]
    if any(not part.isidentifier() for part in module_parts):
        return None
    module = ".".join(["pangolin", *module_parts])
    search_dir = package_dir.joinpath(*module_parts[:-1]) if module_parts else package_dir
    return module, search_dir


def _assert_source_loader(module: str, spec: object, package_dir: Path) -> None:
    loader = getattr(spec, "loader", None) if spec is not None else None
    if not isinstance(loader, importlib.machinery.SourceFileLoader):
        loader_name = type(loader).__name__ if loader is not None else "None"
        raise SourceAuthenticationError(
            f"Pangolin module {module} resolved via non-source loader {loader_name}"
        )
    origin = getattr(spec, "origin", None)
    if type(origin) is not str:
        raise SourceAuthenticationError(f"Pangolin module {module} source origin is missing")
    try:
        Path(origin).resolve().relative_to(package_dir)
    except ValueError as exc:
        raise SourceAuthenticationError(
            f"Pangolin module {module} resolves outside the authenticated package"
        ) from exc


def _assert_pangolin_source_loaders(identity: dict[str, object]) -> None:
    spec = importlib.util.find_spec("pangolin")
    locations = list(spec.submodule_search_locations or ()) if spec is not None else []
    if len(locations) != 1:
        raise SourceAuthenticationError("Pangolin package location unavailable or ambiguous")
    package_dir = Path(locations[0]).resolve()
    _assert_source_loader("pangolin", spec, package_dir)

    modules: dict[str, Path] = {}
    manifest = identity["pangolin_package_manifest"]
    for entry in manifest["files"]:
        candidate = _candidate_module(entry["path"], package_dir)
        if candidate is not None:
            modules[candidate[0]] = candidate[1]
    for module, search_dir in sorted(modules.items()):
        if module == "pangolin":
            continue
        resolved = importlib.machinery.PathFinder.find_spec(module, [str(search_dir)])
        _assert_source_loader(module, resolved, package_dir)

    for module, loaded in sorted(sys.modules.items()):
        if module != "pangolin" and not module.startswith("pangolin."):
            continue
        loaded_spec = getattr(loaded, "__spec__", None) if loaded is not None else None
        _assert_source_loader(module, loaded_spec, package_dir)


def verify_pangolin_package() -> dict[str, object]:
    """Verify installed Pangolin against the sole committed package trust root."""
    if not _SHA256_RE.fullmatch(PINNED_PANGOLIN_PACKAGE_SHA256):
        raise SourceAuthenticationError("committed Pangolin package pin is malformed")
    identity = pangolin_package_identity()
    if identity["pangolin_package_sha256"] != PINNED_PANGOLIN_PACKAGE_SHA256:
        raise SourceAuthenticationError(
            "Pangolin package source identity mismatch with committed pin"
        )
    _assert_pangolin_source_loaders(identity)
    return identity


def probe_pangolin_package_identity(python: Path) -> dict[str, object]:
    """Ask another Python environment for its installed Pangolin source identity."""
    try:
        completed = subprocess.run(
            [str(python), str(Path(__file__).resolve()), "pangolin-package-identity"],
            check=False, capture_output=True, text=True, timeout=180,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise SourceAuthenticationError(f"Pangolin package identity probe failed: {exc}") from exc
    if completed.returncode != 0:
        detail = completed.stderr.strip().splitlines()[-1] if completed.stderr.strip() else "no stderr"
        raise SourceAuthenticationError(f"Pangolin package identity probe failed: {detail}")
    try:
        identity = json.loads(completed.stdout.strip().splitlines()[-1])
    except (IndexError, json.JSONDecodeError) as exc:
        raise SourceAuthenticationError(
            "Pangolin package identity probe returned malformed output"
        ) from exc
    identity = _validated_pangolin_identity(identity)
    if identity["pangolin_package_sha256"] != PINNED_PANGOLIN_PACKAGE_SHA256:
        raise SourceAuthenticationError(
            "Pangolin package identity probe disagrees with committed pin"
        )
    return identity


def source_tree_sha256(entries: Iterable[dict[str, str]]) -> str:
    """Hash sorted ``(repo-relative path, raw blob SHA-256)`` pairs."""
    digest = hashlib.sha256()
    for entry in sorted(entries, key=lambda item: item["path"]):
        digest.update(entry["path"].encode("utf-8"))
        digest.update(b"\0")
        digest.update(entry["git_blob_sha256"].encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def _raw_entry_bytes(path: Path) -> bytes:
    mode = path.lstat().st_mode
    if stat.S_ISREG(mode):
        return path.read_bytes()
    if stat.S_ISLNK(mode):
        return os.readlink(path).encode("utf-8")
    raise SourceAuthenticationError(f"unsupported source entry type: {path}")


def _materialized_entries(source_root: Path) -> dict[str, str]:
    source_dir = source_root / "feng" / "src"
    if not source_dir.is_dir():
        raise SourceAuthenticationError(f"materialized feng/src missing: {source_dir}")
    entries: dict[str, str] = {}
    for directory, dirnames, filenames in os.walk(source_dir, followlinks=False):
        directory_path = Path(directory)
        symlink_dirs = [name for name in dirnames if (directory_path / name).is_symlink()]
        dirnames[:] = [name for name in dirnames if name not in symlink_dirs]
        for name in sorted([*symlink_dirs, *filenames]):
            path = directory_path / name
            relative = path.relative_to(source_root).as_posix()
            entries[relative] = hashlib.sha256(_raw_entry_bytes(path)).hexdigest()
    return entries


def _validated_manifest(
    manifest_path: Path, *, expected_commit: str, expected_tree_sha256: str,
) -> list[dict[str, str]]:
    if not _COMMIT_RE.fullmatch(expected_commit):
        raise SourceAuthenticationError("expected commit must be a full lowercase Git SHA")
    if not _SHA256_RE.fullmatch(expected_tree_sha256):
        raise SourceAuthenticationError("expected source-tree SHA-256 is malformed")
    try:
        payload = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise SourceAuthenticationError(f"source manifest unavailable or invalid: {exc}") from exc
    if type(payload) is not dict or payload.get("schema_version") != SOURCE_MANIFEST_SCHEMA_VERSION:
        raise SourceAuthenticationError("source manifest schema mismatch")
    if payload.get("executed_git_commit") != expected_commit:
        raise SourceAuthenticationError("source manifest commit mismatch")
    if payload.get("source_tree_sha256") != expected_tree_sha256:
        raise SourceAuthenticationError("source manifest tree identity mismatch")
    entries = payload.get("entries")
    if type(entries) is not list or not entries:
        raise SourceAuthenticationError("source manifest entries missing")
    previous = None
    seen: set[str] = set()
    for entry in entries:
        if type(entry) is not dict or set(entry) != {"path", "git_blob_sha256"}:
            raise SourceAuthenticationError("source manifest entry malformed")
        path = entry["path"]
        blob_sha = entry["git_blob_sha256"]
        if type(path) is not str or not path.startswith("feng/src/") or path in seen:
            raise SourceAuthenticationError("source manifest path malformed or duplicated")
        if previous is not None and path <= previous:
            raise SourceAuthenticationError("source manifest entries are not strictly sorted")
        if type(blob_sha) is not str or not _SHA256_RE.fullmatch(blob_sha):
            raise SourceAuthenticationError("source manifest blob SHA-256 malformed")
        seen.add(path)
        previous = path
    if source_tree_sha256(entries) != expected_tree_sha256:
        raise SourceAuthenticationError("source manifest pair digest mismatch")
    return entries


def verify_source_tree(
    *, source_root: Path, manifest_path: Path, expected_commit: str,
    expected_tree_sha256: str,
) -> dict[str, str]:
    """Re-hash the complete staged tree and compare it with the Git manifest."""
    source_root = Path(source_root).resolve()
    entries = _validated_manifest(
        Path(manifest_path), expected_commit=expected_commit,
        expected_tree_sha256=expected_tree_sha256,
    )
    expected = {entry["path"]: entry["git_blob_sha256"] for entry in entries}
    actual = _materialized_entries(source_root)
    if set(actual) != set(expected):
        missing = sorted(set(expected) - set(actual))
        extra = sorted(set(actual) - set(expected))
        raise SourceAuthenticationError(
            f"source file set mismatch (missing={missing[:3]}, extra={extra[:3]})"
        )
    mismatched = sorted(path for path in expected if actual[path] != expected[path])
    if mismatched:
        raise SourceAuthenticationError(f"source bytes mismatch: {mismatched[:3]}")
    return {
        "source_tree_sha256": expected_tree_sha256,
        "executed_git_commit": expected_commit,
    }


def verify_source_tree_from_env() -> dict[str, str]:
    required = {
        "source_root": os.environ.get("FENG_ARME_SOURCE_ROOT", ""),
        "manifest_path": os.environ.get("FENG_ARME_SOURCE_MANIFEST", ""),
        "expected_commit": os.environ.get("FENG_ARME_SEALED_COMMIT", ""),
        "expected_tree_sha256": os.environ.get("FENG_ARME_SOURCE_TREE_SHA256", ""),
    }
    missing = [name for name, value in required.items() if not value]
    if missing:
        raise SourceAuthenticationError(
            f"source authentication environment missing: {', '.join(missing)}"
        )
    return verify_source_tree(
        source_root=Path(required["source_root"]),
        manifest_path=Path(required["manifest_path"]),
        expected_commit=required["expected_commit"],
        expected_tree_sha256=required["expected_tree_sha256"],
    )


def create_manifest_from_git(
    *, git_repo: Path, commit: str, source_root: Path, manifest_path: Path,
) -> dict[str, str]:
    """Build the raw-byte manifest from the sealed commit object store."""
    if not _COMMIT_RE.fullmatch(commit):
        raise SourceAuthenticationError("commit must be a full lowercase Git SHA")
    resolved = subprocess.run(
        ["git", "-C", str(git_repo), "rev-parse", "--verify", f"{commit}^{{commit}}"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    if resolved != commit:
        raise SourceAuthenticationError("sealed commit does not resolve exactly")
    tree = subprocess.run(
        ["git", "-C", str(git_repo), "ls-tree", "-rz", "--full-tree", commit, "--", "feng/src"],
        check=True, capture_output=True,
    ).stdout
    entries = []
    for record in tree.split(b"\0"):
        if not record:
            continue
        try:
            metadata, path_raw = record.split(b"\t", 1)
            _mode, object_type, object_id = metadata.split(b" ", 2)
            path = path_raw.decode("utf-8")
        except (UnicodeDecodeError, ValueError) as exc:
            raise SourceAuthenticationError("unparseable Git tree entry") from exc
        if object_type != b"blob" or not path.startswith("feng/src/"):
            raise SourceAuthenticationError(f"unsupported Git tree entry: {path}")
        raw = subprocess.run(
            ["git", "-C", str(git_repo), "cat-file", "blob", object_id.decode("ascii")],
            check=True, capture_output=True,
        ).stdout
        entries.append({
            "path": path,
            "git_blob_sha256": hashlib.sha256(raw).hexdigest(),
        })
    entries.sort(key=lambda item: item["path"])
    if not entries:
        raise SourceAuthenticationError("sealed commit contains no feng/src files")
    tree_sha = source_tree_sha256(entries)
    payload = {
        "schema_version": SOURCE_MANIFEST_SCHEMA_VERSION,
        "executed_git_commit": commit,
        "source_tree_sha256": tree_sha,
        "entries": entries,
    }
    manifest_path = Path(manifest_path)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{manifest_path.name}.", dir=manifest_path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(payload, handle, sort_keys=True, separators=(",", ":"))
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, manifest_path)
    except BaseException:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise
    return verify_source_tree(
        source_root=source_root, manifest_path=manifest_path,
        expected_commit=commit, expected_tree_sha256=tree_sha,
    )


def run_authenticated_module(module: str, arguments: list[str]) -> None:
    """Verify in this process, then import and execute ``module``."""
    verify_source_tree_from_env()
    source_root = Path(os.environ["FENG_ARME_SOURCE_ROOT"]).resolve()
    source_dir = str(source_root / "feng" / "src")
    sys.dont_write_bytecode = True
    sys.path[:] = [source_dir, *[
        item for item in sys.path
        if item and str(Path(item).resolve()) != source_dir
    ]]
    sys.argv = [module, *arguments]
    runpy.run_module(module, run_name="__main__", alter_sys=False)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    create = subparsers.add_parser("create-manifest")
    create.add_argument("--git-repo", type=Path, required=True)
    create.add_argument("--commit", required=True)
    create.add_argument("--source-root", type=Path, required=True)
    create.add_argument("--manifest", type=Path, required=True)
    run = subparsers.add_parser("run-module")
    run.add_argument("module")
    run.add_argument("arguments", nargs=argparse.REMAINDER)
    subparsers.add_parser("pangolin-package-identity")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "pangolin-package-identity":
            print(json.dumps(verify_pangolin_package(), sort_keys=True, separators=(",", ":")))
            return 0
        if args.command == "create-manifest":
            identity = create_manifest_from_git(
                git_repo=args.git_repo, commit=args.commit,
                source_root=args.source_root, manifest_path=args.manifest,
            )
            print(identity["source_tree_sha256"])
            return 0
        arguments = args.arguments[1:] if args.arguments[:1] == ["--"] else args.arguments
        run_authenticated_module(args.module, arguments)
        return 0
    except SourceAuthenticationError as exc:
        print(f"ERROR: source authentication failed: {exc}", file=sys.stderr)
        return 70


if __name__ == "__main__":
    raise SystemExit(main())
