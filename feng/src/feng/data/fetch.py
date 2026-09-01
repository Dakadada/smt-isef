"""Hash-verifying fetch boundary for the clean-room Feng arm."""

from __future__ import annotations

from pathlib import Path
from typing import Callable

from .hashes import verify_hash
from .manifest import ManifestEntry

Downloader = Callable[[str, Path], None]


def verify_download(path: str | Path, entry: ManifestEntry) -> Path:
    """Verify authoritative MD5 and any declared secondary SHA-256."""
    artifact = Path(path)
    if not artifact.is_file():
        raise FileNotFoundError(artifact)
    verify_hash(artifact, entry.md5, "md5")
    if entry.sha256 is not None:
        verify_hash(artifact, entry.sha256, "sha256")
    return artifact


def fetch(entry: ManifestEntry, destination: str | Path, downloader: Downloader | None = None) -> Path:
    """Run an injected downloader and reject its output unless both hashes match."""
    if downloader is None:
        raise NotImplementedError("network transport is intentionally not implemented in Phase 1")
    artifact = Path(destination)
    downloader(entry.url, artifact)
    return verify_download(artifact, entry)
