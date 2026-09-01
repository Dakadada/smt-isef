"""Streaming MD5 and SHA-256 integrity helpers for the clean-room Feng arm."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Literal

HashAlgorithm = Literal["md5", "sha256"]


def file_hash(path: str | Path, algorithm: HashAlgorithm, chunk_size: int = 1024 * 1024) -> str:
    """Return a lowercase hexadecimal digest without loading the whole file."""
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    digest = hashlib.new(algorithm)
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b""):
            digest.update(chunk)
    return digest.hexdigest()


def md5_file(path: str | Path) -> str:
    """Return the file's MD5 digest for source-manifest compatibility."""
    return file_hash(path, "md5")


def sha256_file(path: str | Path) -> str:
    """Return the file's SHA-256 digest for strong integrity checking."""
    return file_hash(path, "sha256")


def verify_hash(path: str | Path, expected: str, algorithm: HashAlgorithm) -> None:
    """Raise ValueError when a file does not match its expected digest."""
    actual = file_hash(path, algorithm)
    if actual.casefold() != expected.casefold():
        raise ValueError(f"{algorithm} mismatch for {Path(path).name}: expected {expected}, got {actual}")

