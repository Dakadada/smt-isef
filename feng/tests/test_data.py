"""Verify manifest integrity primitives for the clean-room Feng arm."""

import hashlib

import pytest

from feng.data.fetch import fetch, verify_download
from feng.data.hashes import md5_file, sha256_file, verify_hash
from feng.data.manifest import (
    GENCODE_V50_PRIMARY_ASSEMBLY_FASTA,
    GENCODE_V50_PRIMARY_ASSEMBLY_GTF,
    Manifest,
    ManifestEntry,
)


PAYLOAD = b"feng manifest fixture\n"


def _entry(name: str = "fixture") -> ManifestEntry:
    return ManifestEntry(
        name=name,
        url="https://example.org/fixture.dat",
        md5=hashlib.md5(PAYLOAD).hexdigest(),
        sha256=hashlib.sha256(PAYLOAD).hexdigest(),
    )


def test_hash_helpers_and_download_verification(tmp_path) -> None:
    artifact = tmp_path / "fixture.dat"
    artifact.write_bytes(PAYLOAD)
    assert md5_file(artifact) == hashlib.md5(PAYLOAD).hexdigest()
    assert sha256_file(artifact) == hashlib.sha256(PAYLOAD).hexdigest()
    assert verify_download(artifact, _entry()) == artifact


def test_hash_mismatch_is_rejected(tmp_path) -> None:
    artifact = tmp_path / "fixture.dat"
    artifact.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="mismatch"):
        verify_hash(artifact, _entry().sha256, "sha256")


def test_manifest_schema_validates_urls_hashes_and_names() -> None:
    entry = _entry()
    assert Manifest((entry,)).entries == (entry,)
    with pytest.raises(ValueError, match="URL"):
        ManifestEntry("bad", "relative.dat", "0" * 32, "0" * 64)
    with pytest.raises(ValueError, match="md5"):
        ManifestEntry("bad", "https://example.org/x", "xyz", "0" * 64)
    with pytest.raises(ValueError, match="all-zero"):
        ManifestEntry("bad", "https://example.org/x", "0" * 32, None)
    with pytest.raises(ValueError, match="all-zero"):
        ManifestEntry("bad", "https://example.org/x", "1" * 32, "0" * 64)
    with pytest.raises(ValueError, match="unique"):
        Manifest((_entry(), _entry()))


def test_gencode_entries_pin_official_md5_without_unverified_sha256() -> None:
    assert GENCODE_V50_PRIMARY_ASSEMBLY_FASTA.md5 == "da1a11258be075cfa7af718162c894e7"
    assert GENCODE_V50_PRIMARY_ASSEMBLY_GTF.md5 == "289b91e5e95e8b0450d223246f10a12e"
    assert GENCODE_V50_PRIMARY_ASSEMBLY_FASTA.sha256 is None
    assert GENCODE_V50_PRIMARY_ASSEMBLY_GTF.sha256 is None


def test_md5_only_manifest_entry_is_verifiable(tmp_path) -> None:
    artifact = tmp_path / "fixture.dat"
    artifact.write_bytes(PAYLOAD)
    entry = ManifestEntry(
        name="md5-only",
        url="https://example.org/fixture.dat",
        md5=hashlib.md5(PAYLOAD).hexdigest(),
        sha256=None,
    )
    assert verify_download(artifact, entry) == artifact


def test_fetch_requires_injected_transport_and_verifies_output(tmp_path) -> None:
    entry = _entry()
    destination = tmp_path / "download.dat"
    with pytest.raises(NotImplementedError):
        fetch(entry, destination)

    def downloader(url, path) -> None:
        assert url == entry.url
        path.write_bytes(PAYLOAD)

    assert fetch(entry, destination, downloader) == destination
