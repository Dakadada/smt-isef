"""Typed download-manifest schema for the clean-room Feng arm."""

from __future__ import annotations

from dataclasses import dataclass
import re
from urllib.parse import urlparse


def _require_hex(value: str, length: int, label: str) -> str:
    normalized = value.casefold()
    if not re.fullmatch(rf"[0-9a-f]{{{length}}}", normalized):
        raise ValueError(f"{label} must be exactly {length} hexadecimal characters")
    if normalized == "0" * length:
        raise ValueError(f"{label} must be recorded; all-zero digests are forbidden")
    return normalized


@dataclass(frozen=True, slots=True)
class ManifestEntry:
    """One artifact with authoritative MD5 and optional secondary SHA-256."""

    name: str
    url: str
    md5: str
    sha256: str | None

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("manifest entry name must not be empty")
        parsed = urlparse(self.url)
        if parsed.scheme not in {"https", "http", "ftp"} or not parsed.netloc:
            raise ValueError("manifest URL must be an absolute HTTP(S) or FTP URL")
        object.__setattr__(self, "md5", _require_hex(self.md5, 32, "md5"))
        if self.sha256 is not None:
            object.__setattr__(self, "sha256", _require_hex(self.sha256, 64, "sha256"))


@dataclass(frozen=True, slots=True)
class Manifest:
    """A uniquely named collection of clean-room source artifacts."""

    entries: tuple[ManifestEntry, ...]

    def __post_init__(self) -> None:
        names = [entry.name for entry in self.entries]
        if len(names) != len(set(names)):
            raise ValueError("manifest entry names must be unique")


# GENCODE publishes these MD5 values in the release 50 checksum manifest. SHA-256
# is optional secondary metadata and may be recorded only after MD5 verification.
GENCODE_V50_PRIMARY_ASSEMBLY_FASTA = ManifestEntry(
    name="gencode_v50_primary_assembly_fasta",
    url=(
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/"
        "GRCh38.primary_assembly.genome.fa.gz"
    ),
    md5="da1a11258be075cfa7af718162c894e7",
    sha256=None,
)

GENCODE_V50_PRIMARY_ASSEMBLY_GTF = ManifestEntry(
    name="gencode_v50_primary_assembly_gtf",
    url=(
        "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/"
        "gencode.v50.primary_assembly.annotation.gtf.gz"
    ),
    md5="289b91e5e95e8b0450d223246f10a12e",
    sha256=None,
)

GENCODE_V50_MANIFEST = Manifest(
    entries=(
        GENCODE_V50_PRIMARY_ASSEMBLY_FASTA,
        GENCODE_V50_PRIMARY_ASSEMBLY_GTF,
    )
)
