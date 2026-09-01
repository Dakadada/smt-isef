"""Regression tests for GENCODE provenance and STAR completion contracts."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
import subprocess

import pytest

from feng.data.manifest import (
    GENCODE_V50_PRIMARY_ASSEMBLY_FASTA,
    GENCODE_V50_PRIMARY_ASSEMBLY_GTF,
)


FENG_ROOT = Path(__file__).resolve().parents[1]
FETCH_SCRIPT = FENG_ROOT / "scripts" / "fetch_gencode.sh"
STAR_SCRIPT = FENG_ROOT / "scripts" / "star_index.sbatch"


def _decompress(archive: Path, output: Path, archive_md5: str, format_name: str):
    return subprocess.run(
        [
            "bash",
            "-c",
            'source "$1"; decompress_atomic "$2" "$3" "$4" "$5"',
            "bash",
            str(FETCH_SCRIPT),
            str(archive),
            str(output),
            archive_md5,
            format_name,
        ],
        capture_output=True,
        text=True,
        check=False,
    )


@pytest.mark.parametrize(
    ("format_name", "payload"),
    [
        ("fasta", b">chrToy\nACGT\n"),
        ("gtf", b'chrToy\ttest\tgene\t1\t4\t.\t+\t.\tgene_id "toy";\n'),
    ],
)
def test_decompressed_artifact_is_bound_to_archive_and_content(
    tmp_path: Path,
    format_name: str,
    payload: bytes,
) -> None:
    archive = tmp_path / f"input.{format_name}.gz"
    output = tmp_path / f"input.{format_name}"
    archive.write_bytes(gzip.compress(payload, mtime=0))
    archive_md5 = hashlib.md5(archive.read_bytes()).hexdigest()

    first = _decompress(archive, output, archive_md5, format_name)
    assert first.returncode == 0, first.stderr
    assert output.read_bytes() == payload
    provenance = output.with_name(output.name + ".provenance")
    assert f"archive_md5={archive_md5}" in provenance.read_text(encoding="utf-8")

    output.write_bytes(payload + b"tampered\n")
    provenance.write_text(
        provenance.read_text(encoding="utf-8").replace(archive_md5, "f" * 32),
        encoding="utf-8",
    )
    second = _decompress(archive, output, archive_md5, format_name)
    assert second.returncode == 0, second.stderr
    assert output.read_bytes() == payload
    assert f"archive_md5={archive_md5}" in provenance.read_text(encoding="utf-8")


def test_invalid_decompressed_format_is_never_published(tmp_path: Path) -> None:
    archive = tmp_path / "invalid.fa.gz"
    output = tmp_path / "invalid.fa"
    archive.write_bytes(gzip.compress(b"not a FASTA\n", mtime=0))
    result = _decompress(archive, output, hashlib.md5(archive.read_bytes()).hexdigest(), "fasta")
    assert result.returncode != 0
    assert not output.exists()


def test_fetch_script_and_manifest_share_authoritative_md5_values() -> None:
    script = FETCH_SCRIPT.read_text(encoding="utf-8")
    assert GENCODE_V50_PRIMARY_ASSEMBLY_FASTA.md5 in script
    assert GENCODE_V50_PRIMARY_ASSEMBLY_GTF.md5 in script
    assert "/home/dashii/" not in script


def test_star_completion_contract_is_full_and_input_bound() -> None:
    script = STAR_SCRIPT.read_text(encoding="utf-8")
    required_files = {
        "Genome",
        "SA",
        "SAindex",
        "chrLength.txt",
        "chrName.txt",
        "chrNameLength.txt",
        "chrStart.txt",
        "exonGeTrInfo.tab",
        "exonInfo.tab",
        "geneInfo.tab",
        "genomeParameters.txt",
        "sjdbInfo.txt",
        "sjdbList.fromGTF.out.tab",
        "sjdbList.out.tab",
        "transcriptInfo.tab",
    }
    assert all(name in script for name in required_files)
    assert all(
        field in script
        for field in (
            "feng_index_manifest.json",
            "fasta_md5",
            "gtf_md5",
            "star_version",
            "run_threads",
            "sjdb_overhang",
            "limit_genome_generate_ram",
        )
    )
    assert 'index_files_complete "${FENG_STAR_INDEX_DIR}" && marker_matches' in script
    assert "/home/dashii/" not in script
