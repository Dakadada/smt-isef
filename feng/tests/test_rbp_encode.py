"""Unit tests for the Arm B ENCODE RBP-layer loaders.

Synthetic fixtures only; no network and no staged multi-GB data required. The
real staged files live on scratch and are indexed by the committed manifest,
which is parsed here for schema and honesty invariants (roles, liftover flags).
"""

from __future__ import annotations

import gzip
import hashlib
import io
from pathlib import Path
import tarfile

import pandas as pd
import pytest

from feng.rbp import encode


REPO_ROOT = Path(__file__).resolve().parents[2]
MANIFEST = REPO_ROOT / "feng" / "artifacts" / "rbp" / "encode_manifest.tsv"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _source(
    role: str,
    *,
    assembly: str = encode.SOURCE_ASSEMBLY,
    sha256: str = "1" * 64,
) -> encode.EncodeFile:
    return encode.EncodeFile(
        collection="ENCOLL",
        accession="ENCFFTEST",
        role=role,
        target_rbp="HNRNPA1",
        cell_line="K562",
        file_format="bed" if role == "eclip_peaks" else "tar",
        output_type="peaks" if role == "eclip_peaks" else "rMATS",
        assembly=assembly,
        needs_liftover=True,
        md5sum="0" * 32,
        sha256=sha256,
        file_size=1,
        experiment="ENCSRTEST",
        url="https://example.test/file",
        curated_subset=True,
    )


def _write_eclip(path: Path) -> None:
    # two hg19 narrowPeak rows in eCLIP CLIPper form
    rows = [
        ["chr5", "70076000", "70076040", "HNRNPA1_peak_1", "200", "+", "3.5", "9.1", "-1", "20"],
        ["chr5", "70077000", "70077055", "HNRNPA1_peak_2", "150", "-", "2.1", "4.4", "-1", "30"],
    ]
    body = "\n".join("\t".join(r) for r in rows) + "\n"
    with gzip.open(path, "wt") as handle:
        handle.write(body)


def _write_rmats_tar(path: Path) -> None:
    header = [
        "ID", "GeneID", "geneSymbol", "chr", "strand",
        "PValue", "FDR", "IncLevel1", "IncLevel2", "IncLevelDifference",
    ]
    se = [
        ["1", "ENSG1", "SMN2", "chr5", "+", "0.001", "0.01", "0.8", "0.5", "0.30"],
        ["2", "ENSG2", "GENEB", "chr7", "-", "0.20", "0.40", "0.4", "0.5", "-0.10"],
    ]
    a3 = [["1", "ENSG3", "GENEC", "chr1", "+", "0.05", "0.10", "0.6", "0.7", "-0.10"]]
    def _member(tar: tarfile.TarFile, name: str, rows: list[list[str]]) -> None:
        text = "\t".join(header) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n"
        data = text.encode()
        info = tarfile.TarInfo(name)
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    with tarfile.open(path, "w:gz") as tar:
        _member(tar, "knockdown/SE.MATS.JC.txt", se)
        _member(tar, "knockdown/A3SS.MATS.JC.txt", a3)
        _member(tar, "knockdown/README.txt", [["ignore"]])


def test_read_eclip_peaks_schema(tmp_path: Path) -> None:
    bed = tmp_path / "ENCFFTEST.bed.gz"
    _write_eclip(bed)
    frame = encode.read_eclip_peaks(
        bed, source=_source("eclip_peaks", sha256=_sha256(bed))
    )
    assert list(frame.columns[: len(encode.ECLIP_PEAK_COLUMNS)]) == list(encode.ECLIP_PEAK_COLUMNS)
    assert len(frame) == 2
    assert frame["signal_value"].iloc[0] == pytest.approx(3.5)
    assert str(frame["start"].dtype) == "Int64"
    assert set(frame["source_assembly"]) == {encode.SOURCE_ASSEMBLY}
    assert set(frame["assembly"]) == {encode.SOURCE_ASSEMBLY}
    assert "target_assembly" not in frame.columns
    assert frame.attrs["assembly"] == encode.SOURCE_ASSEMBLY


def test_read_kd_splicing_rmats(tmp_path: Path) -> None:
    tar = tmp_path / "ENCFFRMATS.tar.gz"
    _write_rmats_tar(tar)
    tables = encode.read_kd_splicing_rmats(
        tar, source=_source("kd_splicing_rmats", sha256=_sha256(tar))
    )
    assert set(tables) == {"SE", "A3SS"}  # README ignored, only MATS.JC parsed
    se = tables["SE"]
    for col in encode._RMATS_REQUIRED_COLUMNS:
        assert col in se.columns
    # delta-PSI sign is preserved and numeric (direction only; no stat computed)
    assert se["IncLevelDifference"].iloc[0] == pytest.approx(0.30)
    assert se["IncLevelDifference"].iloc[1] == pytest.approx(-0.10)
    assert set(se["source_assembly"]) == {encode.SOURCE_ASSEMBLY}


def test_rmats_missing_columns_raise(tmp_path: Path) -> None:
    bad = tmp_path / "bad.tar.gz"
    with tarfile.open(bad, "w:gz") as tar:
        data = b"ID\tgeneSymbol\n1\tX\n"
        info = tarfile.TarInfo("SE.MATS.JC.txt")
        info.size = len(data)
        tar.addfile(info, io.BytesIO(data))
    with pytest.raises(ValueError, match="missing required columns"):
        encode.read_kd_splicing_rmats(
            bad, source=_source("kd_splicing_rmats", sha256=_sha256(bad))
        )


def test_encode_file_hard_rejects_non_hg19_assembly() -> None:
    with pytest.raises(ValueError, match="expected source assembly"):
        _source("eclip_peaks", assembly=encode.TARGET_ASSEMBLY)


def test_encode_file_validates_sha256() -> None:
    source = _source("eclip_peaks")
    values = {field: getattr(source, field) for field in source.__dataclass_fields__}
    values["sha256"] = "not-a-sha"
    with pytest.raises(ValueError, match="invalid sha256"):
        encode.EncodeFile(**values)


def test_coordinate_join_rejects_mixed_or_missing_assembly() -> None:
    hg19 = pd.DataFrame({"chrom": ["chr1"]})
    hg19.attrs["assembly"] = "hg19"
    grch38 = pd.DataFrame({"chrom": ["chr1"]})
    grch38.attrs["assembly"] = "GRCh38"
    with pytest.raises(ValueError, match="across assemblies"):
        encode.assert_same_assembly(hg19, grch38)
    with pytest.raises(ValueError, match="lacks assembly provenance"):
        encode.assert_same_assembly(hg19, pd.DataFrame())


def test_loaded_frame_merge_hard_rejects_mismatched_assembly(tmp_path: Path) -> None:
    bed = tmp_path / "ENCFFTEST.bed.gz"
    _write_eclip(bed)
    hg19 = encode.read_eclip_peaks(
        bed, source=_source("eclip_peaks", sha256=_sha256(bed))
    )
    grch38 = encode.coordinate_frame(
        pd.DataFrame({"chrom": ["chr5"], "annotation": ["v50"]}),
        assembly=encode.TARGET_ASSEMBLY,
    )
    with pytest.raises(ValueError, match="across assemblies"):
        hg19.merge(grch38, on="chrom")
    with pytest.raises(ValueError, match="across assemblies"):
        grch38.merge(hg19, on="chrom")


def test_reader_rejects_sha256_mismatch(tmp_path: Path) -> None:
    bed = tmp_path / "ENCFFTEST.bed.gz"
    _write_eclip(bed)
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        encode.read_eclip_peaks(bed, source=_source("eclip_peaks"))


@pytest.mark.skipif(not MANIFEST.exists(), reason="manifest artifact not present")
def test_manifest_parses_and_invariants() -> None:
    rows = encode.load_manifest(MANIFEST)
    assert rows
    roles = {r.role for r in rows}
    assert roles == {"eclip_peaks", "kd_splicing_rmats"}
    # every staged ENCODE file for this layer is hg19 -> liftover required, flagged
    assert all(r.assembly == encode.SOURCE_ASSEMBLY for r in rows)
    assert all(r.needs_liftover for r in rows)
    assert all(r.sha256 is None or len(r.sha256) == 64 for r in rows)
    curated = [r for r in rows if r.curated_subset]
    assert all(r.sha256 is not None for r in curated)
    combos = {(r.target_rbp, r.cell_line) for r in curated}
    # curated subset carries both modalities for each RBP x cell combo
    for rbp, cell in combos:
        modes = {r.role for r in curated if (r.target_rbp, r.cell_line) == (rbp, cell)}
        assert modes == {"eclip_peaks", "kd_splicing_rmats"}
