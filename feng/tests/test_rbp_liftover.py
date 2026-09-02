"""Unit tests for the Arm B hg19->GRCh38 liftover bridge.

Synthetic fixtures only; no network and no staged chain required. A minimal but
format-valid UCSC chain drives real coordinate arithmetic through the same parser
that would consume the production ``hg19ToHg38.over.chain.gz`` -- coordinates are
computed, never hard-coded into the frame.
"""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path

import pandas as pd
import pytest

from feng.rbp import encode, liftover


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# --- chain fixtures ---------------------------------------------------------
# A single ungapped block: hg19 chr1 [1000,2000) -> hg38 chr1 [1100,2100).
# Query strand '+', so a target position p maps to p + 100.
_PLUS_CHAIN = "chain 4900 chr1 249250621 + 1000 2000 chr1 248956422 + 1100 2100 1\n1000\n\n"

# A minus-strand chain: hg19 chr3 [1000,2000) -> hg38 chr3 on the '-' strand.
_MINUS_CHAIN = "chain 4900 chr3 198022430 + 1000 2000 chr3 198295559 - 1000 2000 2\n1000\n\n"

# Two chains covering the SAME hg19 target region -> multi-mapping.
_MULTI_CHAIN = (
    "chain 4900 chr1 249250621 + 1000 2000 chr1 248956422 + 1100 2100 1\n1000\n\n"
    "chain 4800 chr1 249250621 + 1000 2000 chr2 242193529 + 5000 6000 2\n1000\n\n"
)


def _write_chain(path: Path, text: str, *, gz: bool = False) -> Path:
    if gz:
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)
    return path


def _write_eclip(path: Path, rows: list[tuple]) -> None:
    lines = [
        "\t".join(str(v) for v in row)
        for row in rows
    ]
    path.write_text("\n".join(lines) + "\n")


def _eclip_source(bed: Path) -> encode.EncodeFile:
    return encode.EncodeFile(
        collection="ENCOLL",
        accession="ENCFFLIFT",
        role="eclip_peaks",
        target_rbp="HNRNPA1",
        cell_line="HepG2",
        file_format="bed",
        output_type="peaks",
        assembly=encode.SOURCE_ASSEMBLY,
        needs_liftover=True,
        md5sum="0" * 32,
        sha256=_sha256(bed),
        file_size=bed.stat().st_size,
        experiment="ENCSR456FVU",
        url="https://example.invalid/peaks.bed",
        curated_subset=True,
    )


def _hg19_eclip_frame(tmp_path: Path, rows: list[tuple]) -> pd.DataFrame:
    bed = tmp_path / "peaks.bed"
    _write_eclip(bed, rows)
    return encode.read_eclip_peaks(bed, source=_eclip_source(bed))


# --- chain parsing / interval math -----------------------------------------
def test_chain_parses_and_lifts_known_interval(tmp_path: Path) -> None:
    chain = liftover.load_chain_file(_write_chain(tmp_path / "c.chain", _PLUS_CHAIN))
    assert chain.n_chains == 1
    # hg19 chr1:1500-1520 -> hg38 chr1:1600-1620 (offset +100), strand preserved.
    results = chain.lift_interval("chr1", 1500, 1520)
    assert results == [("chr1", 1600, 1620, "+")]


def test_gzipped_chain_loads(tmp_path: Path) -> None:
    chain = liftover.load_chain_file(
        _write_chain(tmp_path / "c.chain.gz", _PLUS_CHAIN, gz=True)
    )
    assert chain.lift_interval("chr1", 1500, 1520) == [("chr1", 1600, 1620, "+")]


def test_minus_strand_chain_reorders_and_flips(tmp_path: Path) -> None:
    chain = liftover.load_chain_file(_write_chain(tmp_path / "m.chain", _MINUS_CHAIN))
    results = chain.lift_interval("chr3", 1500, 1520)
    assert len(results) == 1
    chrom, start, end, strand = results[0]
    assert chrom == "chr3"
    assert end - start == 20  # length preserved across the strand flip
    assert start < end  # endpoints reordered into a valid half-open interval
    assert strand == "-"


def test_unmapped_interval_returns_empty(tmp_path: Path) -> None:
    chain = liftover.load_chain_file(_write_chain(tmp_path / "c.chain", _PLUS_CHAIN))
    assert chain.lift_interval("chr1", 5000, 5020) == []  # outside target span
    assert chain.lift_interval("chrX", 1500, 1520) == []  # unknown contig


def test_missing_chain_file_raises(tmp_path: Path) -> None:
    with pytest.raises(liftover.LiftoverError, match="chain file not found"):
        liftover.load_chain_file(tmp_path / "nope.chain")


# --- frame-level lift + provenance -----------------------------------------
def test_lift_eclip_peaks_stamps_grch38_and_counts(tmp_path: Path) -> None:
    frame = _hg19_eclip_frame(
        tmp_path,
        [
            ("chr1", 1500, 1520, "p_ok", 1, "+", 2.0, 5.0, -1, 10),      # lifts
            ("chr1", 5000, 5020, "p_unmapped", 1, "+", 2.0, 5.0, -1, 10),  # unmapped
        ],
    )
    chain = liftover.load_chain_file(_write_chain(tmp_path / "c.chain", _PLUS_CHAIN))
    result = liftover.lift_eclip_peaks(frame, chain)

    assert result.n_input == 2
    assert result.n_lifted == 1
    assert result.n_unmapped == 1
    assert result.n_multimapped == 0

    lifted = result.frame
    assert set(lifted["assembly"]) == {encode.TARGET_ASSEMBLY}
    assert set(lifted["source_assembly"]) == {encode.SOURCE_ASSEMBLY}
    assert set(lifted["target_assembly"]) == {encode.TARGET_ASSEMBLY}
    assert not lifted["needs_liftover"].any()
    assert lifted.attrs["assembly"] == encode.TARGET_ASSEMBLY
    # coordinates actually moved
    assert lifted.iloc[0]["start"] == 1600
    assert lifted.iloc[0]["end"] == 1620
    # chain provenance recorded
    assert set(lifted["liftover_chain_sha256"]) == {chain.sha256}
    assert set(lifted["liftover_chain_path"]) == {str(chain.path)}


def test_multimapped_interval_dropped_and_counted(tmp_path: Path) -> None:
    frame = _hg19_eclip_frame(
        tmp_path, [("chr1", 1500, 1520, "p", 1, "+", 2.0, 5.0, -1, 10)]
    )
    chain = liftover.load_chain_file(_write_chain(tmp_path / "multi.chain", _MULTI_CHAIN))
    result = liftover.lift_eclip_peaks(frame, chain)
    assert result.n_lifted == 0
    assert result.n_multimapped == 1
    assert result.n_unmapped == 0


def test_lift_refuses_non_hg19_frame(tmp_path: Path) -> None:
    grch38 = encode.coordinate_frame(
        pd.DataFrame({"chrom": ["chr1"], "start": [1500], "end": [1520]}),
        assembly=encode.TARGET_ASSEMBLY,
    )
    chain = liftover.load_chain_file(_write_chain(tmp_path / "c.chain", _PLUS_CHAIN))
    with pytest.raises(liftover.LiftoverError, match="not uniformly hg19|only 'hg19'"):
        liftover.lift_eclip_peaks(grch38, chain)


# --- the fence: lifted frame ACCEPTED, raw hg19 frame still REJECTED --------
def test_lifted_frame_joins_grch38_through_fence(tmp_path: Path) -> None:
    frame = _hg19_eclip_frame(
        tmp_path, [("chr1", 1500, 1520, "p", 1, "+", 2.0, 5.0, -1, 10)]
    )
    chain = liftover.load_chain_file(_write_chain(tmp_path / "c.chain", _PLUS_CHAIN))
    lifted = liftover.lift_eclip_peaks(frame, chain).frame

    universe = encode.coordinate_frame(
        pd.DataFrame({"chrom": ["chr1"], "exon_key": ["chr1:1595-1625:+"]}),
        assembly=encode.TARGET_ASSEMBLY,
    )
    merged = lifted.merge(universe, on="chrom")  # fence must ALLOW this
    assert len(merged) == 1
    assert merged.attrs["assembly"] == encode.TARGET_ASSEMBLY


def test_raw_hg19_frame_still_rejected_by_fence(tmp_path: Path) -> None:
    frame = _hg19_eclip_frame(
        tmp_path, [("chr1", 1500, 1520, "p", 1, "+", 2.0, 5.0, -1, 10)]
    )
    universe = encode.coordinate_frame(
        pd.DataFrame({"chrom": ["chr1"], "exon_key": ["chr1:1595-1625:+"]}),
        assembly=encode.TARGET_ASSEMBLY,
    )
    with pytest.raises(ValueError, match="across assemblies"):
        frame.merge(universe, on="chrom")  # unlifted hg19 must STILL fail
