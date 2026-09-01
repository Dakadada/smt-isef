"""GATE-1 subprocess smoke tests for the clean-room Feng arm toolchain."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess

import pytest


FIXTURES = Path(__file__).resolve().parent / "fixtures"
GENOME = FIXTURES / "tiny_genome.fa"
GTF = FIXTURES / "tiny.gtf"
READS = FIXTURES / "tiny_R1.fastq"


def _locate_tool(name: str) -> str | None:
    prefix = os.environ.get("FENG_ENV_PREFIX")
    if prefix:
        candidate = Path(prefix) / "bin" / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)
    return shutil.which(name)


STAR = _locate_tool("STAR")
INTARNA = _locate_tool("IntaRNA")
RMATS = _locate_tool("rmats.py")


def _run(command: list[str], *, cwd: Path, timeout: int) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["LC_ALL"] = "C"
    return subprocess.run(
        command,
        cwd=cwd,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
        timeout=timeout,
    )


def _assert_success(result: subprocess.CompletedProcess[str], command_name: str) -> None:
    assert result.returncode == 0, (
        f"{command_name} failed with exit {result.returncode}\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )


def _build_tiny_star_alignment(work_dir: Path) -> Path:
    index_dir = work_dir / "index"
    index_dir.mkdir(parents=True)
    index_result = _run(
        [
            str(STAR),
            "--runMode",
            "genomeGenerate",
            "--runThreadN",
            "1",
            "--genomeDir",
            str(index_dir),
            "--genomeFastaFiles",
            str(GENOME),
            "--sjdbGTFfile",
            str(GTF),
            "--sjdbOverhang",
            "49",
            "--genomeSAindexNbases",
            "3",
            "--genomeChrBinNbits",
            "10",
            "--limitGenomeGenerateRAM",
            "1000000000",
        ],
        cwd=work_dir,
        timeout=60,
    )
    _assert_success(index_result, "STAR genomeGenerate")
    for required in ("Genome", "SA", "SAindex", "genomeParameters.txt"):
        assert (index_dir / required).stat().st_size > 0

    output_prefix = work_dir / "align."
    alignment_result = _run(
        [
            str(STAR),
            "--runThreadN",
            "1",
            "--genomeDir",
            str(index_dir),
            "--readFilesIn",
            str(READS),
            "--outFileNamePrefix",
            str(output_prefix),
            "--outSAMtype",
            "BAM",
            "SortedByCoordinate",
            "--alignIntronMin",
            "20",
            "--alignIntronMax",
            "2000",
        ],
        cwd=work_dir,
        timeout=60,
    )
    _assert_success(alignment_result, "STAR alignment")

    bam = work_dir / "align.Aligned.sortedByCoord.out.bam"
    final_log = work_dir / "align.Log.final.out"
    assert bam.stat().st_size > 0
    log_text = final_log.read_text(encoding="utf-8")
    uniquely_mapped = next(
        int(line.split("|", 1)[1].strip())
        for line in log_text.splitlines()
        if line.strip().startswith("Uniquely mapped reads number")
    )
    assert uniquely_mapped >= 1
    return bam


def _read_single_fasta_sequence(path: Path) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines and lines[0].startswith(">")
    return "".join(lines[1:])


@pytest.mark.skipif(STAR is None, reason="STAR is absent; set FENG_ENV_PREFIX or PATH")
def test_star_tiny_index_and_alignment(tmp_path: Path) -> None:
    bam = _build_tiny_star_alignment(tmp_path)
    assert bam.name == "align.Aligned.sortedByCoord.out.bam"


@pytest.mark.skipif(INTARNA is None, reason="IntaRNA is absent; set FENG_ENV_PREFIX or PATH")
def test_intarna_one_shot(tmp_path: Path) -> None:
    target = _read_single_fasta_sequence(FIXTURES / "intarna_target.fa")
    query = _read_single_fasta_sequence(FIXTURES / "intarna_query.fa")
    output = tmp_path / "interaction.csv"
    result = _run(
        [
            str(INTARNA),
            f"--target={target}",
            f"--query={query}",
            "--outMode=C",
            f"--out={output}",
        ],
        cwd=tmp_path,
        timeout=30,
    )
    _assert_success(result, "IntaRNA")
    lines = [line for line in output.read_text(encoding="utf-8").splitlines() if line]
    assert len(lines) >= 2
    assert ";" in lines[0]
    assert ";" in lines[1]


@pytest.mark.skipif(
    STAR is None or RMATS is None,
    reason="STAR or rmats.py is absent; set FENG_ENV_PREFIX or PATH",
)
def test_rmats_one_shot_from_tiny_star_bams(tmp_path: Path) -> None:
    source_bam = _build_tiny_star_alignment(tmp_path / "alignment")
    group1_bam = tmp_path / "group1.bam"
    group2_bam = tmp_path / "group2.bam"
    shutil.copyfile(source_bam, group1_bam)
    shutil.copyfile(source_bam, group2_bam)

    b1 = tmp_path / "b1.txt"
    b2 = tmp_path / "b2.txt"
    b1.write_text(str(group1_bam) + "\n", encoding="utf-8")
    b2.write_text(str(group2_bam) + "\n", encoding="utf-8")
    output_dir = tmp_path / "rmats-output"
    temporary_dir = tmp_path / "rmats-tmp"
    output_dir.mkdir()
    temporary_dir.mkdir()

    result = _run(
        [
            str(RMATS),
            "--b1",
            str(b1),
            "--b2",
            str(b2),
            "--gtf",
            str(GTF),
            "--od",
            str(output_dir),
            "--tmp",
            str(temporary_dir),
            "--readLength",
            "50",
            "--nthread",
            "1",
            "-t",
            "single",
            "--task",
            "both",
            "--statoff",
        ],
        cwd=tmp_path,
        timeout=120,
    )
    _assert_success(result, "rmats.py")

    event_table = output_dir / "fromGTF.SE.txt"
    rows = event_table.read_text(encoding="utf-8").splitlines()
    assert len(rows) >= 2
    header = rows[0].split("\t")
    assert {"ID", "GeneID", "chr", "strand"}.issubset(header)

