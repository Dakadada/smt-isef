"""Streaming alignment parsing and strand-aware GENCODE gene overlap."""

from __future__ import annotations

from bisect import bisect_left
from dataclasses import dataclass
import gzip
from math import comb
from pathlib import Path
import re
from typing import Iterable, TextIO


@dataclass(frozen=True, slots=True)
class Hit:
    query_id: str
    chrom: str
    start_0base: int
    end_0base: int
    query_strand: str


@dataclass(frozen=True, slots=True)
class GeneInterval:
    start_0base: int
    end_0base: int
    strand: str
    gene_id: str


def expected_random_hits(length: int, mismatches: int, genome_size: int) -> float:
    """IID equal-base expectation across both genomic orientations."""
    if length < 1 or mismatches < 0 or mismatches > length or genome_size < length:
        raise ValueError("invalid length, mismatch count, or genome size")
    one_orientation = sum(
        comb(length, mismatch_count) * 3**mismatch_count
        for mismatch_count in range(mismatches + 1)
    ) / 4**length
    return 2.0 * (genome_size - length + 1) * one_orientation


def parse_bowtie_line(line: str, query_lengths: dict[str, int]) -> Hit:
    """Parse Bowtie's four-column output after `--suppress 5,6,7,8`."""
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 4 or fields[1] not in {"+", "-"}:
        raise ValueError(f"invalid Bowtie record: {line!r}")
    query_id, strand, chrom, offset = fields
    if query_id not in query_lengths:
        raise ValueError(f"unknown Bowtie query: {query_id}")
    start = int(offset)
    return Hit(query_id, chrom, start, start + query_lengths[query_id], strand)


def parse_blast_line(line: str) -> Hit:
    """Parse `qseqid sseqid sstart send` BLAST tabular output."""
    fields = line.rstrip("\n").split("\t")
    if len(fields) != 4:
        raise ValueError(f"invalid BLAST record: {line!r}")
    query_id, chrom, raw_start, raw_end = fields
    sstart, send = int(raw_start), int(raw_end)
    strand = "+" if sstart <= send else "-"
    start_1based, end_1based = sorted((sstart, send))
    return Hit(query_id, chrom, start_1based - 1, end_1based, strand)


class GeneIndex:
    """Small in-memory interval index over GTF `gene` spans."""

    def __init__(self, intervals: Iterable[tuple[str, GeneInterval]]):
        grouped: dict[str, list[GeneInterval]] = {}
        for chrom, interval in intervals:
            grouped.setdefault(chrom, []).append(interval)
        self._data = {}
        for chrom, values in grouped.items():
            values.sort(key=lambda interval: interval.start_0base)
            starts = [interval.start_0base for interval in values]
            prefix_max: list[int] = []
            maximum = -1
            for interval in values:
                maximum = max(maximum, interval.end_0base)
                prefix_max.append(maximum)
            self._data[chrom] = (values, starts, prefix_max)

    @classmethod
    def from_gtf(cls, path: Path) -> "GeneIndex":
        path = Path(path)
        opener = gzip.open if path.suffix == ".gz" else open
        gene_id_pattern = re.compile(r'(?:^|;\s*)gene_id "([^"]+)"')
        intervals = []
        with opener(path, "rt") as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9 or fields[2] != "gene":
                    continue
                match = gene_id_pattern.search(fields[8])
                if match is None:
                    raise ValueError(f"gene feature lacks gene_id: {line[:200]!r}")
                intervals.append(
                    (
                        fields[0],
                        GeneInterval(
                            start_0base=int(fields[3]) - 1,
                            end_0base=int(fields[4]),
                            strand=fields[6],
                            gene_id=match.group(1),
                        ),
                    )
                )
        if not intervals:
            raise ValueError(f"no gene features found in {path}")
        return cls(intervals)

    @classmethod
    def from_table(cls, path: Path) -> "GeneIndex":
        intervals = []
        with Path(path).open() as handle:
            header = next(handle).rstrip("\n").split("\t")
            expected = ["chrom", "start_0base", "end_0base", "strand", "gene_id"]
            if header != expected:
                raise ValueError(f"unexpected gene-table header: {header}")
            for line in handle:
                chrom, start, end, strand, gene_id = line.rstrip("\n").split("\t")
                intervals.append(
                    (chrom, GeneInterval(int(start), int(end), strand, gene_id))
                )
        return cls(intervals)

    def write_table(self, path: Path) -> None:
        path = Path(path)
        temporary = path.with_suffix(path.suffix + ".tmp")
        with temporary.open("w") as handle:
            handle.write("chrom\tstart_0base\tend_0base\tstrand\tgene_id\n")
            for chrom in sorted(self._data):
                for interval in self._data[chrom][0]:
                    handle.write(
                        f"{chrom}\t{interval.start_0base}\t{interval.end_0base}\t"
                        f"{interval.strand}\t{interval.gene_id}\n"
                    )
        temporary.replace(path)

    def overlaps(self, hit: Hit) -> tuple[set[str], set[str]]:
        data = self._data.get(hit.chrom)
        if data is None:
            return set(), set()
        intervals, starts, prefix_max = data
        index = bisect_left(starts, hit.end_0base) - 1
        any_strand: set[str] = set()
        compatible: set[str] = set()
        required_gene_strand = "+" if hit.query_strand == "-" else "-"
        while index >= 0 and prefix_max[index] > hit.start_0base:
            interval = intervals[index]
            if interval.end_0base > hit.start_0base:
                any_strand.add(interval.gene_id)
                if interval.strand == required_gene_strand:
                    compatible.add(interval.gene_id)
            index -= 1
        return any_strand, compatible


def _write_raw_counts(path: Path, raw_counts: dict[str, int]) -> None:
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w") as handle:
        handle.write("unique_query_id\traw_alignment_count\n")
        for query_id in sorted(raw_counts):
            handle.write(f"{query_id}\t{raw_counts[query_id]}\n")
    temporary.replace(path)


def map_alignment_gene_keys(
    lines: Iterable[str],
    *,
    output: TextIO,
    raw_counts_path: Path,
    aligner: str,
    gene_index: GeneIndex,
    query_lengths: dict[str, int],
) -> None:
    """Stream hits to external-sort keys using memory independent of hit count.

    The emitted keys are ``query_id, A|C, gene_id`` for any-strand and
    transcript-compatible overlaps. They intentionally retain duplicates;
    GNU sort performs disk-backed ``-u`` downstream. Only one integer raw hit
    counter per query is held in memory.
    """
    if aligner not in {"bowtie1-v3", "blastn-short"}:
        raise ValueError(f"unsupported aligner: {aligner}")
    parser = (
        (lambda line: parse_bowtie_line(line, query_lengths))
        if aligner == "bowtie1-v3"
        else parse_blast_line
    )
    raw_counts = {query_id: 0 for query_id in query_lengths}
    for line in lines:
        if not line.strip():
            continue
        hit = parser(line)
        if hit.query_id not in raw_counts:
            raise ValueError(f"alignment contains unknown query: {hit.query_id}")
        raw_counts[hit.query_id] += 1
        any_overlap, compatible_overlap = gene_index.overlaps(hit)
        for gene_id in any_overlap:
            output.write(f"{hit.query_id}\tA\t{gene_id}\n")
        for gene_id in compatible_overlap:
            output.write(f"{hit.query_id}\tC\t{gene_id}\n")
    _write_raw_counts(raw_counts_path, raw_counts)


def reduce_sorted_gene_keys(
    lines: Iterable[str],
    *,
    raw_counts_path: Path,
    query_lengths: dict[str, int],
) -> dict[str, dict[str, int]]:
    """Count externally sorted unique keys with O(number of queries) memory."""
    any_counts = {query_id: 0 for query_id in query_lengths}
    compatible_counts = {query_id: 0 for query_id in query_lengths}
    previous = None
    for line in lines:
        key = line.rstrip("\n")
        if not key or key == previous:
            continue
        previous = key
        fields = key.split("\t")
        if len(fields) != 3:
            raise ValueError(f"invalid gene key: {line!r}")
        query_id, overlap_class, _gene_id = fields
        if query_id not in query_lengths:
            raise ValueError(f"gene key contains unknown query: {query_id}")
        if overlap_class == "A":
            any_counts[query_id] += 1
        elif overlap_class == "C":
            compatible_counts[query_id] += 1
        else:
            raise ValueError(f"invalid overlap class: {overlap_class}")

    # In the live pipeline the mapper publishes this atomic sidecar before
    # closing its stdout. Reading it only after the sorted-key stream reaches
    # EOF avoids a reducer-startup race while retaining streaming behavior.
    raw_counts: dict[str, int] = {}
    with Path(raw_counts_path).open() as handle:
        header = next(handle).rstrip("\n").split("\t")
        if header != ["unique_query_id", "raw_alignment_count"]:
            raise ValueError(f"unexpected raw-count header: {header}")
        for line in handle:
            query_id, count = line.rstrip("\n").split("\t")
            raw_counts[query_id] = int(count)
    if set(raw_counts) != set(query_lengths):
        raise ValueError("raw-count queries differ from query FASTA")
    return {
        query_id: {
            "raw_alignment_count": raw_counts[query_id],
            "distinct_genes_any_strand": any_counts[query_id],
            "distinct_genes_transcript_compatible": compatible_counts[query_id],
        }
        for query_id in query_lengths
    }
