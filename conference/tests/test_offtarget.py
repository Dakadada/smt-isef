from pathlib import Path
import io
import tracemalloc

import pytest

from conference.aso_pipeline.offtarget import (
    GeneIndex,
    expected_random_hits,
    map_alignment_gene_keys,
    parse_blast_line,
    parse_bowtie_line,
    reduce_sorted_gene_keys,
)


def test_expected_random_hit_estimate_exposes_short_mer_saturation():
    genome_size = 3_088_269_832
    assert expected_random_hits(8, 3, genome_size) == pytest.approx(
        168_607_016.51, rel=1e-9
    )
    assert expected_random_hits(10, 3, genome_size) == pytest.approx(
        21_653_136.96, rel=1e-9
    )
    assert expected_random_hits(18, 3, genome_size) == pytest.approx(
        2_108.96, abs=0.01
    )


def test_alignment_parsers_convert_to_same_zero_based_reverse_hit():
    bowtie = parse_bowtie_line("q1\t-\tchr1\t10\n", {"q1": 10})
    blast = parse_blast_line("q1\tchr1\t20\t11\n")
    assert (bowtie.start_0base, bowtie.end_0base, bowtie.query_strand) == (
        10,
        20,
        "-",
    )
    assert blast == bowtie


def test_gene_index_reports_any_and_transcript_compatible_overlaps(tmp_path: Path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        "chr1\ttest\tgene\t11\t30\t.\t+\t.\tgene_id \"PLUS\"; gene_name \"P\";\n"
        "chr1\ttest\tgene\t16\t40\t.\t-\t.\tgene_id \"MINUS\"; gene_name \"M\";\n"
    )
    index = GeneIndex.from_gtf(gtf)
    reverse_hit = parse_bowtie_line("q1\t-\tchr1\t10\n", {"q1": 10})
    forward_hit = parse_bowtie_line("q2\t+\tchr1\t15\n", {"q2": 10})
    assert index.overlaps(reverse_hit) == ({"PLUS", "MINUS"}, {"PLUS"})
    assert index.overlaps(forward_hit) == ({"PLUS", "MINUS"}, {"MINUS"})


def test_external_gene_key_pipeline_preserves_raw_and_distinct_counts(tmp_path: Path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        "chr1\ttest\tgene\t11\t30\t.\t+\t.\tgene_id \"PLUS\";\n"
        "chr1\ttest\tgene\t16\t40\t.\t-\t.\tgene_id \"MINUS\";\n"
    )
    raw_counts = tmp_path / "raw.tsv"
    keys = io.StringIO()
    map_alignment_gene_keys(
        ["q1\t-\tchr1\t10\n", "q1\t-\tchr1\t10\n"],
        output=keys,
        raw_counts_path=raw_counts,
        aligner="bowtie1-v3",
        gene_index=GeneIndex.from_gtf(gtf),
        query_lengths={"q1": 10},
    )
    counts = reduce_sorted_gene_keys(
        sorted(set(keys.getvalue().splitlines())),
        raw_counts_path=raw_counts,
        query_lengths={"q1": 10},
    )
    assert counts == {
        "q1": {
            "raw_alignment_count": 2,
            "distinct_genes_any_strand": 2,
            "distinct_genes_transcript_compatible": 1,
        }
    }


def test_gene_key_mapper_memory_does_not_grow_with_duplicate_hits(tmp_path: Path):
    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        "chr1\ttest\tgene\t1\t100\t.\t+\t.\tgene_id \"PLUS\";\n"
    )

    class NullWriter:
        def write(self, value: str) -> int:
            return len(value)

    raw_counts = tmp_path / "raw.tsv"
    lines = ("q1\t-\tchr1\t10\n" for _ in range(100_000))
    tracemalloc.start()
    map_alignment_gene_keys(
        lines,
        output=NullWriter(),
        raw_counts_path=raw_counts,
        aligner="bowtie1-v3",
        gene_index=GeneIndex.from_gtf(gtf),
        query_lengths={"q1": 10},
    )
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert peak < 5_000_000
    assert raw_counts.read_text().splitlines()[-1] == "q1\t100000"


def test_reducer_waits_for_key_stream_before_opening_mapper_sidecar(tmp_path: Path):
    raw_counts = tmp_path / "raw.tsv"

    def keys_then_sidecar():
        raw_counts.write_text("unique_query_id\traw_alignment_count\nq1\t2\n")
        yield "q1\tA\tGENE1\n"

    counts = reduce_sorted_gene_keys(
        keys_then_sidecar(), raw_counts_path=raw_counts, query_lengths={"q1": 10}
    )
    assert counts["q1"]["raw_alignment_count"] == 2
