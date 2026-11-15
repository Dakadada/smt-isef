from io import StringIO
from pathlib import Path
import subprocess

import pytest

from conference.aso_pipeline.offtarget import (
    GeneIndex,
    map_alignment_gene_keys,
    reduce_sorted_gene_keys,
)


TOOL_DIR = Path("/orcd/scratch/orcd/011/dashii/conference_smn2_aso/tools/bin")
BOWTIE = TOOL_DIR / "bowtie"
BOWTIE_BUILD = TOOL_DIR / "bowtie-build"
STREAM_WRAPPER = Path("conference/scripts/bowtie_stream.sh")


@pytest.mark.skipif(
    not BOWTIE.is_file() or not BOWTIE_BUILD.is_file(), reason="Bowtie 1 unavailable"
)
def test_real_bowtie_stdout_reaches_mapper_and_known_genes(tmp_path: Path):
    query = "TCACTTTCATAATGCTGG"
    reverse_complement = "CCAGCATTATGAAAGTGA"
    genome_sequence = (
        "N" * 10
        + query
        + "N" * 10
        + reverse_complement
        + "N" * 10
        + query
        + "N" * 10
    )
    genome = tmp_path / "genome.fa"
    genome.write_text(f">chr1\n{genome_sequence}\n")
    index = tmp_path / "tiny"
    subprocess.run(
        [str(BOWTIE_BUILD), str(genome), str(index)],
        check=True,
        capture_output=True,
        text=True,
    )
    query_fasta = tmp_path / "query.fa"
    query_fasta.write_text(f">known_query\n{query}\n")

    aligned = subprocess.run(
        [
            "bash",
            str(STREAM_WRAPPER),
            str(BOWTIE),
            str(index),
            str(query_fasta),
            "1",
            "3",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    records = aligned.stdout.splitlines()
    assert len(records) == 3
    assert all(len(record.split("\t")) == 4 for record in records)

    gtf = tmp_path / "genes.gtf"
    gtf.write_text(
        "chr1\ttest\tgene\t11\t28\t.\t-\t.\tgene_id \"MINUS_A\";\n"
        "chr1\ttest\tgene\t39\t56\t.\t+\t.\tgene_id \"PLUS\";\n"
        "chr1\ttest\tgene\t67\t84\t.\t-\t.\tgene_id \"MINUS_B\";\n"
    )
    raw_counts = tmp_path / "raw.tsv"
    keys = StringIO()
    map_alignment_gene_keys(
        [record + "\n" for record in records],
        output=keys,
        raw_counts_path=raw_counts,
        aligner="bowtie1-v3",
        gene_index=GeneIndex.from_gtf(gtf),
        query_lengths={"known_query": len(query)},
    )
    counts = reduce_sorted_gene_keys(
        sorted(set(keys.getvalue().splitlines())),
        raw_counts_path=raw_counts,
        query_lengths={"known_query": len(query)},
    )
    assert counts["known_query"] == {
        "raw_alignment_count": 3,
        "distinct_genes_any_strand": 3,
        "distinct_genes_transcript_compatible": 3,
    }

