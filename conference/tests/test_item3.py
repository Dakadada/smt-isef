import csv
from pathlib import Path
import subprocess
import sys

from conference.aso_pipeline.item3 import classify_query, finalize, prepare_manifest


def _write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_8mer_is_saturated_and_receives_no_numeric_count():
    result = classify_query("3UP8", "GCUGGCAG", mandatory_family=False)
    assert result.status == "saturated"
    assert result.run_alignment is False
    assert result.gene_count is None
    assert "uninformative" in result.reason


def test_f10_family_override_runs_despite_borderline_saturation_estimate():
    result = classify_query("F10", "AUAAUGCUGG", mandatory_family=True)
    assert result.expected_hits > 10_000_000
    assert result.status == "scheduled_mandatory"
    assert result.run_alignment is True
    assert result.gene_count is None


def test_f14_runs_below_predeclared_saturation_threshold():
    result = classify_query("F14", "UUUCAUAAUGCUGG", mandatory_family=False)
    assert result.expected_hits < 10_000_000
    assert result.status == "scheduled"
    assert result.run_alignment is True


def test_manifest_writes_family_fastas_and_one_f10_fasta_per_unique_query(
    tmp_path: Path,
):
    query_table = tmp_path / "queries.tsv"
    rows = [
        ("F18__real", "F18", "", "True", "UCACUUUCAUAAUGCUGG"),
        ("F14__real", "F14", "", "True", "UUUCAUAAUGCUGG"),
        ("F10__real", "F10", "", "True", "AUAAUGCUGG"),
        ("F10__shuffle_001", "F10", "1", "False", "AUA AUGG CUG".replace(" ", "")),
        ("3UP8__real", "3UP8", "", "True", "GCUGGCAG"),
    ]
    with query_table.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["state", "query_id", "aso_id", "draw_id", "is_real", "sequence_rna_5to3"]
        )
        for query_id, aso_id, draw_id, is_real, sequence in rows:
            writer.writerow(["[candidate]", query_id, aso_id, draw_id, is_real, sequence])

    prepare_manifest(query_table, tmp_path)

    assert (tmp_path / "item3_unique_queries_F18.fa").is_file()
    assert (tmp_path / "item3_unique_queries_F14.fa").is_file()
    assert (tmp_path / "item3_unique_queries_F10.fa").is_file()
    f10_singletons = sorted((tmp_path / "item3_unique_queries_F10").glob("*.fa"))
    assert len(f10_singletons) == 2
    assert not (tmp_path / "item3_unique_queries_3UP8.fa").exists()


def test_live_mapper_sort_reducer_pipe_is_bounded_and_race_free(tmp_path: Path):
    gene_table = tmp_path / "genes.tsv"
    gene_table.write_text(
        "chrom\tstart_0base\tend_0base\tstrand\tgene_id\n"
        "chr1\t0\t100\t+\tPLUS\n"
    )
    query_fasta = tmp_path / "queries.fa"
    query_fasta.write_text(">q1\nAAAAAAAAAA\n")
    raw_counts = tmp_path / "raw.tsv"
    output = tmp_path / "counts.tsv"

    mapper = subprocess.Popen(
        [
            sys.executable,
            "-m",
            "conference.aso_pipeline.item3",
            "map",
            "--aligner",
            "bowtie1-v3",
            "--gene-table",
            str(gene_table),
            "--query-fasta",
            str(query_fasta),
            "--raw-counts",
            str(raw_counts),
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
    )
    sorter = subprocess.Popen(
        ["sort", "-u"], stdin=mapper.stdout, stdout=subprocess.PIPE, text=True
    )
    assert mapper.stdout is not None
    mapper.stdout.close()
    reducer = subprocess.Popen(
        [
            sys.executable,
            "-m",
            "conference.aso_pipeline.item3",
            "reduce",
            "--aligner",
            "bowtie1-v3",
            "--query-fasta",
            str(query_fasta),
            "--raw-counts",
            str(raw_counts),
            "--output",
            str(output),
        ],
        stdin=sorter.stdout,
        text=True,
    )
    assert sorter.stdout is not None
    sorter.stdout.close()
    assert mapper.stdin is not None
    mapper.stdin.write("q1\t-\tchr1\t10\n" * 10_000)
    mapper.stdin.close()

    assert mapper.wait(timeout=20) == 0
    assert sorter.wait(timeout=20) == 0
    assert reducer.wait(timeout=20) == 0
    with output.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["raw_alignment_count"] == "10000"
    assert rows[0]["distinct_genes_transcript_compatible"] == "1"


def test_finalize_emits_persisted_identity_real_headline_and_scramble_tables(
    tmp_path: Path,
):
    manifest_rows = []
    for aso_id in ("F18", "F14", "F10"):
        for is_real, draw_id, suffix in ((True, "", "real"), (False, "1", "null")):
            manifest_rows.append(
                {
                    "state": "[candidate]",
                    "query_id": f"{aso_id}__{suffix}",
                    "aso_id": aso_id,
                    "draw_id": draw_id,
                    "is_real": str(is_real),
                    "sequence_rna_5to3": "A" * int(aso_id[1:]),
                    "unique_query_id": f"Q_{aso_id}_{suffix}",
                    "length_nt": int(aso_id[1:]),
                    "expected_random_both_strand_hits_v3": 1.0,
                    "preflight_status": "scheduled",
                    "preflight_reason": "fixture",
                    "run_alignment": True,
                }
            )
    manifest_rows.extend(
        [
            {
                "state": "[candidate]",
                "query_id": f"3UP8__{suffix}",
                "aso_id": "3UP8",
                "draw_id": draw_id,
                "is_real": str(is_real),
                "sequence_rna_5to3": "GCTGGCAG",
                "unique_query_id": f"Q_3UP8_{suffix}",
                "length_nt": 8,
                "expected_random_both_strand_hits_v3": 100_000_000,
                "preflight_status": "saturated",
                "preflight_reason": "uninformative",
                "run_alignment": False,
            }
            for is_real, draw_id, suffix in ((True, "", "real"), (False, "1", "null"))
        ]
    )
    manifest = tmp_path / "manifest.tsv"
    _write_tsv(manifest, manifest_rows)

    bowtie_values = {"F18": (11, 5), "F14": (12, 6), "F10": (13, 9)}
    blast_values = {"F18": (21, 7), "F14": (22, 8), "F10": (23, 10)}
    count_paths = []
    for aligner, values in (("bowtie1-v3", bowtie_values), ("blastn-short", blast_values)):
        rows = []
        for aso_id, (real_count, null_count) in values.items():
            for suffix, count in (("real", real_count), ("null", null_count)):
                rows.append(
                    {
                        "state": "[candidate]",
                        "unique_query_id": f"Q_{aso_id}_{suffix}",
                        "aligner": aligner,
                        "raw_alignment_count": count * 10,
                        "distinct_genes_any_strand": count + 1,
                        "distinct_genes_transcript_compatible": count,
                    }
                )
        path = tmp_path / f"{aligner}.tsv"
        _write_tsv(path, rows)
        count_paths.append(path)

    finalize(
        manifest=manifest,
        bowtie_counts=count_paths[0],
        blast_counts=count_paths[1],
        output_dir=tmp_path,
    )

    identity = list(
        csv.DictReader((tmp_path / "item3_query_identity_map.tsv").open(), delimiter="\t")
    )
    assert len(identity) == 8
    assert {row["is_real"] for row in identity} == {"true", "false"}
    assert next(row for row in identity if row["query_id"] == "F18__null")[
        "scramble_index"
    ] == "1"

    real = list(
        csv.DictReader(
            (tmp_path / "item3_real_aso_gene_counts.tsv").open(), delimiter="\t"
        )
    )
    assert [row["aso_id"] for row in real] == ["F18", "F14", "F10"]
    assert all(row["is_real"] == "true" for row in real)
    assert real[0]["bowtie1_v3_distinct_genes"] == "11"
    assert real[0]["blastn_short_distinct_genes"] == "21"

    headline = list(
        csv.DictReader(
            (tmp_path / "item3_headline_f18_vs_f10.tsv").open(), delimiter="\t"
        )
    )
    assert headline == [
        {
            "state": "[candidate]",
            "metric": "distinct_genes_transcript_compatible",
            "f18_bowtie1_v3": "11",
            "f10_bowtie1_v3": "13",
            "f18_blastn_short": "21",
            "f10_blastn_short": "23",
        }
    ]

    null_summary = list(
        csv.DictReader(
            (tmp_path / "item6_alignment_scramble_summary.tsv").open(), delimiter="\t"
        )
    )
    assert len(null_summary) == 6
    f18_bowtie = next(
        row
        for row in null_summary
        if row["aso_id"] == "F18" and row["aligner"] == "bowtie1-v3"
    )
    assert f18_bowtie["scramble_n"] == "1"
    assert f18_bowtie["scramble_mean"] == "5.0"
    assert f18_bowtie["scramble_q05"] == "5.0"
    assert f18_bowtie["scramble_median"] == "5.0"
    assert f18_bowtie["scramble_q95"] == "5.0"
