"""CPU-only tests for the J4 scale-out config builder. Hand-checked tiny loci.

Toy coordinates sit >5050 nt from contig start so the Pangolin context window fits
(except the deliberate chrM / unknown-scaffold cases, which must be skipped).
"""
import hashlib
import json
from pathlib import Path

import pytest

from conference.model_v2 import gpu_perturbation_atlas as atlas
from conference.model_v2 import j4_build_configs as jb
from feng.titration import gpu_scores as gs

# locus_id, chrom, start1, end1, strand  (padded +-100, start clipped at 1)
TOY_ROWS = [
    ("TX1:ex1", "chr1", 1000000, 1000500, "+"),   # exon 1000100-1000400, donor only
    ("TX1:ex2", "chr1", 1002000, 1002400, "+"),   # exon 1002100-1002300, acceptor+donor
    ("TX1:ex3", "chr1", 1003000, 1003300, "+"),   # exon 1003100-1003200, acceptor only
    ("TXM:ex1", "chr1", 1009000, 1009500, "-"),   # exon 1009100-1009400, donor = 1009100
    ("TXM:ex2", "chr1", 1008000, 1008400, "-"),   # exon 1008100-1008300, acceptor = 1008300
    ("TXCLIP:ex1", "chr1", 50, 400, "+"),         # left pad clipped -> ambiguous
    ("TXCLIP:ex2", "chr1", 1000600, 1000900, "+"),  # exon 1000700-1000800, acceptor 1000700
    ("TXGAP:ex1", "chr1", 1011000, 1011400, "+"),  # indices {1,3} -> whole transcript skipped
    ("TXGAP:ex3", "chr1", 1012000, 1012400, "+"),
    ("TXSINGLE:ex1", "chr1", 1013000, 1013500, "+"),  # one exon -> no splice site
    ("TXMT:ex1", "chrM", 1000, 1500, "+"),        # donor 1400 -> context underruns chrM
    ("TXMT:ex2", "chrM", 5000, 5400, "+"),        # acceptor 5100 -> fits chrM
    ("TXALT:ex1", "chrUn_toy", 1000000, 1000500, "+"),  # unbounded scaffold -> skipped
    ("TXALT:ex2", "chrUn_toy", 1002000, 1002400, "+"),
]


@pytest.fixture
def loci_file(tmp_path):
    lines = ["\t".join(str(f) for f in row) for row in TOY_ROWS]
    lines.append("#locus_id\tchrom\tstart1\tend1\tstrand")  # trailing header, as in J2 output
    lines.append("")
    path = tmp_path / "toy_loci.tsv"
    path.write_text("\n".join(lines))
    return path


@pytest.fixture
def built(loci_file, tmp_path):
    out = tmp_path / "configs"
    summary = jb.build(loci_file, out, configs_per_shard=3)
    return summary, out


def manifest_rows(out_dir: Path):
    rows = []
    for line in (out_dir / "j4_configs_manifest.tsv").read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        rows.append(dict(zip(jb.MANIFEST_COLUMNS, line.split("\t"))))
    return rows


def shard_configs(out_dir: Path):
    configs = []
    for path in sorted((out_dir / "shards").glob("shard_*.jsonl")):
        for line in path.read_text().splitlines():
            configs.append((path.name, line, json.loads(line)))
    return configs


def by_id(out_dir: Path):
    return {cfg["config_id"]: cfg for _, _, cfg in shard_configs(out_dir)}


# --------------------------------------------------------------- parsing / bounds
def test_read_loci_skips_comments_and_sorts(loci_file):
    loci = jb.read_loci(loci_file)
    assert len(loci) == len(TOY_ROWS)
    assert [(l.transcript_id, l.exon_index) for l in loci] == sorted(
        (l.transcript_id, l.exon_index) for l in loci)


@pytest.mark.parametrize("line,message", [
    ("TX9:ex1\tchr1\t100\t200", "5 tab-separated"),
    ("TX9_ex1\tchr1\t1000\t2000\t+", "not '<transcript>:ex"),
    ("TX9:ex1\tchr1\t2000\t1000\t+", "start1 <= end1"),
    ("TX9:ex1\tchr1\t1000\t2000\t*", "strand must be"),
    ("TX9:ex1\tchr1\tx\t2000\t+", "non-integer"),
])
def test_read_loci_fails_closed(tmp_path, line, message):
    path = tmp_path / "bad.tsv"
    path.write_text(line + "\n")
    with pytest.raises(ValueError, match=message):
        jb.read_loci(path)


def test_read_loci_rejects_duplicates(tmp_path):
    path = tmp_path / "dup.tsv"
    path.write_text("TX9:ex1\tchr1\t1000\t2000\t+\nTX9:ex1\tchr1\t1000\t2000\t+\n")
    with pytest.raises(ValueError, match="duplicate locus_id"):
        jb.read_loci(path)


def test_exon_bounds_trims_pad():
    locus = jb.Locus("TX:ex1", "TX", 1, "chr1", 1000000, 1000500, "+")
    assert jb.exon_bounds(locus) == (1000100, 1000400, None)


def test_exon_bounds_flags_clipped_and_empty():
    clipped = jb.Locus("TX:ex1", "TX", 1, "chr1", 100, 400, "+")
    assert jb.exon_bounds(clipped) == (None, None, "padding_clipped_ambiguous")
    empty = jb.Locus("TX:ex1", "TX", 1, "chr1", 1000000, 1000150, "+")
    assert jb.exon_bounds(empty) == (None, None, "padding_trim_empty")


# ------------------------------------------------------------------- site geometry
def test_plus_strand_sites(built):
    _, out = built
    configs = by_id(out)
    assert configs["TX1_ex1__donor"]["windows"][0]["sites"] == [
        {"position": 1000400, "kind": "donor"}]
    assert configs["TX1_ex2__acceptor"]["windows"][0]["sites"] == [
        {"position": 1002100, "kind": "acceptor"}]
    assert configs["TX1_ex2__donor"]["windows"][0]["sites"] == [
        {"position": 1002300, "kind": "donor"}]
    assert configs["TX1_ex3__acceptor"]["windows"][0]["sites"] == [
        {"position": 1003100, "kind": "acceptor"}]
    assert "TX1_ex1__acceptor" not in configs  # transcript 5' end is not an acceptor
    assert "TX1_ex3__donor" not in configs     # transcript 3' end is not a donor


def test_minus_strand_sites_are_transcript_oriented(built):
    _, out = built
    configs = by_id(out)
    donor = configs["TXM_ex1__donor"]["windows"][0]
    acceptor = configs["TXM_ex2__acceptor"]["windows"][0]
    assert donor["strand"] == "-" and donor["sites"] == [{"position": 1009100, "kind": "donor"}]
    assert acceptor["sites"] == [{"position": 1008300, "kind": "acceptor"}]
    assert donor["center"] == 1009100 and donor["distance"] == jb.DEFAULT_HALF_SPAN


def test_smn2_site_convention_matches_gpu_scores():
    # donor = last exonic base (70076574 +), not intron first base 70076575
    assert jb.site_position(70076474, gs.EXON7_LAST_BASE_1B, "+", "donor") == gs.EXON7_LAST_BASE_1B
    assert jb.site_position(gs.EXON7_LAST_BASE_1B, 70076674, "-", "donor") == gs.EXON7_LAST_BASE_1B


def test_skipped_rows_recorded_with_reasons(built):
    _, out = built
    rows = {(r["locus_id"], r["reason"]) for r in manifest_rows(out) if r["status"] == "skipped"}
    assert ("TXCLIP:ex1", "padding_clipped_ambiguous") in rows
    assert ("TXGAP:ex1", "exon_index_noncontiguous") in rows
    assert ("TXGAP:ex3", "exon_index_noncontiguous") in rows
    assert ("TXSINGLE:ex1", "single_exon_transcript_no_splice_site") in rows
    assert ("TXMT:ex1", "context_out_of_chrom") in rows
    assert ("TXALT:ex1", "context_out_of_chrom") in rows
    assert ("TXALT:ex2", "context_out_of_chrom") in rows
    assert not any(c["config_id"].startswith("TXGAP") for _, _, c in shard_configs(out))


def test_clipped_transcript_other_exon_still_emitted(built):
    _, out = built
    assert by_id(out)["TXCLIP_ex2__acceptor"]["windows"][0]["sites"][0]["position"] == 1000700


def test_chrom_bounded_site_emitted(built):
    _, out = built
    assert by_id(out)["TXMT_ex2__acceptor"]["windows"][0]["center"] == 5100


def test_expected_config_and_skip_counts(built):
    summary, out = built
    # TX1: 4, TXM: 2, TXCLIP: 1, TXMT: 1, plus the curated SMN2 config
    assert summary["configs"] == 9
    assert summary["skipped"] == 7
    assert len(manifest_rows(out)) == summary["configs"] + summary["skipped"]


# ---------------------------------------------------------------- grid arithmetic
def test_default_grid_interval_and_counts(built):
    _, out = built
    grid = by_id(out)["TX1_ex2__donor"]["grid"]
    assert (grid["start"], grid["end"]) == (1002300 - 50, 1002300 + 50)
    assert grid["lengths"] == [10, 14, 18] and grid["stride"] == 5
    counts = [len(jb.footprint_starts(grid["start"], grid["end"], n, grid["stride"]))
              for n in grid["lengths"]]
    assert counts == [19, 18, 17]                  # hand-checked: 101 nt interval, stride 5
    assert jb.footprint_count(grid["start"], grid["end"], grid["lengths"], grid["stride"]) == 54


def test_every_footprint_inside_interval(built):
    _, out = built
    for _, _, config in shard_configs(out):
        grid = config["grid"]
        for length in grid["lengths"]:
            starts = jb.footprint_starts(grid["start"], grid["end"], length, grid["stride"])
            assert starts, "empty footprint set"
            assert starts[0] == grid["start"]
            assert starts[-1] + length - 1 <= grid["end"]
            assert starts[-1] + grid["stride"] + length - 1 > grid["end"]


def test_tiny_grid_hand_checked(loci_file, tmp_path):
    out = tmp_path / "tiny"
    jb.build(loci_file, out, half_span=10, lengths=(5, 7), stride=5, curated=False)
    grid = by_id(out)["TX1_ex2__donor"]["grid"]
    site = 1002300
    assert (grid["start"], grid["end"]) == (site - 10, site + 10)
    assert list(jb.footprint_starts(grid["start"], grid["end"], 5, 5)) == [
        site - 10, site - 5, site, site + 5]
    assert list(jb.footprint_starts(grid["start"], grid["end"], 7, 5)) == [
        site - 10, site - 5, site]
    assert jb.footprint_count(grid["start"], grid["end"], (5, 7), 5) == 7


def test_footprint_enumeration_matches_engine_tasks(built):
    _, out = built
    for _, _, config in shard_configs(out):
        window = atlas.Window(**dict(config["windows"][0],
                                     sites=tuple(atlas.Site(**s) for s in config["windows"][0]["sites"])))
        grid = atlas.Grid(**dict(config["grid"], lengths=tuple(config["grid"]["lengths"])))
        engine = sorted((start, length) for _, start, length in atlas.tasks([window], grid))
        mine = sorted((start, length) for length in grid.lengths
                      for start in jb.footprint_starts(grid.start, grid.end, length, grid.stride))
        assert engine == mine


def test_configs_validate_against_engine_dataclasses(built):
    _, out = built
    for _, _, config in shard_configs(out):
        assert config["assembly"] == "GRCh38" and config["annotation"] == "GENCODE v50"
        window = atlas.Window(**dict(config["windows"][0],
                                     sites=tuple(atlas.Site(**s) for s in config["windows"][0]["sites"])))
        grid = atlas.Grid(**dict(config["grid"], lengths=tuple(config["grid"]["lengths"])))
        low, high = gs.window_bounds_1b(window.center, window.distance, gs.PANGOLIN_CONTEXT)
        assert low >= 1 and low <= grid.start <= grid.end <= high
        assert low <= atlas.GRCH38_PRIMARY_LENGTHS[window.chrom]
        assert high <= atlas.GRCH38_PRIMARY_LENGTHS[window.chrom]


# ------------------------------------------------------------------ curated config
def test_curated_smn2_config_is_last_and_fine_grained(built):
    _, out = built
    name, _, config = shard_configs(out)[-1]
    assert config["config_id"] == jb.CURATED_ID
    window, grid = config["windows"][0], config["grid"]
    assert window["chrom"] == "chr5" and window["strand"] == "+"
    assert window["center"] == gs.EXON7_LAST_BASE_1B == 70076574
    assert window["sites"] == [{"position": 70076574, "kind": "donor"}]
    assert (grid["start"], grid["end"]) == (70076584 - 30, 70076601 + 30) == (70076554, 70076631)
    assert grid["stride"] == 1 and grid["lengths"] == [8, 10, 14, 18]
    counts = [len(jb.footprint_starts(grid["start"], grid["end"], n, 1)) for n in grid["lengths"]]
    assert counts == [71, 69, 65, 61]              # 78 nt interval, stride 1 -> 79 - length
    assert jb.footprint_count(grid["start"], grid["end"], grid["lengths"], 1) == 266


def test_curated_can_be_disabled(loci_file, tmp_path):
    summary = jb.build(loci_file, tmp_path / "nocur", configs_per_shard=3, curated=False)
    assert summary["configs"] == 8
    assert jb.CURATED_ID not in by_id(tmp_path / "nocur")


# ------------------------------------------------------ accounting / sharding / sha
def test_row_and_inference_accounting(built):
    summary, out = built
    rows = [r for r in manifest_rows(out) if r["status"] == "emitted"]
    for row in rows:
        footprints = int(row["footprints"])
        assert int(row["rows"]) == footprints * len(atlas.TISSUES)
        assert int(row["inferences"]) == footprints + 1
    assert summary["total_footprints"] == sum(int(r["footprints"]) for r in rows)
    assert summary["total_rows"] == sum(int(r["rows"]) for r in rows)
    assert summary["total_inferences"] == sum(int(r["inferences"]) for r in rows)
    assert summary["total_footprints"] == 8 * 54 + 266


def test_shard_partition_is_deterministic_and_dense(built):
    summary, out = built
    shards = sorted((out / "shards").glob("shard_*.jsonl"))
    assert [p.name for p in shards] == [f"shard_{i:04d}.jsonl" for i in range(summary["shards"])]
    sizes = [len(p.read_text().splitlines()) for p in shards]
    assert sizes == [3, 3, 3]
    emitted = [r for r in manifest_rows(out) if r["status"] == "emitted"]
    assert [int(r["shard"]) for r in emitted] == [0, 0, 0, 1, 1, 1, 2, 2, 2]


def test_manifest_sha256_matches_shard_lines(built):
    _, out = built
    digests = {json.loads(line)["config_id"]: hashlib.sha256(line.encode()).hexdigest()
               for _, line, _ in shard_configs(out)}
    for row in manifest_rows(out):
        if row["status"] == "emitted":
            assert row["config_sha256"] == digests[row["config_id"]]


def test_shard_table_and_provenance_hashes(built, loci_file):
    summary, out = built
    table = [line.split("\t") for line in (out / "j4_config_shards.tsv").read_text().splitlines()
             if not line.startswith("#")]
    assert len(table) == summary["shards"]
    for row in table:
        record = dict(zip(jb.SHARD_COLUMNS, row))
        assert record["file_sha256"] == jb.sha256_file(out / record["shard_file"])
        assert summary["shard_sha256"][record["shard_file"]] == record["file_sha256"]
    provenance = json.loads((out / "j4_build_provenance.json").read_text())
    assert provenance["loci_sha256"] == jb.sha256_file(loci_file)
    assert provenance["configs"] == summary["configs"]
    assert provenance["parameters"]["exon_pad"] == 100
    assert provenance["parameters"]["pangolin_context"] == gs.PANGOLIN_CONTEXT


def test_build_is_byte_identical_across_runs(loci_file, tmp_path):
    first, second = tmp_path / "a", tmp_path / "b"
    jb.build(loci_file, first, configs_per_shard=3)
    jb.build(loci_file, second, configs_per_shard=3)
    names = sorted(p.relative_to(first).as_posix() for p in first.rglob("*") if p.is_file())
    assert names == sorted(p.relative_to(second).as_posix() for p in second.rglob("*") if p.is_file())
    for name in names:
        assert (first / name).read_bytes() == (second / name).read_bytes(), name


def test_rebuild_over_existing_output_is_stable(loci_file, tmp_path):
    out = tmp_path / "same"
    first = jb.build(loci_file, out, configs_per_shard=3)
    snapshot = {p.name: p.read_bytes() for p in (out / "shards").iterdir()}
    second = jb.build(loci_file, out, configs_per_shard=3)
    assert first == second
    assert {p.name: p.read_bytes() for p in (out / "shards").iterdir()} == snapshot


def test_build_rejects_lengths_longer_than_interval(loci_file, tmp_path):
    with pytest.raises(ValueError, match="tiling interval"):
        jb.build(loci_file, tmp_path / "bad", half_span=5, lengths=(18,))


def test_cli_prints_summary(loci_file, tmp_path, capsys):
    assert jb.main(["--loci", str(loci_file), "--out-dir", str(tmp_path / "cli"),
                    "--configs-per-shard", "4"]) == 0
    printed = json.loads(capsys.readouterr().out)
    assert printed["configs"] == 9 and printed["total_footprints"] == 8 * 54 + 266
