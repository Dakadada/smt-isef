"""J6: independent Hamming oracle, annotation, coordinates, and crash recovery."""
import dataclasses
import hashlib
import json
import random
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "feng/src"))
from conference.model_v2 import trans_intervention as ti


def module():
    from conference.model_v2 import j6_site_universe
    return j6_site_universe


def test_module_exposes_site_enumerator():
    assert callable(module().enumerate_sites)


def rc(seq):
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def oracle(seqs, genes, asos, k):
    rows = set()
    for gene in genes:
        for aso in asos:
            query = rc(aso.sequence_5to3)
            for start in range(gene.start1, gene.end1 - len(query) + 2):
                target = seqs[gene.chrom][start - 1:start + len(query) - 1]
                if gene.strand == "-":
                    target = rc(target)
                positions = tuple(sorted(len(query) - i for i, (a, b) in
                                         enumerate(zip(query, target)) if a != b or b == "N"))
                if len(positions) <= k:
                    rows.add((aso.aso_id, gene.chrom, start, start + len(query) - 1,
                              gene.strand, gene.gene_id, target, positions))
    return rows


def key(site):
    return (site.aso_id, site.chrom, site.start1, site.end1, site.strand,
            site.gene_id, site.target_sequence, site.mismatch_positions)


@pytest.mark.parametrize("k", range(4))
@pytest.mark.parametrize("seed", [0, 17, 98])
def test_seed_extend_equals_independent_oracle_with_edges_N_overlaps(seed, k):
    j = module()
    rng = random.Random(seed)
    seq = list("".join(rng.choice("ACGTN") for _ in range(190)))
    target = "AACGTT"
    for pos, planted in [(0, target), (15, "AATGTT"), (40, rc(target)),
                         (65, "AACNTT"), (100, "AAAAAAA"), (184, target)]:
        seq[pos:pos + len(planted)] = planted
    seqs = {"chrT": "".join(seq), "chrTiny": "AC"}
    genes = [j.Gene("g+", "chrT", 1, 190, "+"),
             j.Gene("g-", "chrT", 1, 190, "-"),
             j.Gene("partial", "chrT", 3, 187, "+"),
             j.Gene("tiny", "chrTiny", 1, 2, "+")]
    asos = [ti.ASO("q", rc(target)), ti.ASO("overlap", "TTTT"), ti.ASO("short", "A")]
    genome = ti.InMemoryGenome(seqs)
    expected = oracle(seqs, genes, asos, k)
    unsharded = j.enumerate_sites(genome, genes, asos, max_mismatch=k)
    assert {key(s) for s in unsharded} == expected
    sharded = [s for shard in j.make_shards(genes, window_size=13)
               for s in j.enumerate_sites(genome, genes, asos, max_mismatch=k, shard=shard)]
    assert {key(s) for s in sharded} == expected
    assert len(sharded) == len(expected)
    assert {s.site_id for s in sharded} == {s.site_id for s in unsharded}
    assert all(s.mismatch_count == len(s.mismatch_positions) and s.evidence_class == "derived"
               for s in sharded)


@pytest.mark.parametrize("strand,sequence,want", [("+", "TACGTA", (4,)), ("-", "TACGTA", (4,))])
def test_mismatch_positions_are_one_based_ASO_five_prime(strand, sequence, want):
    j = module()
    sense = sequence
    genome = ti.InMemoryGenome({"chrT": sense if strand == "+" else rc(sense)})
    sites = j.enumerate_sites(genome, [j.Gene("g", "chrT", 1, 6, strand)],
                              [ti.ASO("a", rc("TATGTA"))], max_mismatch=1)
    assert len(sites) == 1
    assert sites[0].mismatch_positions == want
    assert sites[0].target_sequence == sense


def test_annotation_collapses_transcripts_but_not_overlapping_genes():
    j = module()
    transcripts = (j.Transcript("t2", 1, 16, ((1, 8), (13, 16)), ((1, 4),)),
                   j.Transcript("t1", 1, 16, ((1, 8), (13, 16)), ((1, 4),)))
    genes = [j.Gene("g", "chrT", 1, 16, "+", transcripts=transcripts),
             j.Gene("other", "chrT", 1, 16, "+")]
    sites = j.enumerate_sites(ti.InMemoryGenome({"chrT": "A" * 16}), genes,
                              [ti.ASO("a", "TTTT")], max_mismatch=0)
    assert len(sites) == 26
    by_start = {s.start1: s for s in sites if s.gene_id == "g"}
    assert by_start[1].region_class == "UTR"
    assert by_start[5].region_class == "exon"
    assert by_start[7].region_class == "junction-spanning"
    assert by_start[9].region_class == "intron"
    assert by_start[1].transcript_ids == ("t1", "t2")
    assert all(s.region_class == "not_available" and s.transcript_ids is None
               for s in sites if s.gene_id == "other")


def test_alternative_transcripts_preserve_multiple_region_classes():
    j = module()
    txs = (j.Transcript("exonic", 1, 8, ((1, 8),)),
           j.Transcript("intronic", 1, 8, ((1, 2), (7, 8))))
    gene = j.Gene("g", "chrT", 1, 8, "+", transcripts=txs)
    sites = j.enumerate_sites(ti.InMemoryGenome({"chrT": "AACCGGTT"}), [gene],
                              [ti.ASO("a", "CCGG")], max_mismatch=0)
    assert len(sites) == 1
    assert sites[0].region_class == "transcript_dependent"
    assert sites[0].region_classes == ("exon", "intron")


def test_gene_gtf_reader_is_gene_only_region_bounded_and_assembly_fenced(tmp_path):
    j = module()
    path = tmp_path / "genes.gtf"
    path.write_text('##description: genome (GRCh38), version 50\n'
                    'chrT\tt\tgene\t1\t20\t.\t+\t.\tgene_id "g1"; gene_name "SMN2";\n'
                    'chrT\tt\texon\tBROKEN\tIGNORED\t.\t+\t.\t.\n'
                    'chrT\tt\tgene\t30\t40\t.\t-\t.\tgene_id "g2";\n')
    genes = j.read_gene_bodies(path, region=("chrT", 10, 25))
    assert [(g.gene_id, g.start1, g.end1, g.gene_name) for g in genes] == [("g1", 1, 20, "SMN2")]
    for header in ["GRCh37, version 50", "GRCh38, version 49", "unknown"]:
        path.write_text("##description: " + header + "\n")
        with pytest.raises(ValueError):
            j.read_gene_bodies(path)


def test_region_limits_ownership_without_cutting_boundary_hit():
    j = module()
    gene = j.Gene("g", "chrT", 1, 12, "+")
    genome = ti.InMemoryGenome({"chrT": "AAAAAAAAAAAA"})
    shards = j.make_shards([gene], window_size=4, region=("chrT", 4, 6))
    sites = [s for shard in shards for s in j.enumerate_sites(
        genome, [gene], [ti.ASO("a", "TTTT")], max_mismatch=0, shard=shard)]
    assert [s.start1 for s in sites] == [4, 5, 6]
    assert len({s.site_id for s in sites}) == 3
    assert j.make_shards([gene], window_size=4) == j.make_shards([gene, gene], window_size=4)


def tiny_job(tmp_path):
    j = module()
    genome = ti.InMemoryGenome({"chrT": "AACCGGTT"})
    genes = [j.Gene("g", "chrT", 1, 8, "+")]
    asos = [ti.ASO("a", "CCGG")]
    shard = j.make_shards(genes, window_size=8)[0]
    return j, dict(genome=genome, genes=genes, asos=asos, shard=shard,
                   output_dir=tmp_path, max_mismatch=0,
                   identities={"fasta": {"test": "synthetic"}, "gtf": {"test": "synthetic"}})


def test_atomic_parquet_manifest_and_restart_idempotence(tmp_path):
    import pyarrow.parquet as pq
    j, kwargs = tiny_job(tmp_path)
    path = j.run_shard(**kwargs)
    done = path.with_suffix(".done")
    manifest = json.loads(done.read_text())
    assert manifest["rows"] == 1
    assert manifest["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert manifest["params"]["max_mismatch"] == 0
    assert manifest["identities"] == kwargs["identities"]
    table = pq.read_table(path)
    assert table.to_pylist()[0]["mismatch_positions"] == []
    assert table.to_pylist()[0]["transcript_ids"] is None
    mtimes = path.stat().st_mtime_ns, done.stat().st_mtime_ns
    j.run_shard(**kwargs)
    assert (path.stat().st_mtime_ns, done.stat().st_mtime_ns) == mtimes
    assert not list(tmp_path.glob("*.tmp"))


@pytest.mark.parametrize("damage", ["parquet", "manifest", "missing_done", "params", "identity", "genes", "rows"])
def test_restart_rebuilds_incomplete_corrupt_or_stale_shards(tmp_path, damage):
    j, kwargs = tiny_job(tmp_path)
    path = j.run_shard(**kwargs)
    done = path.with_suffix(".done")
    if damage == "parquet":
        path.write_bytes(b"corrupt")
    elif damage == "manifest":
        done.write_text("{")
    elif damage == "missing_done":
        done.unlink()
    elif damage == "params":
        kwargs["max_mismatch"] = 1
    elif damage == "identity":
        kwargs["identities"]["fasta"]["test"] = "changed"
    elif damage == "genes":
        kwargs["genes"] = [dataclasses.replace(kwargs["genes"][0], gene_id="new")]
    else:
        manifest = json.loads(done.read_text())
        manifest["rows"] = 100
        done.write_text(json.dumps(manifest))
    j.run_shard(**kwargs)
    manifest = json.loads(done.read_text())
    assert manifest["rows"] == 1
    assert manifest["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert manifest["params"]["max_mismatch"] == kwargs["max_mismatch"]
    assert manifest["identities"] == kwargs["identities"]
    if damage == "genes":
        import pyarrow.parquet as pq
        assert pq.read_table(path).to_pylist()[0]["gene_id"] == "new"


def test_empty_shard_has_typed_schema(tmp_path):
    import pyarrow.parquet as pq
    j, kwargs = tiny_job(tmp_path)
    kwargs["asos"] = [ti.ASO("nohit", "AAAA")]
    path = j.run_shard(**kwargs)
    table = pq.read_table(path)
    assert table.num_rows == 0
    assert str(table.schema.field("start1").type) == "int64"
    assert "mismatch_positions" in table.column_names


@pytest.mark.parametrize("bad", [-1, 1.5, True])
def test_invalid_budgets_fail_closed(bad):
    j, kwargs = tiny_job(Path("unused"))
    with pytest.raises((TypeError, ValueError)):
        j.enumerate_sites(kwargs["genome"], kwargs["genes"], kwargs["asos"], max_mismatch=bad)


def test_duplicate_ASO_ids_and_conflicting_gene_ids_fail_closed():
    j, kw = tiny_job(Path("unused"))
    with pytest.raises(ValueError):
        j.enumerate_sites(kw["genome"], kw["genes"], kw["asos"] * 2)
    with pytest.raises(ValueError):
        j.enumerate_sites(kw["genome"], kw["genes"] + [
            dataclasses.replace(kw["genes"][0], strand="-")], kw["asos"])


def test_real_GRCh38_SMN1_SMN2_exact_hits_and_paralog_ambiguity():
    j = module()
    try:
        genome = ti.FaiFastaReader()
    except ti.GenomeAccessError as exc:
        pytest.skip(str(exc))
    genes = [j.Gene("SMN2", "chr5", 70076521, 70076674, "+", "SMN2"),
             j.Gene("SMN1", "chr5", 70951941, 70952094, "+", "SMN1")]
    asos = [ti.ASOS[a] for a in ("F18", "F14", "F10")]
    sites = j.enumerate_sites(genome, genes, asos, max_mismatch=0)
    for aso_id, length in [("F18", 18), ("F14", 14), ("F10", 10)]:
        selected = [s for s in sites if s.aso_id == aso_id and s.start1 in (70076584, 70952004)]
        assert {(s.gene_id, s.start1, s.end1) for s in selected} == {
            ("SMN2", 70076584, 70076584 + length - 1),
            ("SMN1", 70952004, 70952004 + length - 1)}
        assert len({s.site_id for s in selected}) == 2
        assert all(s.paralog_ambiguity for s in selected)
    sharded = [s for shard in j.make_shards(genes, window_size=100)
               for s in j.enumerate_sites(genome, genes, asos, max_mismatch=0, shard=shard)]
    assert sharded == sorted(sites, key=j.site_sort_key)
