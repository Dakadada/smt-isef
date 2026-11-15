"""J6: independent Hamming oracle, annotation, coordinates, and crash recovery."""
import dataclasses
import hashlib
import json
import random
import subprocess
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


def enumerate_synthetic(j, *args, **kwargs):
    return j.enumerate_sites(*args, synthetic_reference=True, **kwargs)


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
    unsharded = enumerate_synthetic(j, genome, genes, asos, max_mismatch=k)
    assert {key(s) for s in unsharded} == expected
    sharded = [s for shard in j.make_shards(genes, window_size=13)
               for s in enumerate_synthetic(j, genome, genes, asos, max_mismatch=k, shard=shard)]
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
    sites = enumerate_synthetic(j, genome, [j.Gene("g", "chrT", 1, 6, strand)],
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
    sites = enumerate_synthetic(j, ti.InMemoryGenome({"chrT": "A" * 16}), genes,
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
    sites = enumerate_synthetic(j, ti.InMemoryGenome({"chrT": "AACCGGTT"}), [gene],
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
    sites = [s for shard in shards for s in enumerate_synthetic(j,
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
                   output_dir=tmp_path, max_mismatch=0, synthetic_reference=True,
                   identities={"fasta": {"test": "synthetic"}, "gtf": {"test": "synthetic"}})


def test_atomic_parquet_manifest_and_restart_idempotence(tmp_path):
    import pyarrow.parquet as pq
    j, kwargs = tiny_job(tmp_path)
    path = j.run_shard(**kwargs)
    done = path.with_suffix(".done")
    manifest = json.loads(done.read_text())
    assert manifest["rows"] == 1 and type(manifest["rows"]) is int
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


@pytest.mark.parametrize("damage", ["parquet", "manifest", "missing_done", "params", "identity", "genes", "rows", "list", "null", "bad_rows", "bad_sha"])
def test_restart_rebuilds_incomplete_corrupt_or_stale_shards(tmp_path, damage):
    j, kwargs = tiny_job(tmp_path)
    path = j.run_shard(**kwargs)
    done = path.with_suffix(".done")
    if damage == "parquet":
        path.write_bytes(b"corrupt")
    elif damage == "manifest":
        done.write_text("{")
    elif damage in ("list", "null"):
        done.write_text("[]" if damage == "list" else "null")
    elif damage in ("bad_rows", "bad_sha"):
        manifest = json.loads(done.read_text())
        manifest["rows" if damage == "bad_rows" else "sha256"] = True if damage == "bad_rows" else []
        done.write_text(json.dumps(manifest))
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
    if damage in ("params", "identity", "genes"):
        with pytest.raises(ValueError, match="incompatible run"):
            j.run_shard(**kwargs)
        assert done.is_file()  # The prior generation is preserved intact.
        return
    j.run_shard(**kwargs)
    manifest = json.loads(done.read_text())
    assert manifest["rows"] == 1 and type(manifest["rows"]) is int
    assert manifest["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert manifest["params"]["max_mismatch"] == kwargs["max_mismatch"]
    assert manifest["identities"] == kwargs["identities"]


def test_empty_shard_has_typed_schema(tmp_path):
    import pyarrow.parquet as pq
    j, kwargs = tiny_job(tmp_path)
    kwargs["asos"] = [ti.ASO("nohit", "AAAA")]
    path = j.run_shard(**kwargs)
    table = pq.read_table(path)
    assert table.num_rows == 0
    assert str(table.schema.field("start1").type) == "int64"
    assert "mismatch_positions" in table.column_names


def test_crash_between_parquet_and_marker_leaves_no_valid_commit(tmp_path, monkeypatch):
    j, kwargs = tiny_job(tmp_path)
    replace = j.os.replace

    def crash_at_done(source, dest):
        if str(dest).endswith(".done"):
            raise OSError("simulated interruption")
        return replace(source, dest)

    monkeypatch.setattr(j.os, "replace", crash_at_done)
    with pytest.raises(OSError, match="simulated interruption"):
        j.run_shard(**kwargs)
    assert len(list(tmp_path.glob("*.parquet"))) == 1
    assert not list(tmp_path.glob("*.done"))
    assert not list(tmp_path.glob("*.tmp"))
    monkeypatch.setattr(j.os, "replace", replace)
    path = j.run_shard(**kwargs)
    assert path.with_suffix(".done").is_file()


def test_stale_real_input_identity_cannot_resume(tmp_path):
    j, kwargs = tiny_job(tmp_path / "output")
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chrT\nAACCGGTT\n")
    kwargs["identities"]["fasta"] = j.file_identity(fasta)
    j.run_shard(**kwargs)
    fasta.write_text(">chrT\nAACCGGTTAAAA\n")
    with pytest.raises(ValueError, match="identity changed"):
        j.run_shard(**kwargs)


def test_reference_mutation_during_enumeration_cannot_commit(tmp_path):
    j, kwargs = tiny_job(tmp_path / "output")
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chrT\nAACCGGTT\n")
    kwargs["identities"]["fasta"] = j.file_identity(fasta)

    class ChangingGenome(ti.InMemoryGenome):
        def fetch(self, chrom, start, end):
            fasta.write_text(">chrT\nAACCGGTTAAAA\n")
            return super().fetch(chrom, start, end)

    kwargs["genome"] = ChangingGenome({"chrT": "AACCGGTT"})
    with pytest.raises(ValueError, match="identity changed"):
        j.run_shard(**kwargs)
    assert not list((tmp_path / "output").glob("*.done"))


def test_cli_with_synthetic_fai_gtf_and_idempotent_restart(tmp_path):
    import os
    import pyarrow.parquet as pq
    fasta = tmp_path / "tiny.fa"
    fasta.write_text(">chrT\nAACCGGTT\n")
    fasta.with_suffix(".fa.fai").write_text("chrT\t8\t6\t8\t9\n")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text('##description: genome (GRCh38), version 50\n'
                   'chrT\tt\tgene\t1\t8\t.\t+\t.\tgene_id "g";\n')
    out = tmp_path / "output"
    cmd = [sys.executable, "-m", "conference.model_v2.j6_site_universe",
           "--fasta", str(fasta), "--gtf", str(gtf), "--output-dir", str(out),
           "--synthetic-reference", "--region", "chrT:1-8", "--max-mismatch", "0", "--window-size", "4"]
    env = dict(os.environ, PYTHONPATH=str(ROOT / "feng/src"))
    subprocess.run(cmd, check=True, env=env, cwd=ROOT, capture_output=True, text=True)
    files = list(out.glob("*.done"))
    assert len(files) == 2
    mtimes = {p.name: p.stat().st_mtime_ns for p in files}
    subprocess.run(cmd, check=True, env=env, cwd=ROOT, capture_output=True, text=True)
    assert {p.name: p.stat().st_mtime_ns for p in files} == mtimes
    assert sum(pq.read_metadata(p).num_rows for p in out.glob("*.parquet")) == 0


@pytest.mark.parametrize("bad", [-1, 1.5, True])
def test_invalid_budgets_fail_closed(bad):
    j, kwargs = tiny_job(Path("unused"))
    with pytest.raises((TypeError, ValueError)):
        enumerate_synthetic(j, kwargs["genome"], kwargs["genes"], kwargs["asos"], max_mismatch=bad)


def test_duplicate_ASO_ids_and_conflicting_gene_ids_fail_closed():
    j, kw = tiny_job(Path("unused"))
    with pytest.raises(ValueError):
        enumerate_synthetic(j, kw["genome"], kw["genes"], kw["asos"] * 2)
    with pytest.raises(ValueError):
        enumerate_synthetic(j, kw["genome"], kw["genes"] + [
            dataclasses.replace(kw["genes"][0], strand="-")], kw["asos"])


def test_paralog_flag_requires_identical_sequence_on_correct_gene_strand():
    j = module()
    genes = [j.Gene("g1", "chr1", 1, 6, "+", "SMN1"),
             j.Gene("g2", "chr2", 1, 6, "-", "SMN2")]
    genome = ti.InMemoryGenome({"chr1": "AACCGG", "chr2": "CCGGTT"})
    asos = [ti.ASO("q", "CCGGTT")]
    sites = enumerate_synthetic(j, genome, genes, asos, max_mismatch=0)
    assert len(sites) == 2 and all(s.paralog_ambiguity for s in sites)
    genome = ti.InMemoryGenome({"chr1": "AACCGG", "chr2": "CCGGTA"})
    sites = enumerate_synthetic(j, genome, genes, asos, max_mismatch=1)
    assert len(sites) == 2 and not any(s.paralog_ambiguity for s in sites)


@pytest.mark.parametrize("damage", [None, "missing", "target", "duplicate_id", "checksum"])
def test_sentinel_checker_requires_six_committed_exact_anchor_records(tmp_path, damage):
    import pyarrow as pa
    import pyarrow.parquet as pq
    j = module()
    for name, start in [("SMN2", 70076584), ("SMN1", 70952004)]:
        rows = []
        for aso, seq in [("F18", "CCAGCATTATGAAAGTGA"), ("F14", "CCAGCATTATGAAA"),
                         ("F10", "CCAGCATTAT")]:
            rows.append(dict(site_id=name + aso, aso_id=aso, chrom="chr5", start1=start,
                             end1=start + len(seq) - 1, strand="+", gene_name=name,
                             gene_id=name, target_sequence=seq, mismatch_count=0,
                             mismatch_positions=[], paralog_ambiguity=True,
                             evidence_class="derived"))
        if name == "SMN1":
            if damage == "missing":
                rows.pop()
            elif damage == "target":
                rows[0]["target_sequence"] = "A" * 18
            elif damage == "duplicate_id":
                rows[0]["site_id"] = "SMN2F18"
        path = tmp_path / f"{name}.parquet"
        pq.write_table(pa.Table.from_pylist(rows), path)
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        path.with_suffix(".done").write_text(json.dumps(dict(rows=len(rows), sha256=sha)))
    if damage == "checksum":
        (tmp_path / "SMN1.parquet").write_bytes(b"invalid")
    if damage is None:
        assert j.check_smn_anchors(tmp_path) == 6
        # Old output files cannot satisfy acceptance for an empty current run.
        with pytest.raises(ValueError, match="missing sentinel anchors"):
            j.check_smn_anchors(tmp_path, paths=[])
    else:
        with pytest.raises(ValueError):
            j.check_smn_anchors(tmp_path)


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


@pytest.mark.parametrize("change", ["window", "region"])
def test_cli_rejects_changed_run_configuration(tmp_path, change):
    fasta, gtf = synthetic_inputs(tmp_path, sequence="A" * 12)
    out = tmp_path / "out"
    args = ["--fasta", str(fasta), "--gtf", str(gtf), "--output-dir", str(out),
            "--synthetic-reference", "--region", "chrT:1-12", "--window-size", "4"]
    j = module()
    assert j.main(args) == 0
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert len(manifest["expected_shards"]) == 3
    before = {p.name: p.read_bytes() for p in out.glob("*.done")}
    args[-1 if change == "window" else -3] = "6" if change == "window" else "chrT:2-12"
    with pytest.raises(ValueError, match="incompatible run"):
        j.main(args)
    assert {p.name: p.read_bytes() for p in out.glob("*.done")} == before
    assert json.loads((out / "run_manifest.json").read_text()) == manifest


def synthetic_inputs(tmp_path, *, sequence="AACCGGTT"):
    fasta = tmp_path / "tiny.fa"
    fasta.write_text(">chrT\n" + sequence + "\n")
    Path(str(fasta) + ".fai").write_text(f"chrT\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence)+1}\n")
    gtf = tmp_path / "genes.gtf"
    gtf.write_text('##description: genome (GRCh38), version 50\n'
                   f'chrT\tt\tgene\t1\t{len(sequence)}\t.\t+\t.\tgene_id "g";\n')
    return fasta, gtf


def test_mismatched_assembly_rejected_before_enumeration(tmp_path):
    fasta, gtf = synthetic_inputs(tmp_path)
    # GRCh37 chr5 is 180,915,260 bp; even an in-bounds request must fail.
    Path(str(fasta) + ".fai").write_text("chr5\t180915260\t6\t8\t9\n")
    gtf.write_text('##description: genome (GRCh38), version 50\n'
                   'chr5\tt\tgene\t1\t8\t.\t+\t.\tgene_id "g";\n')
    with pytest.raises(ValueError, match="GRCh38"):
        module().main(["--fasta", str(fasta), "--gtf", str(gtf),
                       "--output-dir", str(tmp_path / "out"), "--region", "chr5:1-8"])
    assert not (tmp_path / "out").exists()


def test_single_paralog_evidence_is_unknown_not_false():
    j = module()
    gene = j.Gene("g1", "chrT", 1, 6, "+", "SMN1")
    sites = enumerate_synthetic(j, ti.InMemoryGenome({"chrT": "AACCGG"}), [gene],
                              [ti.ASO("q", "CCGGTT")], max_mismatch=0)
    assert len(sites) == 1
    assert sites[0].paralog_ambiguity is None


def test_region_filtered_cli_retains_other_paralog_evidence(tmp_path):
    import pyarrow.parquet as pq
    target = "CCAGCATTAT"
    fasta, gtf = synthetic_inputs(tmp_path, sequence=target + "NNNN" + target)
    gtf.write_text('##description: genome (GRCh38), version 50\n'
                   'chrT\tt\tgene\t1\t10\t.\t+\t.\tgene_id "g1"; gene_name "SMN1";\n'
                   'chrT\tt\tgene\t15\t24\t.\t+\t.\tgene_id "g2"; gene_name "SMN2";\n')
    j = module()
    j.main(["--fasta", str(fasta), "--gtf", str(gtf), "--output-dir", str(tmp_path / "out"),
            "--synthetic-reference", "--region", "chrT:1-10", "--asos", "F10", "--max-mismatch", "0"])
    rows = [r for p in (tmp_path / "out").glob("*.parquet") for r in pq.read_table(p).to_pylist()]
    assert len(rows) == 1 and rows[0]["gene_id"] == "g1"
    assert rows[0]["paralog_ambiguity"] is True
    assert rows[0]["assembly"] == "synthetic"


def test_streaming_unsharded_output_is_lazy_and_bounded():
    j = module()
    class CountingGenome(ti.InMemoryGenome):
        calls = 0
        def fetch(self, *args):
            self.calls += 1
            return super().fetch(*args)
    genome = CountingGenome({"chrT": "A" * 120000})
    genes = [j.Gene("g", "chrT", 1, 120000, "+")]
    stream = j.iter_sites(genome, genes, [ti.ASO("a", "TTTT")], max_mismatch=0,
                          synthetic_reference=True)
    assert genome.calls == 0
    assert next(stream).start1 == 1
    assert genome.calls == 1
    assert sum(1 for _ in stream) == 119996
    assert genome.calls == 3


def test_worker_caches_universe_and_interval_lookup(tmp_path, monkeypatch):
    j = module()
    genes = [j.Gene(f"g{i}", f"chr{i}", 1, 12, "+") for i in range(2000)]
    genome = ti.InMemoryGenome({g.chrom: "AACCGGTTAACC" for g in genes})
    worker = j.SiteWorker(genome, genes, [ti.ASO("a", "CCGG")], synthetic_reference=True)
    def no_rebuild(*args, **kwargs):
        raise AssertionError("universe rebuilt within worker")
    monkeypatch.setattr(j, "_genes", no_rebuild)
    monkeypatch.setattr(j, "asdict", no_rebuild)
    for i in range(20):
        sites = list(worker.iter_sites(shard=j.Shard(f"chr{i}", 1, 12), max_mismatch=0))
        assert len(sites) == 1 and sites[0].gene_id == f"g{i}"
    # Intervals nested in one chromosome must not lose the long outer gene.
    nested = [j.Gene("long", "chrT", 1, 100, "+"), j.Gene("short", "chrT", 50, 60, "+")]
    monkeypatch.undo()
    worker = j.SiteWorker(ti.InMemoryGenome({"chrT": "A" * 100}), nested,
                          [ti.ASO("a", "TTTT")], synthetic_reference=True)
    assert {s.gene_id for s in worker.iter_sites(shard=j.Shard("chrT", 55, 55), max_mismatch=0)} == {"long", "short"}


def test_synthetic_mode_is_explicit_and_cannot_share_production_run(tmp_path):
    j, kw = tiny_job(tmp_path)
    with pytest.raises(ValueError, match="GRCh38"):
        j.enumerate_sites(kw["genome"], kw["genes"], kw["asos"])
    path = j.run_shard(**kw)
    import pyarrow.parquet as pq
    assert pq.read_table(path).to_pylist()[0]["assembly"] == "synthetic"
    manifest = json.loads((tmp_path / "run_manifest.json").read_text())
    assert manifest["config"]["worker"]["reference"]["mode"] == "synthetic"
    # Mode changes are incompatible even when every other parameter is identical.
    worker = j.SiteWorker(kw["genome"], kw["genes"], kw["asos"], synthetic_reference=True)
    changed = json.loads(worker.params_json)
    changed["reference"]["mode"] = "GRCh38_primary"
    # Check the persisted fence independently of reference construction.
    manifest["config"]["worker"] = changed
    (tmp_path / "run_manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="incompatible run"):
        j.run_shard(**kw)


def test_multiworker_manifest_lists_all_expected_shards_and_reuses_context(tmp_path, monkeypatch):
    import pyarrow.parquet as pq
    j, kw = tiny_job(tmp_path)
    worker = j.SiteWorker(kw["genome"], kw["genes"], kw["asos"], synthetic_reference=True)
    shards = j.make_shards(kw["genes"], window_size=4)
    run = j.prepare_run(tmp_path, worker, kw["identities"], shards, max_mismatch=0,
                        window_size=4, shard_count=2)
    with pytest.raises(dataclasses.FrozenInstanceError):
        run.max_mismatch = 1
    before = (tmp_path / "run_manifest.json").stat().st_mtime_ns
    other = j.prepare_run(tmp_path, worker, kw["identities"], shards, max_mismatch=0,
                          window_size=4, shard_count=2)
    assert before == (tmp_path / "run_manifest.json").stat().st_mtime_ns
    assert json.loads(run.manifest_json)["expected_shards"] == [
        "chrT.000000000001-000000000004", "chrT.000000000005-000000000008"]
    def no_rebuild(*args, **kwargs):
        raise AssertionError("repeated universe preparation")
    monkeypatch.setattr(j, "_genes", no_rebuild)
    for context, shard in zip((run, other), shards):
        path = j.run_shard(worker=worker, run=context, shard=shard, output_dir=tmp_path,
                           identities=kw["identities"], max_mismatch=0)
        assert path.with_suffix(".done").is_file()
    assert sum(pq.read_metadata(p).num_rows for p in tmp_path.glob("*.parquet")) == 1
    with pytest.raises(ValueError, match="incompatible run/shard"):
        j.run_shard(worker=worker, run=run, shard=j.Shard("chrT", 1, 8), output_dir=tmp_path,
                    identities=kw["identities"], max_mismatch=0)


def test_legacy_and_unlisted_shards_rejected_without_deletion(tmp_path):
    j, kw = tiny_job(tmp_path)
    legacy = tmp_path / "old.done"
    legacy.write_text("{}")
    with pytest.raises(ValueError, match="legacy"):
        j.run_shard(**kw)
    assert legacy.read_text() == "{}"
    legacy.unlink()
    j.run_shard(**kw)
    legacy.write_text("{}")
    with pytest.raises(ValueError, match="unexpected shard"):
        j.run_shard(**kw)
    assert legacy.read_text() == "{}"


def test_interval_index_prunes_large_disjoint_universe(monkeypatch):
    j = module()
    genes = [j.Gene(f"g{i}", "chrT", 20*i + 1, 20*i + 12, "+") for i in range(8000)]
    worker = j.SiteWorker(ti.InMemoryGenome({"chrT": "A" * 160000}), genes,
                          [ti.ASO("a", "TTTT")], synthetic_reference=True)
    original = j._overlapping
    calls = 0
    def counted(*args):
        nonlocal calls
        calls += 1
        yield from original(*args)
    monkeypatch.setattr(j, "_overlapping", counted)
    # Select the middle gene, independent of how many other genes exist.
    rows = list(worker.iter_sites(shard=j.Shard("chrT", 80001, 80001), max_mismatch=0))
    assert len(rows) == 1 and rows[0].gene_id == "g4000"
    assert calls < 50


def test_synthetic_cli_cannot_pass_production_anchor_check(tmp_path):
    with pytest.raises(SystemExit) as error:
        module().main(["--output-dir", str(tmp_path), "--synthetic-reference", "--check-smn-anchors"])
    assert error.value.code == 2


def test_gtf_reader_releases_consumed_cache_in_bounded_chunks(tmp_path, monkeypatch):
    j = module()
    if not hasattr(j.os, "posix_fadvise"):
        pytest.skip("POSIX cache advice unavailable")
    path = tmp_path / "large.gtf"
    with path.open("w") as handle:
        handle.write('##description: genome (GRCh38), version 50\n')
        for _ in range(4000):
            handle.write('#' + 'x' * 100 + '\n')
        handle.write('chrT\tt\tgene\t1\t8\t.\t+\t.\tgene_id "g";\n')
    monkeypatch.setattr(j, "_GTF_CACHE_CHUNK", 65536)
    advice = j.os.posix_fadvise
    calls = []
    def record(fd, start, length, kind):
        calls.append((start, length, kind))
        return advice(fd, start, length, kind)
    monkeypatch.setattr(j.os, "posix_fadvise", record)
    assert len(j.read_gene_bodies(path)) == 1
    releases = [(start, length) for start, length, kind in calls if kind == j.os.POSIX_FADV_DONTNEED]
    assert releases, "a sequential multi-GB GTF must not leave its full scan in the job's cache"
    assert len(releases) > 1
    assert all(0 < length < 65700 for _, length in releases)  # zero means the entire file
    assert sum(length for _, length in releases) == path.stat().st_size


def test_run_manifest_is_cached_across_shards(tmp_path, monkeypatch):
    j, kw = tiny_job(tmp_path)
    worker = j.SiteWorker(kw["genome"], kw["genes"], kw["asos"], synthetic_reference=True)
    run = j.prepare_run(tmp_path, worker, kw["identities"], [kw["shard"]], max_mismatch=0)
    read_text = Path.read_text
    def read_once(path, *args, **kwargs):
        if path.name == "run_manifest.json":
            raise AssertionError("whole expected-shard manifest re-read per shard")
        return read_text(path, *args, **kwargs)
    monkeypatch.setattr(Path, "read_text", read_once)
    path = j.run_shard(worker=worker, run=run, shard=kw["shard"], output_dir=tmp_path,
                       identities=kw["identities"], max_mismatch=0)
    assert path.with_suffix(".done").is_file()


def test_run_manifest_mutation_during_scan_cannot_commit(tmp_path):
    j, kw = tiny_job(tmp_path)
    class ChangingGenome(ti.InMemoryGenome):
        def fetch(self, *args):
            (tmp_path / "run_manifest.json").write_text("{}")
            return super().fetch(*args)
    worker = j.SiteWorker(ChangingGenome({"chrT": "AACCGGTT"}), kw["genes"], kw["asos"], synthetic_reference=True)
    run = j.prepare_run(tmp_path, worker, kw["identities"], [kw["shard"]], max_mismatch=0)
    with pytest.raises(ValueError, match="immutable manifest changed"):
        j.run_shard(worker=worker, run=run, shard=kw["shard"], output_dir=tmp_path,
                    identities=kw["identities"], max_mismatch=0)
    assert not list(tmp_path.glob("*.done"))
