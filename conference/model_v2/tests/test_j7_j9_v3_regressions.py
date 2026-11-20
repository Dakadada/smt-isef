"""Builder authentication and interrupted/concurrent generation publication."""
import json
from pathlib import Path
import pytest
from conference.model_v2 import j7_context_features as j7
from conference.model_v2 import j9_smn2_prediction as j9
from conference.model_v2.tests.test_j7_j9_j10_regressions import j3_source
from conference.model_v2.tests.test_j7_context_features import tx


def exon_map(tmp_path):
    gtf = tmp_path/'input.gtf'; manifest = tmp_path/'loci.tsv'; out = tmp_path/'exons.json'
    gtf.write_text('##description: GRCh38 version 50\nchr1\tx\texon\t101\t110\t.\t+\t.\tgene_id "g"; transcript_id "t"; exon_number "1"; tag "Ensembl_canonical"; gene_type "protein_coding";\n')
    manifest.write_text('t:ex1\tchr1\t1\t210\t+\n')
    j7.build_exons(gtf, [manifest], out)
    return out


def test_builder_records_implementation(tmp_path):
    assert json.loads(exon_map(tmp_path).read_text())['provenance']['implementation_sha256'] == j7.implementation_hash()


@pytest.mark.parametrize('bad', [None, '0'*64])
@pytest.mark.parametrize('stage', ['shard', 'merge'])
def test_rejects_unauthenticated_builder(tmp_path, bad, stage):
    ex = exon_map(tmp_path)
    source = tmp_path/'j3_shard0000_of0001.parquet'; j3_source(source)
    out = tmp_path/'out'
    data = json.loads(ex.read_text())
    data['provenance']['implementation_sha256'] = bad
    if stage == 'merge':
        # Simulate a self-consistent shard from an unauthenticated map.
        original = j7.read_exons
        j7.read_exons = lambda p: data
        try:
            j7.atomic_json(ex, data)
            j7.run_shard(source, ex, out, 0, 1, {'status': 'BLOCKED'})
        finally:
            j7.read_exons = original
    else:
        j7.atomic_json(ex, data)
    with pytest.raises(ValueError, match='implementation'):
        if stage == 'shard':
            j7.run_shard(source, ex, out, 0, 1, {'status': 'BLOCKED'})
        else:
            j7.merge_gate(out, 1, 1, input_dir=tmp_path, exons=ex)


@pytest.fixture
def output(tmp_path, monkeypatch):
    roots = j9._allowed_roots()
    monkeypatch.setattr(j9, '_allowed_roots', lambda: roots + (tmp_path,))
    out = tmp_path/'dataset'
    j9.write_dataset(out, [], [], {'generation': 1})
    return out


def test_interrupted_write_preserves_published_generation(output, monkeypatch):
    before = {n: (output/n).read_bytes() for n in ('manifest.json', 'features.parquet', 'ranks.parquet')}
    original = j9._write_parquet
    def interrupt(path, *args):
        if path.name == 'ranks.parquet':
            raise RuntimeError('interrupted')
        original(path, *args)
    monkeypatch.setattr(j9, '_write_parquet', interrupt)
    with pytest.raises(RuntimeError, match='interrupted'):
        j9.write_dataset(output, [{'aso_id': 'changed'}], [], {'generation': 2})
    assert {n: (output/n).read_bytes() for n in before} == before


def test_snapshot_survives_publish_and_checks_hashes(output):
    snapshot = j9.read_dataset(output)
    j9.write_dataset(output, [{'aso_id': 'new'}], [], {'generation': 2})
    assert snapshot['manifest']['provenance']['generation'] == 1
    assert snapshot['features'].num_rows == 0
    current = j9.read_dataset(output)
    assert current['manifest']['provenance']['generation'] == 2
    assert current['features'].num_rows == 1
    (output/'features.parquet').write_bytes(b'corrupt')
    with pytest.raises(ValueError, match='hash|integrity'):
        j9.read_dataset(output)


def test_reader_during_publication_sees_old_generation(output, monkeypatch):
    original = j9._write_parquet
    observations = []
    def observe(path, *args):
        original(path, *args)
        observations.append(j9.read_dataset(output)['manifest']['provenance']['generation'])
    monkeypatch.setattr(j9, '_write_parquet', observe)
    j9.write_dataset(output, [{'aso_id': 'new'}], [], {'generation': 2})
    assert observations == [1, 1]
    assert j9.read_dataset(output)['manifest']['provenance']['generation'] == 2


def test_fresh_v3_root_allowed():
    path = Path.home()/'orcd/scratch/feng/model_v2/j9_smn2_v3'
    assert j9.check_source_allowed(path) == path.resolve()


def test_reader_pins_generation_across_concurrent_publish(output, monkeypatch):
    original = j9._verify_output
    published = False
    def publish_during_read(path, metadata):
        nonlocal published
        if not published:
            published = True
            j9.write_dataset(output, [{'aso_id': 'new'}], [], {'generation': 2})
        original(path, metadata)
    monkeypatch.setattr(j9, '_verify_output', publish_during_read)
    snapshot = j9.read_dataset(output)
    assert snapshot['manifest']['provenance']['generation'] == 1
    assert snapshot['features'].num_rows == snapshot['ranks'].num_rows == 0
    assert j9.read_dataset(output)['manifest']['provenance']['generation'] == 2


def test_interrupted_pointer_switch_preserves_generation(output, monkeypatch):
    original = j9.os.replace
    def interrupt(src, dst):
        if Path(dst).name == '.current':
            raise OSError('publish interrupted')
        return original(src, dst)
    monkeypatch.setattr(j9.os, 'replace', interrupt)
    with pytest.raises(OSError, match='publish interrupted'):
        j9.write_dataset(output, [{'aso_id': 'new'}], [], {'generation': 2})
    snapshot = j9.read_dataset(output)
    assert snapshot['manifest']['provenance']['generation'] == 1
    assert snapshot['features'].num_rows == 0
