"""Synthetic reproductions for FIX-J7-J9-J10; never read production results."""
import hashlib
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from conference.model_v2 import j7_context_features as j7
from conference.model_v2 import j9_smn2_prediction as j9
from conference.model_v2 import j10_offtarget_profile as j10
from conference.model_v2.tests.test_j7_context_features import row, tx
from conference.model_v2.tests.test_j9_smn2_prediction import SparseGenome, _patches, SMN1


def j3_source(path, index=0, count=1, reference=None, inputs=None):
    provenance = dict(shard_index=index, shard_count=count,
                      reference=reference or {'path': '/synthetic/ref.fa', 'fai_sha256': 'a'*64},
                      inputs=inputs or [{'path': '/synthetic/j6.parquet', 'sha256': 'b'*64}],
                      params={'window': 120, 'span': 60, 'flank': 80})
    table = pa.Table.from_pylist([dict(row(150), site_id=f's{index}')]).replace_schema_metadata(
        {b'j3_provenance': json.dumps(provenance).encode()})
    pq.write_table(table, path)
    j7.atomic_json(path.with_suffix('.done'), dict(provenance, sha256=j7.sha256(path), rows=1))
    return provenance


def j7_fixture(tmp_path):
    src = tmp_path/'src'; src.mkdir()
    ex = tmp_path/'ex.json'
    j7.atomic_json(ex, dict(transcripts={'g': tx()}, provenance={'test': True, 'implementation_sha256':j7.implementation_hash()}))
    return src, ex, tmp_path/'out', {'status': 'BLOCKED'}


def test_j7_rejects_producer_one_of_two_as_zero_of_one(tmp_path):
    src, ex, out, audit = j7_fixture(tmp_path)
    path = src/'j3_shard0000_of0001.parquet'
    j3_source(path, 1, 2)
    with pytest.raises(ValueError, match='producer|source|shard'):
        j7.run_shard(path, ex, out, 0, 1, audit)
    assert not list(out.glob('*.done'))


@pytest.mark.parametrize('field', ['reference', 'inputs'])
def test_j7_merge_rejects_mixed_producer_identity(tmp_path, field):
    src, ex, out, audit = j7_fixture(tmp_path)
    for i in range(2):
        path = src/f'j3_shard{i:04d}_of0002.parquet'
        changes = {field: ({'path': '/other/ref.fa', 'fai_sha256': 'c'*64} if field == 'reference'
                           else [{'path': '/other/j6.parquet', 'sha256': 'd'*64}])} if i else {}
        j3_source(path, i, 2, **changes)
        j7.run_shard(path, ex, out, i, 2, audit)
    with pytest.raises(ValueError, match='provenance|identity'):
        j7.merge_gate(out, 2, 2, input_dir=src, exons=ex)
    assert not (out/'j7_context.done').exists()


def test_j7_implementation_change_recomputes_cache_and_binds_merge(tmp_path, monkeypatch):
    src, ex, out, audit = j7_fixture(tmp_path)
    path = src/'j3_shard0000_of0001.parquet'; j3_source(path)
    monkeypatch.setattr(j7, 'implementation_hash', lambda: '1'*64, raising=False)
    data = json.loads(ex.read_text())
    data['provenance']['implementation_sha256'] = '1'*64
    j7.atomic_json(ex, data)
    result = j7.run_shard(path, ex, out, 0, 1, audit)
    original = j7.geometry
    def changed(site, transcript, *args, **kwargs):
        return dict(original(site, transcript, *args, **kwargs), donor_distance=999,
                    donor_bin='201-1000')
    monkeypatch.setattr(j7, 'geometry', changed)
    monkeypatch.setattr(j7, 'implementation_hash', lambda: '2'*64, raising=False)
    data = json.loads(ex.read_text())
    data['provenance']['implementation_sha256'] = '2'*64
    j7.atomic_json(ex, data)
    j7.run_shard(path, ex, out, 0, 1, audit)
    assert pq.read_table(result)['donor_distance'][0].as_py() == 999
    merged = j7.merge_gate(out, 1, 1, input_dir=src, exons=ex)
    assert j7.verify_done(merged)['provenance']['common']['implementation_sha256'] == '2'*64
    monkeypatch.setattr(j7, 'implementation_hash', lambda: '3'*64)
    with pytest.raises(ValueError, match='implementation|provenance'):
        j7.merge_gate(out, 1, 1, input_dir=src, exons=ex)


@pytest.mark.parametrize('relative', [
    'data/wetlab/SMN2_Exon7_7_27+28_Complete_Data.xlsx',
    'data/outcomes/values.parquet', 'outcomes/values.tsv', 'ops/WETLAB_MANIFEST.json',
    'tests/not_a_source.py', 'jobs/unlisted_source.py',
    'j9_smn2_v2/outcomes/values.tsv'])
def test_j9_in_tree_forbidden_sources_rejected_before_io(relative, monkeypatch):
    path = j9.REPO_ROOT/'conference/model_v2'/relative
    def no_open(*args, **kwargs):
        pytest.fail('forbidden source reached I/O')
    monkeypatch.setattr('builtins.open', no_open)
    with pytest.raises(j9.ForbiddenSourceError):
        j9._sha256(path)


def atlas_rows(placement, kinds=('donor',)):
    result = []
    for kind in kinds:
        position = getattr(placement.geometry, kind+'_last_exonic_base' if kind == 'donor'
                           else 'acceptor_first_exonic_base')
        for tissue in (0, 2, 4, 6):
            result.append(dict(row_id=f'{kind}-{tissue}', locus='SMN2_exon7', chrom=placement.chrom,
                strand=placement.strand, site_kind=kind, site_position=position, tissue=tissue,
                window_start=70071000, window_end=70082000,
                footprint_start=placement.start, footprint_end=placement.end,
                length=placement.end-placement.start+1, reference_score=.5,
                masked_score=.6 if kind == 'donor' else .1,
                raw_model_delta=.1 if kind == 'donor' else -.4,
                status='ok', error=None, evidence_class='sequence_ablation_proxy',
                evidence_category='proxy', model_evidence=True))
    return result


def committed_atlas(directory, rows, *, genuine=True):
    from conference.model_v2.gpu_perturbation_atlas import row_schema
    directory.mkdir(parents=True, exist_ok=True)
    provenance = dict(model_evidence=genuine, model='Pangolin-1.0.2-splice-probability-four-tissue-three-replicate',
        model_sha256='a'*64, weights_sha256='b'*64, source_sha256='c'*64,
        output_channels=[1, 4, 7, 10])
    plan = {'version': 1, 'provenance': provenance}
    path = directory/'ablation.parquet'
    pq.write_table(pa.Table.from_pylist(rows, schema=row_schema(provenance)), path)
    marker = dict(run_id=j9._digest(plan), plan=plan, provenance=provenance,
                  expected_rows=len(rows), status='complete',
                  ablation={'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                            'rows': len(rows)})
    marker['integrity_sha256'] = j9._digest(marker)
    (directory/'shard_0.done').write_text(json.dumps(marker))
    return marker


@pytest.fixture
def placed():
    genome = SparseGenome(SMN1.search_end+1000, _patches())
    aso = j9.aso_registry()['F18']
    return aso, j9.place_asos(genome, {'F18': aso})['F18']['placement']


def m5_from_committed(tmp_path, monkeypatch, placed, rows, genuine=True):
    monkeypatch.setattr(j9, '_allowed_roots', lambda: (tmp_path,))
    directory = tmp_path/'atlas'; committed_atlas(directory, rows, genuine=genuine)
    atlas, identities = j9.read_atlas(directory)
    return j9.m5_rows(*placed, atlas, atlas_source='synthetic', atlas_inputs=identities)


def test_j9_m5_false_model_evidence_is_unavailable(tmp_path, monkeypatch, placed):
    rows = atlas_rows(placed[1])
    for r in rows:
        r.update(model_evidence=False, evidence_category='not_available')
    result = m5_from_committed(tmp_path, monkeypatch, placed, rows)
    assert all(r['status'] == 'not_available' and r['value_double'] is None for r in result)
    assert j9.build_ranks(result) == []


@pytest.mark.parametrize('corruption', ['integrity', 'run_id', 'schema_provenance', 'stub'])
def test_j9_atlas_requires_valid_j4_provenance(tmp_path, monkeypatch, placed, corruption):
    monkeypatch.setattr(j9, '_allowed_roots', lambda: (tmp_path,))
    marker = committed_atlas(tmp_path, atlas_rows(placed[1]), genuine=corruption != 'stub')
    if corruption == 'integrity': marker['expected_rows'] += 1
    if corruption == 'run_id':
        marker['run_id'] = 'e'*64
        marker['integrity_sha256'] = j9._digest({k:v for k,v in marker.items() if k != 'integrity_sha256'})
    if corruption == 'schema_provenance':
        path = tmp_path/'ablation.parquet'
        table = pq.read_table(path).replace_schema_metadata({b'provenance': b'{}'})
        pq.write_table(table, path)
        marker['ablation']['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
        marker['integrity_sha256'] = j9._digest({k:v for k,v in marker.items() if k != 'integrity_sha256'})
    (tmp_path/'shard_0.done').write_text(json.dumps(marker))
    with pytest.raises(j9.AtlasCoverageError):
        j9.read_atlas(tmp_path)


def test_j9_m5_retains_separate_endpoints(tmp_path, monkeypatch, placed):
    result = m5_from_committed(tmp_path, monkeypatch, placed, atlas_rows(placed[1], ('donor','acceptor')))
    values = {r['feature_name']: r['value_double'] for r in result}
    assert len(values) == len(result)
    assert values['m5_donor_raw_model_delta_mean'] == pytest.approx(.1)
    assert values['m5_acceptor_raw_model_delta_mean'] == pytest.approx(-.4)
    assert values['m5_raw_model_delta_mean'] == pytest.approx(.1)
    primary = next(r for r in result if r['feature_name'] == 'm5_raw_model_delta_mean')
    assert 'primary endpoint=donor' in primary['reason']


@pytest.mark.parametrize('changes', [dict(chrom='chr1'), dict(strand='-'),
    dict(site_position=70076575), dict(site_kind='other'), dict(tissue=1)])
def test_j9_m5_rejects_wrong_endpoint_identity(tmp_path, monkeypatch, placed, changes):
    rows = atlas_rows(placed[1]); rows[0].update(changes)
    with pytest.raises(j9.AtlasCoverageError):
        m5_from_committed(tmp_path, monkeypatch, placed, rows)


def test_j9_m5_rejects_duplicate_tissue_key(tmp_path, monkeypatch, placed):
    rows = atlas_rows(placed[1]); rows.append(dict(rows[0]))
    with pytest.raises(j9.AtlasCoverageError, match='duplicate'):
        m5_from_committed(tmp_path, monkeypatch, placed, rows)


def test_j9_m5_incomplete_tissues_unavailable(tmp_path, monkeypatch, placed):
    result = m5_from_committed(tmp_path, monkeypatch, placed, atlas_rows(placed[1])[:1])
    assert all(r['status'] == 'not_available' for r in result)
    assert j9.build_ranks(result) == []


@pytest.mark.parametrize('baseline,feature', [('M2','m2_duplex_energy_kcal_mol'),
                                           ('M3','m3_accessibility_p_unpaired')])
@pytest.mark.parametrize('case', ['one_variant', 'membership', 'complete'])
def test_j9_stability_requires_full_comparable_alternatives(baseline, feature, case):
    features = []
    for variant, *_ in j9.VIENNA_ALTERNATIVES:
        if case == 'one_variant' and variant != 'w80s40': continue
        for aso_id in j9.RANKED_ASO_IDS:
            if case == 'membership' and variant == 'w200s100' and aso_id == 'F10': continue
            features.append(j9._row(aso=j9.aso_registry()[aso_id], baseline=baseline,
                feature_name=feature, variant_id=variant, value_double=1., evidence_class='proxy',
                status='ok', reason='synthetic', source='synthetic'))
    stability = [r for r in j9.build_ranks(features) if r['section'] == 'stability']
    assert stability
    if case == 'complete':
        assert all(r['stable'] is True and r['n_variants'] == 3 for r in stability)
    else:
        assert all(r['stable'] is None and r['status'] == 'not_available'
                   and r['rank_range'] is None for r in stability)


@pytest.mark.parametrize('mutate', [False, True])
def test_j9_build_records_and_rechecks_reference_identity(tmp_path, monkeypatch, mutate):
    monkeypatch.setattr(j9, '_allowed_roots', lambda: (tmp_path,))
    path = tmp_path/'ref.fa'; path.write_text('>chr1\nACGT\n')
    index = tmp_path/'ref.fa.fai'; index.write_text('chr1\t4\t6\t4\t5\n')
    genome = j9.FaiFastaReader(path)
    original = j9.place_asos
    def placement(*args, **kwargs):
        if mutate: index.write_text('chr1\t3\t6\t4\t5\n')
        return original(*args, **kwargs)
    monkeypatch.setattr(j9, 'place_asos', placement)
    if mutate:
        with pytest.raises(ValueError, match='reference.*changed'):
            j9.build_features(genome=genome, backend=j9.NullBackend(), registry={},
                              cache_path=tmp_path, atlas_dir=tmp_path/'missing')
    else:
        _, provenance = j9.build_features(genome=genome, backend=j9.NullBackend(), registry={},
                                         cache_path=tmp_path, atlas_dir=tmp_path/'missing')
        assert provenance['reference']['path'] == str(path)
        assert provenance['reference']['fai_sha256'] == hashlib.sha256(index.read_bytes()).hexdigest()
        assert provenance['reference']['size'] == 11


@pytest.mark.parametrize('shards,expected', [(0,0), (-1,0), (1,-1)])
def test_j10_invalid_counts_do_not_touch_outputs(tmp_path, shards, expected):
    out = tmp_path/'existing'; out.mkdir()
    marker = out/'j10_profile.done'; marker.write_text('prior commitment')
    with pytest.raises(ValueError, match='shard|expected'):
        j10.merge(out, shard_count=shards, expected_rows=expected)
    assert marker.read_text() == 'prior commitment'
    assert sorted(p.name for p in out.iterdir()) == ['j10_profile.done']


def test_j10_production_cli_rejects_non256_before_io(tmp_path):
    import subprocess, sys
    out = tmp_path/'absent'
    result = subprocess.run([sys.executable, '-m', 'conference.model_v2.j10_offtarget_profile',
        'merge', '--shard-count', '1', '--expected-rows', '0', '--output-dir', str(out)],
        capture_output=True, text=True)
    assert result.returncode != 0 and '256' in result.stderr
    assert not out.exists()

@pytest.mark.parametrize('excluded', ['wetlab', 'outcomes'])
def test_j9_excluded_root_rejected_before_filesystem_resolution(excluded, monkeypatch):
    path = j9.REPO_ROOT/'conference/model_v2/data'/excluded/'never_read'
    def forbidden_io(*args, **kwargs):
        pytest.fail('excluded root reached filesystem metadata I/O')
    monkeypatch.setattr(Path, 'resolve', forbidden_io)
    monkeypatch.setattr(Path, 'stat', forbidden_io)
    with pytest.raises(j9.ForbiddenSourceError):
        j9.check_source_allowed(path)


def test_j9_symlink_into_excluded_root_is_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(j9, '_allowed_roots', lambda: (tmp_path,))
    target = tmp_path/'wetlab'/'never_read'
    link = tmp_path/'alias'; link.symlink_to(target)
    with pytest.raises(j9.ForbiddenSourceError):
        j9._sha256(link)


@pytest.mark.parametrize('field', ['reference', 'inputs'])
def test_j7_merge_rechecks_expected_source_identity(tmp_path, field):
    src, ex, out, audit = j7_fixture(tmp_path)
    path = src/'j3_shard0000_of0001.parquet'; j3_source(path)
    j7.run_shard(path, ex, out, 0, 1, audit)
    j7.merge_gate(out, 1, 1, input_dir=src, exons=ex)
    change = {field: {'path': '/changed/ref.fa', 'fai_sha256': 'e'*64} if field == 'reference'
              else [{'path': '/changed/j6.parquet', 'sha256': 'f'*64}]}
    j3_source(path, **change)
    with pytest.raises(ValueError, match='provenance'):
        j7.merge_gate(out, 1, 1, input_dir=src, exons=ex)
    assert not (out/'j7_context.done').exists()


def test_j9_outcome_named_file_in_allowed_root_is_refused(tmp_path, monkeypatch):
    monkeypatch.setattr(j9, '_allowed_roots', lambda: (tmp_path,))
    with pytest.raises(j9.ForbiddenSourceError):
        j9._sha256(tmp_path/'outcomes.parquet')


def test_j9_primary_donor_never_falls_back_to_acceptor(tmp_path, monkeypatch, placed):
    result = m5_from_committed(tmp_path, monkeypatch, placed, atlas_rows(placed[1], ('acceptor',)))
    assert not any(r['feature_name'] in ('m5_raw_model_delta_mean', 'm5_raw_model_delta_max')
                   and r['status'] == 'ok' for r in result)
