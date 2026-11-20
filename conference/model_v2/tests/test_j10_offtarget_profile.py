"""Hand-counted, synthetic J3 shards only."""
import importlib
from pathlib import Path
import sys

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))


def module():
    return importlib.import_module('conference.model_v2.j10_offtarget_profile')


def row(site='a', **changes):
    r = dict(site_id=site, aso_id='F18', assembly='GRCh38', annotation='GENCODE v50',
             chrom='chr1', start1=10, end1=27, strand='+', gene_id='g', gene_name='G',
             target_sequence='CCAGCATTATGAAAGTGA', mismatch_count=0,
             region_class='not_available', complementarity_status='ok',
             complementarity_mismatch_count=0, duplex_status='ok',
             duplex_energy_kcal_mol=-10., accessibility_status='ok', accessibility_p_unpaired=.1)
    return dict(r, **changes)


def totals(summary, scope='other_sites', mismatch=0, aso='F18'):
    return next(r for r in summary['counts'] if r['aso_id'] == aso and r['scope'] == scope
                and r['mismatch_count'] == mismatch and r['dimension'] == 'all'
                and r['energy_threshold'] == -10 and r['accessibility_threshold'] == .1)


def test_hand_counts_failures_duplicates_and_separate_genes():
    m = module()
    rows = [row(), row(), row('b', duplex_status='failed', duplex_energy_kcal_mol=None),
            row('c', accessibility_status=None, accessibility_p_unpaired=None),
            row('d', gene_id='g2', gene_name='OTHER', duplex_energy_kcal_mol=-9)]
    s = m.summarize(rows, energies=(-10,), accessibilities=(.1,))
    t = totals(s)
    assert (s['input_rows'], s['unique_rows'], s['duplicate_rows']) == (5, 4, 1)
    assert (t['sites'], t['duplex_evaluated'], t['duplex_pass'], t['accessibility_evaluated'],
            t['accessibility_pass'], t['joint_pass']) == (4, 3, 2, 3, 3, 1)
    assert t['accessibility_fraction'] == 1
    assert totals(s, aso='F14', mismatch=3)['accessibility_fraction'] is None
    assert sum(r['sites'] for r in s['counts'] if r['dimension'] == 'channel_status'
               and r['category'] == 'duplex:failed') == 1
    assert sum(r['sites'] for r in s['counts'] if r['dimension'] == 'channel_status'
               and r['category'] == 'accessibility:missing') == 1


def test_threshold_boundaries():
    m = module()
    for inclusive, expected in [(True, 1), (False, 0)]:
        t = totals(m.summarize([row()], energies=(-10,), accessibilities=(.1,), inclusive=inclusive))
        assert (t['duplex_pass'], t['accessibility_pass']) == (expected, expected)


def test_anchors_and_smn_gene_rows():
    m = module()
    rows = [row('on', chrom='chr5', start1=70076584, end1=70076601, gene_name='SMN2'),
            row('para', chrom='chr5', start1=70952004, end1=70952021, gene_name='SMN1'),
            row('near', gene_name='SMN2')]
    s = m.summarize(rows, energies=(-10,), accessibilities=(.1,))
    assert totals(s, 'SMN2_on_target')['sites'] == 1
    assert totals(s, 'SMN1_paralog_anchor')['sites'] == 1
    assert totals(s, 'SMN2_other_sites')['sites'] == 1
    assert len(s['smn_rows']) == 3


@pytest.mark.parametrize('changes', [dict(duplex_energy_kcal_mol=None),
    dict(accessibility_p_unpaired=float('nan')), dict(mismatch_count=4),
    dict(complementarity_mismatch_count=1), dict(assembly='hg19')])
def test_invalid_ok_values_fail_closed(changes):
    with pytest.raises(ValueError):
        module().summarize([row(**changes)])


def test_conflicting_duplicate_rejected():
    with pytest.raises(ValueError, match='duplicate'):
        module().summarize([row(), row(duplex_energy_kcal_mol=-11)])


def write_source(path, rows, index, count):
    m = module()
    pq.write_table(pa.Table.from_pylist(rows), path)
    m.atomic_json(path.with_suffix('.done'), dict(sha256=m.sha256(path), rows=len(rows),
        shard_index=index, shard_count=count, params=dict(window=120, span=60, flank=80)))


@pytest.mark.parametrize('reader', ['load', 'merge'])
def test_reader_keeps_committed_generation_during_replacement(tmp_path, monkeypatch, reader):
    """A writer scheduled after hashing must not replace the payload before reading."""
    import fcntl
    import json

    m = module()
    src, out = tmp_path / 'src', tmp_path / 'out'
    src.mkdir()
    write_source(src / 'j3_shard0000_of0001.parquet', [row()], 0, 1)
    original = m.partial(src, out, 0, shard_count=1)
    base = out / 'j10_partial0000'
    identity = json.loads(base.with_suffix('.done').read_text())['identity']
    replacement = dict(original, counts=m.summarize([row(accessibility_p_unpaired=0)])['counts'])
    replacement_identity = dict(identity, source=dict(identity['source'], inputs_sha256='generation-B'))
    read_metadata = pq.read_metadata

    def replace_between_hash_and_payload(path, *args, **kwargs):
        if Path(path) == base.with_suffix('.parquet'):
            with base.with_suffix('.lock').open('a') as lock:
                try:
                    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                except BlockingIOError:
                    pass  # The reader holds the writer's lock: replacement must wait.
                else:
                    m.commit(base, replacement, replacement_identity)
        return read_metadata(path, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(pq, 'read_metadata', replace_between_hash_and_payload)
        if reader == 'load':
            result, marker = m.load(base, identity)
            assert marker['identity'] == identity
        else:
            result = m.merge(out, shard_count=1, expected_rows=1)
            marker = json.loads((out / 'j10_profile.done').read_text())
            assert marker['identity']['partial_markers'][0]['identity'] == identity
        assert totals(result)['accessibility_pass'] == 1

    # Replacement can commit once the reader has released its snapshot lock.
    with m.locked(base):
        m.commit(base, replacement, replacement_identity)
    assert totals(m.load(base, replacement_identity)[0])['accessibility_pass'] == 0
    with pytest.raises(ValueError, match='incompatible partial identity'):
        m.load(base, identity)


def test_stage_merge_single_pass_recovery_and_completeness(tmp_path):
    m = module()
    src, out = tmp_path / 'src', tmp_path / 'out'
    src.mkdir()
    buckets = [[], []]
    rows = [row(str(i), mismatch_count=i % 4, complementarity_mismatch_count=i % 4) for i in range(12)]
    for r in rows:
        buckets[m.partition(r['site_id'], 2)].append(r)
    for i, rs in enumerate(buckets):
        write_source(src / f'j3_shard{i:04d}_of0002.parquet', rs, i, 2)
    m.partial(src, out, 0, shard_count=2)
    with pytest.raises(ValueError, match='partial'):
        m.merge(out, shard_count=2, expected_rows=12)
    m.partial(src, out, 1, shard_count=2)
    merged = m.merge(out, shard_count=2, expected_rows=12)
    assert merged['counts'] == m.summarize(rows)['counts']
    done = out / 'j10_partial0000.done'
    before = done.read_bytes()
    m.partial(src, out, 0, shard_count=2)
    assert done.read_bytes() == before
    (out / 'j10_partial0000.parquet').write_bytes(b'broken')
    with pytest.raises(ValueError):
        m.merge(out, shard_count=2, expected_rows=12)
    m.partial(src, out, 0, shard_count=2)
    assert m.merge(out, shard_count=2, expected_rows=12)['unique_rows'] == 12
    with pytest.raises(ValueError, match='row'):
        m.merge(out, shard_count=2, expected_rows=13)


def test_default_gate_requires_all_256(tmp_path):
    with pytest.raises(ValueError, match='partial completeness'):
        module().merge(tmp_path)


def test_wrong_partition_and_source_corruption(tmp_path):
    m = module()
    src, out = tmp_path / 'src', tmp_path / 'out'
    src.mkdir()
    r = row()
    wrong = 1 - m.partition(r['site_id'], 2)
    path = src / f'j3_shard{wrong:04d}_of0002.parquet'
    write_source(path, [r], wrong, 2)
    with pytest.raises(ValueError, match='ownership'):
        m.partial(src, out, wrong, shard_count=2)
    path.write_bytes(b'corrupt')
    with pytest.raises(ValueError):
        m.partial(src, out, wrong, shard_count=2)


def test_failed_channels_do_not_use_stale_numeric_values():
    m = module()
    r = row(duplex_status='failed', accessibility_status='not_applicable',
            complementarity_status='failed', complementarity_mismatch_count=None)
    t = totals(m.summarize([r], energies=(-10,), accessibilities=(.1,)))
    assert t['sites'] == 1
    assert t['duplex_pass'] == t['accessibility_pass'] == t['joint_pass'] == 0
    assert t['accessibility_fraction'] is None


def test_different_genes_at_same_coordinates_remain_distinct():
    m = module()
    t = totals(m.summarize([row(), row('b', gene_id='other')], energies=(-10,), accessibilities=(.1,)))
    assert t['sites'] == 2


def test_mixed_j3_configuration_rejected(tmp_path):
    m = module()
    src, out = tmp_path / 'src', tmp_path / 'out'
    src.mkdir()
    for i in range(2):
        sid = next(str(j) for j in range(20) if m.partition(str(j), 2) == i)
        path = src / f'j3_shard{i:04d}_of0002.parquet'
        write_source(path, [row(sid)], i, 2)
        marker = __import__('json').loads(path.with_suffix('.done').read_text())
        marker['params']['window'] += i
        m.atomic_json(path.with_suffix('.done'), marker)
        m.partial(src, out, i, shard_count=2)
    with pytest.raises(ValueError, match='configurations differ'):
        m.merge(out, shard_count=2, expected_rows=2)


def test_incomplete_rerun_invalidates_prior_merge_marker(tmp_path):
    m = module()
    m.atomic_json(tmp_path / 'j10_profile.done', {'old': True})
    with pytest.raises(ValueError):
        m.merge(tmp_path)
    assert not (tmp_path / 'j10_profile.done').exists()


def test_annotation_fence():
    with pytest.raises(ValueError, match='GENCODE'):
        module().summarize([row(annotation='GENCODE v19')])
