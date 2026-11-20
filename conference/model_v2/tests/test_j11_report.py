"""Synthetic-only report contract; no campaign or wet-lab reads."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest


def module():
    path = Path(__file__).parents[1] / 'j11_report.py'
    assert path.exists(), 'J11 generator must exist'
    spec = importlib.util.spec_from_file_location('j11', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def campaign(tmp_path):
    paths = {}

    def artifact(role, data, suffix='.json'):
        path = tmp_path / (role + suffix)
        if suffix == '.parquet':
            pq.write_table(pa.Table.from_pylist(data), path)
        elif suffix == '.json':
            path.write_text(json.dumps(data))
        else:
            path.write_text(data)
        paths[role] = str(path)
        return path

    artifact('jobs', '\n'.join(json.dumps(x) for x in [
        dict(job='J2', kind='phase1_result', slurm_id='12', rows=10),
        dict(job='J3', slurm_id='13', rows=3),
        dict(job='J10', slurm_id='14', state='v1 stale'),
    ]), '.jsonl')
    artifact('j6_gate', dict(rows=3, bad=0, expected=1, missing=[]))
    artifact('j3_gate', dict(rows=3, expected_rows=3, rows_match=True, bad=0,
                           channel_counts={'duplex': {'ok': 2, 'failed': 1}}))
    artifact('j4_gate', dict(complete=1, emitted=1, gaps=0, deep=True))
    counts = [dict(aso_id=aso, scope=scope, mismatch_count=mm,
                   dimension='all', category='all', energy_threshold=e,
                   accessibility_threshold=.1, sites=mm + 1)
              for aso in ['F18', 'F14', 'F10'] for scope in ['other_sites', 'SMN2_on_target']
              for mm in range(4) for e in [-5., -10.]]
    artifact('j10_profile', dict(counts=counts, completeness='passed', input_rows=60,
                               expected_rows=60, config={'implementation_sha256': 'a' * 64}))
    rank = artifact('j9_ranks', [dict(section='rank', baseline='M1', feature_name='length',
        variant_id='v1', aso_id=aso, rank=i + 1, evidence_class='derived', status='ok')
        for i, aso in enumerate(['F18', 'F14', 'F10'])], '.parquet')
    artifact('j9_manifest', dict(outputs=[dict(name=rank.name, rows=3,
        sha256=hashlib.sha256(rank.read_bytes()).hexdigest())]))
    ref = artifact('j4_reference', [dict(locus='SMN2_exon7', site_kind='donor',
        site_position=70076574, tissue=i, reference_score=.1 * (i + 1), status='ok')
        for i in range(4)], '.parquet')
    artifact('j4_done', dict(provenance=dict(model_sha256='a' * 64,
        weights_sha256='b' * 64, weight_files={'w': 'c' * 64}), references=[dict(
        path=str(ref), rows=4, sha256=hashlib.sha256(ref.read_bytes()).hexdigest())],
        plan={'windows': [{'locus': 'SMN2_exon7', 'sites': [dict(kind='donor', position=70076574)]}]}))
    artifact('j7_done', dict(rows=3, sha256='d' * 64, provenance={'common': {'rbp': {'status': 'BLOCKED'}}}))
    for name in ['j3', 'j4', 'j6', 'j7', 'j9', 'j10']:
        artifact('source_' + name, '# synthetic source\n', '.py')
    mod = module()
    lock = mod.capture(paths, git_head='0' * 40)
    return mod, paths, lock, tmp_path / 'report'


def test_missing_artifact_raises(campaign):
    mod, paths, lock, out = campaign
    Path(paths['j3_gate']).unlink()
    with pytest.raises(FileNotFoundError):
        mod.generate(lock, out)
    assert not out.exists()


def test_hash_drift_raises_before_publishing(campaign):
    mod, paths, lock, out = campaign
    Path(paths['j9_ranks']).write_bytes(b'drift')
    with pytest.raises(ValueError, match='hash drift'):
        mod.generate(lock, out)
    assert not out.exists()


@pytest.mark.parametrize('claim', ['This is prospective validation.', 'clinical evidence',
                                  'AUPRC is excellent', 'AUC = 0.9', 'p-value = 0.01'])
def test_forbidden_wording(claim):
    with pytest.raises(ValueError, match='forbidden'):
        module().guard_wording(claim)


def test_required_negative_boundary_allowed():
    mod = module()
    mod.guard_wording(mod.COMPARISON_BOUNDARY)


def test_reject_wetlab_before_opening(tmp_path):
    path = tmp_path / 'data/wetlab/outcomes.xlsx'
    with pytest.raises(ValueError, match='forbidden input'):
        module().capture({'outcome': str(path)}, git_head='0' * 40)


def test_large_read_requires_slurm(tmp_path):
    path = tmp_path / 'large.json'
    with path.open('wb') as stream:
        stream.truncate(3 * 1024 * 1024)
    with pytest.raises(ValueError, match='Slurm'):
        module().capture({'large': str(path)}, git_head='0' * 40)


def test_threshold_repetitions_are_not_added(campaign):
    mod, paths, lock, out = campaign
    counts = mod.site_counts(json.loads(Path(paths['j10_profile']).read_text()))
    assert counts['F18'] == [2, 4, 6, 8]
    profile = json.loads(Path(paths['j10_profile']).read_text())
    profile['counts'][1]['sites'] = 999
    with pytest.raises(ValueError, match='inconsistent'):
        mod.site_counts(profile)


def test_package_is_deterministic_and_complete(campaign):
    mod, paths, lock, out = campaign
    mod.generate(lock, out)
    first = {p.name: p.read_bytes() for p in out.iterdir()}
    mod.generate(lock, out)
    assert first == {p.name: p.read_bytes() for p in out.iterdir()}
    assert len(list(out.glob('*.png'))) == 4
    assert len(list(out.glob('MODEL_CARD_*.md'))) == 4
    text = (out / 'REPORT.md').read_text()
    for boundary in ['RBP join', 'J5 embeddings', 'no verified negatives', '3UP8',
                     'chemistry unknown', 'GTF-confirmed', 'exon 3', 'stale', 'failed']:
        assert boundary in text
    mod.guard_wording(text)
    manifest = json.loads((out / 'MANIFEST.json').read_text())
    for item in manifest['files']:
        assert hashlib.sha256(Path(item['path']).read_bytes()).hexdigest() == item['sha256']
    prediction = json.loads((out / 'PREDICTION_MANIFEST.json').read_text())
    assert prediction['j4_recorded_provenance']['weights_sha256'] == 'b' * 64
    assert prediction['snapshot_sha256'] == mod.digest(prediction['state'])


def test_embedded_j9_hash_is_checked(campaign):
    mod, paths, lock, out = campaign
    path = Path(paths['j9_manifest'])
    value = json.loads(path.read_text())
    value['outputs'][0]['sha256'] = 'f' * 64
    path.write_text(json.dumps(value))
    lock = mod.capture(paths, git_head='0' * 40)
    with pytest.raises(ValueError, match='J9'):
        mod.generate(lock, out)


def test_manifest_rejects_generated_figure_drift(campaign):
    mod, paths, lock, out = campaign
    mod.generate(lock, out)
    (out / 'j10_sites.png').write_bytes(b'corrupt figure')
    with pytest.raises(ValueError, match='hash drift'):
        mod.verify_manifest(out / 'MANIFEST.json')


def test_runtime_version_drift_raises(campaign):
    mod, paths, lock, out = campaign
    lock['runtime_versions']['matplotlib'] = 'wrong-version'
    with pytest.raises(ValueError, match='runtime version drift'):
        mod.generate(lock, out)


def test_failed_reference_is_not_plotted_as_zero(campaign):
    mod, paths, lock, out = campaign
    ref = Path(paths['j4_reference'])
    rows = pq.read_table(ref).to_pylist()
    rows[0]['reference_score'] = None
    rows[0]['status'] = 'failed'
    pq.write_table(pa.Table.from_pylist(rows), ref)
    done = Path(paths['j4_done'])
    marker = json.loads(done.read_text())
    marker['references'][0]['sha256'] = hashlib.sha256(ref.read_bytes()).hexdigest()
    done.write_text(json.dumps(marker))
    lock = mod.capture(paths, git_head='0' * 40)
    with pytest.raises(ValueError, match='unavailable/failed'):
        mod.generate(lock, out)


def test_slurm_ids_include_named_rebuild_stages():
    mod = module()
    records = [dict(slurm_id='29266568/29266569'),
               dict(sentinel='29280507', index='29280508',
                    benchmark='29280509 (afterok index, top-12 density)'),
               dict(slurm_id='mv2-j10-offtarget 256 tasks + merge 29254991', rows=11002701)]
    assert mod.slurm_ids(records) == ['29254991', '29266568', '29266569',
                                      '29280507', '29280508', '29280509']

@pytest.mark.parametrize('text,family,expected', [
    ('', 'J7', 'PROVISIONAL / UNREVIEWED'),
    ('J7: DO-NOT-ADVANCE. J9: ADVANCE.', 'J7', 'PROVISIONAL / UNREVIEWED'),
    ('J7: DO-NOT-ADVANCE. J9: ADVANCE.', 'J9', 'REVIEWED / ADVANCE'),
    ('J7: ADVANCE. J7: DO-NOT-ADVANCE.', 'J7', 'PROVISIONAL / UNREVIEWED'),
    ('J9 v2: ADVANCE', 'J9', 'PROVISIONAL / UNREVIEWED'),
])
def test_version_specific_review_label(text, family, expected):
    assert module().review_label(text, family) == expected


def test_current_support_and_missing_review(campaign):
    mod, paths, lock, out = campaign
    mod.generate(lock, out)
    report = (out / 'REPORT.md').read_text()
    assert 'PROVISIONAL / UNREVIEWED' in report
    assert 'M6 undefined' in report
    assert 'Historical stale J10 v1' not in report
    assert 'PROVISIONAL / UNREVIEWED' in (out / 'MODEL_CARD_J9.md').read_text()


def test_slurm_plural_field():
    assert module().slurm_ids([{'slurm_ids': '29383227/29383235'}]) == ['29383227', '29383235']


def test_bounded_production_metadata_only(tmp_path):
    path = tmp_path / 'j7_context_v3/shards/j7_context.done'
    path.parent.mkdir(parents=True)
    with path.open('wb') as f:
        f.truncate(11 * 1024**2)
    module().bounded([path])
    with path.open('wb') as f:
        f.truncate(17 * 1024**2)
    with pytest.raises(ValueError, match='Slurm'):
        module().bounded([path])


def test_j9_uses_snapshot_reader_and_rejects_generation_drift(campaign, monkeypatch):
    import sys
    import types
    mod, paths, lock, out = campaign
    manifest = json.loads(Path(paths['j9_manifest']).read_text())
    calls = []
    def reader(root):
        calls.append(root)
        return dict(manifest=manifest, ranks=pa.table({'rank': [1]}), features=pa.table({'x': [1]}))
    monkeypatch.setitem(sys.modules, 'conference.model_v2.j9_smn2_prediction',
                        types.SimpleNamespace(read_dataset=reader))
    result = mod.read_j9({k: Path(v) for k, v in paths.items()})
    assert result['features'].num_rows == 1
    assert calls == [Path(paths['j9_manifest']).parent]
    manifest['outputs'] = []
    with pytest.raises(ValueError, match='J9 generation drift'):
        mod.read_j9({k: Path(v) for k, v in paths.items()})
