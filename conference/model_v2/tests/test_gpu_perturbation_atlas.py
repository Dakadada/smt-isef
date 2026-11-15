"""CPU infrastructure tests; StubScorer outputs are NEVER model evidence."""
import importlib.util
import json
import os
from dataclasses import replace
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
import pytest

MODULE = "conference.model_v2.gpu_perturbation_atlas"


def test_engine_exists():
    assert importlib.util.find_spec(MODULE) is not None, "atlas engine is missing"


@pytest.fixture
def atlas():
    return __import__(MODULE, fromlist=["*"])


class StubScorer:
    """Position-sensitive deterministic test double, not Pangolin evidence."""
    context = 2
    provenance = {
        "model": "StubScorer_test_only", "model_sha256": "a" * 64,
        "weights_sha256": "b" * 64, "source_sha256": "c" * 64,
        "model_evidence": False,
    }

    def __init__(self, fail=False, oom_limit=100):
        self.fail, self.oom_limit, self.calls = fail, oom_limit, []

    def predict(self, sequences):
        self.calls.append(list(sequences))
        if len(sequences) > self.oom_limit:
            raise RuntimeError("CUDA out of memory")
        if self.fail:
            raise RuntimeError("test scoring failure")
        return np.array([[[sum((i + 1) * (b != "N") for i, b in enumerate(s)) / 1000
                           + j / 100 + t / 100
                           for j in range(len(s) - 4)] for t in range(4)] for s in sequences])

    def clear_oom(self):
        pass


def tiny(atlas, strand="+"):
    from conference.model_v2.trans_intervention import InMemoryGenome
    return (InMemoryGenome({"chrT": "AACCGGTTAACCGGTTAACCGGTT"}),
            atlas.Window("toy", "chrT", strand, 12, 2,
                         (atlas.Site(11, "donor"), atlas.Site(13, "acceptor"))),
            atlas.Grid((2, 3), 10, 14, 2))


@pytest.mark.parametrize("strand,expected,index,masked", [
    ("+", "TAACCGGTT", 1, "TANNCGGTT"),
    ("-", "AACCGGTTA", 3, "AACCGNNTA"),
])
def test_coordinates_mask_and_strand(atlas, strand, expected, index, masked):
    genome, window, _ = tiny(atlas, strand)
    assert atlas.half_open(8, 16) == (7, 16)
    assert atlas.inclusive(7, 16) == (8, 16)
    assert atlas.oriented_sequence(genome, window, 2) == expected
    assert atlas.output_index(window, 11) == index
    assert atlas.oriented_sequence(genome, window, 2, (10, 11)) == masked


def test_fai_boundary_conversion_once(atlas, tmp_path):
    from conference.model_v2.trans_intervention import FaiFastaReader
    p = tmp_path / "tiny.fa"
    p.write_text(">chrT\nACGT\nTGCA\n")
    Path(str(p) + ".fai").write_text("chrT\t8\t6\t4\t5\n")
    reader = FaiFastaReader(str(p))
    assert reader.fetch("chrT", 4, 5) == "TT"
    assert reader.fetch("chrT", 1, 1) == "A"
    with pytest.raises(ValueError):
        atlas.half_open(0, 2)


def test_grid_and_shards_no_loss_or_overlap(atlas):
    _, window, grid = tiny(atlas)
    tasks = list(atlas.tasks([window], grid))
    assert [(x[1], x[2]) for x in tasks] == [(10, 2), (12, 2), (10, 3), (12, 3)]
    shards = [list(atlas.tasks([window], grid, i, 3)) for i in range(3)]
    assert sorted(sum(shards, [])) == sorted(tasks)
    assert len(set(sum(shards, []))) == 4


def test_batched_oom_backoff_never_drops_and_failures_are_null(atlas):
    results = atlas.score_batch(StubScorer(oom_limit=1), ["AACCGGTTA"] * 5)
    assert len(results) == 5 and all(r[0] is not None for r in results)
    failed = atlas.score_batch(StubScorer(fail=True), ["AACCGGTTA"] * 3)
    assert len(failed) == 3 and all(r[0] is None and "failure" in r[1] for r in failed)


def test_reference_cache_deltas_metadata_and_idempotence(atlas, tmp_path):
    genome, window, grid = tiny(atlas)
    scorer = StubScorer()
    done = atlas.run_atlas([window], grid, genome, scorer, tmp_path, batch_size=2)
    assert done["status"] == "complete" and done["failures"] == 0
    rows = pq.read_table(tmp_path / done["ablation"]["path"]).to_pylist()
    assert len(rows) == 32
    assert all(r["evidence_class"] == "sequence_ablation_proxy" and not r["model_evidence"] for r in rows)
    assert all(r["raw_model_delta"] == pytest.approx(r["masked_score"] - r["reference_score"]) for r in rows)
    assert rows[0]["raw_model_delta"] == pytest.approx(-.007)
    metadata = pq.read_metadata(tmp_path / done["ablation"]["path"]).metadata
    assert json.loads(metadata[b"provenance"])["weights_sha256"] == "b" * 64
    scorer.calls.clear()
    assert atlas.run_atlas([window], grid, genome, scorer, tmp_path, batch_size=2) == done
    assert not scorer.calls
    atlas.run_atlas([window], atlas.Grid((2,), 10, 14, 2), genome, scorer, tmp_path / "other", reference_dir=tmp_path / "reference_scores")
    assert all("N" in s for batch in scorer.calls for s in batch)
    (tmp_path / done["ablation"]["path"]).write_bytes(b"broken")
    with pytest.raises(ValueError, match="integrity"):
        atlas.run_atlas([window], grid, genome, scorer, tmp_path, batch_size=2)


@pytest.mark.parametrize("change", ["weights", "sequence", "strand"])
def test_cache_identity_invalidates(atlas, tmp_path, change):
    from conference.model_v2.trans_intervention import InMemoryGenome
    genome, window, grid = tiny(atlas)
    scorer = StubScorer()
    atlas.run_atlas([window], grid, genome, scorer, tmp_path / "first", reference_dir=tmp_path / "cache")
    scorer.calls.clear()
    if change == "weights":
        scorer.provenance = dict(scorer.provenance, weights_sha256="d" * 64)
    elif change == "sequence":
        genome = InMemoryGenome({"chrT": "T" * 24})
    else:
        window = replace(window, strand="-")
    atlas.run_atlas([window], grid, genome, scorer, tmp_path / "second", reference_dir=tmp_path / "cache")
    assert any("N" not in s for batch in scorer.calls for s in batch)


def test_failed_rows_do_not_claim_success(atlas, tmp_path):
    genome, window, grid = tiny(atlas)
    done = atlas.run_atlas([window], grid, genome, StubScorer(fail=True), tmp_path)
    rows = pq.read_table(tmp_path / done["ablation"]["path"]).to_pylist()
    assert len(rows) == 32 and done["failures"] == 32
    assert done["status"] == "complete_with_failures"
    assert all(r["raw_model_delta"] is None and r["status"] == "failed" and r["error"] for r in rows)


def test_checkpoint_resume_after_interruption(atlas, tmp_path):
    genome, window, grid = tiny(atlas)
    class Interrupted(StubScorer):
        def predict(self, sequences):
            if len(self.calls) == 2:
                raise KeyboardInterrupt
            return super().predict(sequences)
    with pytest.raises(KeyboardInterrupt):
        atlas.run_atlas([window], grid, genome, Interrupted(), tmp_path, batch_size=1)
    scorer = StubScorer()
    done = atlas.run_atlas([window], grid, genome, scorer, tmp_path, batch_size=1)
    assert sum(map(len, scorer.calls)) == 3
    rows = pq.read_table(tmp_path / done["ablation"]["path"]).to_pylist()
    assert len(rows) == len({r["row_id"] for r in rows}) == 32


def test_cpu_gpu_tolerance_hook(atlas):
    if not os.environ.get("SLURM_JOB_ID"):
        pytest.skip("real Pangolin CPU/GPU tolerance requires a Slurm GPU allocation; no login-node inference")
    torch = pytest.importorskip("torch", reason="torch unavailable")
    if not torch.cuda.is_available():
        pytest.skip("CUDA unavailable in this allocation")
    report = atlas.cpu_gpu_tolerance(["ACGT" * 2525 + "A"], atol=1e-5, rtol=1e-4)
    assert report["passed"], report


def test_done_rejects_forged_success_counts(atlas, tmp_path):
    genome, window, grid = tiny(atlas)
    atlas.run_atlas([window], grid, genome, StubScorer(fail=True), tmp_path)
    marker = tmp_path / 'shard_0.done'
    data = json.loads(marker.read_text())
    data['failures'], data['status'] = 0, 'complete'
    marker.write_text(json.dumps(data))
    with pytest.raises(ValueError, match='integrity'):
        atlas.run_atlas([window], grid, genome, StubScorer(), tmp_path)


def test_resume_rejects_duplicate_checkpoint_parts(atlas, tmp_path):
    genome, window, grid = tiny(atlas)
    atlas.run_atlas([window], grid, genome, StubScorer(), tmp_path, batch_size=1)
    (tmp_path / 'shard_0.done').unlink()
    checkpoint = tmp_path / 'shard_0.checkpoint.json'
    state = json.loads(checkpoint.read_text())
    state['parts'][1] = state['parts'][0]
    checkpoint.write_text(json.dumps(state))
    with pytest.raises(ValueError, match='integrity'):
        atlas.run_atlas([window], grid, genome, StubScorer(), tmp_path, batch_size=1)


def test_empty_shard_has_typed_empty_parquet(atlas, tmp_path):
    genome, window, grid = tiny(atlas)
    empty = next(i for i in range(20) if not list(atlas.tasks([window], grid, i, 20)))
    done = atlas.run_atlas([window], grid, genome, StubScorer(), tmp_path, shard=empty, shards=20)
    table = pq.read_table(done['ablation']['path'])
    assert len(table) == done['expected_rows'] == 0
    assert str(table.schema.field('raw_model_delta').type) == 'double'


def test_nan_output_is_failure_not_zero(atlas):
    class Nonfinite(StubScorer):
        def predict(self, seqs):
            values = super().predict(seqs)
            values[0, 0, 0] = np.nan
            return values
    results = atlas.score_batch(Nonfinite(), ['AACCGGTTA'] * 2)
    assert results[0][0] is None and 'nonfinite' in results[0][1]
    assert results[1][0] is not None


def test_real_scorer_refuses_login_node(atlas, monkeypatch):
    monkeypatch.delenv('SLURM_JOB_ID', raising=False)
    with pytest.raises(RuntimeError, match='Slurm'):
        atlas.PangolinScorer()


def test_sigterm_flushes_resume_without_requeue(atlas, tmp_path):
    import subprocess
    import sys
    # Only test-double inference in child. Actual SIGTERM exercises the flush hook.
    child = '''
import os, signal
from pathlib import Path
from conference.model_v2.tests.test_gpu_perturbation_atlas import StubScorer, tiny
from conference.model_v2 import gpu_perturbation_atlas as a
class Terminated(StubScorer):
    def predict(self, seqs):
        if len(self.calls) == 2:
            os.kill(os.getpid(), signal.SIGTERM)
        return super().predict(seqs)
genome, window, grid = tiny(a)
a.run_atlas([window], grid, genome, Terminated(), Path(os.environ['ATLAS_TEST_OUT']), batch_size=1)
'''
    env = dict(os.environ, ATLAS_TEST_OUT=str(tmp_path), OPENBLAS_NUM_THREADS='1',
               OMP_NUM_THREADS='1', MKL_NUM_THREADS='1', NUMEXPR_NUM_THREADS='1')
    result = subprocess.run([sys.executable, '-c', child], env=env, capture_output=True, text=True, timeout=20)
    assert result.returncode == 0, result.stderr
    state = json.loads((tmp_path / 'shard_0.checkpoint.json').read_text())
    assert state['next'] == 1 and not (tmp_path / 'shard_0.done').exists()
    genome, window, grid = tiny(atlas)
    done = atlas.run_atlas([window], grid, genome, StubScorer(), tmp_path, batch_size=1)
    assert done['expected_rows'] == 32 and done['failures'] == 0


@pytest.mark.parametrize('phase', ['ready', 'flush_failure', 'importing', 'readiness_race'])
def test_sentinel_shell_is_single_requeue_owner(tmp_path, phase):
    import signal
    import subprocess
    import sys
    import time
    job = Path(__file__).resolve().parents[1] / 'jobs/j4_sentinel.sbatch'
    fake_env = tmp_path / 'env'
    fake_bin = fake_env / 'bin'
    fake_bin.mkdir(parents=True)
    (fake_bin / 'python').write_text('''#!/bin/bash
trap 'if [[ -n "${J4_READY_FILE:-}" ]]; then : > "$J4_READY_FILE"; fi; echo flushed > "$TEST_DIR/flushed"; exit "$TEST_FLUSH_EXIT"' TERM
if [[ "$TEST_PHASE" != readiness_race && -n "${J4_READY_FILE:-}" ]]; then : > "$J4_READY_FILE"; fi
: > "$TEST_DIR/started"
while true; do sleep 0.1; done
''')
    if phase == 'importing':
        # Execute the real module, holding its first heavyweight import before
        # Python has installed any checkpoint handler or written any output.
        hook = tmp_path / 'hook'
        hook.mkdir()
        (hook / 'sitecustomize.py').write_text('''
import importlib.abc, os, pathlib, sys, time
class DelayNumpy(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'numpy':
            pathlib.Path(os.environ['TEST_DIR'], 'started').touch()
            time.sleep(30)
sys.meta_path.insert(0, DelayNumpy())
''')
        import shlex
        (fake_bin / 'python').write_text(
            '#!/bin/bash\nexport PYTHONPATH=' + shlex.quote(str(hook)) + ':"$PYTHONPATH"\n'
            'exec ' + shlex.quote(sys.executable) + ' "$@"\n')
    (fake_bin / 'python').chmod(0o755)
    (fake_bin / 'scontrol').write_text('''#!/bin/bash
if [[ "$TEST_PHASE" != importing ]]; then test -f "$TEST_DIR/flushed" || exit 1; fi
echo "$*" >> "$TEST_DIR/requeues"
''')
    (fake_bin / 'scontrol').chmod(0o755)
    startup = tmp_path / 'startup'
    startup.write_text('module() { :; }; export -f module\n')
    env = dict(os.environ, BASH_ENV=str(startup), PATH=str(fake_bin) + ':' + os.environ['PATH'],
               J4_GPU_ENV=str(fake_env), J4_OUT_DIR=str(tmp_path / 'out'), TEST_DIR=str(tmp_path), SLURM_JOB_ID='123',
               TEST_PHASE=phase, TEST_FLUSH_EXIT='1' if phase in ('flush_failure', 'readiness_race') else '0')
    child = subprocess.Popen(['bash', str(job)], env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        deadline = time.monotonic() + 5
        while not (tmp_path / 'started').exists() and child.poll() is None and time.monotonic() < deadline:
            time.sleep(.05)
        assert (tmp_path / 'started').exists()
        child.send_signal(signal.SIGTERM)
        _, err = child.communicate(timeout=5)
        if phase in ('flush_failure', 'readiness_race'):
            assert child.returncode == 1, err
            assert not (tmp_path / 'requeues').exists()
            assert b'checkpoint flush failed' in err
        else:
            assert child.returncode == 0, err
            assert (tmp_path / 'requeues').read_text().splitlines() == ['requeue 123']
        if phase == 'importing':
            assert not (tmp_path / 'flushed').exists()
            assert not list((tmp_path / 'out').glob('*.checkpoint.json'))
    finally:
        if child.poll() is None:
            child.kill()
            child.wait()


def test_adapter_matches_existing_channels_and_zero_encodes_n(atlas):
    from contextlib import nullcontext
    class StubTensor:
        def __init__(self, array):
            self.array = array
        def to(self, device):
            return self
        def cpu(self):
            return self
        def numpy(self):
            return self.array
        def __getitem__(self, index):
            return StubTensor(self.array[index])
    class StubTorch:
        from_numpy = staticmethod(StubTensor)
        inference_mode = staticmethod(nullcontext)
    captured = []
    class StubModel:
        def __init__(self, index):
            self.index = index
        def to(self, device):
            return self
        def cpu(self):
            return self
        def __call__(self, value):
            captured.append(value.array)
            size = value.array.shape[2] - 4
            output = np.empty((len(value.array), 12, size))
            for channel in range(12):
                output[:, channel, :] = self.index + channel
            return StubTensor(output)
    scorer = atlas.PangolinScorer.__new__(atlas.PangolinScorer)
    scorer.torch, scorer.device, scorer.context = StubTorch(), 'cpu', 2
    scorer.models = [StubModel(i) for i in range(12)]
    result = scorer.predict(['ACGTNACGT', 'NNNNNNNNN'])
    assert result.shape == (2, 4, 5)
    np.testing.assert_allclose(result[0, :, 0], [2, 8, 14, 20])
    np.testing.assert_array_equal(captured[0][0, :, :5], [
        [1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0], [0, 0, 0, 1, 0]])
    assert not captured[0][1].any()


def test_reference_failure_is_not_model_evidence(atlas, tmp_path):
    genome, window, grid = tiny(atlas)
    scorer = StubScorer(fail=True)
    # Exercise real-scorer provenance semantics without claiming a model run.
    scorer.provenance = dict(scorer.provenance, model_evidence=True)
    done = atlas.run_atlas([window], grid, genome, scorer, tmp_path)
    rows = pq.read_table(done['ablation']['path']).to_pylist()
    assert all(not r['model_evidence'] and r['evidence_category'] == 'not_available' for r in rows)


def test_grid_config_invalid_intervals_fail_closed(atlas):
    for args in [((0,), 10, 14, 1), ((2, 2), 10, 14, 1), ((6,), 10, 14, 1), ((2,), 10, 14, 0)]:
        with pytest.raises(ValueError):
            atlas.Grid(*args)
    _, window, grid = tiny(atlas)
    with pytest.raises(ValueError):
        list(atlas.tasks([window], grid, 1, 1))
    with pytest.raises(ValueError):
        atlas.Window('bad', 'chrT', '-', 12, 2, (atlas.Site(15, 'donor'),))


def test_atomic_publish_keeps_previous_parquet_on_write_error(atlas, tmp_path, monkeypatch):
    target = tmp_path / 'prior.parquet'
    original = b'previous committed file'
    target.write_bytes(original)
    def broken_write(table, path):
        Path(path).write_bytes(b'incomplete replacement')
        raise OSError('simulated disk failure')
    monkeypatch.setattr(atlas.pq, 'write_table', broken_write)
    with pytest.raises(OSError, match='disk failure'):
        atlas.atomic_parquet(target, [], StubScorer.provenance)
    assert target.read_bytes() == original
    assert sorted(p.name for p in tmp_path.iterdir()) == ['prior.parquet']


def test_minus_strand_deltas_and_all_shards_equal_unsharded(atlas, tmp_path):
    genome, window, grid = tiny(atlas, '-')
    done = atlas.run_atlas([window], grid, genome, StubScorer(), tmp_path / 'full')
    rows = pq.read_table(done['ablation']['path']).to_pylist()
    assert rows[0]['raw_model_delta'] == pytest.approx(-.013)
    sharded = []
    for shard in range(3):
        done = atlas.run_atlas([window], grid, genome, StubScorer(), tmp_path / 'sharded', shard=shard, shards=3)
        sharded.extend(pq.read_table(done['ablation']['path']).to_pylist())
    assert sorted(rows, key=lambda r: r['row_id']) == sorted(sharded, key=lambda r: r['row_id'])


def test_resume_preserves_failed_reference_snapshot(atlas, tmp_path, monkeypatch):
    genome, window, grid = tiny(atlas)
    original = atlas.atomic_parquet
    calls = 0
    def interrupted(path, rows, provenance):
        nonlocal calls
        calls += 1
        if calls == 3:
            raise KeyboardInterrupt
        return original(path, rows, provenance)
    monkeypatch.setattr(atlas, 'atomic_parquet', interrupted)
    with pytest.raises(KeyboardInterrupt):
        atlas.run_atlas([window], grid, genome, StubScorer(fail=True), tmp_path, batch_size=1)
    monkeypatch.setattr(atlas, 'atomic_parquet', original)
    scorer = StubScorer()
    done = atlas.run_atlas([window], grid, genome, scorer, tmp_path, batch_size=1)
    assert not scorer.calls
    assert done['failures'] == done['expected_rows'] == 32


def test_resume_rejects_changed_batch_plan(atlas, tmp_path):
    genome, window, grid = tiny(atlas)
    atlas.run_atlas([window], grid, genome, StubScorer(), tmp_path, batch_size=2)
    with pytest.raises(ValueError, match='mismatch'):
        atlas.run_atlas([window], grid, genome, StubScorer(), tmp_path, batch_size=1)


def test_scorer_source_hashes_with_symlinked_repo_root(atlas, tmp_path, monkeypatch):
    import sys
    from types import ModuleType
    from unittest.mock import MagicMock
    from feng import source_auth
    root = Path(atlas.__file__).resolve().parents[2]
    alias = tmp_path / 'repo-link'
    alias.symlink_to(root, target_is_directory=True)
    package = ModuleType('pangolin')
    model = ModuleType('pangolin.model')
    model.__file__ = str(tmp_path / 'model.py')
    Path(model.__file__).write_text('# model fixture; no inference\n')
    package.model = model
    torch = MagicMock(__version__='test')
    monkeypatch.setitem(sys.modules, 'torch', torch)
    monkeypatch.setitem(sys.modules, 'pangolin', package)
    monkeypatch.setitem(sys.modules, 'pangolin.model', model)
    monkeypatch.setattr(source_auth, 'verify_pangolin_package', lambda: {})
    monkeypatch.setattr(atlas.gs, 'PANGOLIN_WEIGHT_FILENAMES', ())
    monkeypatch.setenv('SLURM_JOB_ID', 'test-only-no-inference')
    monkeypatch.setattr(atlas, '__file__', str(alias / 'conference/model_v2/gpu_perturbation_atlas.py'))
    linked = atlas.PangolinScorer(device='cpu').provenance
    monkeypatch.setattr(atlas, '__file__', str(root / 'conference/model_v2/gpu_perturbation_atlas.py'))
    canonical = atlas.PangolinScorer(device='cpu').provenance
    assert linked['source_files'] == canonical['source_files']
    assert linked['source_sha256'] == canonical['source_sha256']
    assert 'feng/src/feng/titration/gpu_scores.py' in linked['source_files']


@pytest.fixture
def indexed_reference(tmp_path):
    # Index-only fixture for assembly validation, never model evidence.
    lengths = [248956422, 242193529, 198295559, 190214555, 181538259,
               170805979, 159345973, 145138636, 138394717, 133797422,
               135086622, 133275309, 114364328, 107043718, 101991189,
               90338345, 83257441, 80373285, 58617616, 64444167,
               46709983, 50818468, 156040895, 57227415, 16569]
    names = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
    fasta = tmp_path / 'reference.fa'
    fasta.write_text('>chr5\nACGT\n')
    Path(str(fasta) + '.fai').write_text(''.join(
        f'{name}\t{length}\t6\t4\t5\n' for name, length in zip(names, lengths)))
    return fasta


@pytest.mark.parametrize('smoke', [False, True])
@pytest.mark.parametrize('override', [False, True])
@pytest.mark.parametrize('defect', ['grch37', 'missing', 'duplicate', 'malformed'])
def test_cli_rejects_wrong_reference_before_model(atlas, indexed_reference, tmp_path, monkeypatch, smoke, override, defect):
    index = Path(str(indexed_reference) + '.fai')
    original = index.read_text()
    if defect == 'grch37':
        index.write_text(original.replace('181538259', '180915260'))
    elif defect == 'missing':
        index.write_text(''.join(line for line in original.splitlines(True) if not line.startswith('chr5\t')))
    elif defect == 'duplicate':
        index.write_text(original + 'chr5\t181538259\t6\t4\t5\n')
    else:
        index.write_text(original + 'malformed\n')
    def unexpected_model():
        pytest.fail('model constructed before reference authentication')
    monkeypatch.setattr(atlas, 'PangolinScorer', unexpected_model)
    config = tmp_path / 'config.json'
    from dataclasses import asdict
    windows, grid = atlas.smoke_config()
    config.write_text(json.dumps({'assembly': 'GRCh38', 'annotation': 'GENCODE v50',
                                 'windows': [asdict(w) for w in windows], 'grid': asdict(grid)}))
    mode = ['--smoke'] if smoke else ['--config', str(config)]
    monkeypatch.setattr(atlas.ti, 'DEFAULT_FASTA', str(indexed_reference))
    fasta_args = ['--fasta', str(indexed_reference)] if override else []
    with pytest.raises(ValueError, match='GRCh38'):
        atlas.main(mode + fasta_args + ['--output-dir', str(tmp_path / 'out')])
    assert not (tmp_path / 'out').exists()


def test_cli_reference_identity_bound_to_outputs_and_alias_stable(atlas, indexed_reference, tmp_path, monkeypatch):
    monkeypatch.setattr(atlas, 'PangolinScorer', StubScorer)
    monkeypatch.setattr(atlas.ti.FaiFastaReader, 'fetch', lambda self, chrom, start, end: 'A' * (end - start + 1))
    ready = tmp_path / 'ready'
    monkeypatch.setenv('J4_READY_FILE', str(ready))
    out = tmp_path / 'out'
    assert atlas.main(['--smoke', '--fasta', str(indexed_reference), '--output-dir', str(out)]) == 0
    assert ready.is_file()
    monkeypatch.delenv('J4_READY_FILE')  # Each shell invocation supplies a fresh marker.
    done = atlas.load_manifest(out / 'shard_0.done')
    identity = done['provenance']['reference']
    assert identity['assembly'] == 'GRCh38'
    assert identity['primary_contig_lengths']['chr5'] == 181538259
    assert identity['fai_sha256'] == atlas.file_sha(str(indexed_reference) + '.fai')
    assert identity['fasta_path'] == str(indexed_reference.resolve())
    metadata = pq.read_metadata(done['ablation']['path']).metadata
    assert json.loads(metadata[b'provenance'])['reference'] == identity
    assert done['plan']['provenance']['reference'] == identity
    alias = tmp_path / 'alias'
    alias.symlink_to(tmp_path, target_is_directory=True)
    assert atlas.main(['--smoke', '--fasta', str(alias / indexed_reference.name), '--output-dir', str(out)]) == 0
    assert atlas.load_manifest(out / 'shard_0.done')['run_id'] == done['run_id']
    # Same sequences with a changed valid index must not reuse old manifests.
    with Path(str(indexed_reference) + '.fai').open('a') as stream:
        stream.write('extra_contig\t4\t6\t4\t5\n')
    with pytest.raises(ValueError, match='mismatch'):
        atlas.main(['--smoke', '--fasta', str(indexed_reference), '--output-dir', str(out)])
