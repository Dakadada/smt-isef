"""Layer E: distinct derived/model_assumption/proxy quantities, never a risk score.

Prespecified sensitivity grid: duplex <= {-5,-10,-15,-20} kcal/mol;
whole-footprint unpaired probability >= {.01,.1,.5}. No outcome fitting.
Fractions use channel-evaluable sites, with missing/failed counts beside them.
J6 gene-specific site IDs collapse exact repeated records, never different genes.
Hash ownership confines all copies of a site ID to one partial; conflicts fail.
Only J3's recorded accessibility window is used, with no extrapolated variants.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from contextlib import contextmanager
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

import pyarrow as pa
import pyarrow.parquet as pq

VERSION = '1.0.0'
ASOS = ('F18', 'F14', 'F10')
ENERGIES = (-5., -10., -15., -20.)
ACCESSIBILITIES = (.01, .1, .5)
CHANNELS = ('complementarity', 'duplex', 'accessibility')
METRICS = ('sites', 'duplex_evaluated', 'duplex_pass', 'accessibility_evaluated',
           'accessibility_pass', 'joint_evaluated', 'joint_pass')
SCOPES = ('other_sites', 'SMN2_on_target', 'SMN1_paralog_anchor', 'SMN2_other_sites', 'SMN1_other_sites')
LABELS = dict(sequence='derived', duplex='model_assumption (unmodified RNA, chemistry not modelled)',
              accessibility='proxy')


def canonical(obj):
    return json.dumps(obj, sort_keys=True, separators=(',', ':'), allow_nan=False)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def partition(site_id, count):
    return int(hashlib.sha256(site_id.encode()).hexdigest(), 16) % count


def sync_dir(path):
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def atomic_json(path, obj):
    path = Path(path)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=path.name, suffix='.tmp')
    try:
        with os.fdopen(fd, 'w') as f:
            f.write(canonical(obj) + '\n')
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp, path)
        sync_dir(path.parent)
    finally:
        Path(tmp).unlink(missing_ok=True)


@contextmanager
def locked(base):
    with base.with_suffix('.lock').open('a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        yield


def config(energies=ENERGIES, accessibilities=ACCESSIBILITIES, inclusive=True):
    if not energies or not accessibilities or any(not math.isfinite(x) for x in energies):
        raise ValueError('finite nonempty threshold grids required')
    if any(not math.isfinite(x) or not 0 <= x <= 1 for x in accessibilities):
        raise ValueError('accessibility thresholds must be probabilities')
    return dict(version=VERSION, implementation_sha256=sha256(__file__), energies=list(energies),
                accessibilities=list(accessibilities), inclusive=inclusive, labels=LABELS)


def scope(row):
    name = row.get('gene_name')
    if name not in ('SMN1', 'SMN2'):
        return 'other_sites'
    start = 70076584 if name == 'SMN2' else 70952004
    length = int(row['aso_id'][1:])
    target = 'CCAGCATTATGAAAGTGA'[:length]
    if (row.get('chrom'), row.get('start1'), row.get('end1'), row.get('strand'),
        row.get('target_sequence', '').upper().replace('U', 'T')) == ('chr5', start, start + length - 1, '+', target):
        return 'SMN2_on_target' if name == 'SMN2' else 'SMN1_paralog_anchor'
    return name + '_other_sites'


def status(row, channel):
    s = row.get(channel + '_status')
    return s if isinstance(s, str) and s else 'missing'


def valid_number(row, field, probability=False):
    value = row.get(field)
    if (not isinstance(value, (float, int)) or isinstance(value, bool) or
        not math.isfinite(value) or (probability and not 0 <= value <= 1)):
        raise ValueError(f'ok channel has invalid {field}')
    return value


def render(counters):
    result = []
    for key in sorted(counters):
        aso, sc, mm, dimension, category, energy, access = key
        r = dict(zip(METRICS, counters[key]))
        r.update(aso_id=aso, scope=sc, mismatch_count=mm, dimension=dimension, category=category,
                 energy_threshold=energy, accessibility_threshold=access)
        r['accessibility_fraction'] = (r['accessibility_pass'] / r['accessibility_evaluated']
                                       if r['accessibility_evaluated'] else None)
        result.append(r)
    return result


def summarize(rows, *, energies=ENERGIES, accessibilities=ACCESSIBILITIES, inclusive=True):
    cfg = config(energies, accessibilities, inclusive)
    counts = defaultdict(lambda: [0] * len(METRICS))
    # Explicit zero bins allow an absent ASO/mismatch to be distinguished from no computation.
    for aso in ASOS:
        for sc in SCOPES:
            for mm in range(4):
                for e in energies:
                    for a in accessibilities:
                        counts[aso, sc, mm, 'all', 'all', e, a]
    seen, smn = {}, []
    input_rows = ignored = duplicates = 0
    for r in rows:
        input_rows += 1
        if r.get('assembly') != 'GRCh38':
            raise ValueError('GRCh38 assembly required')
        if r.get('annotation') != 'GENCODE v50':
            raise ValueError('GENCODE v50 annotation required')
        sid = r.get('site_id')
        if not isinstance(sid, str) or not sid:
            raise ValueError('site_id required')
        # Hash complete records, including annotations and failures. Never select a better duplicate.
        digest = hashlib.sha256(canonical(r).encode()).hexdigest()
        if sid in seen:
            if seen[sid] != digest:
                raise ValueError('conflicting duplicate site_id')
            duplicates += 1
            continue
        seen[sid] = digest
        if r.get('aso_id') not in ASOS:
            ignored += 1
            continue
        mm = r.get('mismatch_count')  # J6 sequence compatibility remains available when J3 fails.
        if type(mm) is not int or mm not in range(4):
            raise ValueError('J6 mismatch_count outside 0..3')
        statuses = {c: status(r, c) for c in CHANNELS}
        if statuses['complementarity'] == 'ok' and r.get('complementarity_mismatch_count') != mm:
            raise ValueError('J3/J6 mismatch disagreement')
        d_ok, a_ok = statuses['duplex'] == 'ok', statuses['accessibility'] == 'ok'
        d = valid_number(r, 'duplex_energy_kcal_mol') if d_ok else None
        a = valid_number(r, 'accessibility_p_unpaired', True) if a_ok else None
        sc = scope(r)
        if sc.startswith('SMN'):
            smn.append(dict(r, j10_scope=sc))
        dimensions = [('all', 'all'), ('region_class', r.get('region_class') or 'not_available'),
                      ('gene_biotype', r.get('gene_biotype') or 'not_available')]
        dimensions += [('channel_status', f'{c}:{statuses[c]}') for c in CHANNELS]
        for e in energies:
            for cutoff in accessibilities:
                dp = d_ok and (d <= e if inclusive else d < e)
                ap = a_ok and (a >= cutoff if inclusive else a > cutoff)
                values = (1, int(d_ok), int(dp), int(a_ok), int(ap), int(d_ok and a_ok), int(dp and ap))
                for dimension, category in dimensions:
                    target = counts[r['aso_id'], sc, mm, dimension, category, e, cutoff]
                    for i, value in enumerate(values):
                        target[i] += value
    return dict(config=cfg, input_rows=input_rows, unique_rows=len(seen), duplicate_rows=duplicates,
                ignored_unique_rows=ignored, profiled_rows=len(seen) - ignored,
                counts=render(counts), smn_rows=sorted(smn, key=lambda r: r['site_id']))


def commit(base, payload, identity):
    """Commit marker last; interrupted writes never become reusable output."""
    base.with_suffix('.done').unlink(missing_ok=True)
    atomic_json(base.with_suffix('.json'), payload)
    fd, tmp = tempfile.mkstemp(dir=base.parent, prefix=base.name, suffix='.tmp')
    os.close(fd)
    try:
        table = pa.Table.from_pylist(payload['counts'])
        table = table.replace_schema_metadata({b'j10_labels': canonical(LABELS).encode()})
        pq.write_table(table, tmp, compression='zstd')
        with open(tmp, 'rb') as f:
            os.fsync(f.fileno())
        os.replace(tmp, base.with_suffix('.parquet'))
        sync_dir(base.parent)
        atomic_json(base.with_suffix('.done'), dict(identity=identity,
            sha256=sha256(base.with_suffix('.parquet')), json_sha256=sha256(base.with_suffix('.json')),
            rows=table.num_rows))
    finally:
        Path(tmp).unlink(missing_ok=True)


def load(base, identity=None):
    """Read one committed generation under the same lock used by partial writers."""
    try:
        with locked(base):
            return _load_locked(base, identity)
    except OSError as exc:
        raise ValueError(f'invalid partial {base}: {exc}') from exc


def _load_locked(base, identity=None):
    """Caller must hold locked(base) across marker, hash, and payload reads."""
    try:
        marker = json.loads(base.with_suffix('.done').read_text())
        if identity is not None and marker['identity'] != identity:
            raise ValueError('incompatible partial identity')
        if (marker['sha256'] != sha256(base.with_suffix('.parquet')) or
            marker['json_sha256'] != sha256(base.with_suffix('.json')) or
            marker['rows'] != pq.read_metadata(base.with_suffix('.parquet')).num_rows):
            raise ValueError('partial checksum/row mismatch')
        payload = json.loads(base.with_suffix('.json').read_text())
        if pq.read_table(base.with_suffix('.parquet')).to_pylist() != payload['counts']:
            raise ValueError('partial JSON/parquet mismatch')
        return payload, marker
    except (OSError, KeyError, TypeError, pa.ArrowException) as exc:
        raise ValueError(f'invalid partial {base}: {exc}') from exc


def source_identity(path, index, count):
    try:
        marker = json.loads(path.with_suffix('.done').read_text())
        digest = sha256(path)
        if (marker['sha256'] != digest or marker['rows'] != pq.read_metadata(path).num_rows or
            marker['shard_index'] != index or marker['shard_count'] != count):
            raise ValueError('J3 source commitment mismatch')
        return dict(path=str(path.resolve()), sha256=digest, rows=marker['rows'],
                    params=marker['params'], reference=marker.get('reference'),
                    inputs_sha256=hashlib.sha256(canonical(marker.get('inputs')).encode()).hexdigest(),
                    shard_index=index, shard_count=count)
    except (OSError, KeyError, TypeError, pa.ArrowException) as exc:
        raise ValueError(f'invalid J3 source {path}') from exc


def partial(input_dir, output_dir, index, *, shard_count=256):
    if not 0 <= index < shard_count:
        raise ValueError('invalid shard index')
    path = Path(input_dir).expanduser() / f'j3_shard{index:04d}_of{shard_count:04d}.parquet'
    out = Path(output_dir).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    base = out / f'j10_partial{index:04d}'
    with locked(base):
        source = source_identity(path, index, shard_count)
        identity = dict(config=config(), source=source)
        try:
            return _load_locked(base, identity)[0]
        except ValueError:
            pass
        base.with_suffix('.done').unlink(missing_ok=True)
        sync_dir(out)
        def stream():
            for batch in pq.ParquetFile(path).iter_batches(batch_size=512):
                for r in batch.to_pylist():
                    if partition(r['site_id'], shard_count) != index:
                        raise ValueError('J3 hash ownership mismatch; cross-shard duplicate possible')
                    yield r
        payload = summarize(stream())
        payload['accessibility_configuration'] = {k: source['params'][k] for k in ('window', 'span', 'flank')}
        payload['accessibility_variants'] = 'one recorded J3 configuration only; no window sensitivity available'
        if source_identity(path, index, shard_count) != source:
            raise ValueError('source changed while summarizing')
        commit(base, payload, identity)
        return payload


def merge(output_dir, *, shard_count=256, expected_rows=11002701):
    if type(shard_count) is not int or shard_count <= 0:
        raise ValueError('shard_count must be a positive integer')
    if type(expected_rows) is not int or expected_rows < 0:
        raise ValueError('expected_rows must be a nonnegative integer')
    out = Path(output_dir).expanduser()
    base = out / 'j10_profile'
    with locked(base):
        base.with_suffix('.done').unlink(missing_ok=True)
        sync_dir(out)
        expected = {f'j10_partial{i:04d}.done' for i in range(shard_count)}
        if {p.name for p in out.glob('j10_partial*.done')} != expected:
            raise ValueError('partial completeness gate: expected every shard exactly once')
        counts = defaultdict(lambda: [0] * len(METRICS))
        totals = dict.fromkeys(('input_rows', 'unique_rows', 'duplicate_rows', 'ignored_unique_rows', 'profiled_rows'), 0)
        identities, smn, params = [], [], None
        accessibility = None
        for i in range(shard_count):
            payload, marker = load(out / f'j10_partial{i:04d}')
            identity = marker['identity']
            source = identity['source']
            if (identity['config'] != config() or source['shard_index'] != i or
                source['shard_count'] != shard_count or payload['input_rows'] != source['rows']):
                raise ValueError('partial configuration/index/row mismatch')
            origin = {k: source[k] for k in ('params', 'reference', 'inputs_sha256')}
            if params is not None and params != origin:
                raise ValueError('J3 configurations differ')
            params = origin
            accessibility = payload['accessibility_configuration']
            identities.append(marker)
            for k in totals:
                totals[k] += payload[k]
            smn.extend(payload['smn_rows'])
            for r in payload['counts']:
                key = tuple(r[k] for k in ('aso_id', 'scope', 'mismatch_count', 'dimension', 'category',
                                          'energy_threshold', 'accessibility_threshold'))
                for j, name in enumerate(METRICS):
                    counts[key][j] += r[name]
        if totals['input_rows'] != expected_rows:
            raise ValueError('input row completeness gate failed')
        result = dict(config=config(), **totals, counts=render(counts),
                      smn_rows=sorted(smn, key=lambda r: r['site_id']), shards=shard_count,
                      expected_rows=expected_rows, completeness='passed', accessibility_configuration=accessibility,
                      accessibility_variants='one recorded J3 configuration only; no window sensitivity available',
                      region_biotype_note='J6 gene-only annotation may have region not_available; biotype absent in J6',
                      fraction_denominator='accessibility_evaluated; None when zero; failures retained as status categories',
                      duplex_scope='minimum anywhere in target+flank; not footprint constrained')
        commit(base, result, dict(config=config(), partial_markers=identities, expected_rows=expected_rows))
        return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=('partial', 'merge'))
    parser.add_argument('--input-dir', default='~/orcd/scratch/feng/model_v2/j3_genome_v1')
    parser.add_argument('--output-dir', default='~/orcd/scratch/feng/model_v2/j10_offtarget_v2')
    parser.add_argument('--shard-index', type=int, default=0)
    parser.add_argument('--shard-count', type=int, default=256)
    parser.add_argument('--expected-rows', type=int, default=11002701)
    args = parser.parse_args()
    if args.shard_count != 256:
        parser.error('production CLI requires exactly 256 shards; synthetic counts use the Python API')
    if args.expected_rows < 0:
        parser.error('expected rows must be nonnegative')
    if args.stage == 'partial':
        result = partial(args.input_dir, args.output_dir, args.shard_index, shard_count=args.shard_count)
    else:
        result = merge(args.output_dir, shard_count=args.shard_count, expected_rows=args.expected_rows)
    print(canonical({k: result[k] for k in ('input_rows', 'unique_rows', 'duplicate_rows', 'profiled_rows')}))


if __name__ == '__main__':
    main()
