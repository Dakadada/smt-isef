"""J11: bounded, outcome-free, hash-locked descriptive report generation.

Capture is explicit and refuses to replace an existing lock. Generate never
refreshes hashes. Model/weight hashes are recorded provenance, not a rehash of
large external weight files. No models, workbooks, or genome shards are loaded.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


COMPARISON_BOUNDARY = (
    'Future independent wet-lab comparison must be made against this snapshot. '
    'Already-seen July gels are NOT prospective validation.'
)
BLOCKED = [
    'M6 undefined: no supervised model or efficacy endpoint is available.',
    'RBP join: BLOCKED pending verified GRCh38 peak bytes, liftover chain hash and mapping QC.',
    'J5 embeddings: BLOCKED / not_available; no embedding output is authenticated here.',
    'Supervised training: BLOCKED; no verified negatives. Responders and matched controls do not supply negatives.',
    '3UP8 unplaced: needs human confirmation; excluded from the three-ASO ranks.',
    'ASO chemistry unknown; unmodified-RNA energies do not establish modified-ASO binding or potency.',
    'SMN coordinates are derived, not GTF-confirmed; paralogs must remain separate.',
    'July exon-7 gels were previously seen; no prospective claim; excluded from fitting, selection and joins.',
    'The exon 3 audit is exploratory only; no training, tuning, selection or exon-7 inference.',
]
EVIDENCE = {
    'observed': 'File existence, recorded run metadata, or supplied measurement; no wet-lab measurements enter this package.',
    'derived': 'Deterministic sequence, geometry, count, rank or reference-model computation; no efficacy label.',
    'proxy': 'Accessibility and N-masking sequence_ablation_proxy; neither occupancy nor a measured splicing change.',
    'model_assumption': 'Unmodified-RNA duplex thermodynamics; chemistry and cellular context are not modeled.',
    'unlabeled': 'No verified outcome class, including matched/scrambled controls; not a negative label.',
    'not_available': 'Missing, failed, blocked or undefined quantities stay explicit; never filled with zero.',
}
SUPPORT = {
    'J2': ('derived', 'Locus-window enumeration and accessibility coverage.', 'No whole-genome completeness or modified-ASO efficacy; stage counts have different universes.'),
    'J3': ('derived / proxy / model_assumption', 'Sequence complementarity, RNA duplex energy and accessibility with explicit channel failures.', 'No potency or measured occupancy; failed rows remain failed.'),
    'J4': ('derived / proxy', 'Reference splice-model outputs and N-mask sequence_ablation_proxy.', 'No measured splicing change or ASO mechanism; tissues are model channels, not replicates.'),
    'J5': ('derived / not_available', 'v4 gate records RBP-context output; authenticated embeddings are absent.', 'No usable RBP join or supervised feature claim.'),
    'J6': ('derived', 'GRCh38 / GENCODE v50 gene-body sequence-site enumeration, F18/F14/F10, mismatch 0–3.', 'No intergenic completeness, placed 3UP8, expression or activity label.'),
    'J7': ('derived / unlabeled / not_available', 'Canonical-transcript geometry and unlabeled matched controls; completion-marker metadata only.', 'RBP join is blocked; isoform sensitivity and negative labels are unavailable.'),
    'J8': ('derived / proxy', 'Genome v2 join manifest records full coverage and explicit missing channels.', 'No validated joint score; J2 exact footprint matching only, with 10/14 nt unavailable.'),
    'J9': ('derived / proxy / model_assumption', 'Descriptive M0–M5 per-feature ordering for three placed ASOs and recorded alternatives.', 'No outcome labels, fitted predictor, efficacy test, or placed 3UP8.'),
    'J10': ('derived / proxy / model_assumption', 'Production v2 site counts by mismatch; sequence counts include SMN anchors.', 'v1 is stale and excluded. v2 review is bounded; no risk or safety claim.'),
}
CAPTIONS = {
    'j10_sites.png': 'Evidence class: derived. Production J10 v2 sites by mismatch, per ASO, summed across disjoint scopes once per threshold-independent bin (including SMN anchors). Limitation: gene-specific site identifiers, not unique genome coordinates; no activity or safety inference; reviewed v2 is descriptive only.',
    'j9_ranks.png': 'Evidence class: derived / proxy / model_assumption, inherited per feature. Descriptive ordinal ranks over three ASOs; each feature/alternative stays separate. Limitation: no verified negatives or outcome labels; ranks are not efficacy estimates; blank cells mean not_available.',
    'j4_reference.png': 'Evidence class: derived. Pangolin reference-score distribution at the SMN2 exon-7 donor (last exonic base) across recorded tissue channels. Limitation: tissue channels are not independent experimental replicates; coordinates are derived, not GTF-confirmed; no measured splicing change.',
    'j2_j3_coverage.png': 'Evidence class: derived. J2 ledger stage rows and J3 gate channel counts. Limitation: distinct stage/universe denominators; missing ledger totals are not_available, not zero; channel failures are retained; no cross-stage coverage fraction is implied.',
}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def safe_input(path):
    path = Path(path).expanduser().absolute()
    for candidate in (path, path.resolve()):
        low = str(candidate).lower()
        if (candidate.suffix.lower() in {'.xlsx', '.xls', '.csv'}
                or '/data/wetlab/' in low or '/wetlab/outputs/' in low
                or '/terminal/' in low or 'sealed' in low or 'phase9' in low):
            raise ValueError(f'forbidden input: {path}')
        if '/wetlab/' in low and candidate.name != 'README.md':
            raise ValueError(f'forbidden input: {path}')
    return path


def bounded(paths):
    """Bound named production metadata to 16 MiB/file, 64 MiB total inline.

    All other inputs retain the 2 MiB/file cap. Genome payloads are not read.
    """
    paths = [safe_input(p) for p in paths]
    sizes = [p.stat().st_size for p in paths]
    metadata = ('j7_context_v3/shards/j7_context.done', 'j8_genome_v2/manifest.json')
    limits = [16 * 1024**2 if any(str(p).endswith(n) for n in metadata)
              else 2 * 1024**2 for p in paths]
    if not os.environ.get('SLURM_JOB_ID') and (any(n > cap for n, cap in zip(sizes, limits))
                                               or sum(sizes) > 64 * 1024**2):
        raise ValueError('large reads require Slurm; use report/generate.sbatch')


def sha256(path):
    h = hashlib.sha256()
    with safe_input(path).open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def capture(paths, *, git_head, repo=None):
    bounded(paths.values())
    entries = [dict(role=role, path=str(safe_input(path)), sha256=sha256(path),
                    bytes=Path(path).stat().st_size) for role, path in sorted(paths.items())]
    return dict(version=1, git_head=git_head, repo=str(repo) if repo else None,
                runtime_versions=runtime_versions(), files=entries)


def runtime_versions():
    return dict(python=sys.version, **{name: importlib.metadata.version(name)
                for name in ('matplotlib', 'numpy', 'pyarrow')})


def verify(lock):
    if lock['runtime_versions'] != runtime_versions():
        raise ValueError('runtime version drift')
    bounded([f['path'] for f in lock['files']])
    roles = set()
    for entry in lock['files']:
        if entry['role'] in roles:
            raise ValueError('duplicate input role')
        roles.add(entry['role'])
        if sha256(entry['path']) != entry['sha256']:
            raise ValueError(f"hash drift: {entry['path']}")
    if lock.get('repo'):
        head = git_head(lock['repo'])
        if head != lock['git_head']:
            raise ValueError('git HEAD drift')
    return {f['role']: Path(f['path']) for f in lock['files']}


def verify_manifest(path):
    manifest = load(path)
    bounded([f['path'] for f in manifest['files']])
    for entry in manifest['files']:
        if sha256(entry['path']) != entry['sha256']:
            raise ValueError(f"hash drift: {entry['path']}")
    return len(manifest['files'])


def git_head(repo):
    return subprocess.check_output(['git', '-C', str(repo), 'rev-parse', 'HEAD'], text=True).strip()


def slurm_ids(records):
    # Some rebuild submissions use named stage fields instead of slurm_id.
    return sorted({match for record in records
                   for key in ('slurm_id', 'slurm_ids', 'sentinel', 'index', 'benchmark')
                   for match in re.findall(r'\b\d{7,}\b', str(record.get(key, '')))})


def review_label(text, family):
    # Only the locked v3 review is passed here; fail closed on ambiguous verdicts.
    verdicts = re.findall(r'\b' + family + r'(?: v3)?\s*:\s*(DO-NOT-ADVANCE|ADVANCE)\b',
                          text.replace('**', ''), re.I)
    return ('REVIEWED / ADVANCE' if verdicts and all(v.upper() == 'ADVANCE' for v in verdicts)
            else 'PROVISIONAL / UNREVIEWED')


def read_j9(paths):
    # Resolve once during capture, then pin the same immutable generation.
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from conference.model_v2.j9_smn2_prediction import read_dataset
    directory = paths['j9_manifest'].parent
    root = directory.parent.parent if directory.parent.name == '.generations' else directory
    dataset = read_dataset(root)
    if dataset['manifest'] != load(paths['j9_manifest']):
        raise ValueError('J9 generation drift')
    return dataset


def guard_wording(text):
    text = text.replace(COMPARISON_BOUNDARY, '')
    if re.search(r'prospective\s+validation|\bclinical\b|\bAU(?:PRC|ROC|C)\b|p[- ]values?', text, re.I):
        raise ValueError('forbidden report claim')


def load(path):
    return json.loads(path.read_text())


def site_counts(profile):
    """Collapse repeated threshold grids only after checking equal site counts."""
    bins = {}
    grids = {}
    for row in profile['counts']:
        if row['dimension'] != 'all' or row['category'] != 'all':
            continue
        key = (row['aso_id'], row['scope'], row['mismatch_count'])
        value = row['sites']
        if type(value) is not int or value < 0:
            raise ValueError('invalid site count')
        if key in bins and bins[key] != value:
            raise ValueError('inconsistent threshold-independent site counts')
        grid = (row['energy_threshold'], row['accessibility_threshold'])
        if grid in grids.setdefault(key, set()):
            raise ValueError('duplicate site-count grid bin')
        grids[key].add(grid)
        bins[key] = value
    if not bins or len({frozenset(g) for g in grids.values()}) != 1:
        raise ValueError('incomplete site-count threshold grid')
    asos = sorted({k[0] for k in bins}, key=lambda a: (-int(a[1:]), a))
    scopes = sorted({k[1] for k in bins})
    # Every ASO/mismatch/scope must be explicit, including actual computed zero bins.
    return {aso: [sum(bins[(aso, scope, mm)] for scope in scopes) for mm in range(4)]
            for aso in asos}


def checked_rows(path, entry, label):
    import pyarrow.parquet as pq
    if sha256(path) != entry['sha256']:
        raise ValueError(f'{label} embedded hash mismatch')
    table = pq.read_table(path)
    if table.num_rows != entry['rows']:
        raise ValueError(f'{label} embedded row count mismatch')
    return table.to_pylist()


def read_data(paths):
    jobs = [json.loads(line) for line in paths['jobs'].read_text().splitlines() if line.strip()]
    data = {role: load(paths[role]) for role in ('j6_gate', 'j3_gate', 'j4_gate',
            'j10_profile', 'j9_manifest', 'j4_done', 'j7_done')}
    data['jobs'] = jobs
    for role in ('j5_gate', 'j8_manifest', 'j10_done'):
        if role in paths:
            data[role] = load(paths[role])
    review = paths.get('review_v3')
    data['reviews'] = {family: review_label(review.read_text() if review else '', family)
                       for family in ('J7', 'J9')}
    entries = [e for e in data['j9_manifest']['outputs'] if e['name'] == paths['j9_ranks'].name]
    if len(entries) != 1:
        raise ValueError('J9 ranks manifest entry required')
    if 'j9_features' in paths:
        dataset = read_j9(paths)
        data['ranks'] = dataset['ranks'].to_pylist()
        data['feature_rows'] = dataset['features'].num_rows
    else:  # Synthetic legacy fixture, never selected by campaign_paths.
        data['ranks'] = checked_rows(paths['j9_ranks'], entries[0], 'J9')
    if 'j10_done' in data and sha256(paths['j10_profile']) != data['j10_done']['json_sha256']:
        raise ValueError('J10 embedded hash mismatch')
    entries = [e for e in data['j4_done']['references']
               if Path(e['path']).resolve() == paths['j4_reference'].resolve()]
    if len(entries) != 1:
        raise ValueError('J4 reference manifest entry required')
    reference = checked_rows(paths['j4_reference'], entries[0], 'J4')
    donors = {(w['locus'], s['position']) for w in data['j4_done']['plan']['windows']
              for s in w['sites'] if s['kind'] == 'donor' and w['locus'] == 'SMN2_exon7'}
    data['reference'] = [r for r in reference if (r['locus'], r['site_position']) in donors
                         and r['site_kind'] == 'donor']
    if not data['reference']:
        raise ValueError('SMN2 exon7 donor reference missing')
    for r in data['reference']:
        score = r['reference_score']
        if r['status'] != 'ok' or score is None or not math.isfinite(score) or not 0 <= score <= 1:
            raise ValueError('SMN2 donor reference unavailable/failed')
    provenance = data['j4_done']['provenance']
    for key in ('model_sha256', 'weights_sha256'):
        if not re.fullmatch('[0-9a-f]{64}', provenance[key]):
            raise ValueError('missing J4 model identity')
    if not provenance['weight_files'] or any(not re.fullmatch('[0-9a-f]{64}', v)
                                            for v in provenance['weight_files'].values()):
        raise ValueError('missing J4 weight identities')
    data['sites'] = site_counts(data['j10_profile'])
    data['coverage'] = [(f"J2 {r.get('kind', 'record')}", r['rows']) for r in jobs
                        if r.get('job') == 'J2' and type(r.get('rows')) is int]
    data['coverage'] += [(f'J3 {channel} {status}', value)
                        for channel, counts in sorted(data['j3_gate']['channel_counts'].items())
                        for status, value in sorted(counts.items()) if type(value) is int]
    return data


def figures(data, output):
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/j11-matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'figure.dpi': 100,
                         'savefig.dpi': 120, 'axes.spines.top': False, 'axes.spines.right': False})

    def save(fig, name):
        fig.tight_layout()
        fig.savefig(output / name, metadata={'Software': 'J11 deterministic descriptive report',
                                            'Description': CAPTIONS[name]})
        plt.close(fig)

    fig, axes = plt.subplots(1, len(data['sites']), figsize=(10, 3), squeeze=False)
    for ax, (aso, values) in zip(axes[0], data['sites'].items()):
        ax.bar(range(4), values, color='#31688e')
        ax.set(title=f'{aso} · production v2', xlabel='Mismatches', ylabel='Sites (symlog)')
        ax.set_xticks(range(4)); ax.set_yscale('symlog', linthresh=1)
        for mm, value in enumerate(values):
            ax.annotate(str(value), (mm, value), xytext=(0, 3), textcoords='offset points', ha='center', fontsize=7)
    save(fig, 'j10_sites.png')

    ranks = [r for r in data['ranks'] if r['section'] == 'rank']
    if not ranks:
        raise ValueError('no J9 ranks available')
    keys = sorted({(r['baseline'], r['feature_name'], r['variant_id']) for r in ranks})
    asos = ['F18', 'F14', 'F10']
    matrix = np.full((len(keys), len(asos)), np.nan)
    seen = set()
    for r in ranks:
        cell = ((r['baseline'], r['feature_name'], r['variant_id']), r['aso_id'])
        if cell in seen:
            raise ValueError('duplicate J9 rank cell')
        seen.add(cell)
        if r['rank'] is not None and r['status'] == 'ok':
            matrix[keys.index(cell[0]), asos.index(r['aso_id'])] = r['rank']
    fig, ax = plt.subplots(figsize=(11, max(3, len(keys) * .28)))
    ax.imshow(np.ma.masked_invalid(matrix), cmap='viridis_r', vmin=1, vmax=3, aspect='auto')
    ax.set_xticks(range(3), asos)
    ax.set_yticks(range(len(keys)), [' / '.join(k) for k in keys], fontsize=7)
    ax.set_title('J9 descriptive per-feature rank · N=3 ASOs')
    for (i, j), value in np.ndenumerate(matrix):
        if math.isfinite(value):
            ax.text(j, i, str(int(value)), ha='center', va='center', color='white' if value > 1.5 else 'black')
    save(fig, 'j9_ranks.png')

    scores = [r['reference_score'] for r in data['reference']]
    fig, axes = plt.subplots(1, 2, figsize=(8, 3))
    axes[0].hist(scores, bins=np.linspace(0, 1, 11), color='#31688e')
    axes[0].set(xlabel='Reference score', ylabel='Tissue-channel count', xlim=(0, 1), title='SMN2 exon-7 donor')
    axes[1].scatter([r['tissue'] for r in data['reference']], scores, color='#31688e')
    axes[1].set(xlabel='Recorded tissue channel', ylabel='Reference score', ylim=(0, 1))
    save(fig, 'j4_reference.png')

    fig, ax = plt.subplots(figsize=(10, 5))
    labels, values = zip(*data['coverage'])
    ax.barh(range(len(values)), values, color=['#b63655' if 'failed' in s else '#31688e' for s in labels])
    ax.set_yticks(range(len(labels)), labels, fontsize=7)
    ax.set_xscale('symlog', linthresh=1)
    ax.set(xlabel='Rows (symlog; distinct universes)', title='J2 ledger / J3 gate coverage')
    for i, value in enumerate(values):
        ax.annotate(str(value), (value, i), xytext=(3, 0), textcoords='offset points', va='center', fontsize=7)
    save(fig, 'j2_j3_coverage.png')


def prediction(lock, data):
    sources = {f['path']: f['sha256'] for f in lock['files'] if f['role'].startswith('source_')}
    for name in ('j3', 'j4', 'j6', 'j7', 'j9', 'j10'):
        if 'source_' + name not in {f['role'] for f in lock['files']}:
            raise ValueError(f'missing predictor source: {name}')
    p = data['j4_done']['provenance']
    identity = {k: p[k] for k in ('model_sha256', 'weights_sha256', 'weight_files')}
    state = dict(git_head=lock['git_head'], source_sha256=sources, j4_recorded_identity=identity,
                 inputs_sha256={f['role']: f['sha256'] for f in lock['files']},
                 input_lock_sha256=digest(lock),
                 predictor_kind='unfitted descriptive baselines; no supervised weights',
                 history='dirty working-tree sources individually hashed; artifacts retain their historical producer identities')
    return dict(state=state, snapshot_sha256=digest(state), j4_recorded_provenance=identity,
                comparison_boundary=COMPARISON_BOUNDARY,
                freeze_scope='Report snapshot only; no preregistration seal, terminal unlock or git freeze tag.',
                weights_verification='Hashes copied from locked J4 provenance; external weight bytes were not reopened.')


def documents(lock, data, output, pred):
    entries = {f['role']: f for f in lock['files']}

    def citation(role):
        e = entries[role]
        return f"`{e['path']}` (SHA256 `{e['sha256']}`)"

    text = ['# J11 descriptive report', '', f"Predictor snapshot SHA256: `{pred['snapshot_sha256']}`.",
            '', COMPARISON_BOUNDARY, '', '## What exists', '',
            'Ledger states are recorded history, not fresh Slurm queries. Gates are locked historical reports; '
            'large output bytes are not revalidated. Unknown totals are not_available. '
            'Current inputs: J5 v4, J8 genome v2, J10 off-target v2, J7/J9 v3. J10 v1 is stale and excluded.', '',
            '| Family | Slurm IDs (all recorded stages) | Rows / units | Gate / status | Evidence class |',
            '|---|---|---|---|---|']
    for family, (evidence, supports, limits) in SUPPORT.items():
        records = [r for r in data['jobs'] if r.get('job') == family or (family == 'J7' and r.get('job') in ('J7_v3', 'J7_J9_v3'))]
        ids = slurm_ids(records)
        ledger_rows = [f"{r.get('kind', 'record')}: {r['rows']}" for r in records if type(r.get('rows')) is int]
        rows = '; '.join(ledger_rows) or 'not_available'
        status = 'ledger only; completion not independently authenticated'
        if family in ('J3', 'J6'):
            gate = data[family.lower() + '_gate']
            rows = str(gate['rows'])
            good = gate.get('bad') == 0 and (family != 'J3' or gate.get('rows_match') is True)
            status = 'historical gate passed' if good else 'gate NOT passed / indeterminate'
        elif family == 'J4':
            gate = data['j4_gate']; rows = f"{gate['complete']} complete / {gate['emitted']} emitted configs; {len(data['reference'])} donor reference rows"
            status = 'historical deep gate passed' if gate.get('deep') is True and gate['gaps'] == 0 and gate['complete'] == gate['emitted'] else 'gate NOT passed / indeterminate'
        elif family == 'J5' and 'j5_gate' in data:
            gate = data['j5_gate']; rows = f"{gate['rows']} rows / {gate['valid_windows']} valid windows"
            status = 'j5_j6_v4 gate passed; embeddings BLOCKED' if gate.get('ok') is True else 'gate NOT passed'
            ids = sorted(set(ids) | {str(i) for i in range(29292355, 29292360)})
        elif family == 'J8' and 'j8_manifest' in data:
            rows = str(data['j8_manifest']['rows'])
            status = 'genome v2 manifest; Slurm 29261102 COMPLETED (campaign record); payloads not rehashed'
        elif family == 'J7':
            rows = str(data['j7_done']['rows']); status = data['reviews']['J7'] + '; completion marker only; RBP join BLOCKED; writer overlap unverified'
        elif family == 'J9':
            feature_entries = [e for e in data['j9_manifest']['outputs'] if e['name'] == 'features.parquet']
            feature_count = str(feature_entries[0]['rows']) if len(feature_entries) == 1 else 'not_available'
            rows = f"{len(data['ranks'])} ranks/stability rows; {feature_count} feature rows (verified when read_dataset is used); 3 placed ASOs"
            status = data['reviews']['J9'] + ('; features and ranks verified via read_dataset()' if 'feature_rows' in data else '; ranks verified; feature bytes not reopened')
            ids = sorted(set(ids) | {'29383229'})
        elif family == 'J10':
            rows = f"{data['j10_profile']['input_rows']} input rows; {len(data['j10_profile']['counts'])} aggregate rows"
            status = 'v2 REVIEWED / ADVANCE (bounded review); completeness: ' + data['j10_profile']['completeness']
        text.append(f"| {family} | {', '.join(ids) or 'not_available'} | {rows} | {status} | {evidence} |")
    text += ['', 'J2 phase-2 total: not_available unless explicitly row-counted in the ledger; no total is inferred from task completion.',
             '', '### Artifact citations (exact bytes)', '']
    text += [f"- **{role}**: {citation(role)}" for role in sorted(entries)]
    text += ['', 'J7 completion-marker output SHA256 is a producer-recorded claim, not reverified large parquet bytes: '
             f"`{data['j7_done']['sha256']}`. No large J7 parquet is cited as byte-verified.",
             '', '## What outputs support and do not support', '']
    for family, (evidence, supports, limits) in SUPPORT.items():
        text += [f'- **{family}** ({evidence}): supports {supports} Does not support: {limits}']
    text += ['', 'Review limits: J7 shard writers are not locked through merge consumption; production writer overlap was not independently verified. J9 FASTA identity uses path/size/mtime and FAI hash, not a full FASTA hash. Review verdicts are not scientific certification.', '', 'J3 channel status counts (failed computations remain explicit):', '',
             '```json', json.dumps(data['j3_gate']['channel_counts'], sort_keys=True, indent=2), '```',
             '', '## Evidence boundary', '']
    text += [f'- **{key}**: {meaning}' for key, meaning in EVIDENCE.items()]
    text += ['', '## Blocked list', ''] + [f'- {item}' for item in BLOCKED]
    text += ['', '## Descriptive figures', '']
    for name, caption in CAPTIONS.items():
        text += [f'![{name}]({name})', '', caption, '']
    text += ['## J9 baseline rank table', '',
             '| Baseline | Feature | Alternative | ASO | Rank | Evidence class | Status |',
             '|---|---|---|---|---|---|---|']
    for r in data['ranks']:
        if r['section'] == 'rank':
            values = [r.get(k) for k in ('baseline', 'feature_name', 'variant_id', 'aso_id', 'rank', 'evidence_class', 'status')]
            text.append('| ' + ' | '.join('not_available' if v is None else str(v) for v in values) + ' |')
    text += ['', '## Regeneration', '',
             'Use the locked generator command in REPRODUCE.md. Missing artifacts, source drift, or git HEAD drift '
             'fail before publication. External weights and genome data are not read. This snapshot does not '
             'certify that all historical producers match current source versions. Wet-lab audit documentation '
             'is cited by path/hash only; no wet-lab numeric data are reproduced.']
    write_text(output / 'REPORT.md', '\n'.join(text) + '\n')
    input_roles = {'J3': ['j3_gate', 'source_j3'], 'J4': ['j4_gate', 'j4_done', 'j4_reference', 'source_j4'],
                   'J9': ['j9_manifest', 'j9_ranks', 'source_j9'], 'J10': ['j10_profile', 'source_j10']}
    for family in ('J3', 'J4', 'J9', 'J10'):
        evidence, supports, limits = SUPPORT[family]
        card = [f'# {family} model card', '', f'Inputs: {"; ".join(citation(r) for r in input_roles[family])}.',
                '', f'Evidence class: {evidence}.', '', f'Output: {supports}',
                '', 'Labels: unlabeled; no outcome labels or verified negatives are read, joined, fitted or selected.',
                '', 'Assumptions: GRCh38 / GENCODE v50; sequence-derived SMN coordinates remain unconfirmed against GTF. '
                'Chemistry unknown; any duplex term assumes unmodified RNA. Canonical context and model tissue channels '
                'do not establish cell-specific ASO response.', '', f'Limitations: {limits}', '', COMPARISON_BOUNDARY]
        if family == 'J9':
            card += ['', 'Review: ' + data['reviews']['J9'] + '; see locked J7J9_V3_REVIEW_RESULT.md. FASTA identity is bounded, not a content hash.']
        if family == 'J4':
            card += ['', 'J4 recorded model identity:', '```json', json.dumps(pred['j4_recorded_provenance'], sort_keys=True, indent=2), '```']
        write_text(output / f'MODEL_CARD_{family}.md', '\n'.join(card) + '\n')
    repo = lock.get('repo') or str(Path(__file__).resolve().parents[2])
    # The portable command is relative to the recorded checkout, not the staging directory.
    script = 'conference/model_v2/j11_report.py'
    report = 'conference/model_v2/report'
    py = '/home/dashii/orcd/pool/smt-isef-env/bin/python'
    env = 'PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1'
    write_text(output / 'REPRODUCE.md', f'''# Reproduce J11

Run from `{repo}` with the runtime versions pinned in INPUT_LOCK.json.
The lock is an explicit report snapshot. Capture refuses to replace it; generate
and verify reject drift. Missing input files raise; hashes are never refreshed.
Initial capture (already performed for this package):

```bash
{env} {py} {script} capture
```

Regenerate and verify:

```bash
{env} {py} {script} generate
{env} {py} {script} verify
```

Inline reads are capped at 2 MiB/file, except the named J7 v3 and J8 v2 metadata
markers (16 MiB/file), and 64 MiB total. Run with a 6 GiB memory ceiling.
For larger inputs, use `sbatch {report}/generate.sbatch`; it is supplied and has
not been submitted by J11. It requeues on preemption and restarts from the
unchanged input lock; a final manifest is the publication checkpoint. No model
weights or genome shards are loaded. Byte-identical PNGs require the pinned
Python, NumPy, PyArrow and Matplotlib versions and the same font environment.
The manifest excludes its own hash; its files include all input citations,
figures, model cards, reproduction instructions and the Slurm template.
Wet-lab files are never opened. Audit documentation is hashed without parsing.
J9 uses read_dataset() and requires the live pointer to still select the locked
generation; pointer drift fails closed and requires an explicit new snapshot.
''')
    write_text(output / 'generate.sbatch', f'''#!/bin/bash
#SBATCH --job-name=mv2-j11-report
#SBATCH --partition=mit_preemptable
#SBATCH --requeue
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output={repo}/{report}/j11-%j.out
set -euo pipefail
cd {repo}
export PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=feng/src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ARROW_NUM_THREADS=1
{py} {script} generate
{py} {script} verify
''')


def write_text(path, text):
    guard_wording(text)
    path.write_text(text)


def write_json(path, value):
    path.write_text(json.dumps(value, sort_keys=True, indent=2, allow_nan=False) + '\n')


def generate(lock, output):
    paths = verify(lock)
    data = read_data(paths)
    pred = prediction(lock, data)
    output = Path(output).absolute()
    # Stage everything first; missing/corrupt sources cannot leave a partial package.
    output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='.j11-', dir=output) as tmp:
        stage = Path(tmp)
        figures(data, stage)
        documents(lock, data, stage, pred)
        write_json(stage / 'INPUT_LOCK.json', lock)
        write_json(stage / 'PREDICTION_MANIFEST.json', pred)
        verify(lock)  # Detect sources changed while rendering.
        files = list(lock['files']) + [dict(role='report_output', path=str(output / p.name),
            sha256=sha256(p), bytes=p.stat().st_size) for p in sorted(stage.iterdir())]
        write_json(stage / 'MANIFEST.json', dict(version=1, files=files,
            self_exclusion='MANIFEST.json excludes its own hash to avoid recursion.',
            generator_versions=runtime_versions()))
        for p in sorted(stage.iterdir(), key=lambda p: (p.name == 'MANIFEST.json', p.name)):
            os.replace(p, output / p.name)
    return pred['snapshot_sha256']


def campaign_paths(repo, scratch):
    base = repo / 'conference/model_v2'
    paths = {role: base / path for role, path in {
        'final_audit': 'ops/FINAL_LAUNCH_AUDIT.md', 'jobs': 'ops/JOBS.jsonl', 'usage_ledger': 'ops/USAGE_LEDGER.jsonl',
        'next_action': 'ops/NEXT_ACTION.json', 'wetlab_manifest': 'ops/WETLAB_MANIFEST.json',
        'scientific_lock': 'SCIENTIFIC_LOCK.md', 'j4_gate': 'ops/briefs/j4_gate.json',
        'source_j3': 'j3_scoring.py', 'source_j4': 'gpu_perturbation_atlas.py',
        'source_j4_configs': 'j4_build_configs.py', 'source_j6': 'j6_site_universe.py',
        'source_j7': 'j7_context_features.py', 'source_j9': 'j9_smn2_prediction.py',
        'source_j10': 'j10_offtarget_profile.py', 'source_j10_robustness': 'j10_robustness.py',
        'source_report': 'j11_report.py', 'wetlab_readme_metadata_only': 'wetlab/README.md',
    }.items()}
    for p in sorted((base / 'ops/briefs').glob('*_RESULT.md')):
        if p.name.startswith('WETLAB'):
            continue
        paths['brief_' + p.stem] = p
    paths.update({role: scratch / path for role, path in {
        'j6_gate': 'j6_genome_v1/gate_report.json', 'j3_gate': 'j3_genome_v1/gate_report.json',
        'j10_profile': 'j10_offtarget_v2/j10_profile.json', 'j10_done': 'j10_offtarget_v2/j10_profile.done',
        'j5_gate': 'j5_gate_v4.json', 'j8_manifest': 'j8_genome_v2/manifest.json',
        'j9_manifest': 'j9_smn2_v3/manifest.json',
        'j9_features': 'j9_smn2_v3/features.parquet', 'j9_readme': 'j9_smn2_v3/README.md', 'j9_ranks': 'j9_smn2_v3/ranks.parquet',
        'j7_done': 'j7_context_v3/shards/j7_context.done',
        'j4_done': 'j4_atlas_v1/SMN2_exon7_issn1_fine/shard_0.done',
    }.items()})
    review = base / 'ops/briefs/J7J9_V3_REVIEW_RESULT.md'
    if review.exists():
        paths['review_v3'] = review
    for role in ('j9_manifest', 'j9_features', 'j9_ranks', 'j9_readme'):
        paths[role] = paths[role].resolve()
    bounded(paths.values())
    refs = load(paths['j4_done'])['references']
    if len(refs) != 1:
        raise ValueError('expected one curated SMN2 reference artifact')
    paths['j4_reference'] = safe_input(refs[0]['path'])
    return paths


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('command', choices=['capture', 'generate', 'verify'])
    parser.add_argument('--repo', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--scratch', type=lambda p: Path(p).expanduser(), default=Path.home() / 'orcd/scratch/feng/model_v2')
    parser.add_argument('--out', type=Path, default=Path(__file__).parent / 'report')
    args = parser.parse_args()
    lock_path = args.out / 'INPUT_LOCK.json'
    if args.command == 'capture':
        if lock_path.exists():
            raise ValueError('refusing to replace existing input lock; choose a fresh output directory')
        lock = capture(campaign_paths(args.repo, args.scratch), git_head=git_head(args.repo), repo=args.repo)
        args.out.mkdir(parents=True, exist_ok=True)
        write_json(lock_path, lock)
        print(f'captured {len(lock["files"])} files / {sum(f["bytes"] for f in lock["files"])} bytes')
    elif args.command == 'generate':
        if (args.out / 'MANIFEST.json').exists():
            verify_manifest(args.out / 'MANIFEST.json')
        print(generate(load(lock_path), args.out))
    else:
        lock = load(lock_path)
        verify(lock)
        print(f'verified {verify_manifest(args.out / "MANIFEST.json")} files')


if __name__ == '__main__':
    main()
