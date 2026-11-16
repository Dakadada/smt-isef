"""J5: observed RBP context, with an accession-bound GRCh38 assembly fence.

RBP knockdown and local steric ASO inhibition are NOT interchangeable.
KD events always carry evidence_class='knockdown_not_steric_block', with a
separate measurement_evidence_class='observed'. No potency or outcome joins.

Design: disk-backed SQLite RTree partitions by RBP/cell/chromosome; files commit
individually, query Parquet streams in bounded batches, hits have explicit caps.
BED is 0-based half-open; public queries/J6/J8 footprints are 1-based inclusive.
Overlap uses the query strand only for +/- peaks. Signed nearest distance is
peak minus query (nearest BASES, hence adjacent=1), reversed on minus-strand
peaks; unstranded peaks use genomic orientation. Within each RBP/cell, files
form a union for presence, but signal measurements retain accession/replicate
identity. No signal averaging, replicate consensus, or cross-cell pooling.
A missing local file does not establish that an assay was never performed.
Refused/partial evidence cannot produce a False presence flag after a miss.

Metadata input: JSON mapping accession -> original ENCODE File object. Its
accession, @id and md5sum must match the local file. Assembly may instead come
from an explicit assembly=/genome= header; conflicting assertions are refused.
No inferred assemblies and no implicit liftover. Optional sample1_role,
sample2_role, sample_order_evidence document rMATS sample order; absent these,
raw sample1-minus-sample2 is observed, but KD direction remains unavailable.

Decoys are same-gene, length/GC matched coordinate shifts, labelled proxy;
they are descriptive background intervals, never verified negative labels.
Default fixed shifts are finite candidates, not a statistical sampling universe.
No p-values. Missing reference/gene bounds yields unavailable controls.

CLI: --encode-root ... --out ... [--metadata encode_file_objects.json]
     [--queries shard.parquet ... | --query-list ordered_paths.txt | --sentinel]
     [--fasta GRCh38.fa --genes gene_bounds.tsv]
Gene TSV: gene_id,chrom,start1,end1,strand,assembly,annotation. It is indexed on
disk; optional FASTA uses J6's pinned primary-contig length fence. A sentinel
uses three known F18/F14/F10 footprints plus the exon-7 anchor +/-500 kb window.
Large production ingestion belongs in Slurm. Requeue restarts at committed files
and output shards. Changed inputs/code require a fresh output directory.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sqlite3
import tarfile
import tempfile

import pyarrow as pa
import pyarrow.parquet as pq
from conference.model_v2.evaluation_core import check_source_not_forbidden

VERSION = '1.0.0'
RBPS = tuple('FUS HNRNPA1 HNRNPC HNRNPL HNRNPM HNRNPU KHSRP PCBP1 PCBP2 PTBP1 RBFOX2 RBM22 SF3B4 SRSF1 SRSF7 SRSF9 TARDBP TIA1 TIAL1 U2AF1 U2AF2'.split())
CELLS = ('HepG2', 'K562')
OUTCOME_COLUMNS = {'label', 'label_state', 'outcome', 'psi', 'delta_psi', 'dpsi', 'efficacy', 'potency'}
QUERY_COLUMNS = ('site_id', 'row_id', 'aso_id', 'chrom', 'start1', 'end1', 'footprint_start',
                 'footprint_end', 'strand', 'gene_id', 'assembly', 'annotation', 'gene_start1', 'gene_end1')


def safe_path(path):
    path = Path(check_source_not_forbidden(str(path))).expanduser().resolve()
    if any(p.lower() in {'wetlab', 'terminal', 'psi25'} for p in path.parts):
        raise ValueError('outcome source forbidden')
    return path


def dumps(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def digests(path):
    sha, md5 = hashlib.sha256(), hashlib.md5()
    with safe_path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            sha.update(block); md5.update(block)
    return sha.hexdigest(), md5.hexdigest()


def atomic_json(path, value):
    path = safe_path(path)
    fd, temp = tempfile.mkstemp(dir=path.parent, suffix='.tmp')
    try:
        with os.fdopen(fd, 'w') as f:
            f.write(dumps(value)+'\n'); f.flush(); os.fsync(f.fileno())
        os.replace(temp, path)
    finally:
        if os.path.exists(temp): os.unlink(temp)


def fence(row):
    if row.get('assembly') != 'GRCh38' or row.get('annotation', 'GENCODE v50') != 'GENCODE v50':
        raise ValueError('assembly/annotation fence requires GRCh38 / GENCODE v50')
    if (row.get('strand') not in ('+', '-') or not isinstance(row.get('chrom'), str)
            or not row['chrom'] or type(row.get('start1')) is not int
            or type(row.get('end1')) is not int or not 1 <= row['start1'] <= row['end1']):
        raise ValueError('invalid 1-based inclusive query')


def normalize_query(row):
    row = dict(row)
    if 'start1' not in row: row['start1'] = row.get('footprint_start')
    if 'end1' not in row: row['end1'] = row.get('footprint_end')
    for key, alias in (('start1', 'footprint_start'), ('end1', 'footprint_end')):
        if alias in row and row[alias] != row[key]: raise ValueError('conflicting footprint coordinates')
    if not (row.get('site_id') or row.get('row_id')): raise ValueError('query identity required')
    fence(row)
    return row


def _number(value, *, bed_missing=False):
    if value in (None, '', 'NA', 'nan', 'NaN') or (bed_missing and value == '-1'): return None
    result = float(value)
    if not math.isfinite(result): raise ValueError('nonfinite measurement')
    return result


def _text(path):
    return gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()


def _headers(path, kind):
    if kind != 'peaks': return []
    headers = []
    with _text(path) as f:
        for line in f:
            if line.startswith(('#', 'track', 'browser')):
                if len(headers) >= 128 or len(line) > 65536: raise ValueError('oversize BED header')
                headers.append(line.strip())
            elif line.strip(): break
    return headers


def verify_assembly(accession, headers, metadata, md5):
    assertions = []
    for line in headers:
        assertions.extend(re.findall(r'(?:assembly|genome)\s*[=:]\s*["\']?([A-Za-z0-9_.-]+)', line, re.I))
    if metadata:
        if (metadata.get('accession') != accession or metadata.get('@id') != f'/files/{accession}/'
                or metadata.get('md5sum') != md5):
            return None, None, 'metadata accession/@id/md5 mismatch'
        if metadata.get('assembly'): assertions.append(metadata['assembly'])
    canonical = ['GRCh38' if a == 'hg38' else a for a in assertions]
    assembly = canonical[0] if canonical and len(set(canonical)) == 1 else None
    source = 'header_and_encode' if headers and metadata else 'encode_accession_md5' if metadata else 'header'
    if assembly != 'GRCh38':
        return assembly, source if assertions else None, 'unverified, conflicting, or non-GRCh38 assembly; no implicit liftover'
    return assembly, source, None


def _peaks(path):
    with _text(path) as f:
        for line in f:
            if not line.strip() or line.startswith(('#', 'track', 'browser')): continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) < 3: raise ValueError('BED needs at least 3 columns')
            chrom, lo, hi = fields[0], int(fields[1]), int(fields[2])
            strand = fields[5] if len(fields) > 5 else '.'
            if not chrom or not 0 <= lo < hi or strand not in ('+', '-', '.'):
                raise ValueError('invalid BED coordinates/strand')
            yield [(chrom, lo+1, hi, strand)], dict(
                name=fields[3] if len(fields)>3 else None,
                score=dict(value=_number(fields[4], bed_missing=True) if len(fields)>4 else None, field='BED_column5', evidence_class='observed' if len(fields)>4 and _number(fields[4], bed_missing=True) is not None else 'not_available'),
                signal=dict(value=_number(fields[6], bed_missing=True) if len(fields)>6 else None, field='BED_column7_as_given', evidence_class='observed' if len(fields)>6 and _number(fields[6], bed_missing=True) is not None else 'not_available'))


def _kd(path, metadata):
    # No extractall: never materialize or trust tar member paths.
    with tarfile.open(path, 'r|gz') as tar:
        for member in tar:
            if not member.isfile() or not member.name.endswith('.txt') or '.MATS.' not in member.name: continue
            with tar.extractfile(member) as text:
                # TextIOWrapper calls seekable(), absent on Python 3.11's tar _Stream.
                reader = csv.DictReader((line.decode('utf-8') for line in text), delimiter='\t')
                fields = set(reader.fieldnames or ())
                if not {'chr', 'strand', 'IncLevelDifference'} <= fields:
                    raise ValueError('unsupported rMATS schema')
                pairs = [(a,b) for a,b in [('exonStart_0base','exonEnd'), ('longExonStart_0base','longExonEnd'),
                         ('riExonStart_0base','riExonEnd'), ('1stExonStart_0base','1stExonEnd'), ('2ndExonStart_0base','2ndExonEnd')] if {a,b} <= fields]
                if not pairs: raise ValueError('unsupported rMATS event coordinates')
                for event in reader:
                    blocks = [(event['chr'], int(event[a])+1, int(event[b]), event['strand']) for a,b in pairs]
                    if any(not 1 <= lo <= hi or strand not in ('+', '-') for _,lo,hi,strand in blocks):
                        raise ValueError('invalid rMATS event coordinates')
                    raw = _number(event['IncLevelDifference'])
                    if raw is not None and not -1 <= raw <= 1: raise ValueError('rMATS difference outside [-1,1]')
                    roles = (metadata.get('sample1_role'), metadata.get('sample2_role'))
                    sign = (1 if roles == ('knockdown', 'control') else -1 if roles == ('control', 'knockdown') else None) if metadata.get('sample_order_evidence') else None
                    delta = sign*raw if sign is not None and raw is not None else None
                    yield blocks, dict(event_id=event.get('ID'), member=member.name,
                        event_type=Path(member.name).name.split('.')[0], raw_sample1_minus_sample2=raw,
                        event_blocks=[dict(chrom=c,start1=a,end1=b,strand=s) for c,a,b,s in blocks],
                        kd_minus_control=delta, magnitude=abs(delta) if delta is not None else None,
                        direction=('increased_inclusion' if delta>0 else 'decreased_inclusion' if delta<0 else 'unchanged') if delta is not None else 'not_available',
                        direction_evidence_class='observed' if delta is not None else 'not_available',
                        units='inclusion_fraction',
                        sample_order_evidence=metadata.get('sample_order_evidence'),
                        evidence_class='knockdown_not_steric_block', measurement_evidence_class='observed' if raw is not None else 'not_available')


class ContextIndex:
    """Disk-backed, one writer; the SQLite page cache is bounded to 16 MiB."""
    def __init__(self, path, *, readonly=False):
        self.path = safe_path(path)
        self.code_sha = digests(__file__)[0]
        self._cached_contract = None
        self._availability_cache = {}
        self._cache_version = None
        if not readonly: self.path.parent.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(self.path.as_uri()+'?mode=ro',uri=True) if readonly else sqlite3.connect(self.path)
        self.db.row_factory = sqlite3.Row
        self.db.executescript('PRAGMA cache_size=-16384; PRAGMA temp_store=FILE;')
        if not readonly: self.db.executescript('''
            PRAGMA cache_size=-16384; PRAGMA temp_store=FILE; PRAGMA foreign_keys=ON;
            CREATE TABLE IF NOT EXISTS config(version TEXT NOT NULL, code_sha TEXT NOT NULL);
            CREATE TABLE IF NOT EXISTS files(id INTEGER PRIMARY KEY, path TEXT UNIQUE, rbp TEXT, cell TEXT, kind TEXT, manifest TEXT);
            CREATE TABLE IF NOT EXISTS chroms(id INTEGER PRIMARY KEY, name TEXT UNIQUE);
            CREATE TABLE IF NOT EXISTS intervals(id INTEGER PRIMARY KEY, file_id INTEGER REFERENCES files(id),
                group_id INT, chrom_id INT, lo INT, hi INT, strand TEXT, kind TEXT, source_key TEXT, payload TEXT);
            CREATE INDEX IF NOT EXISTS near_lo ON intervals(group_id,chrom_id,kind,strand,lo);
            CREATE INDEX IF NOT EXISTS near_hi ON intervals(group_id,chrom_id,kind,strand,hi);
            CREATE VIRTUAL TABLE IF NOT EXISTS regions USING rtree(id,lo,hi,g0,g1,c0,c1);
            CREATE TRIGGER IF NOT EXISTS remove_region AFTER DELETE ON intervals BEGIN DELETE FROM regions WHERE id=old.id; END;
            CREATE TABLE IF NOT EXISTS genes(gene_id TEXT PRIMARY KEY, chrom TEXT, lo INT, hi INT, strand TEXT);
            CREATE TABLE IF NOT EXISTS sources(name TEXT PRIMARY KEY, sha256 TEXT);
        ''')
        config = self.db.execute('select * from config').fetchone()
        if config is None:
            if readonly: raise ValueError('missing index contract')
            self.db.execute('insert into config values (?,?)', (VERSION, self.code_sha)); self.db.commit()
        elif tuple(config) != (VERSION, self.code_sha):
            self.db.close(); raise ValueError('index contract changed; use a fresh output directory')

    def __enter__(self): return self
    def __exit__(self, *args): self.db.close()

    def add_file(self, path, rbp, cell, *, kind='peaks', metadata=None):
        if rbp not in RBPS or cell not in CELLS or kind not in ('peaks','kd'): raise ValueError('invalid assay identity')
        path, metadata = safe_path(path), metadata or {}
        stat = path.stat()
        sha, md5 = digests(path)
        accession = path.name.split('.')[0]
        headers = _headers(path, kind)
        assembly, source, reason = verify_assembly(accession, headers, metadata, md5)
        callers = [value for line in headers for value in re.findall(r'(?:peak_caller|peak_calling_method)\s*[=:]\s*([^\s]+)',line)]
        manifest = dict(path=str(path), accession=accession, rbp=rbp, cell_line=cell, kind=kind,
                        assembly=assembly, assembly_source=source, sha256=sha, md5sum=md5,
                        row_count=0, interval_count=0, peak_calling_method=metadata.get('peak_calling_method') or (callers[0] if len(set(callers))==1 else None),
                        strand_counts={'+':0,'-':0,'.':0},
                        replicate_info=metadata.get('biological_replicates'), technical_replicates=metadata.get('technical_replicates'),
                        analysis_step_version=metadata.get('analysis_step_version'), output_type=metadata.get('output_type'),
                        file_format_type=metadata.get('file_format_type'), metadata=metadata, headers=headers,
                        evidence_class='observed', status='refused' if reason else 'indexed', refusal_reason=reason)
        old = self.db.execute('select id,manifest from files where path=?', (str(path),)).fetchone()
        if old:
            previous = json.loads(old['manifest'])
            if (previous['rbp'],previous['cell_line'],previous['kind']) != (rbp,cell,kind): raise ValueError('file assay identity changed')
            if previous['sha256'] == sha and previous['metadata'] == metadata and previous['headers'] == headers:
                return previous
        gid = RBPS.index(rbp)*2 + CELLS.index(cell)
        replicates = set()
        with self.db:
            if old:
                fid = old['id']
                self.db.execute('delete from intervals where file_id=?', (fid,))
                self.db.execute('update files set rbp=?,cell=?,kind=? where id=?', (rbp,cell,kind,fid))
            else:
                fid = self.db.execute('insert into files(path,rbp,cell,kind) values (?,?,?,?)', (str(path),rbp,cell,kind)).lastrowid
            self.db.execute('SAVEPOINT content')
            try:
                rows = _peaks(path) if kind == 'peaks' else _kd(path, metadata)
                for blocks, payload in rows:
                    manifest['row_count'] += 1
                    for _,_,_,strand in blocks: manifest['strand_counts'][strand] += 1
                    if kind == 'peaks': replicates.update(re.findall(r'rep\d+', payload.get('name') or ''))
                    if reason: continue
                    payload.update(accession=accession, source_row=manifest['row_count'])
                    source_key=str(path)+'\n'+str(manifest['row_count']).zfill(12)
                    for chrom, lo, hi, strand in blocks:
                        self.db.execute('insert or ignore into chroms(name) values (?)', (chrom,))
                        cid = self.db.execute('select id from chroms where name=?', (chrom,)).fetchone()[0]
                        ident = self.db.execute('insert into intervals(file_id,group_id,chrom_id,lo,hi,strand,kind,source_key,payload) values (?,?,?,?,?,?,?,?,?)',
                            (fid,gid,cid,lo,hi,strand,kind,source_key,dumps(payload))).lastrowid
                        self.db.execute('insert into regions values (?,?,?,?,?,?,?)', (ident,lo,hi,gid,gid,cid,cid))
                        manifest['interval_count'] += 1
                after = path.stat()
                if (stat.st_size,stat.st_mtime_ns,stat.st_ino) != (after.st_size,after.st_mtime_ns,after.st_ino):
                    raise ValueError('input changed during ingestion')
            except (ValueError, KeyError, OSError, tarfile.TarError, csv.Error) as exc:
                self.db.execute('ROLLBACK TO content')
                manifest.update(status='refused', refusal_reason=str(exc), row_count=None, interval_count=0)
            self.db.execute('RELEASE content')
            if manifest['replicate_info'] is None: manifest['replicate_info'] = sorted(replicates) or None
            self.db.execute('update files set manifest=? where id=?', (dumps(manifest),fid))
        return manifest

    def manifests(self):
        return [json.loads(r[0]) for r in self.db.execute('select manifest from files order by path')]

    def state(self):
        return self.db.total_changes,self.db.execute('PRAGMA data_version').fetchone()[0]

    def _invalidate(self):
        current=self.state()
        if current!=self._cache_version:
            self._cached_contract=None; self._availability_cache.clear(); self._cache_version=current

    def contract(self):
        self._invalidate()
        if self._cached_contract is not None: return self._cached_contract
        digest = hashlib.sha256(dumps(self.manifests()).encode())
        for row in self.db.execute('select * from sources order by name'):
            digest.update(dumps(tuple(row)).encode())
        for row in self.db.execute('select * from genes order by gene_id'):
            digest.update(dumps(tuple(row)).encode()); digest.update(b'\n')
        self._cached_contract=digest.hexdigest()
        return self._cached_contract

    def _availability(self, rbp, cell, kind):
        self._invalidate()
        key=(rbp,cell,kind)
        if key in self._availability_cache: return dict(self._availability_cache[key])
        files = [json.loads(r[0]) for r in self.db.execute('select manifest from files where rbp=? and cell=? and kind=? order by path', (rbp,cell,kind))]
        valid = [f for f in files if f['status'] == 'indexed']
        status = ('no_local_assay' if not files else 'unverified_or_refused' if not valid
                  else 'partial' if len(valid) != len(files) else 'ok')
        result = dict(status=status, missing=status != 'ok', evidence_class='derived' if valid else 'not_available',
                    files=[f['accession'] for f in files], verified_files=[f['accession'] for f in valid])
        self._availability_cache[key]=result
        return dict(result)

    def _hits(self, row, gid, cid, kind, cap):
        # RTree float32 bounds are outward-rounded: exact integer recheck is vital.
        return self.db.execute('''select i.* from regions r join intervals i on i.id=r.id
            where r.lo<=? and r.hi>=? and r.g0<=? and r.g1>=? and r.c0<=? and r.c1>=?
            and i.group_id=? and i.chrom_id=? and i.kind=? and i.lo<=? and i.hi>=?
            and i.strand in (?,'.') order by i.source_key,i.lo,i.hi limit ?''',
            (row['end1'],row['start1'],gid,gid,cid,cid,gid,cid,kind,row['end1'],row['start1'],row['strand'],cap+1)).fetchall()

    @staticmethod
    def _measurement(hit, row):
        value = json.loads(hit['payload'])
        distance = (hit['lo']-row['end1'] if hit['lo']>row['end1'] else hit['hi']-row['start1'] if hit['hi']<row['start1'] else 0)
        value.update(peak_start1=hit['lo'], peak_end1=hit['hi'], peak_strand=hit['strand'],
                     signed_distance=distance*(-1 if hit['strand']=='-' else 1),
                     distance_orientation='peak_strand' if hit['strand'] in ('+','-') else 'genomic', evidence_class='derived')
        return value

    def context(self, row, rbp, cell, *, max_hits=128):
        row = normalize_query(row)
        if max_hits < 1 or max_hits > 10000: raise ValueError('invalid hit cap')
        gid = RBPS.index(rbp)*2 + CELLS.index(cell)
        item = self.db.execute('select id from chroms where name=?', (row['chrom'],)).fetchone()
        cid = item[0] if item else None
        result = self._availability(rbp,cell,'peaks') | dict(rbp=rbp,cell_line=cell,presence=None,nearest=None,
                      overlaps=[],overlaps_truncated=False,presence_evidence_class='not_available')
        hits = self._hits(row,gid,cid,'peaks',max_hits) if cid is not None else []
        result['presence'] = True if hits else False if not result['missing'] else None
        if result['presence'] is not None: result['presence_evidence_class'] = 'derived'
        result['overlaps'] = [self._measurement(h,row) for h in hits[:max_hits]]
        result['overlaps_truncated'] = len(hits)>max_hits
        nearest = hits[:1]
        if not nearest and cid is not None:
            candidates = []
            for strand in (row['strand'], '.'):
                for column, bound, op, direction in [('lo',row['end1'],'>','asc'), ('hi',row['start1'],'<','desc')]:
                    h = self.db.execute(f'''select * from intervals where group_id=? and chrom_id=? and kind='peaks'
                        and strand=? and {column}{op}? order by {column} {direction}, source_key,lo,hi limit 1''', (gid,cid,strand,bound)).fetchone()
                    if h: candidates.append(h)
            if candidates:
                nearest = [min(candidates, key=lambda h:(abs(self._measurement(h,row)['signed_distance']),h['source_key'],h['lo'],h['hi']))]
        if nearest: result['nearest'] = self._measurement(nearest[0],row)
        kd = self._availability(rbp,cell,'kd')
        kh = self._hits(row,gid,cid,'kd',max_hits) if cid is not None else []
        # MXE has two blocks: deduplicate same event payload without collapsing readouts.
        seen = set(); events = []
        for h in kh[:max_hits]:
            if h['payload'] not in seen: events.append(json.loads(h['payload'])); seen.add(h['payload'])
        result['knockdown'] = kd | dict(events=events,events_truncated=len(kh)>max_hits,
                                       evidence_class='knockdown_not_steric_block', measurement_evidence_class='observed' if kd['verified_files'] else 'not_available')
        return result

    def load_genes(self, path):
        with self.db, safe_path(path).open() as f:
            self.db.execute('delete from genes')
            for row in csv.DictReader(f, delimiter='\t'):
                row['start1'],row['end1'] = int(row['start1']),int(row['end1']); fence(row)
                self.db.execute('insert into genes values (?,?,?,?,?)', (row['gene_id'],row['chrom'],row['start1'],row['end1'],row['strand']))
            self.db.execute('insert or replace into sources values (?,?)',('genes',digests(path)[0]))


def matched_decoys(row, fetch, *, shifts=(-1000,-500,-100,100,500,1000), gc_tolerance=0.05):
    """Finite deterministic same-gene GC-matched proxy controls, no null p-values."""
    row = normalize_query(row)
    if not 0 <= gc_tolerance <= 1: raise ValueError('invalid GC tolerance')
    lo,hi = row.get('gene_start1'),row.get('gene_end1')
    if fetch is None or lo is None or hi is None: return []
    if not row.get('gene_id') or not 1 <= lo <= row['start1'] <= row['end1'] <= hi:
        raise ValueError('query outside declared gene bounds')
    length = row['end1']-row['start1']+1
    if length>1_100_000: raise ValueError('control interval exceeds bounded reference fetch')
    def gc(a,b):
        seq = fetch(row['chrom'],a,b).upper()
        if len(seq)!=length: raise ValueError('short reference fetch')
        return (seq.count('G')+seq.count('C'))/length if set(seq)<=set('ACGT') else None
    target = gc(row['start1'],row['end1'])
    if target is None: return []
    result = []
    for shift in sorted(set(shifts)):
        a,b = row['start1']+shift,row['end1']+shift
        if a<lo or b>hi or not (b<row['start1'] or a>row['end1']): continue
        value = gc(a,b)
        if value is not None and abs(value-target)<=gc_tolerance:
            result.append(row | dict(site_id=str(row.get('site_id',row.get('row_id')))+f':shift{shift}',start1=a,end1=b,
                                     footprint_start=a,footprint_end=b,coordinate_shift=shift,gc_fraction=value,
                                     query_gc_fraction=target,evidence_class='proxy'))
    return result


def _schema():
    return pa.schema([(n,pa.int64() if n in ('start1','end1') else pa.string()) for n in
        ('site_id','aso_id','chrom','start1','end1','strand','gene_id','assembly','annotation','evidence_class','rbp_context','matched_controls')])


def _valid_done(out, contract):
    marker = out.with_suffix('.done')
    if not marker.exists(): return None
    try: old = json.loads(marker.read_text())
    except (ValueError,OSError): return None
    if not isinstance(old,dict): return None
    if old.get('contract') is not None and old['contract'] != contract: raise ValueError('output contract changed; use a fresh output')
    try:
        if (old.get('contract') == contract and type(old.get('rows')) is int and old['rows']>=0
                and digests(out)[0] == old.get('sha256') and pq.read_metadata(out).num_rows == old['rows']): return old
    except (OSError,pa.ArrowException): pass
    return None


def _atomic_parquet(out, rows, schema, contract, batch_size):
    out = safe_path(out); out.parent.mkdir(parents=True,exist_ok=True)
    previous = _valid_done(out,contract)
    if previous is not None: return previous
    out.with_suffix('.done').unlink(missing_ok=True)
    fd,temp = tempfile.mkstemp(dir=out.parent,suffix='.parquet.tmp'); os.close(fd)
    count, batch = 0, []
    try:
        with pq.ParquetWriter(temp,schema) as writer:
            for row in rows:
                batch.append(row); count+=1
                if len(batch)>=batch_size:
                    writer.write_table(pa.Table.from_pylist(batch,schema=schema)); batch=[]
            if batch: writer.write_table(pa.Table.from_pylist(batch,schema=schema))
        with open(temp,'rb') as f: os.fsync(f.fileno())
        os.replace(temp,out)
        result = dict(contract=contract,rows=count,sha256=digests(out)[0])
        atomic_json(out.with_suffix('.done'),result)
        return result
    finally:
        if os.path.exists(temp): os.unlink(temp)


def write_features(source, out, index, *, batch_size=128, fetch=None, reference_contract=None, require_query_done=False):
    if not 1 <= batch_size <= 4096: raise ValueError('invalid batch size')
    if fetch is not None and reference_contract is None: raise ValueError('reference provenance required for GC controls')
    def check_reference():
        if reference_contract and reference_contract.get('fasta'):
            from conference.model_v2.j6_site_universe import file_identity
            identity=reference_contract['fasta']
            if file_identity(safe_path(identity['path']))!=identity:
                raise ValueError('reference changed during output')
    check_reference()
    source,out = safe_path(source),safe_path(out)
    if source == out: raise ValueError('output would overwrite query')
    before = source.stat()
    pf = pq.ParquetFile(source)
    if require_query_done:
        try: marker=json.loads(source.with_suffix('.done').read_text())
        except (OSError,ValueError) as exc: raise ValueError('query completion marker required') from exc
        if (not isinstance(marker,dict) or type(marker.get('rows')) is not int
                or marker['rows']!=pf.metadata.num_rows or marker.get('sha256')!=digests(source)[0]):
            raise ValueError('query completion marker count/hash mismatch')
    if OUTCOME_COLUMNS & {c.lower() for c in pf.schema_arrow.names}: raise ValueError('outcome columns forbidden')
    columns = [c for c in QUERY_COLUMNS if c in pf.schema_arrow.names]
    contract = dict(version=VERSION,code_sha256=digests(__file__)[0],query_sha256=digests(source)[0],
                    index_sha256=index.contract(),reference=reference_contract,assembly='GRCh38',annotation='GENCODE v50',
                    cell_pooling='none',replicate_rule='within_cell_union_presence_separate_signal',max_hits=128,
                    decoy_shifts=[-1000,-500,-100,100,500,1000],gc_tolerance=0.05)
    changes = index.state()
    def rows():
        for batch in pf.iter_batches(batch_size=batch_size,columns=columns):
            for original in batch.to_pylist():
                row = normalize_query(original)
                gene = index.db.execute('select * from genes where gene_id=?',(row.get('gene_id'),)).fetchone()
                if gene:
                    if (gene['chrom'],gene['strand']) != (row['chrom'],row['strand']): raise ValueError('gene/query geometry mismatch')
                    if (row.get('gene_start1',gene['lo']),row.get('gene_end1',gene['hi'])) != (gene['lo'],gene['hi']): raise ValueError('gene bounds disagree')
                    row.update(gene_start1=gene['lo'],gene_end1=gene['hi'])
                contexts = [index.context(row,rbp,cell) for rbp in RBPS for cell in CELLS]
                decoys = matched_decoys(row,fetch)
                controls = dict(evidence_class='proxy' if decoys else 'not_available',
                    status='descriptive' if decoys else 'no_reference_gene_bounds_or_matching_candidates',n=len(decoys),
                    intervals=[{k:d[k] for k in ('start1','end1','gene_id','coordinate_shift','gc_fraction','query_gc_fraction','evidence_class')} for d in decoys],
                    baseline_overlap_rate=[])
                for context in contexts:
                    flags=[index.context(d,context['rbp'],context['cell_line'])['presence'] for d in decoys]
                    available=[flag for flag in flags if flag is not None]
                    controls['baseline_overlap_rate'].append(dict(rbp=context['rbp'],cell_line=context['cell_line'],n=len(available),
                        overlap_rate=sum(available)/len(available) if available else None,evidence_class='proxy' if available else 'not_available'))
                yield {k:row.get(k) for k in ('aso_id','chrom','start1','end1','strand','gene_id')} | dict(
                    site_id=row.get('site_id') or str(row['row_id']),assembly='GRCh38',annotation='GENCODE v50',evidence_class='derived',
                    rbp_context=dumps(contexts),matched_controls=dumps(controls))
        after=source.stat()
        check_reference()
        if (before.st_size,before.st_mtime_ns,before.st_ino)!=(after.st_size,after.st_mtime_ns,after.st_ino) or index.state() != changes:
            raise ValueError('query/index changed during output')
    return _atomic_parquet(out,rows(),_schema(),contract,batch_size)


def inventory(root):
    """Only explicit RBP/cell folders; never recursively enumerate scratch."""
    root = safe_path(root)
    for kind,folder,suffix in [('peaks','eclip_peaks','.bed.gz'),('kd','kd_splicing','.tar.gz')]:
        for rbp in RBPS:
            for cell in CELLS:
                directory = root/folder/rbp/cell
                if directory.is_dir():
                    for path in sorted(directory.iterdir()):
                        if path.name.endswith(suffix) and path.is_file(): yield path,rbp,cell,kind


def sentinel_rows():
    base = dict(chrom='chr5',strand='+',gene_id='SMN2_sentinel_anchor',assembly='GRCh38',annotation='GENCODE v50')
    for aso,length in [('F18',18),('F14',14),('F10',10)]:
        yield base | dict(site_id='j5-sentinel-'+aso,aso_id=aso,start1=70076584,end1=70076584+length-1)
    yield base | dict(site_id='j5-sentinel-exon7-1Mb',aso_id=None,start1=69576584,end1=70576584)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--encode-root',type=Path)
    parser.add_argument('--index-dir',type=Path,help='Reuse completed index read-only; no ingestion in this job')
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--metadata',type=Path)
    parser.add_argument('--queries',nargs='+',type=Path)
    parser.add_argument('--query-list',type=Path)
    parser.add_argument('--sentinel',action='store_true')
    parser.add_argument('--require-query-done',action='store_true')
    parser.add_argument('--fasta',type=Path); parser.add_argument('--genes',type=Path)
    args=parser.parse_args(argv)
    if bool(args.encode_root)==bool(args.index_dir): parser.error('select --encode-root or --index-dir')
    if args.index_dir and args.metadata: parser.error('--metadata belongs to the index build')
    if sum(map(bool,(args.queries,args.query_list,args.sentinel)))>1: parser.error('select one query source')
    out=safe_path(args.out); out.mkdir(parents=True,exist_ok=True)
    metadata=json.loads(safe_path(args.metadata).read_text()) if args.metadata else {}
    if not isinstance(metadata,dict): raise ValueError('metadata must map accession to ENCODE File object')
    fetch,reference=None,None
    if args.fasta:
        from conference.model_v2.j6_site_universe import FaiFastaReader, reference_identity, file_identity
        genome=FaiFastaReader(safe_path(args.fasta))
        reference=reference_identity(genome) | dict(fasta=file_identity(args.fasta))
        fetch=genome.fetch
    indexdir=safe_path(args.index_dir) if args.index_dir else out
    if args.index_dir and not (indexdir/'index.done').is_file(): raise ValueError('completed index required')
    with ContextIndex(indexdir/'context.sqlite',readonly=bool(args.index_dir)) as index:
        if args.index_dir:
            complete=json.loads((indexdir/'index.done').read_text())
            if (not isinstance(complete,dict) or complete.get('index_sha256')!=index.contract()
                    or complete.get('provenance_sha256')!=digests(indexdir/'provenance.parquet')[0]):
                raise ValueError('completed index contract mismatch')
            if args.genes:
                genes=index.db.execute('select sha256 from sources where name=?',('genes',)).fetchone()
                if not genes or genes[0]!=digests(args.genes)[0]: raise ValueError('gene provenance does not match completed index')
        else:
            _build_index(args,index,out,metadata)
        manifest=index.manifests()
        if args.sentinel:
            source=out/'sentinel_queries.parquet'
            _atomic_parquet(source,sentinel_rows(),pa.Table.from_pylist(list(sentinel_rows())).schema,dict(version=VERSION,kind='sentinel'),128)
            sources=iter([source])
        elif args.query_list:
            def query_paths():
                with safe_path(args.query_list).open() as f:
                    for line in f:
                        if line.strip(): yield safe_path(line.strip())
            sources=query_paths()
        else: sources=iter(args.queries or ())
        count=0
        for source in sources:
            source=safe_path(source)
            # Digest absolute path avoids equal shard basenames from different dirs.
            name=source.stem+'-'+hashlib.sha256(str(source).encode()).hexdigest()[:12]+'.parquet'
            result=write_features(source,out/'features'/name,index,fetch=fetch,reference_contract=reference,require_query_done=args.require_query_done)
            count+=result['rows']
        print(dumps(dict(files=len(manifest),indexed=sum(m['status']=='indexed' for m in manifest),
                         refused=sum(m['status']=='refused' for m in manifest),rows=count,evidence_class='derived')))
    return 0


def _build_index(args,index,out,metadata):
        (out/'index.done').unlink(missing_ok=True)
        discovered=list(inventory(args.encode_root))  # <=21*2*2 assay folders, not genomic data
        paths={str(safe_path(p)) for p,_,_,_ in discovered}
        # A removed source is missing evidence, never an old binding observation.
        with index.db:
            for old in index.db.execute('select id,path from files').fetchall():
                if old['path'] not in paths:
                    index.db.execute('delete from intervals where file_id=?',(old['id'],)); index.db.execute('delete from files where id=?',(old['id'],))
        for path,rbp,cell,kind in discovered:
            index.add_file(path,rbp,cell,kind=kind,metadata=metadata.get(path.name.split('.')[0]))
        if args.genes: index.load_genes(args.genes)
        manifest=index.manifests()
        # Fixed schema keeps wholly unavailable assembly columns typed as strings.
        ms=pa.schema([(k,pa.int64() if k in ('row_count','interval_count') else pa.string()) for k in
            ('path','accession','rbp','cell_line','kind','assembly','assembly_source','sha256','row_count','interval_count','peak_calling_method','replicate_info','status','refusal_reason','evidence_class','provenance')])
        records=({k:(dumps(m[k]) if k=='replicate_info' else m[k]) for k in ms.names if k!='provenance'} | dict(provenance=dumps(m)) for m in manifest)
        result=_atomic_parquet(out/'provenance.parquet',records,ms,dict(version=VERSION,index_sha256=index.contract()),128)
        atomic_json(out/'index.done',dict(index_sha256=index.contract(),provenance_sha256=result['sha256']))


if __name__=='__main__': raise SystemExit(main())
