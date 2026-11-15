"""J8: bounded-batch derived-table join. No outcomes or scoring backends.

Coordinates are 1-based inclusive GRCh38/GENCODE v50. Geometry is relative to
curated J2 exon interiors (the manifests include exactly 100 bases of padding).
Distances use the footprint midpoint and transcript orientation. J4 overlap
summaries are explicitly distinct from exact-footprint predictions. JSON feature
cells retain alternative/context identity, status, evidence and missing values.

The SQLite index checkpoints each verified input atomically. Its cache is 32 MiB;
Arrow batches are bounded, but upstream Parquet footers also consume memory.
Full indexing is a Slurm operation; --bounded selects an explicitly partial atlas.
A measured full-genome index is ~37 GB, so --max-index-gib must be raised past its
20 GiB default (see jobs/j8_full.sbatch for the per-source measurements).

Known upstream coverage gap, not a defect here: J2 folded u in {8,12,18} while the
ASO footprints are 10/14/18 nt, so only F18 sites (~0.01% of J6) can ever match a
J2 grid cell exactly. Every other site's j2_accessibility is explicitly
not_available; no length substitution is ever performed. J3's own per-footprint
accessibility channel is unaffected and remains the usable accessibility source.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3
import struct
import tempfile
import zlib
from bisect import bisect_left, bisect_right
from itertools import groupby

import pyarrow as pa
import pyarrow.parquet as pq
from conference.model_v2.evaluation_core import check_source_not_forbidden

VERSION = '1.0.0'
ASSEMBLY, ANNOTATION = 'GRCh38', 'GENCODE v50'
OUTCOME_COLUMNS = {'label', 'label_state', 'outcome', 'psi', 'delta_psi', 'dpsi', 'efficacy', 'potency'}
J3_DROP_COLUMNS = {'target_with_flank', 'flank_left', 'flank_right', 'duplex_structure'}


def safe_path(path):
    p = Path(check_source_not_forbidden(path)).expanduser().resolve()
    if any(x.lower() in {'wetlab', 'terminal', 'psi25'} for x in p.parts):
        raise ValueError('outcome source forbidden')
    return p


def sha256(path):
    h = hashlib.sha256()
    with safe_path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def dumps(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def atomic_json(path, value):
    path = Path(path)
    fd, temp = tempfile.mkstemp(dir=path.parent, suffix='.tmp')
    try:
        with os.fdopen(fd, 'w') as stream:
            stream.write(dumps(value) + '\n'); stream.flush(); os.fsync(stream.fileno())
        os.replace(temp, path)
    finally:
        if os.path.exists(temp): os.unlink(temp)


def fence(row):
    if row.get('assembly') != ASSEMBLY or row.get('annotation') != ANNOTATION:
        raise ValueError('assembly/annotation fence requires GRCh38 / GENCODE v50')
    if row.get('strand') not in ('+', '-') or not 1 <= row['start1'] <= row['end1']:
        raise ValueError('invalid 1-based inclusive coordinates')


def load_loci(paths):
    result = {}
    for path in paths:
        with safe_path(path).open() as stream:
            for line in stream:
                if not line.strip() or line.startswith('#'): continue
                lid, chrom, start, end, strand = line.strip().split('\t')
                row = dict(locus_id=lid, chrom=chrom, start1=int(start), end1=int(end), strand=strand)
                fence(dict(row, assembly=ASSEMBLY, annotation=ANNOTATION))
                if row['end1'] - row['start1'] < 200: raise ValueError('locus lacks 100nt exon padding')
                if lid in result and result[lid] != row: raise ValueError('conflicting locus geometry')
                result[lid] = row
    # These are full canonical transcript exon lists, not selected interior exons.
    # Terminal transcript ends are not splice junctions. Fail on incomplete lists.
    transcripts={}
    for lid in result:
        if ':ex' in lid:
            tx, number=lid.rsplit(':ex',1)
            transcripts.setdefault(tx,[]).append(int(number))
    for tx, numbers in transcripts.items():
        if sorted(numbers)!=list(range(1,max(numbers)+1)):
            raise ValueError('incomplete canonical exon list: '+tx)
    for lid,row in result.items():
        if row['start1']==1:
            raise ValueError('clipped exon padding needs explicit unpadded geometry')
        if ':ex' in lid:
            tx, number=lid.rsplit(':ex',1); number=int(number)
            row['acceptor_position']=(row['start1']+100 if row['strand']=='+' else row['end1']-100) if number>1 else None
            row['donor_position']=(row['end1']-100 if row['strand']=='+' else row['start1']+100) if number<max(transcripts[tx]) else None
    return sorted(result.values(), key=lambda x:x['locus_id'])


class GeometryIndex:
    """Small curated exon index; nearest boundaries via bisect, overlaps by chromosome."""
    def __init__(self, loci):
        self.groups = {}; self.boundaries = {}; self.intervals={}
        for loc in loci:
            key = (loc['chrom'], loc['strand'])
            self.groups.setdefault(key, []).append(loc)
        for key, group in self.groups.items():
            ordered=sorted(group,key=lambda l:(l['start1'],l['locus_id']))
            ends=[]; maximum=0
            for loc in ordered:
                maximum=max(maximum,loc['end1']); ends.append(maximum)
            self.intervals[key]=([l['start1'] for l in ordered],ends,ordered)
            for kind in ('donor', 'acceptor'):
                positions=[]
                for loc in group:
                    default=loc['end1']-100 if (kind == 'donor') == (key[1] == '+') else loc['start1']+100
                    position=loc.get(kind+'_position',default)
                    if position is not None: positions.append((position,loc['locus_id']))
                self.boundaries[key, kind] = sorted(positions)

    def overlapping(self,row):
        item=self.intervals.get((row['chrom'],row['strand']))
        if item is None: return []
        starts,ends,ordered=item
        lo=bisect_left(ends,row['start1']); hi=bisect_right(starts,row['end1'])
        return [l for l in ordered[lo:hi] if l['end1']>=row['start1']]

    def annotate(self, row):
        fence(row)
        key = (row['chrom'], row['strand']); group = self.groups.get(key, [])
        result = dict(status='ok' if group else 'not_available', evidence_class='derived' if group else 'not_available',
                      anchor='footprint_midpoint', donor_distance=None, acceptor_distance=None,
                      donor_locus=None, acceptor_locus=None, region_class='not_available')
        mid = (row['start1']+row['end1'])/2
        for kind in ('donor','acceptor'):
            boundaries = self.boundaries.get((key,kind), [])
            i = bisect_left(boundaries, (mid, ''))
            options = boundaries[max(0,i-1):i+1]
            if options:
                pos, lid = min(options, key=lambda x:(abs(x[0]-mid),x))
                result[kind+'_distance'] = (mid-pos)*(1 if row['strand']=='+' else -1)
                result[kind+'_locus'] = lid
        local=self.overlapping(row)
        overlaps = [l for l in local if l['start1']+100 <= row['end1'] and l['end1']-100 >= row['start1']]
        if overlaps:
            result['region_class'] = ('exonic' if any(l['start1']+100 <= row['start1'] and row['end1'] <= l['end1']-100 for l in overlaps)
                                      else 'exon_intron_spanning')
        elif local:
            result['region_class'] = 'intronic_flank'
        elif group: result['region_class'] = 'outside_curated_exons'
        return result


def geometry(row, loci):
    return GeometryIndex(loci).annotate(row)


def accessibility(row, loci, grids):
    fence(row); hits = []
    for loc in loci:
        if (row['chrom'],row['strand']) != (loc['chrom'],loc['strand']): continue
        if not loc['start1'] <= row['start1'] <= row['end1'] <= loc['end1']: continue
        end = row['end1']-loc['start1']+1 if row['strand']=='+' else loc['end1']-row['start1']+1
        length = row['end1']-row['start1']+1
        for g in grids:
            if (g['locus_id'],g['strand'],g['pos_end_1b_transcript'],g['u']) != (loc['locus_id'],row['strand'],end,length): continue
            value = g.get('p_unpaired'); ok = value is not None and math.isfinite(value) and 0 <= value <= 1
            hits.append(dict(locus_id=loc['locus_id'], W=g['W'], L=g['L'], u=g['u'],
                             value=value if ok else None, status='ok' if ok else 'failed',
                             evidence_class='proxy' if ok else 'not_available'))
    return sorted(hits,key=lambda x:(x['locus_id'],x['W'],x['L'],x['u'])) or [dict(
        value=None,status='not_available',evidence_class='not_available',reason='no exact locus/position/u grid cell; no length substitution')]


def atlas_summary(row, atlas):
    fence(row); hits=[]
    for a in atlas:
        if (a['chrom'],a['strand']) != (row['chrom'],row['strand']): continue
        if a['footprint_start'] > row['end1'] or a['footprint_end'] < row['start1']: continue
        exact=(a['footprint_start'],a['footprint_end'])==(row['start1'],row['end1'])
        v=a.get('raw_model_delta')
        ok=a.get('status')=='ok' and a.get('model_evidence') is True and v is not None and math.isfinite(v)
        hits.append({k:a.get(k) for k in ('locus','config_id','site_kind','site_position','tissue','footprint_start','footprint_end')} |
                    dict(value=v if ok else None, status='ok' if ok else 'failed', exact=exact,
                         evidence_class='proxy' if ok else 'not_available'))
    hits.sort(key=lambda a:dumps(a))
    exact=[a for a in hits if a['exact']]; valid=[a['value'] for a in hits if a['status']=='ok']
    exact_ok=bool(exact) and all(a['status']=='ok' for a in exact)
    return dict(status='ok' if valid else ('failed' if hits else 'not_available'),
                evidence_class='proxy' if valid else 'not_available', evidence_label='sequence_ablation_proxy',
                exact_status='ok' if exact_ok else ('failed' if exact else 'not_available'),
                exact_mean=sum(a['value'] for a in exact)/len(exact) if exact_ok else None,
                overlap_mean=sum(valid)/len(valid) if valid else None,
                n_failed=sum(a['status']!='ok' for a in hits), n_overlap=len(hits), cells=hits,
                interpretation='raw N-ablation deltas; overlap is descriptive, not a site-specific prediction')


def input_identity(path, kind):
    path=safe_path(path); marker=(path.parent.parent/(path.stem+'.done')) if kind=='j4' else path.with_suffix('.done')
    m=json.loads(marker.read_text())
    entry=m.get('ablation',m)
    digest=sha256(path); rows=pq.read_metadata(path).num_rows
    if entry.get('sha256')!=digest or entry.get('rows')!=rows:
        raise ValueError('invalid input commit marker')
    if kind=='j2':
        coordinates=(pq.read_schema(path).metadata or {}).get(b'coordinates',b'')
        if coordinates:
            valid=b'GRCh38; GENCODE v50' in coordinates
        else:
            fasta=str(m.get('fasta',''))
            valid=('/gencode/v50/GRCh38.primary_assembly.genome.fa' in fasta and
                   Path(str(m.get('loci',''))).name in ('j2_loci_phase1.tsv','j2_loci_phase2.tsv'))
        if not valid: raise ValueError('J2 assembly fence requires coordinate metadata or pinned upstream FASTA/loci provenance')
    return dict(path=str(path),sha256=digest,rows=rows,kind=kind,marker_sha256=sha256(marker))


class FeatureIndex:
    """Single writer; transactional per-shard resume. No directory discovery."""
    def __init__(self,path,max_bytes=20*1024**3):
        self.path=safe_path(path); self.path.parent.mkdir(parents=True,exist_ok=True)
        self.db=sqlite3.connect(self.path); self.max_bytes=max_bytes
        self.db.executescript('''PRAGMA cache_size=-32768; PRAGMA temp_store=FILE;
          CREATE TABLE IF NOT EXISTS inputs(path TEXT PRIMARY KEY, identity TEXT);
          CREATE TABLE IF NOT EXISTS j3(id TEXT PRIMARY KEY,payload TEXT);
          CREATE TABLE IF NOT EXISTS j2(locus TEXT,strand TEXT,u INTEGER,w INTEGER,l INTEGER,first INTEGER,last INTEGER,values_blob BLOB,
            PRIMARY KEY(locus,strand,u,w,l,first));
          CREATE TABLE IF NOT EXISTS j4(id TEXT PRIMARY KEY,chrom TEXT,strand TEXT,start INTEGER,end INTEGER,payload TEXT);
          CREATE INDEX IF NOT EXISTS j4_overlap ON j4(chrom,strand,start,end);''')
    def __enter__(self): return self
    def __exit__(self,*args): self.db.close()
    def identities(self):
        return [json.loads(x[0]) for x in self.db.execute('SELECT identity FROM inputs ORDER BY path')]
    def ingest(self,kind,path,*,config_id=None):
        if kind not in ('j2','j3','j4'): raise ValueError('unknown index kind')
        identity=input_identity(path,kind); identity['config_id']=config_id
        prior=self.db.execute('SELECT identity FROM inputs WHERE path=?',(identity['path'],)).fetchone()
        if prior:
            if prior[0]!=dumps(identity): raise ValueError('index input changed')
            return
        pf=pq.ParquetFile(path)
        if OUTCOME_COLUMNS.intersection(x.lower() for x in pf.schema_arrow.names): raise ValueError('outcome columns forbidden')
        if kind=='j4':
            coords=(pf.schema_arrow.metadata or {}).get(b'coordinates',b'')
            if b'GRCh38; GENCODE v50' not in coords: raise ValueError('atlas assembly fence')
        # Project away the bulk sequence/structure columns at read time: the J3 shards
        # are ~357 MB each and ~95% of that is flank/structure text this join discards.
        columns=[n for n in pf.schema_arrow.names if n not in J3_DROP_COLUMNS] if kind=='j3' else None
        try:
            with self.db:
                for batch in pf.iter_batches(batch_size=512,columns=columns):
                    rows=batch.to_pylist()
                    if kind=='j2':
                        # Bounded compressed contiguous chunks; never one SQLite row per
                        # base/grid cell (the full grid contains billions of cells).
                        key=lambda r:(r['locus_id'],r['strand'],r['u'],r['W'],r['L'])
                        for group,items in groupby(rows,key=key):
                            chunk=list(items); positions=[r['pos_end_1b_transcript'] for r in chunk]
                            if positions!=list(range(positions[0],positions[0]+len(positions))):
                                raise ValueError('J2 positions must be contiguous and ordered within a grid chunk')
                            previous=self.db.execute('SELECT 1 FROM j2 WHERE locus=? AND strand=? AND u=? AND w=? AND l=? AND first<=? AND last>=?',
                                                     (*group,positions[-1],positions[0])).fetchone()
                            if previous: raise ValueError('duplicate J2 grid positions')
                            values=[r['p_unpaired'] if r['p_unpaired'] is not None else float('nan') for r in chunk]
                            blob=zlib.compress(struct.pack('!'+str(len(values))+'d',*values))
                            self.db.execute('INSERT INTO j2 VALUES (?,?,?,?,?,?,?,?)',(*group,positions[0],positions[-1],blob))
                    for row in rows:
                        if kind=='j3':
                            fence(row)
                            # Retain provenance/status and channel fields, omit large redundant flanks.
                            keep={k:v for k,v in row.items() if k not in J3_DROP_COLUMNS}
                            for channel,columns in {'duplex':('duplex_energy_kcal_mol',),
                                    'accessibility':('accessibility_p_unpaired',),
                                    'complementarity':('complementarity_wc_pairs','complementarity_gu_pairs','complementarity_mismatch_count')}.items():
                                if channel+'_status' not in keep: continue
                                invalid=any(keep.get(c) is None or not math.isfinite(keep[c]) for c in columns)
                                if keep[channel+'_status']=='ok' and invalid:
                                    keep[channel+'_status']='failed'
                                    keep[channel+'_reason']='upstream ok status with missing/nonfinite numeric channel'
                                if keep[channel+'_status']!='ok':
                                    for c in columns: keep[c]=None
                                    keep[channel+'_evidence_class']='not_available'
                            self.db.execute('INSERT INTO j3 VALUES (?,?)',(row['site_id'],dumps(keep)))
                        elif kind=='j2':
                            continue
                        else:
                            row=dict(row,config_id=config_id)
                            self.db.execute('INSERT INTO j4 VALUES (?,?,?,?,?,?)',(str(config_id)+':'+row['row_id'],row['chrom'],row['strand'],row['footprint_start'],row['footprint_end'],dumps(row)))
                    size=self.db.execute('PRAGMA page_count').fetchone()[0]*self.db.execute('PRAGMA page_size').fetchone()[0]
                    if size>self.max_bytes: raise ValueError('index disk budget exceeded; raise explicitly in Slurm')
                self.db.execute('INSERT INTO inputs VALUES (?,?)',(identity['path'],dumps(identity)))
        except sqlite3.IntegrityError as exc:
            raise ValueError('duplicate feature keys; index transaction rolled back') from exc
    def annotate(self,row,loci,geo):
        fence(row)
        hit=self.db.execute('SELECT payload FROM j3 WHERE id=?',(row['site_id'],)).fetchone()
        channels=json.loads(hit[0]) if hit else dict(status='not_available',evidence_class='not_available',reason='site_id absent from indexed J3 inputs')
        if hit:
            for key in ('chrom','strand','start1','end1','assembly','annotation','aso_id'):
                if channels.get(key)!=row.get(key): raise ValueError('site identity collision: '+key)
        local=[l for l in geo.overlapping(row) if l['start1']<=row['start1'] and row['end1']<=l['end1']]
        grids=[]
        for l in local:
            pos=row['end1']-l['start1']+1 if row['strand']=='+' else l['end1']-row['start1']+1
            for w,span,u,first,blob in self.db.execute('SELECT w,l,u,first,values_blob FROM j2 WHERE locus=? AND strand=? AND u=? AND first<=? AND last>=?',
                      (l['locus_id'],row['strand'],row['end1']-row['start1']+1,pos,pos)):
                v=struct.unpack_from('!d',zlib.decompress(blob),8*(pos-first))[0]
                grids.append(dict(locus_id=l['locus_id'],strand=row['strand'],pos_end_1b_transcript=pos,W=w,L=span,u=u,p_unpaired=v))
        atlas=[json.loads(x[0]) for x in self.db.execute('SELECT payload FROM j4 WHERE chrom=? AND strand=? AND start<=? AND end>=?',
                                                      (row['chrom'],row['strand'],row['end1'],row['start1']))]
        return dict(row,j3_channels=dumps(channels),curated_geometry=dumps(geo.annotate(row)),
                    j2_accessibility=dumps(accessibility(row,local,grids)),j4_ablation=dumps(atlas_summary(row,atlas)))


def prepare_join(index,locus_paths,output_dir):
    """One shared immutable catalog and geometry index for all output shards."""
    output_dir=safe_path(output_dir); output_dir.mkdir(parents=True,exist_ok=True)
    loci=load_loci(locus_paths)
    locus_inputs=[dict(path=str(safe_path(p)),sha256=sha256(p),kind='loci') for p in locus_paths]
    catalog=dict(version=VERSION,inputs=index.identities(),loci=locus_inputs,
                 assembly=ASSEMBLY,annotation=ANNOTATION,
                 geometry='canonical transcript exons +/-100nt; terminal nonsplice boundaries excluded')
    path=output_dir/'feature_inputs.json'
    if path.exists():
        if json.loads(path.read_text())!=catalog: raise ValueError('join catalog inputs changed; use fresh output')
    else: atomic_json(path,catalog)
    return dict(loci=loci,geo=GeometryIndex(loci),catalog_sha256=sha256(path),
                catalog_path=str(path),locus_inputs=locus_inputs,code_sha256=sha256(__file__),
                index_version=(index.db.total_changes,index.db.execute('PRAGMA data_version').fetchone()[0]))


def join_shard(j6_path,output,index,locus_paths,*,batch_size=512,context=None):
    if batch_size<1: raise ValueError('positive batch size required')
    source=input_identity(j6_path,'j6'); output=safe_path(output); output.parent.mkdir(parents=True,exist_ok=True)
    context=context or prepare_join(index,locus_paths,output.parent)
    if context['index_version']!=(index.db.total_changes,index.db.execute('PRAGMA data_version').fetchone()[0]):
        raise ValueError('feature index changed after catalog snapshot')
    inputs=[source,dict(path=context['catalog_path'],sha256=context['catalog_sha256'],kind='feature_catalog')]+context['locus_inputs']
    contract=dict(version=VERSION,code_sha256=context['code_sha256'],inputs=inputs,
                  catalog_sha256=context['catalog_sha256'],assembly=ASSEMBLY,annotation=ANNOTATION,exon_padding=100)
    marker=output.with_suffix('.done')
    if marker.exists():
        try: old=json.loads(marker.read_text())
        except (ValueError,OSError): old={}
        if not isinstance(old,dict): old={}
        if old.get('contract') is not None and old['contract']!=contract: raise ValueError('join inputs or implementation changed; use a fresh output')
        if output.exists() and sha256(output)==old.get('sha256') and pq.read_metadata(output).num_rows==source['rows']:
            return old
    marker.unlink(missing_ok=True)
    pf=pq.ParquetFile(j6_path)
    if OUTCOME_COLUMNS.intersection(x.lower() for x in pf.schema_arrow.names): raise ValueError('outcome columns forbidden')
    fields=[f.with_metadata({b'evidence_class':b'derived'}) for f in pf.schema_arrow]
    for name in ('j3_channels','curated_geometry','j2_accessibility','j4_ablation'):
        fields.append(pa.field(name,pa.string(),metadata={b'evidence_class':b'derived',b'cell_evidence':b'JSON status/evidence_class; J3 channel-specific companions'}))
    schema=pa.schema(fields,metadata={b'coordinates':b'GRCh38; GENCODE v50; 1-based inclusive',b'j8_version':VERSION.encode()})
    loci=context['loci']; geo=context['geo']; n=0
    fd,temp=tempfile.mkstemp(dir=output.parent,suffix='.parquet.tmp'); os.close(fd)
    try:
        with pq.ParquetWriter(temp,schema,compression='zstd') as writer:
            for batch in pf.iter_batches(batch_size=batch_size):
                records=[index.annotate(r,loci,geo) for r in batch.to_pylist()]
                writer.write_table(pa.Table.from_pylist(records,schema=schema)); n+=len(records)
        if n!=source['rows']: raise ValueError('row-count conservation failed')
        with open(temp,'rb') as stream: os.fsync(stream.fileno())
        os.replace(temp,output)
        result=dict(contract=contract,inputs=inputs,rows=n,sha256=sha256(output))
        atomic_json(marker,result)
        return result
    finally:
        if os.path.exists(temp): os.unlink(temp)


def production_inventory(root,*,bounded=False):
    """Enumerate only names declared by manifests / fixed upstream array sizes."""
    root=safe_path(root)
    manifest=root/'j6_genome_v1/run_manifest.json'
    names=json.loads(manifest.read_text())['expected_shards']
    j6=[root/'j6_genome_v1'/(n+'.parquet') for n in (names[:2] if bounded else names)]
    entries=[('j3',root/'j3_genome_v1'/f'j3_shard{i:04d}_of0256.parquet',None) for i in range(1 if bounded else 256)]
    for phase,count in ((1,32),(2,128)):
        entries += [('j2',root/f'j2_phase{phase}'/f'j2_shard{i:04d}.parquet',None) for i in range(1 if bounded else count)]
    config=root/'j4_configs/j4_configs_manifest.tsv'
    with config.open() as stream:
        rows=csv.DictReader((l.lstrip('#') for l in stream if l.startswith('#config_id') or not l.startswith('#')),delimiter='\t')
        for r in rows:
            if r['status']!='emitted': continue
            cid=r['config_id']; entries.append(('j4',root/'j4_atlas_v1'/cid/'ablation_maps/shard_0.parquet',cid))
            if bounded: break
    return j6,entries,[manifest,config]


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root',default='~/orcd/scratch/feng/model_v2'); p.add_argument('--out',required=True)
    p.add_argument('--bounded',action='store_true'); p.add_argument('--max-index-gib',type=float,default=20)
    p.add_argument('--loci',nargs='+',default=[str(Path(__file__).parent/'manifests'/f'j2_loci_phase{i}.tsv') for i in (1,2)])
    args=p.parse_args(argv); out=safe_path(args.out); out.mkdir(parents=True,exist_ok=True)
    j6,entries,manifests=production_inventory(args.root,bounded=args.bounded)
    locus_paths=args.loci
    with FeatureIndex(out/'index.sqlite',int(args.max_index_gib*1024**3)) as index:
        for kind,path,cid in entries: index.ingest(kind,path,config_id=cid)
        context=prepare_join(index,locus_paths,out)
        count=0
        for path in j6: count+=join_shard(path,out/path.name,index,locus_paths,context=context)['rows']
        if not args.bounded:
            gate=safe_path(args.root)/'j6_genome_v1/gate_report.json'
            expected=json.loads(gate.read_text())['rows']
            if count!=expected: raise ValueError('global J6 row-count conservation failed')
            manifests.append(gate)
        atomic_json(out/'manifest.json',dict(rows=count,n_shards=len(j6),coverage='partial' if args.bounded else 'full',
          inputs=index.identities(),source_manifests=[dict(path=str(p),sha256=sha256(p)) for p in manifests],
          note='J2 exact u only; no outcomes; per-output .done contains J6 and locus input hashes',
          j2_coverage_warning='J2 folded u in {8,12,18}; ASO footprints are 10/14/18 nt, so only '
                              '18-nt (F18) sites can match a J2 grid cell. All other sites report '
                              'j2_accessibility not_available by design; no length substitution.'))


if __name__=='__main__': main()
