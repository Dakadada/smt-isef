"""Derived GRCh38/GENCODE v50 canonical-only context, independent of outcomes.

Distances use the footprint transcript-5' base: downstream positive, donor origin
last exon base, acceptor origin first exon base. Only internal splice boundaries
count; transcript ends are not splice sites. Region uses the entire footprint;
junction-proximal means crossing a boundary or within 10 bp of one. Bins are
absolute distance bins; signed distances are retained. Terminal/unselected genes
are not_available. J2 covers canonical protein-coding transcripts on chr1-22,X;
no isoform robustness claim. Controls are unlabeled background positions matched
within the same exon/intron, footprint length, region and distance bins, never
verified biological negatives. RBP knockdown is not equivalent to steric block.

CLI build/shard/merge performs heavy I/O ONLY in the supplied Slurm jobs.
No knockdown, sealed, outcome, or wetlab files are opened.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile

VERSION = '1.0.0'
ASSEMBLY = 'GRCh38'
ANNOTATION = 'GENCODE v50'


def implementation_hash():
    return sha256(__file__)


def source_identity(path, index, count):
    """Authenticate the producer partition and retain its reference/universe."""
    import pyarrow.parquet as pq
    marker = verify_done(path)
    try:
        producer = {k: marker[k] for k in
                    ('shard_index', 'shard_count', 'reference', 'inputs', 'params')}
        if (type(producer['shard_index']) is not int or
                type(producer['shard_count']) is not int or
                (producer['shard_index'], producer['shard_count']) != (index, count) or
                not isinstance(producer['reference'], dict) or not producer['reference'] or
                not isinstance(producer['inputs'], list) or not producer['inputs']):
            raise ValueError('invalid producer identity')
        metadata = pq.read_schema(path).metadata or {}
        if json.loads(metadata[b'j3_provenance']) != producer:
            raise ValueError('producer schema/marker provenance mismatch')
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError('invalid J3 producer shard provenance') from exc
    return dict(path=str(Path(path).resolve()), sha256=marker['sha256'],
                rows=marker['rows'], **producer)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(dir=path.parent, prefix=path.name+'.')
    try:
        with os.fdopen(fd, 'w') as f:
            json.dump(value, f, sort_keys=True, allow_nan=False)
            f.write('\n')
            f.flush()
            os.fsync(f.fileno())
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)


def build_exons(gtf, manifests, output):
    """One streaming GTF scan; validate every J2 exon and complete exon sets.

    Padding clipped at coordinate 1 cannot be inverted, so compare forward
    padded GTF coordinates instead. Exon coordinates always come from the GTF.
    Ambiguous canonical transcript selection for a gene fails closed.
    """
    expected = {}
    for path in manifests:
        with open(path) as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                key, chrom, lo, hi, strand = line.rstrip().split('\t')
                value = (chrom, int(lo), int(hi), strand)
                if key in expected:
                    raise ValueError('duplicate manifest exon '+key)
                expected[key] = value
    selected = {k.rsplit(':ex', 1)[0] for k in expected}
    found, transcripts = set(), {}
    verified_header = False
    h = hashlib.sha256()
    with open(gtf, 'rb') as f:
        for raw in f:
            h.update(raw)
            if raw.startswith(b'#'):
                if b'GRCh38' in raw and re.search(rb'version 50\b', raw):
                    verified_header = True
                continue
            if not verified_header:
                raise ValueError('GRCh38 GENCODE version 50 description required')
            cols = raw.decode().rstrip().split('\t')
            if len(cols) != 9 or cols[2] != 'exon':
                continue
            attrs = dict(re.findall(r'(\w+)\s+"?([^";]+)"?;', cols[8]))
            tid = attrs.get('transcript_id')
            if tid not in selected:
                continue
            # Repeated tags cannot be represented by dict; inspect original text.
            if 'tag "Ensembl_canonical"' not in cols[8] or 'gene_type "protein_coding"' not in cols[8]:
                raise ValueError('selected exon is not canonical protein-coding')
            number = attrs.get('exon_number', '').strip()
            key = tid+':ex'+number
            lo, hi = int(cols[3]), int(cols[4])
            if key in found or expected.get(key) != (cols[0], max(1, lo-100), hi+100, cols[6]):
                raise ValueError('manifest/GTF exon mismatch '+key)
            found.add(key)
            gene = attrs['gene_id']
            tx = transcripts.setdefault(gene, dict(gene_id=gene, transcript_id=tid,
                                      chrom=cols[0], strand=cols[6], exons=[]))
            if (tx['transcript_id'],tx['chrom'],tx['strand']) != (tid,cols[0],cols[6]):
                raise ValueError('ambiguous canonical gene '+gene)
            tx['exons'].append([lo, hi])
    if not expected or found != set(expected):
        raise ValueError('missing J2 canonical exons in GTF')
    for tx in transcripts.values():
        tx['exons'].sort()
        _validate_tx(tx)
    provenance = dict(version=VERSION, implementation_sha256=implementation_hash(),
                      assembly=ASSEMBLY, annotation=ANNOTATION,
                      gtf_path=str(Path(gtf).resolve()), gtf_sha256=h.hexdigest(),
                      manifests=[dict(path=str(Path(p).resolve()),sha256=sha256(p)) for p in manifests],
                      selection='J2 Ensembl_canonical protein_coding chr1-22,X; canonical only',
                      isoform_sensitivity='not_evaluated', exon_count=len(found))
    atomic_json(output, dict(transcripts=transcripts, provenance=provenance))


def _validate_tx(tx):
    exons = sorted(tx['exons'])
    if not exons or tx['strand'] not in ('+', '-'):
        raise ValueError('invalid transcript')
    if any(not 1 <= lo <= hi for lo,hi in exons) or any(a[1] >= b[0]-1 for a,b in zip(exons,exons[1:])):
        raise ValueError('invalid or overlapping exons')


def read_exons(path):
    data = json.loads(Path(path).read_text())
    if data.get('provenance', {}).get('implementation_sha256') != implementation_hash():
        raise ValueError('exon builder implementation provenance mismatch')
    for gene, tx in data['transcripts'].items():
        _validate_tx(tx)
        if tx['gene_id'] != gene:
            raise ValueError('gene mismatch')
    return data


def distance_bin(value):
    if value is None:
        return 'not_available'
    n = abs(value)
    return next(label for limit,label in [(0,'0'),(10,'1-10'),(50,'11-50'),(200,'51-200'),(1000,'201-1000'),(float('inf'),'>1000')] if n <= limit)


def geometry(site, tx, proximal=10):
    if site['assembly'] != ASSEMBLY or site['annotation'] != ANNOTATION:
        raise ValueError('GRCh38 / GENCODE v50 required')
    lo, hi = site['start1'], site['end1']
    if type(lo) is not int or type(hi) is not int or not 1 <= lo <= hi or site['strand'] not in ('+', '-'):
        raise ValueError('invalid site coordinates')
    result = dict(evidence_class='derived', geometry_status='not_available',
                  transcript_id=None, isoform_sensitivity='not_evaluated', anchor1=None,
                  donor_position1=None, acceptor_position1=None, donor_distance=None,
                  acceptor_distance=None, donor_bin='not_available', acceptor_bin='not_available',
                  region_class='not_available', exon_length=None, intron_length=None)
    if tx is None or any(site[k] != tx[k] for k in ('gene_id','chrom','strand')):
        return result
    exons = sorted(tx['exons'])
    if lo < exons[0][0] or hi > exons[-1][1]:
        return result
    sign = 1 if site['strand']=='+' else -1
    anchor = lo if sign==1 else hi
    donors = [e[1] for e in exons[:-1]] if sign==1 else [e[0] for e in exons[1:]]
    acceptors = [e[0] for e in exons[1:]] if sign==1 else [e[1] for e in exons[:-1]]
    # Deterministic tie: upstream boundary in transcript orientation.
    nearest = lambda xs: min(xs, key=lambda x:(abs(anchor-x),x*sign)) if xs else None
    donor, acceptor = nearest(donors), nearest(acceptors)
    introns = [[a[1]+1,b[0]-1] for a,b in zip(exons,exons[1:])]
    exon = next((e for e in exons if e[0]<=lo<=hi<=e[1]),None)
    intron = next((e for e in introns if e[0]<=lo<=hi<=e[1]),None)
    boundaries = donors+acceptors
    near = any((lo-proximal)<=b<=(hi+proximal) for b in boundaries)
    result.update(geometry_status='ok', transcript_id=tx['transcript_id'], anchor1=anchor,
                  donor_position1=donor, acceptor_position1=acceptor,
                  donor_distance=(anchor-donor)*sign if donor is not None else None,
                  acceptor_distance=(anchor-acceptor)*sign if acceptor is not None else None,
                  region_class='junction-proximal' if near or (exon is None and intron is None) else ('exon' if exon else 'intron'),
                  exon_length=exon[1]-exon[0]+1 if exon else None,
                  intron_length=intron[1]-intron[0]+1 if intron else None)
    for kind in ('donor','acceptor'):
        result[kind+'_bin'] = distance_bin(result[kind+'_distance'])
    return result


def matched_control(site, tx):
    """Bounded deterministic rejection sampling, null if no match in 64 draws.

    No outcome/sequence matching. Same host segment, length, region and absolute
    distance bins. Controls cannot overlap the source footprint; they may overlap
    other sites and are NOT verified negatives or peak-free by assumption.
    """
    base = geometry(site, tx)
    if base['geometry_status'] != 'ok':
        return None
    exons = sorted(tx['exons'])
    segments = exons+[[a[1]+1,b[0]-1] for a,b in zip(exons,exons[1:])]
    segment = next((e for e in segments if e[0]<=site['start1']<=site['end1']<=e[1]),None)
    if segment is None:
        return None
    width = site['end1']-site['start1']+1
    n = segment[1]-segment[0]-width+2
    for i in range(64):
        seed = f"j7-v1:{site['site_id']}:{site['gene_id']}:{i}".encode()
        start = segment[0]+int.from_bytes(hashlib.sha256(seed).digest()[:8],'big')%n
        end = start+width-1
        if start<=site['end1'] and end>=site['start1']:
            continue
        other = geometry(dict(site,start1=start,end1=end),tx)
        if all(other[k]==base[k] for k in ('region_class','donor_bin','acceptor_bin')):
            return start
    return None


def rbp_audit(manifest):
    """Metadata-only audit: never opens KD tables or claims bytes verified."""
    with open(manifest) as f:
        rows = [r for r in csv.DictReader(f,delimiter='\t') if r['role']=='eclip_peaks']
    assemblies = sorted({r['assembly'] for r in rows})
    return dict(status='BLOCKED', manifest_path=str(Path(manifest).resolve()),
                manifest_sha256=sha256(manifest), eclip_files=len(rows), assemblies=assemblies,
                source_cells=sorted({r['cell_line'] for r in rows}),
                recorded_sha256_count=sum(bool(r.get('sha256')) for r in rows),
                reason='No verified GRCh38 peak bytes with accession/cell/assay/source SHA256 and liftover chain SHA256 + mapping/drop QC supplied. Indexed hg19 sources cannot join. Source bytes not reverified; no KD files read.',
                claim_boundary='RBP knockdown is not steric-block equivalence')


def write_done(path, provenance):
    import pyarrow.parquet as pq
    atomic_json(Path(path).with_suffix('.done'),dict(provenance=provenance,
                sha256=sha256(path),rows=pq.read_metadata(path).num_rows))


def verify_done(path):
    import pyarrow.parquet as pq
    try:
        marker=json.loads(Path(path).with_suffix('.done').read_text())
        if not isinstance(marker,dict) or marker['sha256']!=sha256(path) or marker['rows']!=pq.read_metadata(path).num_rows:
            raise ValueError('hash/count mismatch')
        return marker
    except Exception as exc:
        raise ValueError(f'uncommitted/damaged parquet {path}') from exc


def context_schema():
    import pyarrow as pa
    strings=['site_id','gene_id','transcript_id','geometry_status','evidence_class','isoform_sensitivity',
             'region_class','donor_bin','acceptor_bin','rbp_status','rbp_cell_type','rbp_assay',
             'rbp_accession','rbp_source_sha256','rbp_evidence_class','control_status','control_evidence_class']
    ints=['anchor1','donor_position1','acceptor_position1','donor_distance','acceptor_distance',
          'exon_length','intron_length','control_start1','control_end1','rbp_peak_count']
    return pa.schema([(s,pa.string()) for s in strings]+[(s,pa.int64()) for s in ints]+[('rbp_missing',pa.bool_())])


def run_shard(source, exons, output, index, count, audit):
    import pyarrow as pa
    import pyarrow.parquet as pq
    if not 0<=index<count or audit['status']!='BLOCKED':
        raise ValueError('invalid shard or unsupported RBP provenance')
    source, exons, output=Path(source),Path(exons),Path(output)
    origin=source_identity(source,index,count)
    data=read_exons(exons)
    expected=dict(version=VERSION,implementation_sha256=implementation_hash(),
                  index=index,count=count,source=str(source.resolve()),
                  source_sha256=origin['sha256'],producer=origin,exons_sha256=sha256(exons),
                  geometry_provenance=data['provenance'],rbp=audit,
                  geometry_policy='5prime anchor; proximal=10; canonical only',
                  control_policy='same segment/length/region/bins; 64 draws; unlabeled')
    output.mkdir(parents=True,exist_ok=True)
    path=output/f'j7_shard{index:04d}_of{count:04d}.parquet'
    with path.with_suffix('.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        try:
            if verify_done(path)['provenance']==expected:
                return path
        except ValueError:
            pass
        path.with_suffix('.done').unlink(missing_ok=True)
        fd,temp=tempfile.mkstemp(dir=output,prefix=path.name+'.'); os.close(fd)
        schema=context_schema().with_metadata({b'j7_provenance':json.dumps(expected,sort_keys=True).encode()})
        try:
            with pq.ParquetWriter(temp,schema,compression='zstd') as writer:
                for batch in pq.ParquetFile(source).iter_batches(batch_size=512,columns=['site_id','gene_id','chrom','start1','end1','strand','assembly','annotation']):
                    rows=[]
                    for site in batch.to_pylist():
                        tx=data['transcripts'].get(site['gene_id'])
                        g=geometry(site,tx)
                        control=matched_control(site,tx)
                        rows.append(dict(g,site_id=site['site_id'],gene_id=site['gene_id'],
                            rbp_status='BLOCKED',rbp_missing=True,rbp_peak_count=None,
                            rbp_cell_type=None,rbp_assay=None,rbp_accession=None,rbp_source_sha256=None,
                            rbp_evidence_class='not_available',control_start1=control,
                            control_end1=control+site['end1']-site['start1'] if control is not None else None,
                            control_status='ok' if control is not None else 'not_available',control_evidence_class='unlabeled'))
                    writer.write_table(pa.Table.from_pylist(rows,schema=schema))
            if (source_identity(source,index,count)!=origin or
                    sha256(exons)!=expected['exons_sha256'] or
                    implementation_hash()!=expected['implementation_sha256']):
                raise ValueError('input changed during run')
            with open(temp,'rb') as f: os.fsync(f.fileno())
            os.replace(temp,path)
            write_done(path,expected)
        finally:
            Path(temp).unlink(missing_ok=True)
    return path


def _merge_gate(output, count, expected_rows, input_dir, exons):
    """Stream exact shard set; SQLite duplicate gate bounds RAM for 11M sites."""
    import sqlite3
    import pyarrow.parquet as pq
    output=Path(output)
    paths=[output/f'j7_shard{i:04d}_of{count:04d}.parquet' for i in range(count)]
    if count<1 or expected_rows<0 or set(output.glob('j7_shard*.parquet'))!=set(paths):
        raise ValueError('unexpected/missing shard set')
    read_exons(exons)  # Authenticate builder independently of shard commitments.
    merged=output/'j7_context.parquet'
    merged.with_suffix('.done').unlink(missing_ok=True)
    with tempfile.TemporaryDirectory(dir=output,prefix='j7_merge_') as temporary:
        db=sqlite3.connect(str(Path(temporary)/'ids.sqlite'))
        db.execute('CREATE TABLE ids (site_id TEXT PRIMARY KEY)')
        temp=Path(temporary)/'merged.parquet'
        total=0; shared=None; sources=set(); markers=[]
        try:
            with pq.ParquetWriter(temp,context_schema(),compression='zstd') as writer:
                for i,path in enumerate(paths):
                    marker=verify_done(path); p=marker['provenance']
                    source=Path(input_dir)/f'j3_shard{i:04d}_of{count:04d}.parquet'
                    origin=source_identity(source,i,count)
                    if (p.get('producer')!=origin or p['source']!=str(source.resolve()) or
                            p['source_sha256']!=origin['sha256'] or
                            p.get('implementation_sha256')!=implementation_hash() or
                            p['exons_sha256']!=sha256(exons)):
                        raise ValueError('source/implementation/exon provenance changed')
                    common={k:v for k,v in p.items() if k not in ('index','source','source_sha256','producer')}
                    common['producer_identity']={k:origin[k] for k in ('reference','inputs','params','shard_count')}
                    if p['index']!=i or p['count']!=count or (shared is not None and common!=shared) or p['source'] in sources:
                        raise ValueError('inconsistent shard provenance')
                    shared=common; sources.add(p['source']); markers.append(marker['sha256'])
                    if not pq.read_schema(path).equals(context_schema(),check_metadata=False):
                        raise ValueError('unexpected context schema')
                    for batch in pq.ParquetFile(path).iter_batches(batch_size=512):
                        rows=batch.to_pylist()
                        for r in rows:
                            for kind in ('donor','acceptor'):
                                distance=r[kind+'_distance']
                                if (distance is None)!=(r[kind+'_position1'] is None) or distance_bin(distance)!=r[kind+'_bin']:
                                    raise ValueError('inconsistent distance/null/bin')
                            if r['geometry_status']=='not_available' and any(r[k] is not None for k in ('anchor1','donor_distance','acceptor_distance','exon_length','intron_length')):
                                raise ValueError('fabricated unavailable geometry')
                            if not r['site_id'] or r['evidence_class']!='derived' or not r['rbp_missing'] or r['rbp_peak_count'] is not None:
                                raise ValueError('invalid context evidence')
                        try:
                            db.executemany('INSERT INTO ids VALUES (?)',[(r['site_id'],) for r in rows])
                        except sqlite3.IntegrityError as exc:
                            raise ValueError('duplicate site_id') from exc
                        writer.write_batch(batch); total+=len(rows)
                    db.commit()
            if total!=expected_rows:
                raise ValueError(f'expected {expected_rows} rows, got {total}')
            os.replace(temp,merged)
            write_done(merged,dict(shards_sha256=markers,common=shared))
        finally:
            db.close()
    return merged


def merge_gate(output, count, expected_rows, *, input_dir, exons):
    if type(count) is not int or count<1 or type(expected_rows) is not int or expected_rows<0:
        raise ValueError('invalid shard count/expected rows')
    output=Path(output)
    output.mkdir(parents=True,exist_ok=True)
    with (output/'j7_merge.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        return _merge_gate(output,count,expected_rows,input_dir,exons)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command',required=True)
    b=sub.add_parser('build'); b.add_argument('--gtf',required=True); b.add_argument('--manifests',nargs='+',required=True); b.add_argument('--output',required=True)
    s=sub.add_parser('shard'); s.add_argument('--input-dir',required=True); s.add_argument('--exons',required=True); s.add_argument('--output-dir',required=True); s.add_argument('--index',type=int,required=True); s.add_argument('--count',type=int,default=256); s.add_argument('--rbp-manifest',default='feng/artifacts/rbp/encode_manifest.tsv')
    m=sub.add_parser('merge'); m.add_argument('--output-dir',required=True); m.add_argument('--count',type=int,default=256); m.add_argument('--expected-rows',type=int,required=True); m.add_argument('--input-dir',required=True); m.add_argument('--exons',required=True)
    args=parser.parse_args()
    if args.command=='build':
        build_exons(Path(args.gtf).expanduser(),args.manifests,args.output)
    elif args.command=='shard':
        source=Path(args.input_dir)/f'j3_shard{args.index:04d}_of{args.count:04d}.parquet'
        # Production CLI accepts only committed J3 schema/provenance, not raw J6.
        import pyarrow.parquet as pq
        from conference.model_v2.j3_scoring import _score_schema
        if not set(_score_schema().names)<=set(pq.read_schema(source).names):
            raise ValueError('committed J3 scoring schema required')
        run_shard(source,args.exons,args.output_dir,args.index,args.count,rbp_audit(args.rbp_manifest))
    else:
        merge_gate(args.output_dir,args.count,args.expected_rows,input_dir=args.input_dir,exons=args.exons)


if __name__=='__main__':
    main()
