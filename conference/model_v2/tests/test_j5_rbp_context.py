"""Synthetic J5 tests: no ENCODE, genome, Slurm, or outcome dependencies."""
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from conference.model_v2 import j5_rbp_context as j5


def peak(tmp_path, text, name='ENCFF000AAA.bed.gz'):
    p = tmp_path / name
    with gzip.open(p, 'wt') as f:
        f.write(text)
    return p


def row(**kw):
    return dict(site_id='s1', aso_id='F10', chrom='chr5', start1=101, end1=110,
                strand='+', gene_id='g', assembly='GRCh38', annotation='GENCODE v50', **kw)


def query(**kw):
    return row() | kw


def get(index, **kw):
    return index.context(query(**kw), 'FUS', 'HepG2')


@pytest.fixture
def index(tmp_path):
    with SyntheticIndex(tmp_path / 'context.sqlite') as idx:
        yield idx


def test_verified_manifest_and_half_open_boundaries(index, tmp_path):
    p = peak(tmp_path, '# assembly=GRCh38\nchr5\t100\t110\tx_rep01\t7\t+\t3.5\t0\t0\t-1\n')
    m = index.add_file(p, 'FUS', 'HepG2')
    assert m['sha256'] == hashlib.sha256(p.read_bytes()).hexdigest()
    assert (m['accession'], m['assembly'], m['row_count']) == ('ENCFF000AAA', 'GRCh38', 1)
    assert m['replicate_info'] == ['rep01']
    assert m['peak_calling_method'] is None
    assert get(index)['presence'] is True
    assert get(index)['nearest']['signed_distance'] == 0
    assert get(index)['overlaps'][0]['signal']['value'] == 3.5
    assert get(index)['overlaps'][0]['signal']['evidence_class'] == 'observed'
    adjacent = get(index, start1=111, end1=120)
    assert adjacent['presence'] is False
    assert adjacent['nearest']['signed_distance'] == -1
    assert get(index, start1=90, end1=100)['nearest']['signed_distance'] == 1


@pytest.mark.parametrize('assembly', [None, 'hg19', 'GRCh37'])
def test_assembly_refusal_is_not_a_negative(index, tmp_path, assembly):
    p = peak(tmp_path, (f'# assembly={assembly}\n' if assembly else '') + 'chr5\t100\t110\tx\t0\t+\n')
    m = index.add_file(p, 'FUS', 'HepG2')
    assert m['status'] == 'refused'
    assert m['row_count'] == 1
    c = get(index)
    assert c['presence'] is None and c['missing'] is True
    assert c['status'] == 'unverified_or_refused'


def test_metadata_is_accession_and_checksum_bound(index, tmp_path):
    p = peak(tmp_path, 'chr5\t100\t110\tx\t0\t+\n')
    md = {'accession': 'ENCFF000AAA', '@id': '/files/ENCFF000AAA/',
          'assembly': 'GRCh38', 'md5sum': hashlib.md5(p.read_bytes()).hexdigest(),
          'biological_replicates': [1], 'peak_calling_method': 'CLIPper'}
    assert index.add_file(p, 'FUS', 'HepG2', metadata=md)['status'] == 'indexed'
    assert index.add_file(p, 'FUS', 'HepG2', metadata=md)['status'] == 'indexed'
    assert index.db.execute('select count(*) from intervals').fetchone()[0] == 1
    for change in ({'accession': 'ENCFF999ZZZ'}, {'md5sum': '0'*32}, {'assembly': 'GRCh37'}):
        assert index.add_file(p, 'FUS', 'HepG2', metadata=md | change)['status'] == 'refused'
    assert get(index)['presence'] is None  # stale accepted intervals must be removed


def test_conflicting_header_cannot_be_overridden(index, tmp_path):
    p = peak(tmp_path, '# assembly=GRCh37\n# assembly=GRCh38\nchr5\t100\t110\tx\t0\t+\n')
    assert index.add_file(p, 'FUS', 'HepG2')['status'] == 'refused'


def test_stranded_and_unstranded_files_keep_cell_separate(index, tmp_path):
    p = peak(tmp_path, '# assembly=GRCh38\nchr5\t100\t110\tx\t0\t-\t2\n')
    index.add_file(p, 'FUS', 'HepG2')
    assert get(index)['presence'] is False
    assert get(index, strand='-')['presence'] is True
    assert get(index, start1=111, end1=120, strand='-')['nearest']['signed_distance'] == 1
    un = peak(tmp_path, '# assembly=GRCh38\nchr5\t100\t110\tx\t0\t.\t-1\n', 'ENCFF000AAB.bed.gz')
    index.add_file(un, 'FUS', 'K562')
    c = index.context(query(strand='-'), 'FUS', 'K562')
    assert c['presence'] is True
    assert c['nearest']['distance_orientation'] == 'genomic'
    assert c['overlaps'][0]['signal']['value'] is None
    absent = index.context(row(), 'PCBP2', 'K562')
    assert absent['presence'] is None and absent['status'] == 'no_local_assay'


def test_partial_coverage_and_no_chromosome_peak(index, tmp_path):
    index.add_file(peak(tmp_path, '# assembly=GRCh38\nchr1\t100\t110\tx\t0\t+\n'), 'FUS', 'HepG2')
    assert get(index)['presence'] is False
    assert get(index)['nearest'] is None
    p = peak(tmp_path, 'chr5\t100\t110\tx\t0\t+\n', 'ENCFF000AAB.bed.gz')
    index.add_file(p, 'FUS', 'HepG2')
    assert get(index)['presence'] is None and get(index)['missing']


@pytest.mark.parametrize('updates', [{'assembly': 'hg19'}, {'start1': 0}, {'end1': 100}, {'strand': '.'}, {'annotation': 'GENCODE v49'}])
def test_query_fence(index, updates):
    with pytest.raises(ValueError):
        index.context(query(**updates), 'FUS', 'HepG2')


def kd_tar(tmp_path, header, values):
    p = tmp_path / 'ENCFF000AAC.tar.gz'
    content = (header + '\n' + values + '\n').encode()
    with tarfile.open(p, 'w:gz') as tar:
        info = tarfile.TarInfo('x/MATS_Norm_output/SE.MATS.JunctionCountOnly.txt')
        info.size = len(content)
        tar.addfile(info, io.BytesIO(content))
    md = {'accession': 'ENCFF000AAC', '@id': '/files/ENCFF000AAC/', 'assembly': 'GRCh38',
          'md5sum': hashlib.md5(p.read_bytes()).hexdigest()}
    return p, md


def test_kd_observed_evidence_and_documented_sample_order(index, tmp_path):
    p, md = kd_tar(tmp_path, 'ID\tchr\tstrand\texonStart_0base\texonEnd\tIncLevelDifference',
                   'e1\tchr5\t+\t100\t110\t-0.25')
    index.add_file(p, 'FUS', 'HepG2', kind='kd', metadata=md)
    kd = get(index)['knockdown']
    assert kd['events'][0]['evidence_class'] == 'knockdown_not_steric_block'
    assert kd['events'][0]['raw_sample1_minus_sample2'] == -0.25
    assert kd['events'][0]['kd_minus_control'] is None
    md |= {'sample1_role': 'control', 'sample2_role': 'knockdown', 'sample_order_evidence': 'explicit experiment sample sheet'}
    index.add_file(p, 'FUS', 'HepG2', kind='kd', metadata=md)
    event = get(index)['knockdown']['events'][0]
    assert (event['kd_minus_control'], event['magnitude'], event['direction']) == (0.25, 0.25, 'increased_inclusion')
    assert get(index, start1=111, end1=120)['knockdown']['events'] == []


def test_mxe_does_not_fill_intronic_hull(index, tmp_path):
    p, md = kd_tar(tmp_path, 'ID\tchr\tstrand\t1stExonStart_0base\t1stExonEnd\t2ndExonStart_0base\t2ndExonEnd\tIncLevelDifference',
                   'e1\tchr5\t+\t100\t110\t200\t210\t0.2')
    index.add_file(p, 'FUS', 'HepG2', kind='kd', metadata=md)
    assert get(index, start1=150, end1=160)['knockdown']['events'] == []
    assert len(get(index, start1=205, end1=210)['knockdown']['events']) == 1


def test_bad_peak_file_is_atomic_and_refused(index, tmp_path):
    p = peak(tmp_path, '# assembly=GRCh38\nchr5\t100\t110\tx\t0\t+\nchr5\t200\t190\tx\t0\t+\n')
    assert index.add_file(p, 'FUS', 'HepG2')['status'] == 'refused'
    assert index.db.execute('select count(*) from intervals').fetchone()[0] == 0


def test_gc_controls_are_same_gene_deterministic_and_not_negatives():
    sequence = 'ACGT'*100
    fetch = lambda chrom, lo, hi: sequence[lo-1:hi]
    q = query(start1=101, end1=110, gene_start1=1, gene_end1=400)
    decoys = j5.matched_decoys(q, fetch, shifts=(-100, -20, 20, 100))
    assert decoys == j5.matched_decoys(q, fetch, shifts=(-100, -20, 20, 100))
    assert len(decoys) == 4
    assert all(d['gene_id'] == 'g' and d['evidence_class'] == 'proxy' for d in decoys)
    assert all(1 <= d['start1'] <= d['end1'] <= 400 for d in decoys)
    assert j5.matched_decoys(q, lambda *args: 'N'*10) == []
    assert j5.matched_decoys(query(), fetch) == []


def test_bounded_hits_have_explicit_truncation(index, tmp_path):
    p = peak(tmp_path, '# assembly=GRCh38\n' + 'chr5\t100\t110\tx\t0\t+\t2\n'*9)
    index.add_file(p, 'FUS', 'HepG2')
    c = index.context(row(), 'FUS', 'HepG2', max_hits=2)
    assert len(c['overlaps']) == 2 and c['overlaps_truncated'] is True


def test_streaming_parquet_idempotence_and_corruption(index, tmp_path):
    index.add_file(peak(tmp_path, '# assembly=GRCh38\nchr5\t100\t110\tx\t0\t+\n'), 'FUS', 'HepG2')
    src, out = tmp_path/'queries.parquet', tmp_path/'out.parquet'
    pq.write_table(pa.Table.from_pylist([row(), query(site_id='s2', start1=201, end1=210)]), src)
    result = j5.write_features(src, out, index, batch_size=1)
    assert result['rows'] == 2
    assert result['sha256'] == hashlib.sha256(out.read_bytes()).hexdigest()
    data = pq.read_table(out).to_pylist()
    pairs = json.loads(data[0]['rbp_context'])
    assert len(pairs) == 42 and pairs[0]['cell_line'] == 'HepG2'
    assert data[0]['evidence_class'] == 'derived'
    stamp = out.stat().st_mtime_ns
    assert j5.write_features(src, out, index, batch_size=1) == result
    assert out.stat().st_mtime_ns == stamp
    out.write_bytes(b'corrupt')
    assert j5.write_features(src, out, index)['rows'] == 2
    marker = out.with_suffix('.done')
    marker.write_text('[]')
    assert j5.write_features(src, out, index)['rows'] == 2
    pq.write_table(pa.Table.from_pylist([query(site_id='changed')]), src)
    with pytest.raises(ValueError, match='contract'):
        j5.write_features(src, out, index)


def test_query_outcomes_and_forbidden_paths_are_rejected(index, tmp_path):
    p = tmp_path/'q.parquet'
    pq.write_table(pa.Table.from_pylist([query(psi=0.5)]), p)
    with pytest.raises(ValueError, match='outcome'):
        j5.write_features(p, tmp_path/'o.parquet', index)
    with pytest.raises((ValueError, RuntimeError)):
        j5.safe_path(tmp_path/'wetlab'/'anything.parquet')


def test_j8_footprint_aliases_and_empty_output(index, tmp_path):
    p, out = tmp_path/'q.parquet', tmp_path/'o.parquet'
    r = row(); r['footprint_start'] = r.pop('start1'); r['footprint_end'] = r.pop('end1')
    pq.write_table(pa.Table.from_pylist([r]), p)
    assert j5.write_features(p, out, index)['rows'] == 1
    empty = tmp_path/'empty.parquet'
    pq.write_table(pa.Table.from_pylist([row()]).slice(0, 0), empty)
    assert j5.write_features(empty, tmp_path/'empty-out.parquet', index)['rows'] == 0


def test_kd_minus_one_is_valid_not_missing(index, tmp_path):
    p, md = kd_tar(tmp_path, 'ID\tchr\tstrand\texonStart_0base\texonEnd\tIncLevelDifference',
                   'e1\tchr5\t+\t100\t110\t-1')
    md |= {'sample1_role':'knockdown','sample2_role':'control','sample_order_evidence':'sample sheet'}
    index.add_file(p,'FUS','HepG2',kind='kd',metadata=md)
    assert get(index)['knockdown']['events'][0]['kd_minus_control'] == -1.0


def test_rtree_large_integer_boundary_exact_recheck(index, tmp_path):
    index.add_file(peak(tmp_path, '# assembly=GRCh38\nchr5\t70076583\t70076601\tx\t0\t+\n'), 'FUS','HepG2')
    assert get(index,start1=70076602,end1=70076611)['presence'] is False
    assert get(index,start1=70076601,end1=70076611)['presence'] is True


def test_nearest_unstranded_is_genomic_on_minus_query(index, tmp_path):
    index.add_file(peak(tmp_path, '# assembly=GRCh38\nchr5\t100\t110\tx\t0\t.\n'), 'FUS','HepG2')
    c=get(index,start1=111,end1=120,strand='-')
    assert c['nearest']['signed_distance'] == -1
    assert c['nearest']['distance_orientation'] == 'genomic'


def test_cli_checkpoint_index_and_sentinel(tmp_path, capsys):
    root=tmp_path/'encode'; d=root/'eclip_peaks'/'FUS'/'HepG2'; d.mkdir(parents=True)
    peak(d, '# assembly=GRCh38\nchr5\t70076583\t70076601\tx_rep01\t0\t+\n')
    idxout=tmp_path/'index'; out=tmp_path/'sentinel'
    lookup=lookup_for(tmp_path,[d/'ENCFF000AAA.bed.gz'],assembly='GRCh38')
    opts=['--assembly-lookup',str(lookup.path),'--lookup-sha256',lookup.sha256]
    assert j5.main(['--encode-root',str(root),'--out',str(idxout)]+opts) == 0
    assert (idxout/'index.done').exists()
    args=['--index-dir',str(idxout),'--out',str(out),'--sentinel']+opts
    assert j5.main(args) == 0
    outputs=list((out/'features').glob('*.parquet')); assert len(outputs)==1
    data=pq.read_table(outputs[0]).to_pylist()
    assert {r['aso_id'] for r in data}=={'F18','F14','F10',None}
    assert len(data)==4
    assert all(json.loads(r['rbp_context'])[0]['presence'] is True for r in data)
    stamp=outputs[0].stat().st_mtime_ns
    assert j5.main(args)==0 and outputs[0].stat().st_mtime_ns==stamp
    manifest=pq.read_table(idxout/'provenance.parquet').to_pylist()
    assert manifest[0]['assembly']=='GRCh38'
    (idxout/'index.done').unlink()
    with pytest.raises(ValueError,match='completed index'):
        j5.main(args)


def test_index_external_mutation_invalidates_contract(index, tmp_path):
    first=index.contract()
    p=peak(tmp_path,'# assembly=GRCh38\nchr5\t100\t110\tx\t0\t+\n')
    with SyntheticIndex(index.path) as writer:
        writer.add_file(p,'FUS','HepG2')
    assert index.contract()!=first
    assert get(index)['presence'] is True


def test_integrated_controls_keep_cells_separate(index, tmp_path):
    p=peak(tmp_path,'# assembly=GRCh38\nchr5\t200\t210\tx\t0\t+\n')
    index.add_file(p,'FUS','HepG2')
    source=tmp_path/'q.parquet'; out=tmp_path/'o.parquet'
    pq.write_table(pa.Table.from_pylist([query(gene_start1=1,gene_end1=400)]),source)
    seq='ACGT'*100
    result=j5.write_features(source,out,index,fetch=lambda chrom,a,b: seq[a-1:b],reference_contract={'test_reference':'synthetic'})
    controls=json.loads(pq.read_table(out).to_pylist()[0]['matched_controls'])
    assert result['rows']==1 and controls['n']==2
    h,k=controls['baseline_overlap_rate'][:2]
    assert h['overlap_rate']==0.5 and h['n']==2
    assert k['overlap_rate'] is None and k['n']==0
    assert 'p_value' not in controls and controls['evidence_class']=='proxy'


def test_manifest_strand_and_header_caller(index,tmp_path):
    m=index.add_file(peak(tmp_path,'# assembly=GRCh38 peak_caller=CLIPper\nchr5\t100\t110\tx\t0\t-\n'),'FUS','HepG2')
    assert m['peak_calling_method']=='CLIPper'
    assert m['strand_counts']=={'+':0,'-':1,'.':0}


def test_query_marker_required_mode(index,tmp_path):
    source=tmp_path/'q.parquet'; out=tmp_path/'o.parquet'
    pq.write_table(pa.Table.from_pylist([row()]),source)
    with pytest.raises(ValueError,match='query completion'):
        j5.write_features(source,out,index,require_query_done=True)
    source.with_suffix('.done').write_text(json.dumps({'rows':1,'sha256':j5.digests(source)[0]}))
    assert j5.write_features(source,out,index,require_query_done=True)['rows']==1
    source.with_suffix('.done').write_text(json.dumps({'rows':2,'sha256':j5.digests(source)[0]}))
    with pytest.raises(ValueError,match='query completion'):
        j5.write_features(source,tmp_path/'o2.parquet',index,require_query_done=True)


def test_feature_order_independent_of_ingestion_order(tmp_path):
    p=peak(tmp_path,'# assembly=GRCh38\nchr5\t100\t110\tx\t0\t+\t2\n')
    q=peak(tmp_path,'# assembly=GRCh38\nchr5\t100\t110\ty\t0\t+\t3\n','ENCFF000AAB.bed.gz')
    with SyntheticIndex(tmp_path/'first.sqlite') as first, SyntheticIndex(tmp_path/'second.sqlite') as second:
        for f in (p,q): first.add_file(f,'FUS','HepG2')
        for f in (q,p): second.add_file(f,'FUS','HepG2')
        assert get(first)==get(second)


def test_reference_mutation_prevents_commit(index,tmp_path):
    from conference.model_v2.j6_site_universe import file_identity
    reference=tmp_path/'reference.fa'; reference.write_text('test reference')
    identity=file_identity(reference)
    source,out=tmp_path/'q.parquet',tmp_path/'o.parquet'
    pq.write_table(pa.Table.from_pylist([query(gene_start1=1,gene_end1=400)]),source)
    def fetch(chrom,a,b):
        reference.write_text('changed test reference')
        return ('ACGT'*100)[a-1:b]
    with pytest.raises(ValueError,match='reference changed'):
        j5.write_features(source,out,index,fetch=fetch,reference_contract={'fasta':identity})
    assert not out.with_suffix('.done').exists()


def test_missing_kd_measurement_is_explicitly_unavailable(index,tmp_path):
    p,md=kd_tar(tmp_path,'ID\tchr\tstrand\texonStart_0base\texonEnd\tIncLevelDifference',
                'e1\tchr5\t+\t100\t110\tNA')
    index.add_file(p,'FUS','HepG2',kind='kd',metadata=md)
    event=get(index)['knockdown']['events'][0]
    assert event['raw_sample1_minus_sample2'] is None
    assert event['measurement_evidence_class']=='not_available'


def lookup_for(tmp_path, paths, *, assembly='hg19', kind='peaks', rbp='FUS', cell='HepG2'):
    p = tmp_path/'lookup.tsv'
    fields = ['kind','rbp','cell','file','accession','assembly','file_format','format_type','output_type','md5','lookup','dataset','target_label','biosample_term_name','identity_status']
    import csv
    with p.open('w') as f:
        writer = csv.DictWriter(f, fields, delimiter='\t'); writer.writeheader()
        for path in paths:
            writer.writerow(dict(kind='eclip_peaks' if kind=='peaks' else 'kd_splicing', rbp=rbp, cell=cell,
                file=path.name, accession=path.name.split('.')[0], assembly=assembly,
                file_format='bed' if kind=='peaks' else 'tar', format_type='narrowPeak' if kind=='peaks' else 'None',
                output_type='peaks' if kind=='peaks' else 'differential splicing quantifications',
                md5=hashlib.md5(path.read_bytes()).hexdigest(), lookup='ok',
                dataset='/experiments/ENCSR000AAA/',target_label=rbp,biosample_term_name=cell,identity_status='verified'))
    from conference.model_v2.jobs.j5_encode_lookup import write_provenance
    with p.open() as f: write_provenance(p,list(csv.DictReader(f,delimiter='\t')))
    return j5.AssemblyLookup(p, expected_sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def lift_chain(tmp_path):
    from feng.rbp.liftover import load_chain_file
    p = tmp_path/'test.chain'
    # chr5 [100,200) -> [1000,1100), minus chr6 [100,200) -> [800,900).
    # chr7 gap expands target width; chr8 has two mappings; chr9 splits equally.
    p.write_text('chain 1 chr5 2000 + 100 200 chr5 3000 + 1000 1100 1\n100\n\n'
        'chain 1 chr6 2000 + 100 200 chr6 2000 - 1100 1200 2\n100\n\n'
        'chain 1 chr7 2000 + 100 200 chr7 3000 + 1000 1110 3\n40 10 20\n50\n\n'
        'chain 1 chr8 2000 + 100 200 chr8 3000 + 1000 1100 4\n100\n\n'
        'chain 1 chr8 2000 + 100 200 chr8 3000 + 2000 2100 5\n100\n\n'
        'chain 1 chr9 2000 + 100 200 chr9 3000 + 1000 1100 6\n40 10 10\n50\n')
    return load_chain_file(p)


def test_lookup_hash_pin_and_missing_accession_fail_closed(tmp_path):
    p=peak(tmp_path,'# assembly=GRCh38\nchr5\t100\t110\tx\t0\t+\n')
    lookup=lookup_for(tmp_path,[p],assembly='GRCh38')
    with j5.ContextIndex(tmp_path/'p.sqlite', lookup=lookup) as idx:
        assert idx.add_file(p,'FUS','HepG2')['status']=='indexed'
        q=peak(tmp_path,'# assembly=GRCh38\nchr5\t100\t110\tx\t0\t+\n','ENCFF999AAA.bed.gz')
        assert idx.add_file(q,'FUS','HepG2')['status']=='refused'
    with pytest.raises(ValueError,match='SHA256'):
        j5.AssemblyLookup(tmp_path/'lookup.tsv',expected_sha256='0'*64)
    with gzip.open(p,'wt') as f: f.write('# assembly=GRCh38\nchr5\t100\t110\tchanged\t0\t+\n')
    with j5.ContextIndex(tmp_path/'q.sqlite',lookup=lookup) as idx:
        assert idx.add_file(p,'FUS','HepG2')['status']=='refused'
    with j5.ContextIndex(tmp_path/'default.sqlite',lookup=j5.AssemblyLookup(test_fixture=True)) as idx:
        assert idx.add_file(q,'FUS','HepG2')['status']=='refused'  # header is never a fallback


def test_hg19_lift_accounting_coordinates_and_rejection_audit(tmp_path):
    p=peak(tmp_path,'chr5\t100\t101\ta\t0\t+\nchr5\t199\t200\tb\t0\t.\n'
        'chr6\t100\t110\tc\t0\t+\nchr7\t130\t160\td\t0\t+\n'
        'chr8\t100\t110\te\t0\t+\nchr9\t130\t160\tf\t0\t+\nchr0\t100\t110\tg\t0\t+\n')
    lookup=lookup_for(tmp_path,[p]); chain=lift_chain(tmp_path)
    with j5.ContextIndex(tmp_path/'lift.sqlite',lookup=lookup,chain=chain) as idx:
        m=idx.add_file(p,'FUS','HepG2')
        assert m['status']=='partial' and m['assembly']=='GRCh38' and m['source_assembly']=='hg19'
        a=m['liftover_accounting']
        assert (a['input_intervals'],a['lifted'],a['unmapped'],a['multimapped'],a['split'],a['size_changed'],a['strand_flipped'])==(7,3,1,1,2,1,1)
        assert m['liftover_chain_sha256']==hashlib.sha256((tmp_path/'test.chain').read_bytes()).hexdigest()
        hit=get(idx,start1=1001,end1=1001)['overlaps'][0]
        assert (hit['peak_start1'],hit['peak_end1'])==(1001,1001)
        assert hit['original_hg19_blocks']==[dict(chrom='chr5',start0=100,end0=101,strand='+')]
        assert hit['coordinate_evidence_class']=='derived'
        stored=json.loads(idx.db.execute("select payload from intervals where kind='peaks' order by id limit 1").fetchone()[0])
        assert stored['evidence_class']=='derived'
        assert get(idx,start1=1002,end1=1002)['presence'] is None
        minus=get(idx,chrom='chr6',start1=891,end1=900,strand='-')['overlaps'][0]
        assert minus['peak_strand']=='-'
        assert idx.db.execute('select count(*) from lift_rejections').fetchone()[0]==4
        assert get(idx,start1=1100,end1=1100)['overlaps'][0]['peak_strand']=='.'


def test_hg19_kd_lifts_blocks_preserves_knockdown_label(tmp_path):
    p,md=kd_tar(tmp_path,'ID\tchr\tstrand\texonStart_0base\texonEnd\tIncLevelDifference',
        'e1\tchr6\t+\t100\t110\t-0.25')
    with j5.ContextIndex(tmp_path/'kd.sqlite',lookup=lookup_for(tmp_path,[p],kind='kd'),chain=lift_chain(tmp_path)) as idx:
        m=idx.add_file(p,'FUS','HepG2',kind='kd')
        assert m['status']=='indexed' and m['liftover_accounting']['lifted']==1
        event=get(idx,chrom='chr6',start1=891,end1=900,strand='-')['knockdown']['events'][0]
        assert event['event_blocks']==[dict(chrom='chr6',start1=891,end1=900,strand='-')]
        assert event['raw_sample1_minus_sample2']==-0.25 and event['kd_minus_control'] is None
        assert event['evidence_class']=='knockdown_not_steric_block'
        assert event['coordinate_evidence_class']=='derived'


def test_sentinel_window_and_region_filter_account_for_excluded_rows(tmp_path):
    rows=list(j5.sentinel_rows())
    assert (rows[-1]['start1'],rows[-1]['end1'])==(69076584,71076584)
    p=peak(tmp_path,'chr5\t100\t110\ta\t0\t+\nchr5\t160\t170\tb\t0\t+\n')
    with j5.ContextIndex(tmp_path/'region.sqlite',lookup=lookup_for(tmp_path,[p]),chain=lift_chain(tmp_path),
                         target_region=('chr5',1001,1010)) as idx:
        m=idx.add_file(p,'FUS','HepG2')
        assert m['row_count']==2 and m['interval_count']==1
        assert m['liftover_accounting']['outside_region']==1
        assert get(idx,start1=1001,end1=1010)['presence'] is True
        assert get(idx,start1=1061,end1=1070)['presence'] is None


class SyntheticIndex(j5.ContextIndex):
    """Bind synthetic assay files to explicit per-test lookup snapshots."""
    def __init__(self,path):
        super().__init__(path,lookup=lookup_for(Path(path).parent,[]))

    def add_file(self,path,rbp,cell,*,kind='peaks',metadata=None):
        headers=j5._headers(path,kind)
        import re
        declarations=[v for h in headers for v in re.findall(r'assembly=([A-Za-z0-9]+)',h)]
        assembly=(metadata or {}).get('assembly') or (declarations[0] if declarations else '')
        self.lookup=lookup_for(self.path.parent,[path],assembly=assembly,kind=kind,rbp=rbp,cell=cell)
        return super().add_file(path,rbp,cell,kind=kind,metadata=metadata)


def test_lookup_refuses_wrong_cell_and_hg19_without_chain(tmp_path):
    p=peak(tmp_path,'chr5\t100\t110\tx\t0\t+\n')
    lookup=lookup_for(tmp_path,[p])
    with j5.ContextIndex(tmp_path/'identity.sqlite',lookup=lookup) as idx:
        assert 'chain' in idx.add_file(p,'FUS','HepG2')['refusal_reason']
        with pytest.raises(ValueError,match='identity'):
            idx.add_file(p,'FUS','K562')
    with j5.ContextIndex(tmp_path/'wrong-cell.sqlite',lookup=lookup) as idx:
        assert 'identity' in idx.add_file(p,'FUS','K562')['refusal_reason']


def test_kd_rejects_entire_event_if_one_exon_unmapped(tmp_path):
    p,_=kd_tar(tmp_path,'ID\tchr\tstrand\t1stExonStart_0base\t1stExonEnd\t2ndExonStart_0base\t2ndExonEnd\tIncLevelDifference',
        'e1\tchr5\t+\t100\t110\t300\t310\t0.2')
    with j5.ContextIndex(tmp_path/'event.sqlite',lookup=lookup_for(tmp_path,[p],kind='kd'),chain=lift_chain(tmp_path)) as idx:
        m=idx.add_file(p,'FUS','HepG2',kind='kd')
        assert m['interval_count']==0 and m['status']=='partial'
        assert m['liftover_accounting']['discarded_event_blocks']==1
        assert m['liftover_accounting']['rejected_intervals']==1
        assert get(idx,start1=1001,end1=1010)['knockdown']['events']==[]
        assert idx.db.execute('select count(*) from lift_rejections').fetchone()[0]==1


def test_regional_reverse_chain_filter_and_restart_contract(tmp_path):
    p=peak(tmp_path,'chr6\t100\t110\tx\t0\t+\nchr6\t110\t120\ty\t0\t+\n')
    lookup=lookup_for(tmp_path,[p]); chain=lift_chain(tmp_path)
    path=tmp_path/'reverse.sqlite'
    with j5.ContextIndex(path,lookup=lookup,chain=chain,target_region=('chr6',891,900)) as idx:
        m=idx.add_file(p,'FUS','HepG2')
        assert m['interval_count']==1 and m['liftover_accounting']['outside_region']==1
        assert get(idx,chrom='chr6',start1=891,end1=900,strand='-')['nearest_coverage']=='regional'
        assert idx.add_file(p,'FUS','HepG2')==m
    with pytest.raises(ValueError,match='contract'):
        j5.ContextIndex(path,lookup=lookup,chain=chain,target_region=('chr6',801,900))


def test_kd_unknown_archive_format_is_explicitly_unavailable(tmp_path):
    p=tmp_path/'ENCFF000AAC.tar.gz'
    with tarfile.open(p,'w:gz') as t:
        info=tarfile.TarInfo('unrelated.txt'); info.size=4; t.addfile(info,io.BytesIO(b'text'))
    with j5.ContextIndex(tmp_path/'unknown.sqlite',lookup=lookup_for(tmp_path,[p],kind='kd'),chain=lift_chain(tmp_path)) as idx:
        m=idx.add_file(p,'FUS','HepG2',kind='kd')
        assert m['status']=='refused' and 'not_available' in m['refusal_reason']
        assert get(idx)['knockdown']['measurement_evidence_class']=='not_available'


def build_fixture(tmp_path, *, missing=False):
    from types import SimpleNamespace
    root=tmp_path/'encode'; folder=root/'eclip_peaks'/'FUS'/'HepG2'
    folder.mkdir(parents=True)
    a=peak(folder,'chr5\t0\t10\tx\t0\t+\n')
    b=peak(folder,'chr5\t100\t110\tx\t0\t+\n','ENCFF000AAB.bed.gz')
    lookup=lookup_for(tmp_path,[a,b],assembly='GRCh38')
    if missing: b.unlink()
    out=tmp_path/'built'; out.mkdir()
    with j5.ContextIndex(out/'context.sqlite',lookup=lookup) as idx:
        j5._build_index(SimpleNamespace(encode_root=root,genes=None),idx,out,{})
    return out,lookup,b


def test_expected_missing_source_is_partial_even_on_first_build(tmp_path):
    out,lookup,b=build_fixture(tmp_path,missing=True)
    with j5.ContextIndex(out/'context.sqlite',readonly=True,lookup=lookup) as idx:
        c=get(idx)
        assert c['presence'] is None and c['status']=='partial' and c['missing']
        assert get(idx,start1=1,end1=10)['presence'] is True
        m={m['accession']:m for m in idx.manifests()}['ENCFF000AAB']
        assert m['status']=='missing_source' and m['sha256'] is None
    rows=pq.read_table(out/'provenance.parquet').to_pylist()
    assert len(rows)==2 and any(r['status']=='missing_source' for r in rows)


@pytest.mark.parametrize('mutation',[
    'delete from intervals', 'delete from regions',
    "update intervals set payload='{}'", 'update regions set lo=lo+1000,hi=hi+1000',
])
def test_committed_index_corruption_refused_at_reader_entry(tmp_path,mutation):
    import sqlite3
    out,lookup,_=build_fixture(tmp_path)
    with sqlite3.connect(out/'context.sqlite') as db: db.execute(mutation)
    with pytest.raises(ValueError,match='SHA256|digest|integrity'):
        with j5.ContextIndex(out/'context.sqlite',readonly=True,lookup=lookup) as idx: get(idx)


def test_open_reader_refuses_later_mutation_and_output_reuse(tmp_path):
    import sqlite3
    out,lookup,_=build_fixture(tmp_path)
    source=tmp_path/'queries.parquet'; dest=tmp_path/'features.parquet'
    pq.write_table(pa.Table.from_pylist([query()]),source)
    with j5.ContextIndex(out/'context.sqlite',readonly=True,lookup=lookup) as idx:
        j5.write_features(source,dest,idx)
        with sqlite3.connect(out/'context.sqlite') as db: db.execute('delete from intervals')
        with pytest.raises(ValueError,match='SHA256|digest|integrity'):
            j5.write_features(source,dest,idx)
        with pytest.raises(ValueError,match='SHA256|digest|integrity'): get(idx)


@pytest.mark.parametrize('target,biosample', [('SRSF1','HepG2'),('FUS','K562'),('','')])
def test_independent_encode_biology_refuses_folder_contradiction(tmp_path,target,biosample):
    import csv
    p=peak(tmp_path,'chr5\t100\t110\tx\t0\t+\n')
    lookup=lookup_for(tmp_path,[p],assembly='GRCh38')
    records=list(lookup.rows.values())
    records[0].update(dataset='/experiments/ENCSR000AAA/',target_label=target,biosample_term_name=biosample)
    with lookup.path.open('w') as f:
        w=csv.DictWriter(f,records[0].keys(),delimiter='\t'); w.writeheader(); w.writerows(records)
    from conference.model_v2.jobs.j5_encode_lookup import write_provenance
    write_provenance(lookup.path,records)
    lookup=j5.AssemblyLookup(lookup.path,expected_sha256=j5.digests(lookup.path)[0])
    with j5.ContextIndex(tmp_path/'bio.sqlite',lookup=lookup) as idx:
        manifest=idx.add_file(p,'FUS','HepG2')
        assert manifest['status']=='refused'
        assert 'ENCODE' in manifest['refusal_reason']
        assert get(idx)['presence'] is None


def test_bulk_context_matches_scalar_and_avoids_per_site_sql(index,tmp_path):
    p=peak(tmp_path,'# assembly=GRCh38\nchr5\t100\t110\ta\t0\t+\n'
        'chr5\t120\t200\tb\t0\t.\nchr5\t150\t160\tc\t0\t-\n')
    index.add_file(p,'FUS','HepG2')
    rows=[query(site_id=str(i),start1=i,end1=i+9,strand=s) for i in range(1,240,3) for s in ('+','-')]
    expected=[[index.context(r,rbp,cell) for rbp in j5.RBPS for cell in j5.CELLS] for r in rows]
    statements=[]; index.db.set_trace_callback(statements.append)
    actual=index.contexts_bulk(rows)
    index.db.set_trace_callback(None)
    assert actual==expected
    assert len(statements)<250  # bounded per batch, never per site x assay


def test_committed_index_cannot_be_reopened_writable_after_corruption(tmp_path):
    import sqlite3
    out,lookup,_=build_fixture(tmp_path)
    with sqlite3.connect(out/'context.sqlite') as db: db.execute('delete from regions')
    with pytest.raises(ValueError,match='integrity|SHA256'):
        j5.ContextIndex(out/'context.sqlite',lookup=lookup)


def test_missing_source_reconciliation_survives_restart_and_restore(tmp_path):
    from types import SimpleNamespace
    out,lookup,b=build_fixture(tmp_path)
    original=b.read_bytes(); b.unlink()
    args=SimpleNamespace(encode_root=tmp_path/'encode',genes=None)
    for _ in range(2):
        with j5.ContextIndex(out/'context.sqlite',lookup=lookup) as idx:
            j5._build_index(args,idx,out,{})
            assert get(idx)['presence'] is None
            assert len(idx.manifests())==2
    b.write_bytes(original)
    with j5.ContextIndex(out/'context.sqlite',lookup=lookup) as idx:
        j5._build_index(args,idx,out,{})
        assert get(idx)['presence'] is True


def test_all_expected_sources_missing_is_partial(tmp_path):
    from types import SimpleNamespace
    out,lookup,b=build_fixture(tmp_path,missing=True)
    (b.parent/'ENCFF000AAA.bed.gz').unlink()
    with j5.ContextIndex(out/'context.sqlite',lookup=lookup) as idx:
        j5._build_index(SimpleNamespace(encode_root=tmp_path/'encode',genes=None),idx,out,{})
        c=get(idx)
        assert c['status']=='partial' and c['presence'] is None


def test_encode_lookup_uses_experiment_not_folder_and_preserves_missing_inventory():
    from conference.model_v2.jobs import j5_encode_lookup as job
    records=[dict(accession='ENCFF000AAA',kind='eclip_peaks',rbp='FUS',cell='HepG2',file='ENCFF000AAA.bed.gz')]
    responses={
        '/files/ENCFF000AAA/':dict(accession='ENCFF000AAA',assembly='GRCh38',dataset='/experiments/ENCSR000AAA/',md5sum='abc',file_format_type=None),
        '/experiments/ENCSR000AAA/':dict(target='/targets/SRSF1-human/',biosample_term_name='K562',biosample_ontology='/biosample-types/K562/'),
        '/biosample-types/K562/':dict(term_name='K562'),
        '/targets/SRSF1-human/':dict(label='SRSF1'),
    }
    stamp_encode_objects(responses)
    actual=job.refresh(records,fetch=responses.__getitem__)
    assert len(actual)==1 and actual[0]['lookup']=='ok'
    assert actual[0]['target_label']=='SRSF1' and actual[0]['biosample_term_name']=='K562'
    assert actual[0]['identity_status']=='contradiction'
    assert actual[0]['format_type']=='None'  # ENCODE null subtype is the pinned KD representation
    def unavailable(_): raise OSError('offline')
    failed=job.refresh(records,fetch=unavailable)
    assert len(failed)==1 and failed[0]['lookup']=='lookup_failed'
    assert not failed[0].get('target_label') and failed[0]['error']


def test_j6_array_partition_has_exact_coverage_and_detects_manifest_changes(tmp_path):
    shards=[f'chr5.{i:012d}-{i+50:012d}' for i in range(1,21)]
    path=tmp_path/'run_manifest.json'
    path.write_text(json.dumps(dict(expected_shards=shards)))
    digest=j5.digests(path)[0]
    parts=[j5.genome_queries(tmp_path,i,6,expected_windows=20,manifest_sha256=digest) for i in range(6)]
    assert sorted(p.name for part in parts for p in part)==sorted(s+'.parquet' for s in shards)
    assert len(set(p for part in parts for p in part))==20
    assert max(map(len,parts))-min(map(len,parts))<=1
    path.write_text(json.dumps(dict(expected_shards=shards[::-1])))
    with pytest.raises(ValueError,match='manifest'):
        j5.genome_queries(tmp_path,0,6,expected_windows=20,manifest_sha256=digest)


def test_metadata_biology_cannot_contradict_verified_lookup(tmp_path):
    p=peak(tmp_path,'chr5\t100\t110\tx\t0\t+\n')
    lookup=lookup_for(tmp_path,[p],assembly='GRCh38')
    metadata=dict(accession='ENCFF000AAA',**{'@id':'/files/ENCFF000AAA/'},md5sum=j5.digests(p)[1],
        target={'label':'SRSF1'},biosample_term_name='K562')
    with j5.ContextIndex(tmp_path/'assertions.sqlite',lookup=lookup) as idx:
        m=idx.add_file(p,'FUS','HepG2',metadata=metadata)
        assert m['status']=='refused' and 'contradiction' in m['refusal_reason']
        assert m['lookup_row']['target_label']=='FUS'
        assert get(idx)['presence'] is None


def test_bulk_multichrom_kd_and_tied_nearest_preserve_scalar_results(index,tmp_path):
    p=peak(tmp_path,'# assembly=GRCh38\nchr5\t10\t20\tz\t0\t+\nchr5\t5\t20\ta\t0\t+\n'
        'chr5\t40\t50\tc\t0\t.\nchr6\t30\t40\td\t0\t-\n')
    index.add_file(p,'FUS','HepG2')
    kd,md=kd_tar(tmp_path,'ID\tchr\tstrand\t1stExonStart_0base\t1stExonEnd\t2ndExonStart_0base\t2ndExonEnd\tIncLevelDifference',
        'e1\tchr5\t+\t5\t20\t40\t50\t0.2')
    index.add_file(kd,'FUS','HepG2',kind='kd',metadata=md)
    rows=[query(chrom=c,start1=a,end1=b,strand=s) for c in ('chr5','chr6','chr7')
          for a,b in ((21,30),(1,60),(10,19)) for s in ('+','-')]
    expected=[[index.context(r,rbp,cell,max_hits=1) for rbp in j5.RBPS for cell in j5.CELLS] for r in rows]
    assert index.contexts_bulk(rows,max_hits=1)==expected


def test_benchmark_does_not_time_a_reused_output_after_lost_timing_checkpoint(tmp_path):
    import os
    from conference.model_v2.jobs.j5_benchmark import measure
    out,lookup,_=build_fixture(tmp_path)
    source=tmp_path/'queries.parquet'; dest=tmp_path/'features.parquet'
    pq.write_table(pa.Table.from_pylist([query()]),source)
    with j5.ContextIndex(out/'context.sqlite',readonly=True,lookup=lookup) as idx:
        j5.write_features(source,dest,idx)
        os.utime(dest,ns=(1,1))
        result=measure(idx,source,dest,'test')
        assert dest.stat().st_mtime_ns>1
        assert result['sites']==1 and result['seconds']>0
        assert measure(idx,source,dest,'test')['reused_timing']



def stamp_encode_objects(objects):
    for path,obj in objects.items():
        obj.update({'@id':path,'@type':[{'files':'File','experiments':'Experiment','targets':'Target','biosample-types':'BiosampleType'}[path.split('/')[1]]]})


@pytest.fixture
def encode_identity():
    record=dict(accession='ENCFF000AAA',kind='eclip_peaks',rbp='FUS',cell='HepG2',file='ENCFF000AAA.bed.gz')
    objects={
        '/files/ENCFF000AAA/':dict(accession='ENCFF000AAA',dataset='/experiments/ENCSR000AAA/'),
        '/experiments/ENCSR000AAA/':dict(target='/targets/FUS-human/',biosample_term_name='HepG2',biosample_ontology='/biosample-types/hepg2/'),
        '/targets/FUS-human/':dict(label='FUS'),
        '/biosample-types/hepg2/':dict(term_name='HepG2'),
    }
    stamp_encode_objects(objects)
    return record,objects


@pytest.mark.parametrize('embedded',[False,True])
def test_lookup_validates_canonical_and_embedded_objects(encode_identity,embedded):
    from conference.model_v2.jobs.j5_encode_lookup import refresh
    record,objects=encode_identity
    if embedded:
        exp=objects['/experiments/ENCSR000AAA/']
        exp['target']=dict(objects['/targets/FUS-human/'])
        exp['biosample_ontology']=dict(objects['/biosample-types/hepg2/'])
    assert refresh([record],fetch=objects.__getitem__)[0]['identity_status']=='verified'


@pytest.mark.parametrize('change', ['review_repro','target_path','ontology_path','term_conflict',
    'file_id','file_type','experiment_id','experiment_type','target_id','target_type',
    'ontology_id','ontology_type','embedded_target_label','embedded_ontology_term',
    'embedded_target_type','embedded_ontology_id','missing_ontology'])
def test_lookup_rejects_malformed_identity(encode_identity,change):
    from conference.model_v2.jobs.j5_encode_lookup import refresh
    record,objects=encode_identity
    exp=objects['/experiments/ENCSR000AAA/']
    if change=='review_repro':
        exp.update(target='/biosample-types/fake/',biosample_ontology={'term_name':'K562'})
        objects['/biosample-types/fake/']={'label':'FUS'}
    elif change=='target_path':
        exp['target']='/biosample-types/hepg2/'
        objects[exp['target']]['label']='FUS'
    elif change=='ontology_path':
        exp['biosample_ontology']='/targets/FUS-human/'
        objects[exp['biosample_ontology']]['term_name']='HepG2'
    elif change=='term_conflict': exp['biosample_term_name']='K562'
    elif change=='missing_ontology': del exp['biosample_ontology']
    elif change.startswith('embedded_'):
        target='target' if 'target' in change else 'biosample_ontology'
        exp[target]=dict(objects[exp[target]])
        key,value={'embedded_target_label':('label','SRSF1'),'embedded_ontology_term':('term_name','K562'),
                   'embedded_target_type':('@type',['File']),'embedded_ontology_id':('@id','/targets/FUS-human/')}[change]
        exp[target][key]=value
    else:
        name,field=change.split('_')
        path={'file':'/files/ENCFF000AAA/','experiment':'/experiments/ENCSR000AAA/',
              'target':'/targets/FUS-human/','ontology':'/biosample-types/hepg2/'}[name]
        objects[path]['@'+field]='/experiments/ENCSR999ZZZ/' if field=='id' else ['File' if name!='file' else 'Experiment']
    row=refresh([record],fetch=objects.__getitem__)[0]
    assert row['lookup']=='lookup_failed' or row['identity_status']=='contradiction'
    assert row['identity_status']!='verified'


def test_benchmark_selects_density_extremes_from_only_named_manifests(tmp_path):
    from conference.model_v2.jobs.j5_benchmark import select_windows
    paths=[tmp_path/f'chr1.{i}.parquet' for i in range(8)]
    for p,n in zip(paths,[1,100,3,80,5,60,7,40]):
        p.with_suffix('.done').write_text(json.dumps(dict(rows=n,sha256='a'*64)))
    (tmp_path/'unrelated.done').write_text('invalid')
    selected=select_windows(paths,high=2,median=2)
    assert [p.name for p,_,_ in selected]==['chr1.1.parquet','chr1.3.parquet','chr1.6.parquet','chr1.7.parquet']
    assert [n for _,n,_ in selected]==[100,80,7,40]
    assert [group for _,_,group in selected]==['high','high','median','median']
    paths[0].with_suffix('.done').write_text(json.dumps(dict(rows=-1)))
    with pytest.raises(ValueError,match='count'):
        select_windows(paths,high=2,median=2)


def test_benchmark_profiles_gc_controls_and_checkpoints_reference(tmp_path):
    from conference.model_v2.jobs.j5_benchmark import measure
    out,lookup,_=build_fixture(tmp_path)
    source=tmp_path/'q.parquet'; dest=tmp_path/'features.parquet'
    pq.write_table(pa.Table.from_pylist([query(gene_id='g',gene_start1=1,gene_end1=3000)]),source)
    fetch=lambda chrom,a,b:'G'*(b-a+1)
    with j5.ContextIndex(out/'context.sqlite',readonly=True,lookup=lookup) as idx:
        result=measure(idx,source,dest,'gc',fetch=fetch,reference_contract={'test_reference':'G'})
        assert json.loads(pq.read_table(dest)['matched_controls'][0].as_py())['n']>0
        assert result['reference']=={'test_reference':'G'}
        assert result['disk_bytes']>=result['output_bytes']>0
        assert result['maxrss_mib']>0 and result['seconds']>0
        with pytest.raises(ValueError,match='contract'):
            measure(idx,source,dest,'gc',fetch=fetch,reference_contract={'test_reference':'A'})


@pytest.mark.parametrize('field',['@id','@type'])
@pytest.mark.parametrize('path',['/files/ENCFF000AAA/','/experiments/ENCSR000AAA/',
    '/targets/FUS-human/','/biosample-types/hepg2/'])
def test_lookup_requires_object_identity_fields(encode_identity,path,field):
    from conference.model_v2.jobs.j5_encode_lookup import refresh
    record,objects=encode_identity
    del objects[path][field]
    assert refresh([record],fetch=objects.__getitem__)[0]['lookup']=='lookup_failed'


def test_lookup_ingestion_refuses_explicit_unverified_identity(tmp_path):
    p=peak(tmp_path,'chr5\t100\t110\tx\t0\t+\n')
    lookup=lookup_for(tmp_path,[p],assembly='GRCh38')
    lookup.rows['ENCFF000AAA']['identity_status']='unavailable'
    with j5.ContextIndex(tmp_path/'unverified.sqlite',lookup=lookup) as idx:
        result=idx.add_file(p,'FUS','HepG2')
        assert result['status']=='refused'
        assert 'identity' in result['refusal_reason']


def test_benchmark_window_uses_production_fasta_reader_and_gene_contract(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from conference.model_v2.jobs.j5_benchmark import profile_window
    from conference.model_v2 import j6_site_universe as j6
    out,lookup,_=build_fixture(tmp_path)
    fasta=tmp_path/'tiny.fa'
    fasta.write_text('>chr5\n'+'G'*4000+'\n')
    fasta.with_suffix('.fa.fai').write_text('chr5\t4000\t6\t4000\t4001\n')
    # Production's GRCh38 length fence cannot bless a tiny fixture genome.
    monkeypatch.setattr(j6,'reference_identity',lambda genome:dict(test_reference='tiny'))
    source=tmp_path/'window.parquet'; dest=tmp_path/'profile.parquet'
    pq.write_table(pa.Table.from_pylist([query(gene_id='g',gene_start1=1,gene_end1=3000)]),source)
    source.with_suffix('.done').write_text(json.dumps(dict(rows=1,sha256=j5.digests(source)[0])))
    args=SimpleNamespace(assembly_lookup=lookup.path,lookup_sha256=lookup.sha256,
        index_dir=out,fasta=fasta,genes=None)
    result=profile_window(args,source,dest)
    assert result['reference']['fasta']==j6.file_identity(fasta)
    assert json.loads(pq.read_table(dest)['matched_controls'][0].as_py())['n']>0
    args.genes=tmp_path/'wrong_genes.tsv'; args.genes.write_text('wrong genes')
    with pytest.raises(ValueError,match='gene provenance'):
        profile_window(args,source,tmp_path/'bad.parquet')



def test_lookup_refuses_conflicting_concrete_object_types(encode_identity):
    from conference.model_v2.jobs.j5_encode_lookup import refresh
    record,objects=encode_identity
    objects['/targets/FUS-human/']['@type']=['Target','BiosampleType','Item']
    assert refresh([record],fetch=objects.__getitem__)[0]['lookup']=='lookup_failed'


@pytest.mark.parametrize('status', [None, '', 'unavailable', 'contradiction'])
def test_production_requires_verified_identity(tmp_path, status):
    import csv
    from conference.model_v2.jobs.j5_encode_lookup import write_provenance
    p=peak(tmp_path,'chr5\t100\t110\tx\t0\t+\n')
    lookup=lookup_for(tmp_path,[p],assembly='GRCh38')
    records=list(lookup.rows.values())
    if status is None: records[0].pop('identity_status',None)
    else: records[0]['identity_status']=status
    with lookup.path.open('w') as f:
        writer=csv.DictWriter(f,records[0].keys(),delimiter='\t'); writer.writeheader(); writer.writerows(records)
    write_provenance(lookup.path,records)
    with pytest.raises(ValueError,match='identity_status|identity unverified'):
        j5.AssemblyLookup(lookup.path,expected_sha256=j5.digests(lookup.path)[0])
    if status is None:
        fixture=j5.AssemblyLookup(lookup.path,expected_sha256=j5.digests(lookup.path)[0],test_fixture=True)
        assert fixture.verify('ENCFF000AAA',p,'FUS','HepG2','peaks',j5.digests(p)[1])


@pytest.mark.parametrize('field', ['missing','tsv_sha256','script_sha256','rows','failures','contradictions','utc_time','python_version'])
def test_lookup_requires_matching_provenance(tmp_path,field):
    p=peak(tmp_path,'chr5\t100\t110\tx\t0\t+\n')
    lookup=lookup_for(tmp_path,[p],assembly='GRCh38')
    sidecar=lookup.path.with_suffix('.tsv.json')
    if field=='missing': sidecar.unlink()
    else:
        record=json.loads(sidecar.read_text()); record[field]=None
        sidecar.write_text(json.dumps(record))
    with pytest.raises(ValueError,match='provenance'):
        j5.AssemblyLookup(lookup.path,expected_sha256=lookup.sha256)


def test_gate_covers_all_windows_and_fails_corruption(tmp_path):
    from conference.model_v2.jobs import j5_gate as gate
    out,lookup,_=build_fixture(tmp_path)
    queries=tmp_path/'queries'; queries.mkdir()
    names=['chr5.0001','chr5.0002']
    manifest=queries/'run_manifest.json'; manifest.write_text(json.dumps(dict(expected_shards=names)))
    full=tmp_path/'full'
    with j5.ContextIndex(out/'context.sqlite',readonly=True,lookup=lookup) as idx:
        for shard,name in enumerate(names):
            source=queries/(name+'.parquet')
            pq.write_table(pa.Table.from_pylist([query()]),source)
            j5.atomic_json(source.with_suffix('.done'),dict(rows=1,sha256=j5.digests(source)[0]))
            dest=full/f'shard-{shard}'/'features'/j5.feature_name(source)
            j5.write_features(source,dest,idx,require_query_done=True)
        opts=dict(index=idx,j6_genome_dir=queries,manifest_sha256=j5.digests(manifest)[0],
                  full_out=full,shard_count=2,expected_windows=2)
        result=gate.audit(**opts)
        assert result['ok'] and result['valid_windows']==2 and result['rows']==2
        assert sum(result['peak_status_counts'].values())==84
        source=queries/(names[1]+'.parquet')
        dest=full/'shard-1'/'features'/j5.feature_name(source)
        marker=dest.with_suffix('.done'); original=marker.read_text()
        for mutation in ('gap','hash','rows','contract'):
            marker.write_text(original)
            if mutation=='gap': marker.unlink()
            else:
                data=json.loads(original)
                if mutation=='hash': data['sha256']='0'*64
                elif mutation=='rows': data['rows']=3
                else: data['contract']['index_sha256']='0'*64
                marker.write_text(json.dumps(data))
            failed=gate.audit(**opts)
            assert not failed['ok'] and failed['valid_windows']==1 and len(failed['gaps'])==1


def test_persistent_profile_uses_one_worker_and_consecutive_shard(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from conference.model_v2.jobs import j5_benchmark as bench
    out,lookup,_=build_fixture(tmp_path)
    sources=[]
    for i in range(5):
        p=tmp_path/f'w{i}.parquet'
        pq.write_table(pa.Table.from_pylist([query(site_id=f's{i}')]),p)
        j5.atomic_json(p.with_suffix('.done'),dict(rows=1,sha256=j5.digests(p)[0])); sources.append(p)
    args=SimpleNamespace(assembly_lookup=lookup.path,lookup_sha256=lookup.sha256,index_dir=out,
        fasta=None,genes=None,shard_id=0,shard_count=64)
    seen=[]; original=bench.measure
    def track(index,*a,**kw):
        seen.append(id(index)); return original(index,*a,**kw)
    monkeypatch.setattr(bench,'measure',track)
    target=tmp_path/'profile'; target.mkdir()
    result=bench.profile_persistent(args,sources[:3],target,total_windows=130,shard_windows=5)
    assert len(set(seen))==1 and len(seen)==3
    assert [r['label'] for r in result['windows']]==[p.name for p in sources[:3]]
    assert result['genes_supplied'] is False
    assert result['projection']['array_tasks']==64
    assert result['projection']['recommended_seconds']>=1.5*result['projection']['shard_seconds']
    assert result['projection']['recommended_mem_mib']>=1.5*result['projection']['peak_rss_mib']
    assert all(r['output_bytes']>0 and r['maxrss_mib']>0 for r in result['windows'])


def test_lookup_refresh_writes_auditable_sidecar_offline(tmp_path,monkeypatch,encode_identity):
    from conference.model_v2.jobs import j5_encode_lookup as job
    record,objects=encode_identity
    inventory=tmp_path/'inventory.tsv'
    import csv
    with inventory.open('w') as f:
        w=csv.DictWriter(f,record.keys(),delimiter='\t'); w.writeheader(); w.writerow(record)
    original=job.refresh
    monkeypatch.setattr(job,'refresh',lambda records:original(records,fetch=objects.__getitem__))
    out=tmp_path/'refreshed.tsv'
    assert job.main(['--inventory',str(inventory),'--out',str(out)])==0
    lookup=j5.AssemblyLookup(out,expected_sha256=j5.digests(out)[0])
    sidecar=lookup.verify_provenance()
    assert (sidecar['rows'],sidecar['failures'],sidecar['contradictions'])==(1,0,0)
    assert out.with_suffix('.tsv.sha256').read_text().strip()==sidecar['tsv_sha256']


def test_benchmark_refuses_missing_lookup_sidecar(tmp_path):
    from types import SimpleNamespace
    from conference.model_v2.jobs.j5_benchmark import profile_window
    out,lookup,_=build_fixture(tmp_path)
    lookup.path.with_suffix('.tsv.json').unlink()
    args=SimpleNamespace(assembly_lookup=lookup.path,lookup_sha256=lookup.sha256,
        index_dir=out,fasta=None,genes=None)
    with pytest.raises(ValueError,match='provenance'):
        profile_window(args,tmp_path/'absent.parquet',tmp_path/'out.parquet')


def test_gate_writes_hashed_failure_report_and_exits_nonzero(tmp_path):
    from conference.model_v2.jobs import j5_gate as gate
    report=tmp_path/'gate.json'
    result=gate.main(['--index-dir',str(tmp_path/'absent'),
        '--assembly-lookup',str(tmp_path/'absent.tsv'),'--lookup-sha256','0'*64,
        '--j6-genome-dir',str(tmp_path),'--j6-manifest-sha256','0'*64,
        '--full-out',str(tmp_path),'--report',str(report)])
    assert result==1 and json.loads(report.read_text())['ok'] is False
    assert report.with_suffix('.json.sha256').read_text().strip()==j5.digests(report)[0]
    gate.write_report(report,dict(ok=True))
    assert json.loads(report.read_text())==dict(ok=True)
    assert report.with_suffix('.json.sha256').read_text().strip()==j5.digests(report)[0]


def test_full_allocation_requires_profile_margin_and_matching_config(tmp_path):
    from types import SimpleNamespace
    from conference.model_v2.jobs import j5_benchmark as bench
    out,lookup,_=build_fixture(tmp_path)
    args=SimpleNamespace(assembly_lookup=lookup.path,lookup_sha256=lookup.sha256,index_dir=out,
        fasta=None,genes=None,shard_count=64,j6_manifest_sha256='manifest')
    with j5.ContextIndex(out/'context.sqlite',readonly=True,lookup=lookup) as idx:
        contract=idx.contract()
    record=dict(mode='persistent',windows=[{}]*40,index_sha256=contract,reference=None,
        genes_sha256=None,lookup_sha256=lookup.sha256,manifest_sha256='manifest',
        benchmark_sha256=j5.digests(bench.__file__)[0],
        projection=dict(array_tasks=64,recommended_seconds=900,recommended_mem_mib=3000))
    path=tmp_path/'profile.json'; j5.atomic_json(path,record)
    assert bench.check_allocation(args,path,seconds=900,memory_mib=3000)
    for seconds,memory in [(899,3000),(900,2999)]:
        with pytest.raises(ValueError,match='1.5x'):
            bench.check_allocation(args,path,seconds=seconds,memory_mib=memory)
    record['reference']={'other':'reference'}; j5.atomic_json(path,record)
    with pytest.raises(ValueError,match='configuration'):
        bench.check_allocation(args,path,seconds=900,memory_mib=3000)


@pytest.mark.parametrize('value,expected',[('08:00:00',28800),('1-02:03:04',93784),('00:15:00',900)])
def test_slurm_allocation_time_units(value,expected):
    from conference.model_v2.jobs.j5_benchmark import allocation_seconds
    assert allocation_seconds(value)==expected
