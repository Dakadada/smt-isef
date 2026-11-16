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
    with j5.ContextIndex(tmp_path / 'context.sqlite') as idx:
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
    assert j5.main(['--encode-root',str(root),'--out',str(idxout)]) == 0
    assert (idxout/'index.done').exists()
    args=['--index-dir',str(idxout),'--out',str(out),'--sentinel']
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
    with j5.ContextIndex(index.path) as writer:
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
    with j5.ContextIndex(tmp_path/'first.sqlite') as first, j5.ContextIndex(tmp_path/'second.sqlite') as second:
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
