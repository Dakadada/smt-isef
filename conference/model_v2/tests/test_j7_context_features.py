"""Synthetic geometry and atomic pipeline regressions; no outcome inputs."""
import importlib.util
from pathlib import Path
import json
import pytest

PATH = Path(__file__).parents[1] / 'j7_context_features.py'

def module():
    assert PATH.exists(), 'J7 implementation missing'
    spec = importlib.util.spec_from_file_location('j7', PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

def row(pos, strand='+', end=None):
    return dict(site_id='s', gene_id='g', chrom='chr1', strand=strand,
                start1=pos, end1=end or pos, assembly='GRCh38', annotation='GENCODE v50')

def tx(strand='+'):
    return dict(gene_id='g', transcript_id='t', chrom='chr1', strand=strand,
                exons=[[101,110],[201,220],[301,310]])

@pytest.mark.parametrize('strand,pos,donor,acceptor', [
    ('+',110,0,-91),('+',111,1,-90),('+',200,-20,-1),('+',201,-19,0),
    ('-',301,0,-81),('-',300,1,-80),('-',221,-20,-1),('-',220,-19,0)])
def test_signed_boundaries(strand,pos,donor,acceptor):
    got=module().geometry(row(pos,strand),tx(strand))
    assert got['donor_distance']==donor
    assert got['acceptor_distance']==acceptor
    assert got['evidence_class']=='derived'

@pytest.mark.parametrize('pos,end,region,length', [(105,105,'exon',10),(150,150,'intron',90),
                                                   (109,112,'junction-proximal',None)])
def test_region_and_lengths(pos,end,region,length):
    got=module().geometry(row(pos,end=end),tx(),proximal=0)
    assert got['region_class']==region
    assert got['exon_length' if pos==105 else 'intron_length']==length

@pytest.mark.parametrize('annotation', [None,dict(tx(),chrom='chr2'),dict(tx(),strand='-')])
def test_not_available(annotation):
    got=module().geometry(row(105),annotation)
    assert got['geometry_status']=='not_available'
    assert got['donor_distance'] is None
    assert got['region_class']=='not_available'

def test_outside_transcript_and_single_exon():
    m=module()
    assert m.geometry(row(90),tx())['donor_distance'] is None
    got=m.geometry(row(105),dict(tx(),exons=[[101,110]]))
    assert got['donor_distance'] is None and got['acceptor_distance'] is None
    assert got['exon_length']==10

def test_minus_footprint_anchor():
    got=module().geometry(row(295,'-',300),tx('-'))
    assert got['anchor1']==300 and got['donor_distance']==1

def test_assembly_fence():
    with pytest.raises(ValueError,match='GRCh38'):
        module().geometry(dict(row(105),assembly='hg19'),tx())

def test_controls_deterministic_matched():
    m=module(); r=row(150,end=155)
    a=m.matched_control(r,tx()); b=m.matched_control(r,tx())
    assert a==b and a is not None and a!=150
    assert 111<=a<=195
    c=m.geometry(dict(r,start1=a,end1=a+5),tx())
    d=m.geometry(r,tx())
    assert (c['region_class'],c['donor_bin'],c['acceptor_bin'])==(d['region_class'],d['donor_bin'],d['acceptor_bin'])

def test_build_validation(tmp_path):
    m=module(); manifest=tmp_path/'loci.tsv'; gtf=tmp_path/'a.gtf'; out=tmp_path/'exons.json'
    manifest.write_text('t:ex1\tchr1\t1\t210\t+\nt:ex2\tchr1\t101\t320\t+\n')
    gtf.write_text('##description: GRCh38 version 50\n'+''.join(f'chr1\tx\texon\t{lo}\t{hi}\t.\t+\t.\tgene_id "g"; transcript_id "t"; exon_number "{n}"; tag "Ensembl_canonical"; gene_type "protein_coding";\n' for n,lo,hi in [(1,101,110),(2,201,220)]))
    m.build_exons(gtf,[manifest],out)
    data=m.read_exons(out)
    assert data['transcripts']['g']['exons']==[[101,110],[201,220]]
    assert data['provenance']['gtf_sha256']==m.sha256(gtf)
    gtf.write_text(gtf.read_text().replace('201\t220','202\t220'))
    with pytest.raises(ValueError): m.build_exons(gtf,[manifest],out)

def test_blocked_rbp():
    got=module().rbp_audit(Path('feng/artifacts/rbp/encode_manifest.tsv'))
    assert got['status']=='BLOCKED' and got['eclip_files']==669
    assert got['assemblies']==['hg19'] and got['manifest_sha256']

def write_j3_done(m, path, index, count):
    import pyarrow.parquet as pq
    provenance = dict(shard_index=index, shard_count=count,
                      reference={'path': '/synthetic/ref.fa', 'fai_sha256': 'a'*64},
                      inputs=[{'path': '/synthetic/j6.parquet', 'sha256': 'b'*64}], params={})
    table = pq.read_table(path).replace_schema_metadata({b'j3_provenance': json.dumps(provenance).encode()})
    pq.write_table(table, path)
    m.atomic_json(path.with_suffix('.done'), dict(provenance, sha256=m.sha256(path), rows=table.num_rows))

def test_shard_restart_and_gate(tmp_path):
    import pyarrow as pa
    import pyarrow.parquet as pq
    m=module(); source=tmp_path/'j3_shard0000_of0001.parquet'; ex=tmp_path/'exons.json'; out=tmp_path/'out'
    pq.write_table(pa.Table.from_pylist([row(150)]),source)
    write_j3_done(m,source,0,1)
    m.atomic_json(ex,dict(transcripts={'g':tx()},provenance={'test':True, 'implementation_sha256':m.implementation_hash()}))
    audit=m.rbp_audit(Path('feng/artifacts/rbp/encode_manifest.tsv'))
    result=m.run_shard(source,ex,out,0,1,audit)
    assert pq.read_table(result).to_pylist()[0]['rbp_peak_count'] is None
    stamp=result.stat().st_mtime_ns
    m.run_shard(source,ex,out,0,1,audit)
    assert result.stat().st_mtime_ns==stamp
    merged=m.merge_gate(out,1,input_dir=tmp_path,exons=ex,expected_rows=1)
    assert pq.read_metadata(merged).num_rows==1
    result.write_bytes(b'broken')
    with pytest.raises(ValueError): m.merge_gate(out,1,input_dir=tmp_path,exons=ex,expected_rows=1)
    m.run_shard(source,ex,out,0,1,audit)
    assert pq.read_metadata(result).num_rows==1
    with pytest.raises(ValueError): m.merge_gate(out,1,input_dir=tmp_path,exons=ex,expected_rows=2)

@pytest.mark.parametrize('failure',['duplicate','mixed','missing','null'])
def test_gate_rejects_bad_shards(tmp_path,failure):
    import pyarrow as pa
    import pyarrow.parquet as pq
    m=module(); out=tmp_path/'out'; ex=tmp_path/'ex.json'
    m.atomic_json(ex,dict(transcripts={'g':tx()},provenance={'test':True, 'implementation_sha256':m.implementation_hash()}))
    audit=m.rbp_audit(Path('feng/artifacts/rbp/encode_manifest.tsv'))
    paths=[]
    for i in range(2):
        source=tmp_path/f'j3_shard{i:04d}_of0002.parquet'
        r=dict(row(150),site_id=f's{i}')
        pq.write_table(pa.Table.from_pylist([r]),source); write_j3_done(m,source,i,2)
        paths.append(m.run_shard(source,ex,out,i,2,audit))
    p=paths[1]
    if failure=='missing':
        p.unlink()
    elif failure=='mixed':
        marker=json.loads(p.with_suffix('.done').read_text())
        marker['provenance']['exons_sha256']='different'
        m.atomic_json(p.with_suffix('.done'),marker)
    else:
        table=pq.read_table(p); r=table.to_pylist()[0]
        if failure=='duplicate': r['site_id']='s0'
        else: r['donor_distance']=None
        provenance=json.loads(p.with_suffix('.done').read_text())['provenance']
        pq.write_table(pa.Table.from_pylist([r],schema=table.schema),p)
        m.write_done(p,provenance)
    with pytest.raises(ValueError): m.merge_gate(out,2,input_dir=tmp_path,exons=ex,expected_rows=2)
    assert not (out/'j7_context.done').exists()


def test_not_available_control_and_bins():
    m=module()
    assert m.matched_control(row(110),None) is None
    assert m.matched_control(row(101,end=110),dict(tx(),exons=[[101,110]])) is None
    assert [m.distance_bin(x) for x in [None,0,-1,10,11,50,51,200,201,1000,1001]]==[
        'not_available','0','1-10','1-10','11-50','11-50','51-200','51-200','201-1000','201-1000','>1000']


def test_genuine_missing_annotation_stays_null(tmp_path):
    import pyarrow as pa
    import pyarrow.parquet as pq
    m=module(); source=tmp_path/'source.parquet'; ex=tmp_path/'ex.json'
    pq.write_table(pa.Table.from_pylist([row(105)]),source); write_j3_done(m,source,0,1)
    m.atomic_json(ex,dict(transcripts={},provenance={'test':True, 'implementation_sha256':m.implementation_hash()}))
    result=m.run_shard(source,ex,tmp_path/'out',0,1,m.rbp_audit(Path('feng/artifacts/rbp/encode_manifest.tsv')))
    r=pq.read_table(result).to_pylist()[0]
    assert r['donor_distance'] is None and r['exon_length'] is None
    assert r['geometry_status']=='not_available' and r['control_start1'] is None

@pytest.mark.parametrize('description', ['GRCh37 version 50','GRCh38 version 49'])
def test_gtf_assembly_version_fence(tmp_path,description):
    m=module(); gtf=tmp_path/'bad.gtf'; manifest=tmp_path/'m.tsv'
    manifest.write_text('t:ex1\tchr1\t1\t210\t+\n')
    gtf.write_text('##description: '+description+'\nchr1\tx\texon\t101\t110\t.\t+\t.\tgene_id "g"; transcript_id "t"; exon_number "1"; tag "Ensembl_canonical"; gene_type "protein_coding";\n')
    with pytest.raises(ValueError,match='GRCh38.*50'):
        m.build_exons(gtf,[manifest],tmp_path/'out.json')
