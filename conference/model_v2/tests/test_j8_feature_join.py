"""Synthetic contracts: coordinate orientation, coverage, conservation and restart."""
import importlib
import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq
import pytest


def api():
    name = 'conference.model_v2.j8_feature_join'
    assert importlib.util.find_spec(name), 'J8 implementation absent'
    return importlib.import_module(name)


def site(**kw):
    return dict(site_id='s', aso_id='a', chrom='chr1', start1=119, end1=122,
                strand='+', assembly='GRCh38', annotation='GENCODE v50', **kw)


def write(path, rows, metadata=None):
    t = pa.Table.from_pylist(rows).replace_schema_metadata(metadata)
    pq.write_table(t, path)
    import hashlib
    path.with_suffix('.done').write_text(json.dumps(dict(rows=len(rows), sha256=hashlib.sha256(path.read_bytes()).hexdigest())))
    return path


@pytest.mark.parametrize('strand,donor,acceptor', [('+', -9.5, 10.5), ('-', -10.5, 9.5)])
def test_geometry_signed_midpoint_and_overlap(strand, donor, acceptor):
    m = api()
    row = site(); row['strand'] = strand
    loci = [dict(locus_id='e', chrom='chr1', start1=10, end1=230, strand=strand)]
    g = m.geometry(row, loci)
    assert g['donor_distance'] == donor
    assert g['acceptor_distance'] == acceptor
    assert g['region_class'] == 'exonic'
    assert g['anchor'] == 'footprint_midpoint'
    row['start1'], row['end1'] = 108, 113
    assert m.geometry(row, loci)['region_class'] == 'exon_intron_spanning'


def test_geometry_does_not_join_other_strand_or_chrom():
    m = api()
    assert m.geometry(site(), [dict(locus_id='e', chrom='chr2', start1=10, end1=230, strand='+')])['status'] == 'not_available'


def test_j2_exact_u_and_minus_coordinate():
    m = api(); row = site(); row.update(strand='-', start1=119, end1=122)
    locus = dict(locus_id='e', chrom='chr1', start1=10, end1=230, strand='-')
    grids = [dict(locus_id='e', strand='-', pos_end_1b_transcript=112, W=80,L=40,u=4,p_unpaired=.2),
             dict(locus_id='e', strand='-', pos_end_1b_transcript=112, W=80,L=40,u=8,p_unpaired=.9)]
    result = m.accessibility(row, [locus], grids)
    assert [x['value'] for x in result] == [.2]
    assert m.accessibility(row, [locus], grids[1:])[0]['status'] == 'not_available'
    grids[0]['p_unpaired'] = float('nan')
    assert m.accessibility(row, [locus], grids[:1])[0]['value'] is None


def test_atlas_partial_overlap_not_called_exact_and_failed_not_averaged():
    m = api()
    common = dict(locus='e', chrom='chr1', strand='+', site_kind='donor', site_position=130,
                  tissue=0, config_id='c', model_evidence=True, evidence_class='sequence_ablation_proxy')
    rows = [dict(common, footprint_start=120, footprint_end=123, status='ok', raw_model_delta=-.4),
            dict(common, footprint_start=119, footprint_end=122, status='failed', raw_model_delta=None)]
    r = m.atlas_summary(site(), rows)
    assert r['exact_status'] == 'failed'
    assert r['exact_mean'] is None
    assert r['overlap_mean'] == -.4
    assert r['n_failed'] == 1
    assert r['evidence_label'] == 'sequence_ablation_proxy'


def test_stream_join_conservation_missing_failure_evidence_and_idempotence(tmp_path):
    m = api()
    j6 = write(tmp_path/'j6.parquet', [site(), dict(site(), site_id='t', start1=150,end1=153)])
    j3 = write(tmp_path/'j3.parquet', [dict(site(), duplex_status='failed', duplex_energy_kcal_mol=None)])
    loci = tmp_path/'loci.tsv'; loci.write_text('e\tchr1\t10\t230\t+\n')
    with m.FeatureIndex(tmp_path/'index.sqlite') as index:
        index.ingest('j3', j3)
        out = tmp_path/'out.parquet'
        result = m.join_shard(j6, out, index, [loci], batch_size=1)
        rows = pq.read_table(out).to_pylist()
        assert result['rows'] == len(rows) == 2
        assert json.loads(rows[0]['j3_channels'])['duplex_status'] == 'failed'
        assert json.loads(rows[1]['j3_channels'])['status'] == 'not_available'
        assert all(f.metadata.get(b'evidence_class') for f in pq.read_schema(out))
        stamp = out.stat().st_mtime_ns
        assert m.join_shard(j6, out, index, [loci])['rows'] == 2
        assert out.stat().st_mtime_ns == stamp
        assert len(result['inputs']) == 3
        out.write_bytes(b'corrupt')
        assert m.join_shard(j6, out, index, [loci])['rows'] == 2
        loci.write_text('e\tchr1\t11\t231\t+\n')
        with pytest.raises(ValueError, match='changed'):
            m.join_shard(j6, out, index, [loci])


@pytest.mark.parametrize('mutation', [{'assembly':'GRCh37'}, {'annotation':'GENCODE v49'}, {'strand':'?'}, {'start1':0}])
def test_assembly_coordinate_fence(mutation):
    m=api()
    with pytest.raises(ValueError): m.geometry(dict(site(), **mutation), [])


def test_duplicate_index_and_checksum_rejected(tmp_path):
    m=api(); p=write(tmp_path/'j3.parquet', [site(),site()])
    with m.FeatureIndex(tmp_path/'i.sqlite') as index:
        with pytest.raises(ValueError, match='duplicate'): index.ingest('j3',p)
        p.with_suffix('.done').write_text('{}')
        with pytest.raises(ValueError, match='commit'): index.ingest('j3',p)


def test_forbidden_source_and_outcome_columns(tmp_path):
    m=api()
    with pytest.raises(Exception): m.load_loci([tmp_path/'psi25'/'loci.tsv'])
    p=write(tmp_path/'bad.parquet', [dict(site(), delta_psi=.8)])
    with m.FeatureIndex(tmp_path/'i.sqlite') as index:
        with pytest.raises(ValueError, match='outcome'): index.ingest('j3',p)


def test_index_j2_j4_real_layout_and_join(tmp_path):
    m=api(); loci=tmp_path/'loci.tsv'; loci.write_text('e\tchr1\t10\t230\t+\n')
    j6=write(tmp_path/'j6.parquet',[site()])
    j2=write(tmp_path/'j2.parquet',[dict(locus_id='e',strand='+',pos_end_1b_transcript=113,W=80,L=40,u=4,p_unpaired=.3)], {b'coordinates':b'GRCh38; GENCODE v50'})
    atlas=tmp_path/'cfg'/'ablation_maps'; atlas.mkdir(parents=True)
    j4=write(atlas/'shard_0.parquet',[dict(row_id='r',locus='e',chrom='chr1',strand='+',
        site_kind='donor',site_position=130,tissue=0,footprint_start=119,footprint_end=122,
        raw_model_delta=-.2,status='ok',model_evidence=True,evidence_class='sequence_ablation_proxy')],
        {b'coordinates':b'GRCh38; GENCODE v50; 1-based inclusive'})
    j4.with_suffix('.done').rename(atlas.parent/'shard_0.done')
    with m.FeatureIndex(tmp_path/'i.sqlite') as index:
        index.ingest('j2',j2); index.ingest('j4',j4,config_id='cfg')
        index.ingest('j4',j4,config_id='cfg')
        m.join_shard(j6,tmp_path/'out.parquet',index,[loci])
    row=pq.read_table(tmp_path/'out.parquet').to_pylist()[0]
    assert json.loads(row['j2_accessibility'])[0]['value']==.3
    assert json.loads(row['j4_ablation'])['exact_mean']==-.2
    assert json.loads(row['j4_ablation'])['cells'][0]['config_id']=='cfg'


def test_index_disk_cap_rolls_back_and_bad_schema_fences(tmp_path):
    m=api(); p=write(tmp_path/'j3.parquet',[site()])
    with m.FeatureIndex(tmp_path/'i.sqlite',max_bytes=1) as index:
        with pytest.raises(ValueError,match='budget'): index.ingest('j3',p)
        assert index.identities()==[]
        assert index.db.execute('SELECT COUNT(*) FROM j3').fetchone()[0]==0


def test_j3_nonfinite_score_becomes_explicit_failure(tmp_path):
    m=api(); p=write(tmp_path/'j3.parquet',[dict(site(),duplex_status='ok',duplex_energy_kcal_mol=float('nan'))])
    with m.FeatureIndex(tmp_path/'i.sqlite') as index:
        index.ingest('j3',p)
        r=json.loads(index.db.execute('SELECT payload FROM j3').fetchone()[0])
        assert r['duplex_status']=='failed'
        assert r['duplex_energy_kcal_mol'] is None


def test_j3_bulk_sequence_columns_are_projected_away(tmp_path):
    m=api()
    assert m.J3_DROP_COLUMNS == {'target_with_flank','flank_left','flank_right','duplex_structure'}
    p=write(tmp_path/'j3.parquet',[dict(site(),target_with_flank='A'*5000,flank_left='C'*80,
        flank_right='G'*80,duplex_structure='('*5000,duplex_status='ok',duplex_energy_kcal_mol=-7.5)])
    with m.FeatureIndex(tmp_path/'i.sqlite') as index:
        index.ingest('j3',p)
        r=json.loads(index.db.execute('SELECT payload FROM j3').fetchone()[0])
    assert not m.J3_DROP_COLUMNS & set(r)
    assert r['duplex_energy_kcal_mol']==-7.5


def test_geometry_explicit_canonical_sites_excludes_terminal_boundaries():
    m=api(); loci=[dict(locus_id='e',chrom='chr1',start1=10,end1=230,strand='+',donor_position=130,acceptor_position=None)]
    r=m.geometry(site(),loci)
    assert r['acceptor_distance'] is None
    assert r['donor_distance']==-9.5


def test_j2_chunk_lookup_across_batch_boundary_and_invalid_assembly(tmp_path):
    m=api()
    rows=[dict(locus_id='e',strand='+',pos_end_1b_transcript=i,W=80,L=40,u=4,p_unpaired=i/1000) for i in range(1,601)]
    p=write(tmp_path/'grid.parquet',rows,{b'coordinates':b'GRCh38; GENCODE v50'})
    with m.FeatureIndex(tmp_path/'i.sqlite') as index:
        index.ingest('j2',p)
        assert index.db.execute('SELECT COUNT(*) FROM j2').fetchone()[0] < 10
        loc=dict(locus_id='e',chrom='chr1',start1=10,end1=800,strand='+')
        row=dict(site(),start1=519,end1=522)
        result=index.annotate(row,[loc],m.GeometryIndex([loc]))
        assert json.loads(result['j2_accessibility'])[0]['value']==.513
    bad=write(tmp_path/'bad.parquet',rows[:1],{b'coordinates':b'GRCh37; GENCODE v19'})
    with m.FeatureIndex(tmp_path/'i2.sqlite') as index:
        with pytest.raises(ValueError,match='assembly'): index.ingest('j2',bad)


def test_canonical_terminal_exons_and_gap_refusal(tmp_path):
    m=api(); p=tmp_path/'loci.tsv'
    p.write_text('tx:ex1\tchr1\t10\t230\t+\ntx:ex2\tchr1\t300\t530\t+\n')
    loci=m.load_loci([p])
    assert loci[0]['acceptor_position'] is None
    assert loci[1]['donor_position'] is None
    assert loci[1]['acceptor_position']==400
    p.write_text('tx:ex1\tchr1\t10\t230\t+\ntx:ex3\tchr1\t300\t530\t+\n')
    with pytest.raises(ValueError,match='incomplete'): m.load_loci([p])


def test_shared_catalog_reused_across_output_shards(tmp_path):
    m=api(); p=write(tmp_path/'j6.parquet',[site()]); loc=tmp_path/'loc.tsv'; loc.write_text('e\tchr1\t10\t230\t+\n')
    with m.FeatureIndex(tmp_path/'i.sqlite') as index:
        context=m.prepare_join(index,[loc],tmp_path)
        a=m.join_shard(p,tmp_path/'a.parquet',index,[loc],context=context)
        b=m.join_shard(p,tmp_path/'b.parquet',index,[loc],context=context)
    assert a['contract']['catalog_sha256']==b['contract']['catalog_sha256']
    assert (tmp_path/'feature_inputs.json').is_file()


def test_damaged_done_rebuild_and_context_rejects_late_index_mutation(tmp_path):
    m=api(); p=write(tmp_path/'j6.parquet',[site()]); loc=tmp_path/'loc.tsv'; loc.write_text('e\tchr1\t10\t230\t+\n')
    with m.FeatureIndex(tmp_path/'i.sqlite') as index:
        out=tmp_path/'a.parquet'; m.join_shard(p,out,index,[loc]); out.with_suffix('.done').write_text('{')
        assert m.join_shard(p,out,index,[loc])['rows']==1
        context=m.prepare_join(index,[loc],tmp_path)
        j3=write(tmp_path/'j3.parquet',[site()]); index.ingest('j3',j3)
        with pytest.raises(ValueError,match='changed'): m.join_shard(p,tmp_path/'b.parquet',index,[loc],context=context)


def test_cli_bounded_inventory_without_directory_scan(tmp_path):
    m=api(); root=tmp_path/'root'; root.mkdir()
    for name in ['j6_genome_v1','j3_genome_v1','j2_phase1','j2_phase2','j4_configs']:
        (root/name).mkdir()
    (root/'j6_genome_v1/run_manifest.json').write_text(json.dumps({'expected_shards':['one']}))
    write(root/'j6_genome_v1/one.parquet',[site()])
    write(root/'j3_genome_v1/j3_shard0000_of0256.parquet',[site()])
    for phase in (1,2):
        write(root/f'j2_phase{phase}/j2_shard0000.parquet',[
            dict(locus_id=f'e{phase}',strand='+',pos_end_1b_transcript=113,W=80,L=40,u=4,p_unpaired=.2)],
            {b'coordinates':b'GRCh38; GENCODE v50'})
    (root/'j4_configs/j4_configs_manifest.tsv').write_text('#config_id\tstatus\nnone\tskipped\n')
    loci=tmp_path/'loci.tsv'; loci.write_text('e1\tchr1\t10\t230\t+\n')
    out=tmp_path/'out'
    m.main(['--root',str(root),'--out',str(out),'--bounded','--loci',str(loci)])
    result=json.loads((out/'manifest.json').read_text())
    assert result['rows']==1 and result['coverage']=='partial'
    assert len(result['inputs'])==3


@pytest.mark.parametrize('kind', ['j3', 'j4', 'catalog', 'j6_manifest', 'j4_manifest', 'output_marker'])
def test_all_sidecars_reject_forbidden_symlink_before_read(tmp_path, monkeypatch, kind):
    m=api(); forbidden=tmp_path/'wetlab'; forbidden.mkdir()
    target=forbidden/'synthetic.json'; target.write_text('{}')
    original=Path.open
    def guarded_open(path, *args, **kwargs):
        assert forbidden not in path.resolve().parents, 'forbidden target was opened'
        return original(path,*args,**kwargs)
    monkeypatch.setattr(Path,'open',guarded_open)
    with m.FeatureIndex(tmp_path/'i.sqlite') as index:
        if kind in ('j3','j4'):
            p=tmp_path/'cfg'/'ablation_maps'/'s.parquet' if kind=='j4' else tmp_path/'s.parquet'
            p.parent.mkdir(parents=True,exist_ok=True)
            marker=p.parent.parent/'s.done' if kind=='j4' else p.with_suffix('.done')
            marker.symlink_to(target)
            action=lambda: m.input_identity(p,kind)
        elif kind=='catalog':
            (tmp_path/'feature_inputs.json').symlink_to(target)
            action=lambda: m.prepare_join(index,[],tmp_path)
        elif kind=='output_marker':
            p=write(tmp_path/'s.parquet',[site()]); out=tmp_path/'out.parquet'
            out.with_suffix('.done').symlink_to(target)
            action=lambda: m.join_shard(p,out,index,[])
        else:
            root=tmp_path/'root'; (root/'j6_genome_v1').mkdir(parents=True)
            (root/'j4_configs').mkdir()
            manifest=root/'j6_genome_v1/run_manifest.json'
            manifest.write_text('{"expected_shards":[]}')
            link=manifest if kind=='j6_manifest' else root/'j4_configs/j4_configs_manifest.tsv'
            link.unlink(missing_ok=True); link.symlink_to(target)
            action=lambda: m.production_inventory(root,bounded=True)
        with pytest.raises(ValueError,match='forbidden'): action()


def test_missing_j8_contract_rebuilds_even_with_matching_payload(tmp_path):
    m=api(); p=write(tmp_path/'in.parquet',[site()]); out=tmp_path/'out.parquet'
    with m.FeatureIndex(tmp_path/'i.sqlite') as index:
        context=m.prepare_join(index,[],tmp_path)
        old=m.join_shard(p,out,index,[],context=context)
        del old['contract']; out.with_suffix('.done').write_text(json.dumps(old))
        context['code_sha256']='changed-code'
        result=m.join_shard(p,out,index,[],context=context)
        assert result['contract']['code_sha256']=='changed-code'


def test_j4_full_density_benchmark_and_prefix_scaling(tmp_path, capsys):
    """Full per-locus density: 216 cells/config, 1,000 configs, 18nt spans.

    A bounded synthetic slice of the 20,395-config atlas; VM work, not flaky
    wall-clock thresholds, gates scaling. Timings/RSS are sizing evidence only.
    """
    import resource
    import time
    m=api(); geo=m.GeometryIndex([])
    common=dict(locus='e',chrom='chr1',strand='+',site_kind='donor',site_position=1,
                tissue=0,config_id='c',status='ok',model_evidence=True,raw_model_delta=-.2)
    with m.FeatureIndex(tmp_path/'bench.sqlite') as index:
        def add(lo,hi):
            with index.db:
                index.db.executemany('INSERT INTO j4 VALUES (?,?,?,?,?,?)',
                    ((f'{i}:{j}','chr1','+',i*1000+j//12+1,i*1000+j//12+18,
                      m.dumps(dict(common,footprint_start=i*1000+j//12+1,footprint_end=i*1000+j//12+18)))
                     for i in range(lo,hi) for j in range(216)))
        def work(row):
            count=0
            def tick():
                nonlocal count
                count+=1
                return 0
            index.db.set_progress_handler(tick,1)
            try: result=index.annotate(row,[],geo)
            finally: index.db.set_progress_handler(None,0)
            return count,json.loads(result['j4_ablation'])
        started=time.perf_counter(); add(0,10)
        small,_=work(dict(site(),start1=10001,end1=10018))
        add(10,1000); build=time.perf_counter()-started
        miss,summary=work(dict(site(),start1=1000001,end1=1000018))
        assert summary['n_overlap']==0
        assert miss < small*2+100, (small,miss)
        hit,summary=work(dict(site(),start1=999018,end1=999018))
        assert summary['n_overlap']==216  # both inclusive endpoints retained
        assert hit < 12000, hit
        timings={}
        for name,offset in [('dense',18),('miss',500)]:
            started=time.perf_counter()
            for i in range(300):
                pos=(i*37%1000)*1000+offset
                index.annotate(dict(site(),start1=pos,end1=pos+17),[],geo)
            timings[name+'_rows_s']=300/(time.perf_counter()-started)
        index_bytes=index.path.stat().st_size
        width_bytes=index.db.execute("SELECT SUM(pgsize) FROM dbstat WHERE name='j4_width'").fetchone()[0]
        with capsys.disabled():
            print('J4_BENCH',m.dumps(dict(rows=216000,configs=1000,build_s=build,
                index_bytes=index_bytes,width_index_bytes=width_bytes,maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                small_miss_vm=small,full_miss_vm=miss,dense_hit_vm=hit,**timings)))


def test_j4_bound_tracks_wide_new_intervals_and_restart(tmp_path):
    m=api(); path=tmp_path/'i.sqlite'; geo=m.GeometryIndex([])
    with m.FeatureIndex(path) as index:
        index.db.execute('INSERT INTO j4 VALUES (?,?,?,?,?,?)',('short','chr1','+',1,2,'{}'))
        index.db.commit()
        assert json.loads(index.annotate(dict(site(),start1=100,end1=100),[],geo)['j4_ablation'])['n_overlap']==0
        row=dict(chrom='chr1',strand='+',footprint_start=1,footprint_end=100,status='ok',model_evidence=True,raw_model_delta=.4)
        index.db.execute('INSERT INTO j4 VALUES (?,?,?,?,?,?)',('wide','chr1','+',1,100,m.dumps(row)))
        index.db.commit()
        assert json.loads(index.annotate(dict(site(),start1=100,end1=100),[],geo)['j4_ablation'])['overlap_mean']==.4
    with m.FeatureIndex(path) as index:
        assert json.loads(index.annotate(dict(site(),start1=100,end1=100),[],geo)['j4_ablation'])['n_overlap']==1


@pytest.mark.parametrize('sidecar', ['-journal','-wal','-shm'])
def test_sqlite_implicit_sidecars_forbidden(tmp_path,sidecar):
    m=api(); target=tmp_path/'wetlab'/'synthetic'; target.parent.mkdir(); target.write_bytes(b'')
    (tmp_path/('i.sqlite'+sidecar)).symlink_to(target)
    with pytest.raises(ValueError,match='forbidden'): m.FeatureIndex(tmp_path/'i.sqlite')


def test_atomic_manifest_output_forbidden(tmp_path):
    m=api(); target=tmp_path/'wetlab'/'synthetic.json'; target.parent.mkdir(); target.write_text('{}')
    output=tmp_path/'manifest.json'; output.symlink_to(target)
    with pytest.raises(ValueError,match='forbidden'): m.atomic_json(output,{'rows':1})
    assert target.read_text()=='{}'


def test_full_gate_manifest_forbidden_before_read(tmp_path,monkeypatch):
    m=api(); root=tmp_path/'root'; (root/'j6_genome_v1').mkdir(parents=True)
    target=tmp_path/'wetlab'/'synthetic.json'; target.parent.mkdir(); target.write_text('{}')
    (root/'j6_genome_v1/gate_report.json').symlink_to(target)
    loci=tmp_path/'empty.tsv'; loci.write_text('')
    monkeypatch.setattr(m,'production_inventory',lambda *a,**kw:([],[],[]))
    original=Path.open
    def guarded_open(path,*a,**kw):
        assert target.parent not in path.resolve().parents, 'forbidden target was opened'
        return original(path,*a,**kw)
    monkeypatch.setattr(Path,'open',guarded_open)
    with pytest.raises(ValueError,match='forbidden'):
        m.main(['--root',str(root),'--out',str(tmp_path/'out'),'--loci',str(loci)])
