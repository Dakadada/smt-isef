"""Rank/null/ablation contracts with independent hand-checked expectations."""
import importlib
from collections import Counter
import pytest


def api():
    name='conference.model_v2.j10_robustness'
    assert importlib.util.find_spec(name), 'J10 implementation absent'
    return importlib.import_module(name)


def test_rank_reversal_ties_and_missing_are_explicit():
    m=api()
    r=m.rank_stability({'a':{'x':1,'y':2,'z':3},'b':{'x':3,'y':2,'z':1}})
    assert r['comparisons'][0]['spearman'] == pytest.approx(-1)
    assert r['comparisons'][0]['kendall'] == pytest.approx(-1)
    assert r['sign_flip_frequency'] == 1
    assert r['stability_flag'] == 'unstable'
    tied=m.rank_stability({'a':{'x':1,'y':1},'b':{'x':2,'y':2}})
    assert tied['stability_flag'] == 'indeterminate'
    assert tied['comparisons'][0]['kendall'] is None
    missing=m.rank_stability({'a':{'x':1,'y':None},'b':{'x':1,'y':2}})
    assert missing['stability_flag'] == 'indeterminate'
    assert missing['missing']['a'] == ['y']


def test_alternatives_fixed_candidates_and_direction():
    m=api(); rows=[dict(candidate_id='x',variant_id='w80',wc=3,energy=-4,mismatch_count=0),
                  dict(candidate_id='y',variant_id='w80',wc=1,energy=-2,mismatch_count=2)]
    settings=[m.Setting('wc','wc','w80',True,3),m.Setting('energy','energy','w80',False,3),m.Setting('k0','wc','w80',True,0)]
    scores=m.alternative_scores(rows,settings,['x','y'])
    assert scores['wc']['x'] > scores['wc']['y']
    assert scores['energy']['x'] > scores['energy']['y']
    assert scores['k0']['y'] is None
    with pytest.raises(ValueError): m.alternative_scores(rows+rows,settings,['x','y'])


@pytest.mark.parametrize('seq', ['ACGACGTTAC','AAAA','A','ACACACAC','GGCCAATTGGCCAATT'])
def test_dinucleotide_shuffle_preserves_edges_length_gc_and_seed(seq):
    m=api(); shuffled=m.dinucleotide_shuffle(seq,seed=42)
    assert len(shuffled)==len(seq)
    assert Counter(zip(shuffled,shuffled[1:]))==Counter(zip(seq,seq[1:]))
    assert Counter(shuffled)==Counter(seq)
    assert shuffled==m.dinucleotide_shuffle(seq,seed=42)


def test_decoys_same_locus_exclude_real_match_gc_and_deterministic():
    m=api(); r=m.shifted_decoys('ACGTACGTACGT', start=4,length=4,n=20,seed=3,gc_tolerance=0)
    assert r['sampling_universe']['eligible_count']==8
    assert all(d['start'] != 4 and len(d['sequence'])==4 for d in r['decoys'])
    assert len({d['start'] for d in r['decoys']})==8
    assert r==m.shifted_decoys('ACGTACGTACGT',start=4,length=4,n=20,seed=3,gc_tolerance=0)


def test_null_pvalue_blocked_without_universe_or_minimum_or_failed_score():
    m=api()
    for universe,values in [(None,[1]*30),({'definition':'fixed decoys'},[1,2]),({'definition':'fixed decoys'},[None]*30)]:
        with pytest.raises(m.MetricBlocked): m.empirical_null(3,values,universe=universe,min_decoys=20)
    r=m.empirical_null(3,[1,2,3,4],universe={'definition':'four fixed matched decoys'},min_decoys=4)
    assert r['p_value']==pytest.approx(3/5)
    assert r['distribution']==[1,2,3,4]


def test_nulls_do_not_invent_uncomputed_channel_values():
    m=api()
    result=m.null_channels('ACGTACGT',start=2,length=4,n=25,seed=2,
        scorers={'gc':lambda seq:sum(x in 'GC' for x in seq),'duplex':None},min_decoys=20)
    assert result['duplex']['status']=='blocked'
    assert result['gc']['shuffle']['sampling_universe']['method']=='randomized_euler_trail_nonuniform'
    assert result['gc']['shuffle']['p_value'] is None


def test_ablation_calls_core_guards_and_missing_features_stay_blocked():
    m=api()
    features={'a':{'m1_wc_pairs':1},'b':{'m1_wc_pairs':2},'c':{'m1_wc_pairs':3}}
    r=m.baseline_ablation(features, source_paths=['derived.parquet'], seed=8)
    assert set(r['baselines'])=={'M0','M1','M2','M3','M4','M5'}
    assert r['supervised']['status']=='blocked'
    assert 'negative' in r['supervised']['reason']
    assert r['baselines']['M1']['status']=='blocked'
    assert r['mode']=='descriptive_rank_agreement_only'


def test_complete_ablation_descriptive_leave_out():
    m=api()
    features={x:{c:(i+1)*(-1 if 'energy' in c else 1) for c in m.ec.baseline_feature_columns('M5')} for i,x in enumerate('abc')}
    r=m.baseline_ablation(features,source_paths=['derived.parquet'],seed=1)
    assert r['baselines']['M5']['status']=='ok'
    assert len(r['baselines']['M5']['leave_out'])==5
    assert r['baselines']['M1']['leave_out']['complementarity']['status']=='blocked'
    with pytest.raises(m.ec.ForbiddenSourceError): m.baseline_ablation(features, source_paths=['psi25/outcomes.csv'])


def test_j9_adapter_preserves_seed_and_uses_mm3_counts():
    m=api(); rows=[]
    for c,i in zip(['F18','F14','F10'],[3,2,1]):
        for col in m.ec.baseline_feature_columns('M1'):
            rows.append(dict(aso_id=c,feature_name=col,variant_id='mm3' if 'n_sites' in col or 'n_exact_sites' in col else 'default',
                             value_double=i,value_int=None,status='ok'))
    rows.append(dict(aso_id='3UP8',feature_name='forbidden',variant_id='default',status='ok',value_int=10))
    r=m.j9_report(rows,source='derived.parquet',seed=99)
    assert r['ablation']['baselines']['M1']['status']=='ok'
    assert r['ablation']['seed']==99
    assert 'forbidden' not in r['per_feature']


def test_rank_stable_and_missing_alternative_and_invalid_inputs():
    m=api()
    r=m.rank_stability({'a':{'x':1,'y':2,'z':2},'b':{'x':2,'y':4,'z':4}})
    assert r['stability_flag']=='stable'
    assert r['n_tied_comparisons']==1
    assert m.rank_stability({'one':{'a':1}})['stability_flag']=='indeterminate'
    with pytest.raises(ValueError): m.Setting('a','c',max_mismatch=-1)
    with pytest.raises(ValueError): m.dinucleotide_shuffle('AN',seed=1)
    with pytest.raises(ValueError): m.shifted_decoys('ACGT',start=4,length=2,n=3,seed=1)


def test_null_lower_tail_and_failed_callback():
    m=api()
    assert m.empirical_null(2,[1,2,3],universe={'definition':'three matched decoys'},min_decoys=3,higher_is_extreme=False)['p_value']==.75
    r=m.null_channels('ACGTACGT',start=0,length=4,n=4,seed=1,scorers={'missing':lambda s:None},min_decoys=2)
    assert r['missing']['shift']['p_value'] is None


def test_cli_is_bounded_excludes_3up8_and_restarts(tmp_path):
    import json
    import pyarrow as pa
    import pyarrow.parquet as pq
    m=api(); j9=tmp_path/'j9.parquet'; j3=tmp_path/'j3.parquet'; out=tmp_path/'result.json'
    pq.write_table(pa.Table.from_pylist([dict(aso_id=c,feature_name='m1_wc_pairs',variant_id='default',value_int=i,value_double=None,status='ok')
                        for i,c in enumerate(['F18','F14','F10','3UP8'])]),j9)
    base=dict(assembly='GRCh38',annotation='GENCODE v50',chrom='chr1',strand='+',start1=10,end1=13,mismatch_count=0,aso_sequence='ACGT',
              complementarity_wc_pairs=4,complementarity_status='ok',duplex_energy_kcal_mol=-4.,duplex_status='ok',
              accessibility_p_unpaired=.1,accessibility_status='ok',target_with_flank='ACGTACGTACGT',footprint_start1_in_context=5)
    pq.write_table(pa.Table.from_pylist([dict(base,site_id=f's{i}',aso_id=c) for i,c in enumerate(['3UP8','F18','F14','F10'])]),j3)
    j3.with_suffix('.done').write_text(json.dumps(dict(rows=4,params=dict(window=120,span=60))))
    args=['--j9',str(j9),'--j3',str(j3),'--out',str(out),'--limit','2']
    m.main(args); stamp=out.stat().st_mtime_ns; m.main(args)
    r=json.loads(out.read_text())
    assert r['j3_sample']['rows']==2
    assert r['j9']['excluded']==['3UP8']
    assert out.stat().st_mtime_ns==stamp
    assert all(n['site_id']!='s0' for n in r['j3_sample']['nulls'])
    assert all(n['observed_wc_pairs']==4 for n in r['j3_sample']['nulls'])
    assert r['j3_sample']['nulls'][0]['channels']['duplex']['status']=='blocked'
    # The one computed window is labelled from the marker, never assumed; the other
    # two declared J2 windows stay present-but-missing rather than silently reused.
    assert r['contract']['j3_fold_variant']=='w120s60'
    assert r['j3_sample']['j3_fold_variant']=='w120s60'
    absent={s['name'] for s in r['j3_sample']['settings']}&{'w80s40','w200s100'}
    assert absent=={'w80s40','w200s100'}
    assert all(r['j3_sample']['stability']['missing'][n]==['s1','s2'] for n in absent)


def test_complementarity_scorer_matches_j3_and_feeds_a_real_null():
    m=api()
    from conference.model_v2.j3_scoring import complementarity
    aso='ACGUACGU'; target='ACGUACGU'
    scorer=m.complementarity_scorer(aso)
    assert scorer(target)==float(complementarity(aso,target)['wc_pairs'])
    assert scorer('ACGTACGT')==scorer(target)     # DNA/RNA alphabets agree
    with pytest.raises(ValueError): scorer('AAAA')  # unequal lengths refuse, never pad
    with pytest.raises(ValueError): m.complementarity_scorer('')
    # A real (not blocked) shift null: uniform-without-replacement universe, p-value emitted.
    r=m.null_channels('ACGTACGTTGCAACGTACGTTGCA',start=0,length=8,n=12,seed=5,
                      scorers={'complementarity':scorer},min_decoys=3,gc_tolerance=1)
    assert r['complementarity']['shift']['status']=='ok'
    assert 0 < r['complementarity']['shift']['p_value'] <= 1
    assert r['complementarity']['shift']['n_decoys']>=3
    assert r['complementarity']['shuffle']['status']=='blocked'


def test_fold_variant_refuses_unlabelled_or_undeclared_marker(tmp_path):
    import json
    m=api(); p=tmp_path/'j3.parquet'; p.write_bytes(b'x')
    done=p.with_suffix('.done')
    done.write_text(json.dumps({'rows': 1}))
    with pytest.raises(ValueError, match='lacks explicit fold window'): m.j3_fold_variant(p)
    done.write_text(json.dumps({'params': {'window': 70, 'span': 35}}))
    with pytest.raises(ValueError, match='outside the declared J2 variant set'): m.j3_fold_variant(p)
    done.write_text(json.dumps({'params': {'window': 200, 'span': 100}}))
    assert m.j3_fold_variant(p)=='w200s100'
