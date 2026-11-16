"""J10: descriptive rank robustness and declared-universe perturbation nulls.

No fitting, labels, outcomes, Vienna calls or GPU inference. Alternative columns
must already be computed. Missing scores remain missing on a fixed candidate set.
The default descriptive composite is the mean of oriented within-candidate-set
percentile ranks: a model assumption, not an efficacy score. Supervised requests
always pass through evaluation_core and are blocked for these unlabeled tables.
"""
from __future__ import annotations

import argparse
import hashlib
from dataclasses import asdict, dataclass
from itertools import combinations
import json
import math
from pathlib import Path
import random
import warnings

from scipy.stats import kendalltau, spearmanr, rankdata
from conference.model_v2 import evaluation_core as ec
from conference.model_v2.j3_scoring import complementarity as j3_complementarity
from conference.model_v2.j8_feature_join import atomic_json, safe_path, sha256, input_identity, dumps

MetricBlocked=ec.MetricBlocked
# The fold windows J2 actually computed (jobs/j2_fold.py windows=((80,40),(120,60),(200,100))).
FOLD_VARIANTS=('w80s40','w120s60','w200s100')


def j3_fold_variant(path):
    """Read the single J3 fold window from its commit marker; never assume a label."""
    marker=safe_path(safe_path(path).with_suffix('.done'))
    params=json.loads(marker.read_text()).get('params',{})
    window,span=params.get('window'),params.get('span')
    if type(window) is not int or type(span) is not int:
        raise ValueError('J3 marker lacks explicit fold window/span; cannot label the accessibility variant')
    label='w%ds%d'%(window,span)
    if label not in FOLD_VARIANTS: raise ValueError('J3 fold window is outside the declared J2 variant set: '+label)
    return label


def finite(value):
    return isinstance(value,(int,float)) and not isinstance(value,bool) and math.isfinite(value)


@dataclass(frozen=True)
class Setting:
    name: str
    channel: str
    variant_id: str='default'
    higher_is_better: bool=True
    max_mismatch: int | None=None

    def __post_init__(self):
        if not self.name or not self.channel or not self.variant_id: raise ValueError('setting names required')
        if self.max_mismatch is not None and (type(self.max_mismatch) is not int or self.max_mismatch<0):
            raise ValueError('nonnegative max_mismatch required')


def alternative_scores(rows,settings,candidates):
    """Input wide rows keyed by (candidate_id, variant_id); no variant fallback."""
    if len(candidates)!=len(set(candidates)) or not candidates: raise ValueError('unique nonempty candidates required')
    if len(settings)!=len({s.name for s in settings}): raise ValueError('duplicate settings')
    index={}
    for row in rows:
        key=(row['candidate_id'],row.get('variant_id','default'))
        if key in index: raise ValueError('duplicate candidate/variant row')
        index[key]=row
    result={}
    for setting in settings:
        scores={}
        for candidate in candidates:
            r=index.get((candidate,setting.variant_id),{})
            v=r.get(setting.channel)
            status=r.get(setting.channel+'_status',r.get('status','ok'))
            mm=r.get('mismatch_count')
            eligible=setting.max_mismatch is None or (finite(mm) and mm<=setting.max_mismatch)
            scores[candidate]=(v if setting.higher_is_better else -v) if finite(v) and status=='ok' and eligible else None
        result[setting.name]=scores
    return result


def rank_stability(alternatives,*,candidates=None,min_agreement=.8,max_flip_fraction=0):
    """Tau-b / Spearman on pairwise complete cases, with full-set coverage guard.

    Sign flips compare strict pair orderings across every pair of settings. Ties
    are not flips and are reported separately. Incomplete coverage cannot be
    declared stable; an observed reversal can still establish instability.
    """
    if not -1<=min_agreement<=1 or not 0<=max_flip_fraction<=1: raise ValueError('invalid stability thresholds')
    names=sorted(alternatives)
    ids=sorted(candidates if candidates is not None else {c for scores in alternatives.values() for c in scores})
    if len(ids)!=len(set(ids)): raise ValueError('duplicate candidates')
    missing={n:[c for c in ids if not finite(alternatives[n].get(c))] for n in names}
    comparisons=[]; flips=0; strict=0; ties=0
    for a,b in combinations(names,2):
        common=[c for c in ids if c not in missing[a] and c not in missing[b]]
        x=[alternatives[a][c] for c in common]; y=[alternatives[b][c] for c in common]
        tau=rho=None
        if len(common)>=2 and len(set(x))>1 and len(set(y))>1:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                tau=float(kendalltau(x,y).statistic); rho=float(spearmanr(x,y).statistic)
            tau=tau if finite(tau) else None; rho=rho if finite(rho) else None
        f=s=t=0
        for i,j in combinations(range(len(common)),2):
            dx=x[i]-x[j]; dy=y[i]-y[j]
            if dx==0 or dy==0: t+=1
            else: s+=1; f+=int(dx*dy<0)
        flips+=f; strict+=s; ties+=t
        comparisons.append(dict(settings=[a,b],n_common=len(common),kendall=tau,spearman=rho,
                                sign_flips=f,strict_order_comparisons=s,tied_order_comparisons=t,
                                sign_flip_frequency=f/s if s else None,evidence_class='derived'))
    frequency=flips/strict if strict else None
    known=[r for r in comparisons if r['kendall'] is not None and r['spearman'] is not None]
    unstable=any(min(r['kendall'],r['spearman'])<min_agreement for r in known) or (frequency is not None and frequency>max_flip_fraction)
    complete=bool(comparisons) and not any(missing.values()) and len(known)==len(comparisons)
    flag='unstable' if unstable else ('stable' if complete else 'indeterminate')
    return dict(stability_flag=flag,n_candidates=len(ids),missing=missing,comparisons=comparisons,
                sign_flip_frequency=frequency,n_sign_flips=flips,n_strict_comparisons=strict,n_tied_comparisons=ties,
                thresholds=dict(min_agreement=min_agreement,max_flip_fraction=max_flip_fraction),evidence_class='derived')


def dna(sequence):
    value=sequence.upper().replace('U','T')
    if not value or set(value)-set('ACGT'): raise ValueError('nonempty unambiguous sequence required')
    return value


def dinucleotide_shuffle(sequence,*,seed):
    """Randomized Euler trail; exact dinucleotide multiset, NOT a uniform sampler."""
    sequence=dna(sequence); rng=random.Random(seed); edges={b:[] for b in 'ACGT'}
    for a,b in zip(sequence,sequence[1:]): edges[a].append(b)
    for values in edges.values(): rng.shuffle(values)
    stack=[sequence[0]]; trail=[]
    while stack:
        if edges[stack[-1]]: stack.append(edges[stack[-1]].pop())
        else: trail.append(stack.pop())
    return ''.join(reversed(trail))


def shifted_decoys(sequence,*,start,length,n,seed,gc_tolerance=0):
    """Uniform without replacement over eligible 0-based starts in this locus.

    Excludes the true start (overlap with the true footprint is allowed and stated).
    GC tolerance is an absolute fraction difference, with exact length matching.
    """
    sequence=dna(sequence)
    if any(type(x) is not int for x in (start,length,n)) or not 0<=start<=len(sequence)-length or length<1 or n<1:
        raise ValueError('invalid footprint or draw count')
    if not 0<=gc_tolerance<=1: raise ValueError('invalid GC tolerance')
    gc=lambda s:sum(b in 'GC' for b in s)/len(s)
    target=sequence[start:start+length]; target_gc=gc(target)
    eligible=[i for i in range(len(sequence)-length+1) if i!=start and abs(gc(sequence[i:i+length])-target_gc)<=gc_tolerance+1e-12]
    starts=sorted(random.Random(seed).sample(eligible,min(n,len(eligible))))
    return dict(decoys=[dict(start=i,sequence=sequence[i:i+length]) for i in starts],
                sampling_universe=dict(definition='all same-locus starts except observed start, exact length, GC-matched; overlap allowed',
                    method='uniform_without_replacement',locus_sequence_sha256=__import__('hashlib').sha256(sequence.encode()).hexdigest(),
                    locus_length=len(sequence),observed_start=start,length=length,target_gc=target_gc,
                    gc_tolerance=gc_tolerance,eligible_count=len(eligible),eligible_starts=eligible),seed=seed)


def empirical_null(observed,values,*,universe,min_decoys=20,higher_is_extreme=True):
    """One-sided plus-one empirical tail; refuses missing scores or undefined universe."""
    if not isinstance(universe,dict) or not isinstance(universe.get('definition'),str) or not universe['definition'].strip():
        raise MetricBlocked('empirical_null','sampling universe is not defined')
    if type(min_decoys) is not int or min_decoys<1: raise ValueError('positive minimum decoy count required')
    if universe.get('method')=='randomized_euler_trail_nonuniform':
        raise MetricBlocked('empirical_null','Euler-trail draw probabilities are not uniform; descriptive null only')
    if universe.get('method')=='uniform_without_replacement':
        raise MetricBlocked('empirical_null','J6 discovery selection is not modeled: random positions are not '
                            'exchangeable with discovery-selected near-matches; descriptive distribution only, '
                            'no prespecified selection-aware null is justified')
    if not finite(observed) or any(not finite(v) for v in values):
        raise MetricBlocked('empirical_null','missing/failed observed or decoy scores; no selective deletion')
    if len(values)<min_decoys: raise MetricBlocked('empirical_null',f'n_decoys={len(values)} below minimum={min_decoys}')
    extreme=sum(v>=observed if higher_is_extreme else v<=observed for v in values)
    return dict(distribution=list(values),observed=observed,n_decoys=len(values),minimum_decoys=min_decoys,
                p_value=(1+extreme)/(1+len(values)),tail='upper' if higher_is_extreme else 'lower',
                sampling_universe=universe,evidence_class='derived',interpretation='conditional matched-decoy tail; not biological efficacy evidence')


def null_channels(sequence,*,start,length,n,seed,scorers,min_decoys=20,gc_tolerance=0,directions=None):
    """Callbacks are lookups/statistical descriptors only; absent channels are blocked.

    Entire target flank is dinucleotide shuffled. Only draws also matching the
    footprint GC are retained; acceptance counts and the nonuniform law are stated.
    """
    sequence=dna(sequence); target=sequence[start:start+length]
    shifted=shifted_decoys(sequence,start=start,length=length,n=n,seed=seed,gc_tolerance=gc_tolerance)
    target_gc=sum(b in 'GC' for b in target)/length
    shuffled=[]
    for i in range(n):
        s=dinucleotide_shuffle(sequence,seed=seed+i)[start:start+length]
        if abs(sum(b in 'GC' for b in s)/length-target_gc)<=gc_tolerance+1e-12: shuffled.append(s)
    universe=dict(definition='dinucleotide-preserving permutations of the supplied flank, conditioned on footprint GC/length',
                  method='randomized_euler_trail_nonuniform',locus_sequence_sha256=shifted['sampling_universe']['locus_sequence_sha256'],
                  length=length,gc_tolerance=gc_tolerance,n_attempted=n,n_accepted=len(shuffled),n_unique=len(set(shuffled)))
    result={}
    for channel,scorer in sorted(scorers.items()):
        if scorer is None:
            result[channel]=dict(status='blocked',p_value=None,reason='decoy channel not precomputed; no scoring backend in J10',evidence_class='not_available'); continue
        def score(s):
            try:
                v=scorer(s)
            except (KeyError,ValueError,RuntimeError): return None
            return v if finite(v) else None
        observed=score(target); parts={}
        for name,sequences,u in [('shift',[d['sequence'] for d in shifted['decoys']],shifted['sampling_universe']),('shuffle',shuffled,universe)]:
            values=[score(s) for s in sequences]
            try:
                parts[name]=dict(status='ok',**empirical_null(observed,values,universe=u,min_decoys=min_decoys,
                                                   higher_is_extreme=(directions or {}).get(channel,True)))
            except MetricBlocked as exc:
                parts[name]=dict(status='blocked',reason=exc.reason,p_value=None,distribution=values,observed=observed,
                                 n_decoys=len(values),minimum_decoys=min_decoys,sampling_universe=u,evidence_class='derived')
        result[channel]=parts
    return result


def complementarity_scorer(aso_sequence):
    """Reuse the J3 descriptor itself: pure sequence arithmetic, no ViennaRNA/GPU/IO.

    This is the only channel a decoy can be rescored for inside J10. Duplex energy
    and accessibility need a folding backend, so their decoy scorers stay None and
    their nulls are reported blocked rather than approximated.
    """
    if not isinstance(aso_sequence,str) or not aso_sequence.strip(): raise ValueError('ASO sequence required')
    def score(sequence):
        return float(j3_complementarity(aso_sequence,sequence)['wc_pairs'])
    return score


def _direction(column):
    return -1 if any(x in column for x in ('energy','mismatch','abs_','n_sites','n_exact_sites')) else 1


def descriptive_scores(features,columns):
    """Prespecified equal-weight oriented rank composite; all cells must be present."""
    ids=sorted(features)
    if not ids or not columns: raise MetricBlocked('rank_composite','no feature columns / candidates')
    missing=[(c,k) for c in ids for k in columns if not finite(features[c].get(k))]
    if missing: raise MetricBlocked('rank_composite',f'{len(missing)} missing feature cells, e.g. {missing[:3]}')
    scores={c:0. for c in ids}
    for k in columns:
        ranks=rankdata([_direction(k)*features[c][k] for c in ids],method='average')
        for c,r in zip(ids,ranks): scores[c]+=float(r)/len(ids)/len(columns)
    return scores


def baseline_ablation(features,*,source_paths,seed=0):
    """Use evaluation_core registry, lineage guard, M0 and supervised metric guard.

    No fake split groups or fitted ablation are created: these candidates have no
    verified negatives. Leave-outs here compare descriptive fixed composites.
    """
    if not source_paths: raise ValueError('feature provenance required')
    for p in source_paths: ec.check_source_not_forbidden(p)
    rows=[dict(observation_id=c,label_state='unlabeled',assay_id='unassayed',study_id='derived_feature_atlas',
               independent_unit_id=c,provenance={'source':str(source_paths[0]),'source_paths':list(map(str,source_paths))}) for c in sorted(features)]
    table=ec.validate_table(rows)
    random_scores=ec.permuted_reference_scores(table.ids,seed)
    try:
        ec.compute_metric('auroc',table,random_scores,mode='pu')
    except MetricBlocked as exc:
        supervised=dict(status='blocked',reason='No verified negatives or measured labels. '+exc.reason)
    else: raise RuntimeError('unlabeled supervised metric unexpectedly passed evaluation_core')
    baselines={}
    for name in ('M0','M1','M2','M3','M4','M5'):
        columns=ec.baseline_feature_columns(name)
        if name=='M0':
            baselines[name]=dict(status='ok',scores=random_scores,evidence_class='model_assumption',seed=seed); continue
        try:
            full=descriptive_scores(features,columns)
        except MetricBlocked as exc:
            baselines[name]=dict(status='blocked',reason=exc.reason,evidence_class='not_available'); continue
        leave={}
        for family in ec.BASELINES[name].families:
            kept=tuple(c for c in columns if c not in family.columns)
            try:
                scores=descriptive_scores(features,kept)
                leave[family.name]=dict(status='ok',agreement=rank_stability({'full':full,'leave_out':scores}),columns=list(kept))
            except MetricBlocked as exc:
                leave[family.name]=dict(status='blocked',reason=exc.reason)
        baselines[name]=dict(status='ok',scores=full,leave_out=leave,evidence_class='model_assumption')
    return dict(mode='descriptive_rank_agreement_only',supervised=supervised,baselines=baselines,
                composite='equal-weight oriented percentile ranks; fixed assumption, not fitted',
                ablation_fit='not run: no supervised training/splits; evaluation_core registry and metric guards called',seed=seed)


def j9_report(rows,*,source,seed=0,fold_variant='w120s60',mismatch_variant='mm3'):
    candidates=['F18','F14','F10']; scores={}; features={c:{} for c in candidates}
    baseline_variants=('default',fold_variant,mismatch_variant)
    for row in rows:
        if row['aso_id'] not in candidates: continue
        value=row.get('value_double') if row.get('value_double') is not None else row.get('value_int')
        key=(row['feature_name'],row['variant_id'])
        scores.setdefault(key,{})
        if row['aso_id'] in scores[key]: raise ValueError('duplicate J9 feature cell')
        scores[key][row['aso_id']]=value if row['status']=='ok' and finite(value) else None
        if row['variant_id'] in baseline_variants:
            if row['feature_name'] in features[row['aso_id']]: raise ValueError('ambiguous baseline feature alternative')
            features[row['aso_id']][row['feature_name']]=scores[key][row['aso_id']]
    reports={}
    for feature in sorted({k[0] for k in scores}):
        alternatives={v:{c:(_direction(feature)*values[c] if finite(values.get(c)) else None) for c in candidates}
                      for (f,v),values in scores.items() if f==feature}
        reports[feature]=rank_stability(alternatives,candidates=candidates)
    # Compare the prespecified thermodynamic/complementarity/accessibility channels.
    channels={f:{c:(_direction(f)*scores.get((f,v),{}).get(c) if finite(scores.get((f,v),{}).get(c)) else None) for c in candidates}
              for f,v in [('m1_wc_pairs','default'),('m2_duplex_energy_kcal_mol',fold_variant),('m3_accessibility_p_unpaired',fold_variant)]}
    return dict(candidates=candidates,excluded=['3UP8'],per_feature=reports,channel_agreement=rank_stability(channels),
                baseline_variants=list(baseline_variants),
                ablation=baseline_ablation(features,source_paths=[source],seed=seed),evidence_class='derived')


def j9_identity(path):
    """Authenticate the feature output against J9's published commit manifest."""
    import pyarrow.parquet as pq
    path=safe_path(path); marker=safe_path(path.parent/'manifest.json')
    manifest=json.loads(marker.read_text())
    entries=[e for e in manifest.get('outputs',[]) if e.get('name')==path.name]
    digest=sha256(path); rows=pq.read_metadata(path).num_rows
    if len(entries)!=1 or entries[0].get('rows')!=rows or entries[0].get('sha256')!=digest:
        raise ValueError('invalid J9 input commit manifest')
    return dict(path=str(path),sha256=digest,rows=rows,kind='j9',
                marker_path=str(marker),marker_sha256=sha256(marker))


def site_null(row,seed):
    """Keep unavailable sites in rank coverage while blocking their decoy nulls."""
    length=row['end1']-row['start1']+1
    try:
        if row.get('complementarity_status')!='ok' or not finite(row.get('complementarity_wc_pairs')):
            raise ValueError('upstream complementarity unavailable')
        sequence=row.get('target_with_flank'); offset=row.get('footprint_start1_in_context')
        if not isinstance(sequence,str) or type(offset) is not int or not 1<=offset<=len(sequence)-length+1:
            raise ValueError('upstream context/footprint offset unavailable or invalid')
        sequence=dna(sequence); start=offset-1
        scorer=complementarity_scorer(row.get('aso_sequence'))
        observed=scorer(sequence[start:start+length])
    except (ValueError,TypeError) as exc:
        reasons=[str(row[k]) for k in ('context_reason','complementarity_reason','duplex_reason','accessibility_reason') if row.get(k)]
        reason='; '.join([str(exc)]+list(dict.fromkeys(reasons)))
        return dict(site_id=row['site_id'],aso_id=row['aso_id'],observed_wc_pairs=None,
                    channels={c:dict(status='blocked',p_value=None,reason=reason,evidence_class='not_available')
                              for c in ('complementarity','duplex','accessibility')})
    if observed!=float(row['complementarity_wc_pairs']):
        raise ValueError('decoy scorer disagrees with the committed J3 complementarity cell: '+row['site_id'])
    return dict(site_id=row['site_id'],aso_id=row['aso_id'],observed_wc_pairs=observed,
                channels=null_channels(sequence,start=start,length=length,n=50,seed=seed,
                         scorers={'complementarity':scorer,'duplex':None,'accessibility':None}))


def payload_sha256(result):
    return hashlib.sha256(dumps({k:v for k,v in result.items() if k!='payload_sha256'}).encode()).hexdigest()


def valid_checkpoint(result,contract):
    """Require a current complete result schema and an intact canonical payload."""
    def fields(obj,schema):
        return isinstance(obj,dict) and all(k in obj and type(obj[k]) is t for k,t in schema.items())
    def stability(obj):
        return fields(obj,dict(stability_flag=str,n_candidates=int,missing=dict,comparisons=list,
                               n_sign_flips=int,n_strict_comparisons=int,n_tied_comparisons=int,thresholds=dict))
    try:
        if not fields(result,dict(contract=dict,j9=dict,j3_sample=dict,blocked=list,payload_sha256=str)): return False
        if result['contract']!=contract or result['payload_sha256']!=payload_sha256(result): return False
        j9=result['j9']; sample=result['j3_sample']
        if not fields(j9,dict(candidates=list,excluded=list,per_feature=dict,channel_agreement=dict,
                              baseline_variants=list,ablation=dict,evidence_class=str)): return False
        if j9['candidates']!=['F18','F14','F10'] or j9['excluded']!=['3UP8']: return False
        if not stability(j9['channel_agreement']) or not all(stability(x) for x in j9['per_feature'].values()): return False
        if not fields(j9['ablation'],dict(mode=str,supervised=dict,baselines=dict)): return False
        if set(j9['ablation']['baselines'])!={'M0','M1','M2','M3','M4','M5'}: return False
        if not fields(sample,dict(rows=int,sampling_universe=str,j3_fold_variant=str,settings=list,stability=dict,nulls=list)): return False
        if not 1<=sample['rows']<=contract['limit'] or sample['j3_fold_variant']!=contract['j3_fold_variant']: return False
        if not stability(sample['stability']) or sample['stability']['n_candidates']!=sample['rows']: return False
        if len(sample['settings'])!=9 or not all(fields(s,dict(name=str,channel=str,variant_id=str,higher_is_better=bool)) for s in sample['settings']): return False
        if len(sample['nulls'])!=min(3,sample['rows']): return False
        for null in sample['nulls']:
            if not fields(null,dict(site_id=str,aso_id=str,channels=dict)) or 'observed_wc_pairs' not in null: return False
            if set(null['channels'])!={'complementarity','duplex','accessibility'}: return False
            for channel in null['channels'].values():
                if not isinstance(channel,dict): return False
                parts=[channel] if 'status' in channel else [channel.get('shift'),channel.get('shuffle')]
                if not all(fields(p,dict(status=str,reason=str)) and p['status']=='blocked' and 'p_value' in p and p['p_value'] is None for p in parts): return False
                if 'status' not in channel:
                    for part in parts:
                        if not fields(part,dict(distribution=list,n_decoys=int,minimum_decoys=int,sampling_universe=dict,evidence_class=str)): return False
                        if 'observed' not in part or len(part['distribution'])!=part['n_decoys']: return False
        return bool(result['blocked']) and all(isinstance(x,str) for x in result['blocked'])
    except (ValueError,TypeError,KeyError):
        return False


def main(argv=None):
    import pyarrow.parquet as pq
    from itertools import islice
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--j9',required=True); p.add_argument('--j3',required=True); p.add_argument('--out',required=True)
    p.add_argument('--limit',type=int,default=300); p.add_argument('--seed',type=int,default=1729)
    args=p.parse_args(argv)
    if not 2<=args.limit<=1000: p.error('bounded sample limit must be 2..1000')
    paths=[safe_path(args.j9),safe_path(args.j3)]
    identities=[j9_identity(paths[0]),input_identity(paths[1],'j3')]
    variant=j3_fold_variant(paths[1])
    # Bounded output checkpoint; exact source SHA is acquired in Slurm, not login.
    contract=dict(inputs=identities,seed=args.seed,limit=args.limit,
                  j3_fold_variant=variant,code_sha256=sha256(__file__))
    output=safe_path(args.out); output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists():
        try: prior=json.loads(output.read_text())
        except (ValueError,OSError): prior={}
        if valid_checkpoint(prior,contract): return
        if isinstance(prior,dict) and prior.get('contract') is not None and prior['contract']!=contract:
            raise ValueError('changed robustness inputs; choose fresh output')
    j9rows=pq.read_table(paths[0]).to_pylist()
    columns=['site_id','aso_id','aso_sequence','assembly','annotation','chrom','strand','start1','end1','mismatch_count',
             'complementarity_wc_pairs','complementarity_status','duplex_energy_kcal_mol','duplex_status',
             'accessibility_p_unpaired','accessibility_status','target_with_flank','footprint_start1_in_context']
    pf=pq.ParquetFile(paths[1])
    columns += [k for k in ('context_reason','complementarity_reason','duplex_reason','accessibility_reason') if k in pf.schema_arrow.names]
    raw=list(islice((r for b in pf.iter_batches(batch_size=64,columns=columns) for r in b.to_pylist() if r['aso_id']!='3UP8'),args.limit))
    from conference.model_v2.j8_feature_join import fence
    for r in raw: fence(r)
    rows=[dict(candidate_id=r['site_id'],variant_id=variant,mismatch_count=r['mismatch_count'],
               wc=r['complementarity_wc_pairs'],wc_status=r['complementarity_status'],
               energy=r['duplex_energy_kcal_mol'],energy_status=r['duplex_status'],
               access=r['accessibility_p_unpaired'],access_status=r['accessibility_status']) for r in raw]
    settings=[Setting('complementarity','wc',variant),Setting('duplex','energy',variant,False),Setting('accessibility','access',variant)]
    settings += [Setting(f'max_mismatch_{k}','wc',variant,True,k) for k in (0,1,2,3)]
    # Explicit missing alternatives, never silently reuse the one computed J3 window.
    settings += [Setting(v,'access',v) for v in FOLD_VARIANTS if v!=variant]
    ids=[r['candidate_id'] for r in rows]
    alternatives=alternative_scores(rows,settings,ids)
    nulls=[site_null(r,args.seed+i) for i,r in enumerate(raw[:3])]
    result=dict(contract=contract,j9=j9_report(j9rows,source=str(paths[0]),seed=args.seed,fold_variant=variant),j3_sample=dict(
        rows=len(raw),sampling_universe='first N non-3UP8 rows in explicitly named J3 shard; convenience sample, not genome representative',
        j3_fold_variant=variant,settings=[asdict(s) for s in settings],stability=rank_stability(alternatives),nulls=nulls),
        blocked=['Duplex-energy and accessibility decoy nulls blocked: both need a folding backend, '
                 'and J10 runs no scoring backend. Only the pure-sequence complementarity channel is rescored.',
                 'Shift distributions are descriptive: no selection-aware null for J6-discovered near-matches; p_value=None.',
                 'J3 alternative fold windows are absent; coverage explicitly incomplete.',
                 'No verified negatives: supervised metrics blocked by evaluation_core.'])
    result['payload_sha256']=payload_sha256(result)
    if not valid_checkpoint(result,contract): raise ValueError('invalid robustness result schema')
    atomic_json(output,result)


if __name__=='__main__': main()
