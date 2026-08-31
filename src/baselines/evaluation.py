"""Identity-component-disjoint, provenance-locked baseline evaluation."""
import csv
import json
import random
from pathlib import Path
from .features import FEATURE_SETS
from .models import CapacityMatchedLinear

PROTOCOL_MIN_CI_FAMILIES = 8
IDENTITY_FIELDS = ("study", "locus", "intervention", "chemistry", "sample")


class EvaluationConfig:
    def __init__(self, neutral_threshold=5.0, sensitivity_thresholds=(2.0,10.0), min_ci_families=8,
                 max_effective_df=8.0, ridge_grid=(0.1,1.0,10.0,100.0,1000.0), kmer_k_grid=(3,4,5), seed=20260828):
        if int(min_ci_families) < PROTOCOL_MIN_CI_FAMILIES: raise ValueError("min_ci_families cannot be lower than protocol minimum 8")
        self.neutral_threshold=float(neutral_threshold); self.sensitivity_thresholds=tuple(map(float,sensitivity_thresholds))
        self.min_ci_families=int(min_ci_families); self.max_effective_df=float(max_effective_df)
        self.ridge_grid=tuple(map(float,ridge_grid)); self.kmer_k_grid=tuple(map(int,kmer_k_grid)); self.seed=int(seed)
        if not self.ridge_grid or any(x <= 0 for x in self.ridge_grid): raise ValueError("positive pre-registered ridge grid required")

    @classmethod
    def from_path(cls,path):
        with Path(path).open() as h: cfg=json.load(h)
        return cls(cfg["neutral_threshold"],cfg["sensitivity_thresholds"],cfg["min_ci_families"],cfg["max_effective_df"],cfg["ridge_grid"],cfg["kmer_k_grid"],cfg["seed"])


def _truth(value,threshold): return 1 if value>=threshold else (-1 if value<=-threshold else 0)
def _balanced_accuracy(y,p):
    classes=sorted(set(y)); return sum(sum(a==c and b==c for a,b in zip(y,p))/sum(a==c for a in y) for c in classes)/len(classes)
def _metrics(rows,threshold):
    y=[_truth(r["observed"],threshold) for r in rows]; p=[_truth(r["predicted"],threshold) for r in rows]
    return {"directional_balanced_accuracy":_balanced_accuracy(y,p),"directional_concordance":sum(a==b for a,b in zip(y,p))/len(y)}


def _theta_valid(row):
    kind=row["theta_evidence_type"]
    if kind=="none": return not row["theta_input_hash"] and not row["theta_protocol_hash"]
    if kind!="matched-occupancy-evidence/v1": return False
    return all(len(row[x])==64 and all(c in "0123456789abcdef" for c in row[x].lower()) for x in ("theta_input_hash","theta_protocol_hash"))


def load_harmonized(path):
    """Load only schema-proven treated-minus-matched-control contrasts."""
    with Path(path).open(newline="") as h: rows=list(csv.DictReader(h,delimiter="\t"))
    required={"cluster_family","cluster_study","cluster_locus","cluster_biological_sample","cluster_intervention","cluster_chemistry","cluster_replicate","event","primary_event","condition","magnitude","magnitude_eligible","coverage_pass","contrast_type","control_cluster","matched_control_count","sign_convention","units","contrast_schema","theta_evidence_type","theta_input_hash","theta_protocol_hash"}
    if rows and not required<=set(rows[0]): raise ValueError(f"harmonized labels missing columns: {sorted(required-set(rows[0]))}")
    grouped={}
    for row in rows:
        for key in ("cluster_family","cluster_study","cluster_locus","cluster_biological_sample","cluster_intervention","cluster_chemistry","cluster_replicate"):
            if not row[key].strip(): raise ValueError(f"missing identity field: {key}")
        if row["condition"]!="treated" or row["primary_event"].lower()!="true" or row["coverage_pass"].lower()!="true" or not row["magnitude"]: continue
        if row["contrast_schema"]!="harmonized-treated-minus-matched-control/v1" or row["contrast_type"]!="treated-minus-control" or row["sign_convention"]!="treated-minus-control" or not row["units"].strip() or not row["control_cluster"].strip() or row["matched_control_count"]!="1":
            raise ValueError("contrast provenance is not a verified one-control treated-minus-control contrast")
        eligible=row["magnitude_eligible"].lower()=="true"
        if not _theta_valid(row) or (eligible and row["theta_evidence_type"]=="none"): raise ValueError("invalid typed theta-evidence provenance")
        grouped.setdefault(row["cluster_intervention"],[]).append(row)
    out={}
    for key,reps in grouped.items():
        invariant=("cluster_family","cluster_study","cluster_locus","cluster_intervention","cluster_chemistry","event","control_cluster","units","sign_convention")
        for field in invariant:
            if len({r[field] for r in reps})!=1: raise ValueError(f"intervention has inconsistent {field}")
        if len({r["cluster_replicate"] for r in reps})!=len(reps): raise ValueError("duplicate replicate cluster identity")
        first=reps[0]
        out[key]={"family":first["cluster_family"],"study":first["cluster_study"],"locus":first["cluster_locus"],"sample":tuple(sorted({r["cluster_biological_sample"] for r in reps})),"intervention":key,"chemistry":first["cluster_chemistry"],"replicate_clusters":tuple(sorted(r["cluster_replicate"] for r in reps)),"event":first["event"],"observed":sum(float(r["magnitude"]) for r in reps)/len(reps),"magnitude_eligible":all(r["magnitude_eligible"].lower()=="true" for r in reps)}
    return out


def _validate_labels(labels):
    for key,row in labels.items():
        missing=[f for f in ("family","study","locus","sample","intervention","chemistry","replicate_clusters","event","observed","magnitude_eligible") if f not in row]
        if missing: raise ValueError(f"label {key} missing fields: {missing}")
        for field in IDENTITY_FIELDS:
            value=row[field]
            values=value if isinstance(value,(tuple,list,set)) else (value,)
            if not values or any(v is None or not str(v).strip() for v in values): raise ValueError(f"label {key} missing identity {field}")


def _components(records):
    parent=list(range(len(records)))
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    def union(a,b):
        a,b=find(a),find(b)
        if a!=b: parent[b]=a
    seen={}
    for i,row in enumerate(records):
        for field in IDENTITY_FIELDS:
            vals=row[field] if isinstance(row[field],(tuple,list,set)) else (row[field],)
            for value in vals:
                token=(field,str(value)); union(i,seen[token]) if token in seen else seen.setdefault(token,i)
    groups={}
    for i in range(len(records)): groups.setdefault(find(i),[]).append(i)
    return [groups[k] for k in sorted(groups,key=lambda x:min(groups[x]))]


def _assert_identity_disjoint(train,test):
    overlap={}
    for field in IDENTITY_FIELDS:
        def vals(rows):
            answer=set()
            for row in rows:
                value=row[field]; answer.update(value if isinstance(value,(tuple,list,set)) else (value,))
            return answer
        common=vals(train)&vals(test)
        if common: overlap[field]=sorted(map(str,common))
    if overlap: raise AssertionError(f"identity leakage across fold: {overlap}")


def identity_disjoint_component_splits(records):
    """Return leave-component-out folds and authenticate every identity boundary.

    This is the only outer-split primitive used by the decisive benchmark.  A
    component is the transitive closure of shared study, locus, intervention,
    chemistry, or biological-sample identity.
    """
    records = list(records)
    folds = []
    for indices in _components(records):
        held = frozenset(indices)
        train = tuple(row for i, row in enumerate(records) if i not in held)
        test = tuple(row for i, row in enumerate(records) if i in held)
        _assert_identity_disjoint(train, test)
        folds.append((train, test, held))
    return tuple(folds)


def _select_features(row,selectors): return {k:v for k,v in row.items() if any(k==s or (s.startswith("_") and k.endswith(s)) or (s.endswith("_") and k.startswith(s)) for s in selectors)}
def _select_k(row,k):
    if k is None: return row
    return {name:value for name,value in row.items() if not name.startswith("kmer_k") or name.startswith(f"kmer_k{k}_")}
def _weights(records):
    family_interventions={}; int_n={}
    for row in records: family_interventions.setdefault(row["family"],set()).add(row["intervention"]); int_n[row["intervention"]]=int_n.get(row["intervention"],0)+1
    return [1/len(family_interventions[r["family"]])/int_n[r["intervention"]] for r in records]


def _cluster_uncertainty(rows,threshold,seed,draws=10000):
    grouped={}
    for row in rows: grouped.setdefault(row["family"],[]).append(row)
    families=sorted(grouped); rng=random.Random(seed); values={k:[] for k in ("directional_balanced_accuracy","directional_concordance")}; supports=[]; effective=[]
    for _ in range(draws):
        picked=[rng.choice(families) for _ in families]; sampled=[]
        for draw_id,family in enumerate(picked):
            sampled.extend(dict(row,bootstrap_family_instance=f"{draw_id}:{family}") for row in grouped[family])
        metric=_metrics(sampled,threshold)
        for key in values: values[key].append(metric[key])
        counts=[picked.count(f) for f in set(picked)]; supports.append(len(counts)); effective.append(sum(counts)**2/sum(x*x for x in counts))
    intervals={}
    for key,samples in values.items(): samples.sort(); intervals[key]=[samples[int(.025*draws)],samples[min(draws-1,int(.975*draws))]]
    loo={key:[] for key in values}
    for held in families:
        metric=_metrics([r for r in rows if r["family"]!=held],threshold)
        for key in loo: loo[key].append(metric[key])
    return {"intervals":intervals,"bootstrap_unique_family_support":{"min":min(supports),"mean":sum(supports)/draws,"max":max(supports)},"bootstrap_effective_clusters":{"min":min(effective),"mean":sum(effective)/draws,"max":max(effective)},"leave_one_family_out_sensitivity":{k:[min(v),max(v)] for k,v in loo.items()}}


class BenchmarkHarness:
    def __init__(self,config=None): self.config=config or EvaluationConfig()
    def _choose_ridge(self,training,feature_key):
        comps=_components(training); scores=[]
        has_tuned_k=any(any(name.startswith("kmer_k") for name in row[feature_key]) for row in training)
        k_candidates=self.config.kmer_k_grid if has_tuned_k else (None,)
        for k in k_candidates:
          for ridge in self.config.ridge_grid:
            errors=[]; valid=True
            for indices in comps:
                held=set(indices); tr=[r for i,r in enumerate(training) if i not in held]; va=[r for i,r in enumerate(training) if i in held]
                if not tr: continue
                try: model=CapacityMatchedLinear(self.config.max_effective_df,ridge).fit([_select_k(r[feature_key],k) for r in tr],[r["observed"] for r in tr],_weights(tr))
                except ValueError: valid=False; break
                errors.extend(abs(a-b) for a,b in zip(model.predict([_select_k(r[feature_key],k) for r in va]),[r["observed"] for r in va]))
            if valid and errors: scores.append((sum(errors)/len(errors),ridge,k))
        if not scores: raise ValueError("no pre-registered ridge satisfies EDF budget")
        _,ridge,k=min(scores,key=lambda x:(x[0],x[1],-1 if x[2] is None else x[2])); return ridge,k

    def run(self,labels,features,*,physics_effective_degrees_of_freedom):
        _validate_labels(labels)
        physics_edf=float(physics_effective_degrees_of_freedom)
        if not 0 < physics_edf <= self.config.max_effective_df: raise ValueError("physics adapter exceeds matched EDF budget")
        records=[dict(label,all_features=features[key]) for key,label in labels.items() if key in features]
        splits=identity_disjoint_component_splits(records); comps=[sorted(held) for _,_,held in splits]
        if len(comps)<2: raise ValueError("at least two identity-disjoint connected components required")
        output={}
        for name,selectors in FEATURE_SETS.items():
            prepared=[dict(r,selected=_select_features(r["all_features"],selectors)) for r in records]; predictions=[]; fold_meta=[]
            for fold_id,indices in enumerate(comps):
                held=set(indices); tr=[r for i,r in enumerate(prepared) if i not in held]; te=[r for i,r in enumerate(prepared) if i in held]
                _assert_identity_disjoint(tr,te); ridge,selected_k=self._choose_ridge(tr,"selected")
                model=CapacityMatchedLinear(self.config.max_effective_df,ridge).fit([_select_k(r["selected"],selected_k) for r in tr],[r["observed"] for r in tr],_weights(tr))
                predictions.extend(dict(r,predicted=p) for r,p in zip(te,model.predict([_select_k(r["selected"],selected_k) for r in te])))
                fold_meta.append({"held_out_component":fold_id,"ridge":ridge,"selected_k":selected_k,"baseline_effective_degrees_of_freedom":model.effective_degrees_of_freedom,"physics_effective_degrees_of_freedom":physics_edf,"matched_edf_budget":self.config.max_effective_df,"identity_disjoint":True})
            metrics=_metrics(predictions,self.config.neutral_threshold); metrics["sensitivities"]={str(t):_metrics(predictions,t) for t in self.config.sensitivity_thresholds}
            eligible=[r for r in predictions if r["magnitude_eligible"]]; eligible_families={r["family"] for r in eligible}
            if len(eligible_families)>=PROTOCOL_MIN_CI_FAMILIES:
                mae=sum(abs(r["predicted"]-r["observed"]) for r in eligible)/len(eligible)
                xm=sum(r["predicted"] for r in eligible)/len(eligible); ym=sum(r["observed"] for r in eligible)/len(eligible)
                denom=sum((r["predicted"]-xm)**2 for r in eligible)
                calibration=None if denom==0 else sum((r["predicted"]-xm)*(r["observed"]-ym) for r in eligible)/denom
                status="decisive-matched-theta-evidence"
            else: mae=None; calibration=None; status="withheld-fewer-than-8-theta-evidence-families"
            families={r["family"] for r in predictions}; ci_ok=len(families)>=self.config.min_ci_families
            metrics.update({"mae":mae,"calibration_slope":calibration,"magnitude_n":len(eligible),"magnitude_status":status,"n_families":len(families),"n_identity_components":len(comps),"n_interventions":len(predictions),"cluster_uncertainty":_cluster_uncertainty(predictions,self.config.neutral_threshold,self.config.seed+sum(map(ord,name))) if ci_ok else None,"ci_status":"family-bootstrap-plus-leave-one-family-out" if ci_ok else f"withheld-fewer-than-{self.config.min_ci_families}-families","folds":fold_meta,"predictions":[{k:r[k] for k in ("family","study","locus","sample","intervention","chemistry","event","observed","predicted","magnitude_eligible")} for r in predictions]})
            output[name]=metrics
        return output
