import math,re
from collections import defaultdict
from pathlib import Path
from .common import Checkpoints, read_tsv, validate_config, write_tsv

FIELDS=["family","sample","run_accession","intervention","condition","event","psi","inclusion_count","skipping_count","informative_reads","coverage_pass","contrast","pair_set_id","control_mean_psi","dpsi","direction","direction_band_2","direction_band_10","therapeutic_polarity_multiplier","therapeutic_effect","magnitude_eligible","batch_covariates","quantification_tier","approximation_label","confound_flags"]
def psi(up,down,skip):
    denom=(up+down)/2+skip
    return None if denom==0 else 100*((up+down)/2)/denom
def direction(delta,band=5):
    if delta is None or math.isnan(delta): return "MISSING"
    if abs(delta)<band:return "no_change"
    return "increase" if delta>0 else "decrease"
def _sample_map(cfg,family,count):
    study=cfg["studies"][family]; acc=count.get("run_accession","")
    if acc in study.get("sample_map",{}): return study["sample_map"][acc]
    matches=[r for r in study.get("sample_rules",[]) if re.search(r["regex"],count["sample"],re.I)]
    if len(matches)!=1: raise ValueError(f"{family}/{acc}: sample arm mapping matched {len(matches)} rules; require exactly one frozen accession mapping")
    return matches[0]
def run(config_path,counts_path=None):
    cfg=validate_config(config_path); work=Path(cfg["paths"]["work_dir"]); counts_path=Path(counts_path or work/"counts/all_counts.tsv"); output=work/"psi/sample_effects.tsv"; cp=Checkpoints(cfg["paths"]["checkpoints"],"compute",config_path,cfg)
    rows=[]; raw=read_tsv(counts_path)
    if not raw: raise ValueError("empty counts input")
    metadata={f:cfg["studies"][f].get("sample_map",{}) for f in cfg["studies"]}
    if cp.valid("all",[output],[counts_path],metadata): return output
    tiers=defaultdict(set)
    for count in raw: tiers[(count["family"],count["event"])].add(count["quantification_tier"])
    if any(len(v)>1 for v in tiers.values()) and not cfg.get("tier_concordance",{}).get("passed",False): raise ValueError("cross-tier primary-label pooling prohibited until frozen concordance gate passes")
    minimum=cfg["psi"]["minimum_informative_reads"]
    for count in raw:
        mapping=_sample_map(cfg,count["family"],count); informative=int(count["informative_reads"]); passed=informative>=minimum and count.get("count_status") in ("observed","verified_zero")
        inclusion=int(count.get("inclusion_fragments","") or (int(count["inclusion_up"])+int(count["inclusion_down"])))
        skipping=int(count.get("skipping_fragments","") or count["skipping"])
        p=100.0*inclusion/(inclusion+skipping) if passed and inclusion+skipping else None
        if not mapping.get("pair_set_id") or any(str(v).startswith("TO_RESOLVE") for v in mapping.get("batch_covariates",{}).values()): raise ValueError(f"{count['run_accession']}: unresolved pair/covariates")
        study=cfg["studies"][count["family"]]
        rows.append({**count,"intervention":mapping["intervention"],"condition":mapping["condition"],"contrast":mapping["contrast"],"pair_set_id":mapping["pair_set_id"],"psi":"" if p is None else f"{p:.6g}","inclusion_count":inclusion,"skipping_count":skipping,"coverage_pass":str(passed).lower(),"batch_covariates":";".join(f"{k}={v}" for k,v in sorted(mapping.get("batch_covariates",{}).items())),"confound_flags":";".join(study.get("confound_flags",[]))})
    controls=defaultdict(list)
    for row in rows:
        if row["condition"]=="control" and row["coverage_pass"]=="true": controls[(row["family"],row["contrast"],row["event"],row["pair_set_id"],row["quantification_tier"])].append(float(row["psi"]))
    polarity={e["id"]:int(e["therapeutic_polarity_multiplier"]) for e in cfg["events"]}
    for row in rows:
        vals=controls[(row["family"],row["contrast"],row["event"],row["pair_set_id"],row["quantification_tier"])]; control=sum(vals)/len(vals) if vals else None; delta=float(row["psi"])-control if row["psi"]!="" and control is not None else None; mult=polarity[row["event"]]
        row.update({"control_mean_psi":"" if control is None else f"{control:.6g}","dpsi":"" if delta is None else f"{delta:.6g}","direction":direction(delta,5),"direction_band_2":direction(delta,2),"direction_band_10":direction(delta,10),"therapeutic_polarity_multiplier":mult,"therapeutic_effect":"MISSING" if delta is None else f"{delta*mult:.6g}","magnitude_eligible":str(bool(cfg["studies"][row["family"]].get("matched_theta_proxy",False))).lower()})
    write_tsv(output,rows,FIELDS); cp.complete("all",[output],inputs=[counts_path],resolved_metadata=metadata,formula="raw inclusion dPSI=treated-matched-control"); return output
