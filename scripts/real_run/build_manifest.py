#!/usr/bin/env python3
"""Freeze ENA run metadata and accession-level Stage-A arm assignments."""
import argparse
import csv
import hashlib
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

PROJECTS = {
    "GSE307431_PRJNA1321660": ("PRJNA1321660", 20),
    "GSE167762_PRJNA705187": ("PRJNA705187", 12),
    "PRJNA1052769_dCas13d": ("PRJNA1052769", 4),
    "PRJEB13980_U7": ("PRJEB13980", 12),
}
FAMILY_EVENTS = {
    "GSE307431_PRJNA1321660": "SMN2_exon7_inclusion",
    "GSE167762_PRJNA705187": "SMN2_exon7_inclusion",
    "PRJNA1052769_dCas13d": "SMN2_exon7_inclusion",
    "PRJEB13980_U7": "DMD_exon53_inclusion",
}
ENA_FIELDS = ["study_accession", "run_accession", "sample_accession", "secondary_sample_accession",
              "sample_alias", "sample_title", "experiment_alias", "experiment_title", "fastq_ftp",
              "fastq_md5", "fastq_bytes", "library_layout", "library_strategy"]
MAP_FIELDS = ["run_accession", "intervention", "condition", "contrast", "pair_set_id", "covariates"]
MANIFEST_FIELDS = ["family", "project", "run_accession", "sample_accession", "sample_alias",
                   "experiment_title", "fastq_ftp", "fastq_md5", "fastq_bytes", "layout",
                   "library_strategy", "endpoint_eligible", "exclusion_reason", "retrieved_at",
                   "metadata_source", *MAP_FIELDS[1:]]

# Exact GEO sample-to-arm mapping transcribed from docs/dataset-audit.md.  These
# GSMs, not fuzzy treatment words, determine the two GEO study assignments.
GEO_ARMS = {
    **{f"GSM{x}": ("F18OMe", "treated", "OMe") for x in (9224136, 9224142, 9224149)},
    **{f"GSM{x}": ("scramble_OMe", "control", "OMe") for x in (9224135, 9224141, 9224148)},
    **{f"GSM{x}": ("F18MOE", "treated", "MOE") for x in (9224138, 9224144, 9224151)},
    **{f"GSM{x}": ("scramble_MOE", "control", "MOE") for x in (9224137, 9224143, 9224150)},
    **{f"GSM{x}": ("F20PMO", "treated", "PMO") for x in (9224140, 9224146, 9224153)},
    **{f"GSM{x}": ("scramble_PMO", "control", "PMO") for x in (9224139, 9224145, 9224152)},
    **{f"GSM{x}": ("untreated", "control", "untreated_only") for x in (9224134, 9224147)},
    **{f"GSM{x}": ("untreated", "control", "pure_ASO1") for x in range(6063682, 6063685)},
    **{f"GSM{x}": ("VPA_only", "control", "factorial_combo") for x in range(6063685, 6063688)},
    **{f"GSM{x}": ("ASO1", "treated", "pure_ASO1") for x in range(6063688, 6063691)},
    **{f"GSM{x}": ("ASO1_plus_VPA", "treated", "factorial_combo") for x in range(6063691, 6063694)},
}

# Verified from the ENA PRJNA1052769 filereport. Only these four endpoint
# RNA-seq runs belong to the D[333,357] versus non-targeting-guide contrast.
DCAS13D_ARMS = {
    "SRR27213209": ("dCas13d_D[333,357]", "treated"),
    "SRR27213208": ("dCas13d_D[333,357]", "treated"),
    "SRR27213210": ("non_targeting_guide", "control"),
    "SRR27213211": ("non_targeting_guide", "control"),
}

def ena(project, metadata_dir=None):
    if metadata_dir:
        source = Path(metadata_dir) / f"{project}.tsv"
        if not source.is_file(): raise ValueError(f"REQUIRED-INPUT: missing mock/export {source}")
        text = source.read_text(); label = str(source)
    else:
        query = urllib.parse.urlencode({"accession": project, "result": "read_run",
                                        "fields": ",".join(ENA_FIELDS), "format": "tsv", "download": "true"})
        url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + query
        with urllib.request.urlopen(url, timeout=180) as response: text = response.read().decode()
        label = url
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows: raise ValueError(f"REQUIRED-INPUT: ENA returned no runs for {project}")
    return rows, label

def load_overrides(path):
    if not path: return {}
    with Path(path).open(newline="") as handle: rows = list(csv.DictReader((line for line in handle if not line.startswith("#")), delimiter="\t"))
    if not rows or not set(MAP_FIELDS).issubset(rows[0]): raise ValueError("override TSV must contain " + ",".join(MAP_FIELDS))
    if len({r["run_accession"] for r in rows}) != len(rows): raise ValueError("duplicate override run_accession")
    return {r["run_accession"]: r for r in rows}

def text(row): return " ".join(str(row.get(k, "")) for k in ENA_FIELDS)
def gsm(row):
    hits = set(re.findall(r"GSM\d+", text(row), re.I))
    return next(iter(hits)).upper() if len(hits) == 1 else None

def covariates(family, intervention):
    if family == "GSE307431_PRJNA1321660":
        chemistry = {"F18OMe":"2OMe-PS", "scramble_OMe":"2OMe-PS", "F18MOE":"2MOE-PS",
                     "scramble_MOE":"2MOE-PS", "F20PMO":"PMO", "scramble_PMO":"PMO"}.get(intervention,"none")
        return {"chemistry":chemistry,"system":"GM03813","dose":"0" if intervention=="untreated" else "6uM","time":"24h"}
    if family == "GSE167762_PRJNA705187": return {"VPA":"yes" if "VPA" in intervention else "no","system":"HEK293T"}
    if family == "PRJNA1052769_dCas13d": return {"effector":"dCas13d","system":"HEK293T"}
    return {}

def assignment(family, row, override):
    if override:
        try: cov = json.loads(override["covariates"])
        except json.JSONDecodeError as exc: raise ValueError(f"{row['run_accession']}: covariates must be JSON") from exc
        mapped={k:override[k] for k in MAP_FIELDS[1:-1]}; mapped["batch_covariates"]=cov
        return mapped
    if family.startswith("GSE"):
        key = gsm(row)
        if key in GEO_ARMS:
            intervention, condition, contrast = GEO_ARMS[key]
            return {"intervention":intervention,"condition":condition,"contrast":contrast,
                    "pair_set_id":contrast,"batch_covariates":covariates(family,intervention)}
    if family == "PRJNA1052769_dCas13d" and row["run_accession"] in DCAS13D_ARMS:
        intervention, condition = DCAS13D_ARMS[row["run_accession"]]
        contrast = "targeting_vs_nontargeting"
        return {"intervention":intervention,"condition":condition,"contrast":contrast,
                "pair_set_id":contrast,"batch_covariates":covariates(family,intervention)}
    # The audit does not document run-level mappings for these archives. Refuse
    # to infer them from prose-like ENA titles; the operator supplies verified rows.
    return None

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--out",required=True); ap.add_argument("--base-config",default="config/pipeline.yaml")
    ap.add_argument("--metadata-dir",help="offline ENA TSV exports, one PROJECT.tsv each")
    ap.add_argument("--overrides",help="verified accession mapping TSV for non-GEO/ambiguous runs")
    ap.add_argument("--families",help="comma-separated family names to build (default: all families)")
    args=ap.parse_args(); root=Path(__file__).resolve().parents[2]; out=Path(args.out).resolve(); out.mkdir(parents=True,exist_ok=True)
    if args.families:
        selected=list(dict.fromkeys(part.strip() for part in args.families.split(",") if part.strip()))
        unknown=sorted(set(selected)-set(PROJECTS))
        if not selected: raise ValueError("--families must name at least one family")
        if unknown: raise ValueError(f"unknown families: {unknown}")
    else: selected=list(PROJECTS)
    overrides=load_overrides(args.overrides); cfg=json.loads((root/args.base_config).read_text()); cfg["alignment"]["tier"]="targeted"
    cfg["paths"].update({"work_dir":str(out/"work"),"scratch_dir":str(out/"scratch"),"scratch_namespace":str(out/"scratch"),"checkpoints":str(out/"work/checkpoints")})
    cfg["studies"]={family:cfg["studies"][family] for family in selected}
    cfg["events"]=[event for event in cfg["events"] if event["family"] in selected]
    targeted_events={FAMILY_EVENTS[family] for family in selected}
    configured_events={event["id"] for event in cfg["events"]}
    if configured_events != targeted_events:
        raise ValueError(f"family/event mapping disagrees with base config: mapped={sorted(targeted_events)}, configured={sorted(configured_events)}")
    master_events=root/"config/event_junctions.GRCh38.tsv"
    run_events=out/"event_junctions.GRCh38.run.tsv"
    write_event_subset(master_events,run_events,targeted_events)
    cfg["references"]["event_junctions_tsv"]=str(run_events)
    cfg["asset_sha256"]["event_junctions"]=hashlib.sha256(run_events.read_bytes()).hexdigest()
    all_rows=[]; used=set(); now=time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime())
    for family in selected:
        project,expected=PROJECTS[family]
        raw,source=ena(project,args.metadata_dir); cooked=[]; sample_map={}
        for row in raw:
            eligible=row.get("library_strategy","").casefold()=="rna-seq"
            if family == "GSE167762_PRJNA705187": eligible = eligible and gsm(row) in GEO_ARMS
            if family=="PRJNA1052769_dCas13d": eligible = eligible and row["run_accession"] in DCAS13D_ARMS
            mapping=assignment(family,row,overrides.get(row["run_accession"])) if eligible else None
            if eligible and mapping is None: raise ValueError(f"REQUIRED-INPUT: {family}/{row['run_accession']} has no verified arm mapping; add it to --overrides")
            if eligible:
                for key in ("intervention","condition","contrast","pair_set_id"):
                    if not mapping.get(key) or str(mapping[key]).startswith(("REQUIRED","TO_RESOLVE")): raise ValueError(f"{row['run_accession']}: unresolved {key}")
                if not mapping["batch_covariates"] or any(str(v).startswith(("REQUIRED","TO_RESOLVE")) for v in mapping["batch_covariates"].values()): raise ValueError(f"{row['run_accession']}: unresolved covariates")
                sample_map[row["run_accession"]]=mapping; used.add(row["run_accession"])
            fq=row.get("fastq_ftp",""); sums=row.get("fastq_md5","")
            if eligible and (not fq or not sums or len(fq.split(";"))!=len(sums.split(";"))): raise ValueError(f"{row['run_accession']}: missing/inconsistent FASTQ URL/MD5")
            record={"family":family,"project":project,"run_accession":row.get("run_accession",""),"sample_accession":row.get("sample_accession",""),
                    "sample_alias":row.get("sample_alias","") or row.get("secondary_sample_accession",""),"experiment_title":row.get("experiment_title","") or row.get("sample_title","") or row.get("experiment_alias",""),
                    "fastq_ftp":fq,"fastq_md5":sums,"fastq_bytes":row.get("fastq_bytes","") or "","layout":row.get("library_layout","") or "",
                    "library_strategy":row.get("library_strategy","") or "","endpoint_eligible":str(eligible).lower(),"exclusion_reason":"" if eligible else "not verified endpoint RNA-seq",
                    "retrieved_at":now,"metadata_source":source}
            if mapping: record.update({k:(json.dumps(mapping["batch_covariates"],sort_keys=True) if k=="covariates" else mapping[k]) for k in MAP_FIELDS[1:]})
            else: record.update({k:"" for k in MAP_FIELDS[1:]})
            cooked.append(record)
        if len(sample_map)!=expected: raise ValueError(f"REQUIRED-INPUT: {family} expected {expected} endpoint runs, resolved {len(sample_map)}")
        cfg["studies"][family]["sample_map"]=sample_map
        family_path=out/f"{family}.frozen.tsv"; write_tsv(family_path,cooked); all_rows.extend(cooked)
    extra=sorted(set(overrides)-used)
    if extra and not args.families: raise ValueError(f"override accessions were not eligible ENA runs: {extra}")
    write_tsv(out/"all_runs.frozen.tsv",all_rows)
    config_path=out/"stageA.targeted.json"; config_path.write_text(json.dumps(cfg,indent=2,sort_keys=True)+"\n")
    frozen=[out/"all_runs.frozen.tsv",config_path,run_events]+[out/f"{family}.frozen.tsv" for family in selected]
    for path in frozen: (path.with_suffix(path.suffix+".sha256")).write_text(hashlib.sha256(path.read_bytes()).hexdigest()+"  "+path.name+"\n")

def write_event_subset(source,path,required):
    # Keep the verified header and selected definition rows byte-for-byte.
    lines=Path(source).read_bytes().splitlines(keepends=True)
    if not lines: raise ValueError("REQUIRED-INPUT: master event junction table is empty")
    selected=[]; seen=set()
    for line in lines[1:]:
        event=line.split(b"\t",1)[0].decode("utf-8")
        if event in required:
            selected.append(line); seen.add(event)
    if seen != required:
        raise ValueError(f"REQUIRED-INPUT: master event definitions missing {sorted(required-seen)}")
    Path(path).write_bytes(lines[0]+b"".join(selected))

def write_tsv(path,rows):
    with Path(path).open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=MANIFEST_FIELDS,delimiter="\t",extrasaction="ignore"); writer.writeheader(); writer.writerows(rows)
if __name__=="__main__":
    try: main()
    except Exception as exc: print(f"FATAL: {exc}",file=sys.stderr); raise SystemExit(2)
