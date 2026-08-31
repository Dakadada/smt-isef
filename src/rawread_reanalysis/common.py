import csv, hashlib, json, os, re, shutil, tempfile, time
from contextlib import contextmanager
from pathlib import Path

PLACEHOLDER = re.compile(r"^(?:REQUIRED.*|TO_RESOLVE.*|PLACEHOLDER.*|/path/to.*|)$", re.I)
@contextmanager
def file_lock(path):
    """Hold an exclusive inter-process advisory lock for a filesystem transaction."""
    import fcntl
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w") as handle:
        fcntl.flock(handle,fcntl.LOCK_EX)
        yield
def load_config(path):
    text=Path(path).read_text()
    try: return json.loads(text)
    except json.JSONDecodeError:
        try: import yaml
        except ImportError as exc: raise ValueError("config must be JSON-compatible YAML, or install PyYAML") from exc
        return yaml.safe_load(text)
def atomic_json(path,value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True); fd,tmp=tempfile.mkstemp(prefix=f".{path.name}.",dir=path.parent)
    try:
        with os.fdopen(fd,"w") as h: json.dump(value,h,indent=2,sort_keys=True); h.write("\n"); h.flush(); os.fsync(h.fileno())
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)
def sha256(path):
    d=hashlib.sha256()
    with Path(path).open("rb") as h:
        for block in iter(lambda:h.read(1048576),b""): d.update(block)
    return d.hexdigest()
def md5(path):
    d=hashlib.md5() # nosec: archive integrity
    with Path(path).open("rb") as h:
        for block in iter(lambda:h.read(1048576),b""): d.update(block)
    return d.hexdigest()
def read_tsv(path):
    with Path(path).open(newline="") as h: return list(csv.DictReader(h,delimiter="\t"))
def write_tsv(path,rows,fields):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True); fd,tmp=tempfile.mkstemp(prefix=f".{path.name}.",dir=path.parent)
    try:
        with os.fdopen(fd,"w",newline="") as h:
            w=csv.DictWriter(h,fieldnames=fields,delimiter="\t",extrasaction="ignore"); w.writeheader(); w.writerows(rows); h.flush(); os.fsync(h.fileno())
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp): os.unlink(tmp)
def _bundle_hash(path):
    p=Path(path)
    if not p.exists(): raise ValueError(f"REQUIRED-INPUT asset does not exist: {p}")
    if p.is_file(): return sha256(p)
    files=sorted(x for x in p.rglob("*") if x.is_file())
    if not files: raise ValueError(f"reference index bundle is empty: {p}")
    d=hashlib.sha256()
    for f in files: d.update(str(f.relative_to(p)).encode()+b"\0"+bytes.fromhex(sha256(f)))
    return d.hexdigest()
def validate_config(config_path):
    cfg=load_config(config_path); refs=cfg.get("references",{}); hashes=cfg.get("asset_sha256",{})
    tier=cfg.get("alignment",{}).get("tier")
    if tier not in ("targeted","full"): raise ValueError("alignment.tier must be targeted or full")
    anomaly_threshold=cfg.get("psi",{}).get("marker_specificity_anomaly_rate_threshold",0.01)
    if not isinstance(anomaly_threshold,(int,float)) or isinstance(anomaly_threshold,bool) or not 0<=anomaly_threshold<=1: raise ValueError("psi.marker_specificity_anomaly_rate_threshold must be between 0 and 1")
    # Targeted mode scans frozen junction-marker strings in FASTQ and has no
    # genome/index/container dependency. Full mode keeps the all-assets gate.
    assets={"event_junctions":refs.get("event_junctions_tsv")}
    if tier=="full":
        if "GRCh38" not in refs.get("build","") or PLACEHOLDER.match(refs.get("annotation_release","")): raise ValueError("REQUIRED-INPUT: name GRCh38 and a verified Ensembl/GENCODE annotation_release")
        assets.update({"reference_fasta":refs.get("fasta"),"annotation_gtf":refs.get("gtf"),"star_index":refs.get("star_index"),"hisat2_index":refs.get("hisat2_index"),"container":refs.get("container")})
    verified={}
    for name,path in assets.items():
        expected=str(hashes.get(name,""))
        if not re.fullmatch(r"[0-9a-f]{64}",expected) or PLACEHOLDER.match(str(path or "")): raise ValueError(f"REQUIRED-INPUT: non-placeholder SHA-256 and path required for {name}")
        actual=_bundle_hash(path)
        if actual!=expected: raise ValueError(f"SHA-256 mismatch for {name}: expected {expected}, got {actual}")
        verified[name]=actual
    events=read_tsv(assets["event_junctions"]); required={e["id"] for e in cfg.get("events",[])}
    if not events: raise ValueError("REQUIRED-INPUT: event junction table is empty")
    seen=set(); columns=("event","chrom","strand","inclusion_up_start","inclusion_up_end","inclusion_down_start","inclusion_down_end","skipping_start","skipping_end","inclusion_up_marker","inclusion_down_marker","skipping_marker")
    for row in events:
        if any(PLACEHOLDER.match(str(row.get(k,""))) for k in columns): raise ValueError(f"malformed/placeholder event junction row: {row.get('event','<missing>')}")
        if row["event"] in seen: raise ValueError(f"duplicate event definition: {row['event']}")
        seen.add(row["event"])
        if row["strand"] not in ("+","-") or not row["chrom"].startswith("chr"): raise ValueError(f"malformed GRCh38 junction: {row['event']}")
        try: coords=[int(row[k]) for k in columns[3:9]]
        except ValueError as exc: raise ValueError(f"non-integer junction coordinate: {row['event']}") from exc
        if min(coords)<1 or any(coords[i]>=coords[i+1] for i in range(0,6,2)): raise ValueError(f"invalid junction coordinates: {row['event']}")
    if seen!=required: raise ValueError(f"event definitions must exactly cover configured events; missing={sorted(required-seen)}, extra={sorted(seen-required)}")
    scratch=Path(cfg["paths"]["scratch_dir"]).resolve(); namespace=Path(cfg["paths"].get("scratch_namespace",scratch)).resolve()
    if scratch==Path("/") or (scratch!=namespace and namespace not in scratch.parents): raise ValueError("scratch_dir fails namespace containment check")
    cfg["_verified_asset_hashes"]=verified; return cfg
def file_identities(paths): return {str(Path(p)):sha256(p) for p in paths}
class Checkpoints:
    def __init__(self,root,stage,config_path,cfg=None): self.root=Path(root)/stage; self.root.mkdir(parents=True,exist_ok=True); self.config_hash=sha256(config_path); self.assets=(cfg or load_config(config_path)).get("_verified_asset_hashes",{})
    def path(self,item): return self.root/f"{item}.json"
    def valid(self,item,outputs=(),inputs=(),metadata=None):
        try: state=json.loads(self.path(item).read_text())
        except (OSError,json.JSONDecodeError): return False
        if state.get("status")!="complete" or state.get("config_sha256")!=self.config_hash or state.get("verified_asset_hashes")!=self.assets or state.get("inputs")!=file_identities(inputs): return False
        if metadata is not None and state.get("resolved_metadata")!=metadata: return False
        return all(Path(p).is_file() and state.get("outputs",{}).get(str(Path(p)))==sha256(p) for p in outputs)
    def complete(self,item,outputs,inputs=(),resolved_metadata=None,**metadata): atomic_json(self.path(item),{"status":"complete","completed_at":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),"config_sha256":self.config_hash,"verified_asset_hashes":self.assets,"inputs":file_identities(inputs),"outputs":file_identities(outputs),"resolved_metadata":resolved_metadata,**metadata})
def require_space(path,needed_bytes,aggregate_cap=None):
    p=Path(path); p.mkdir(parents=True,exist_ok=True); used=sum(x.stat().st_size for x in p.rglob("*") if x.is_file())
    if aggregate_cap is not None and used+needed_bytes>aggregate_cap: raise RuntimeError(f"scratch aggregate cap exceeded: {used}+{needed_bytes}>{aggregate_cap}")
    free=shutil.disk_usage(p).free
    if free<needed_bytes: raise RuntimeError(f"scratch quota guard: need {needed_bytes} bytes, only {free} free at {p}")
