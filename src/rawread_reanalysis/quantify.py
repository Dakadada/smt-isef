import gzip, json, shlex, shutil, subprocess, time
from pathlib import Path
from .common import Checkpoints, atomic_json, file_lock, read_tsv, validate_config, write_tsv
from .fetch import _accession_lock

COUNT_FIELDS=["family","sample","run_accession","event","inclusion_up","inclusion_down","skipping","inclusion_fragments","skipping_fragments","anomalous_fragments","anomaly_rate","marker_specificity_warning","informative_reads","count_status","quantification_tier","approximation_label"]
def _execute(command,log):
    log.parent.mkdir(parents=True,exist_ok=True)
    with log.open("a") as h: h.write("$ "+shlex.join(command)+"\n"); subprocess.run(command,check=True,stdout=h,stderr=subprocess.STDOUT)
def _definitions(cfg): return {r["event"]:r for r in read_tsv(cfg["references"]["event_junctions_tsv"])}
def _reads(path):
    opener=gzip.open if path.suffix==".gz" else open
    with opener(path,"rt") as h:
        while True:
            name=h.readline().strip()
            if not name: return
            seq=h.readline().strip().upper(); plus=h.readline(); quality=h.readline()
            if not plus or not quality: raise ValueError(f"truncated FASTQ: {path}")
            key=name.split()[0]
            if key.endswith(("/1","/2")): key=key[:-2]
            yield key,seq
def _targeted_fastq(files,event,progress_path=None,heartbeat_seconds=600):
    markers={k:event[f"{k}_marker"].upper() for k in ("inclusion_up","inclusion_down","skipping")}; counts={k:0 for k in (*markers,"inclusion_fragments","skipping_fragments","anomalous_fragments")}
    processed=0
    if progress_path and Path(progress_path).is_file():
        state=json.loads(Path(progress_path).read_text())
        if state.get("event")!=event["event"]: raise ValueError("targeted progress event mismatch")
        # Old checkpoints cannot reconstruct deduplicated fragment counts, so replay them safely.
        if state.get("counting_semantics")==2:
            processed=int(state["processed_fragments"]); counts={k:int(state["counts"].get(k,0)) for k in counts}
    last=time.monotonic(); seen=0
    def consume(sequences):
        nonlocal seen,last
        seen+=1
        if seen<=processed: return
        _count_fragments((sequences,),markers,counts,event)
        if progress_path and time.monotonic()-last>=heartbeat_seconds:
            atomic_json(progress_path,{"event":event["event"],"counting_semantics":2,"processed_fragments":seen,"counts":counts}); last=time.monotonic()
    if len(files)==2:
        import itertools
        for left,right in itertools.zip_longest(_reads(files[0]),_reads(files[1])):
            if left is None or right is None or left[0]!=right[0]: raise ValueError("paired FASTQs are unequal or not name-synchronized")
            consume((left[1],right[1]))
    elif len(files)==1:
        for _,seq in _reads(files[0]): consume((seq,))
    else: raise ValueError("targeted tier requires one single-end FASTQ or one synchronized paired-end pair per run")
    if progress_path:
        try: Path(progress_path).unlink()
        except FileNotFoundError: pass
    return counts
def _count_fragments(fragments,markers,counts,event):
    for sequences in fragments:
        hits={k for seq in sequences for k,m in markers.items() if m in seq}
        inclusion_hits=hits & {"inclusion_up","inclusion_down"}
        if inclusion_hits and "skipping" in hits:
            counts["anomalous_fragments"]+=1
            continue
        for hit in hits: counts[hit]+=1
        if inclusion_hits: counts["inclusion_fragments"]+=1
        elif "skipping" in hits: counts["skipping_fragments"]+=1
def _full(files,run_dir,cfg,log):
    tool=cfg["alignment"]["full_tool"]; bam=run_dir/"aligned.bam"
    if tool=="STAR": cmd=["STAR","--runThreadN",str(cfg["alignment"]["threads"]),"--genomeDir",cfg["references"]["star_index"],"--readFilesIn",*map(str,files),"--outSAMtype","BAM","SortedByCoordinate","--outFileNamePrefix",str(run_dir)+"/"]; bam=run_dir/"Aligned.sortedByCoord.out.bam"; _execute(cmd,log)
    elif tool=="HISAT2":
        sam=run_dir/"aligned.sam"; cmd=["hisat2","-p",str(cfg["alignment"]["threads"]),"-x",cfg["references"]["hisat2_index"],"-U",str(files[0]),"-S",str(sam)]; _execute(cmd,log); _execute(["samtools","sort","-o",str(bam),str(sam)],log); sam.unlink()
    else: raise ValueError(f"unsupported full_tool: {tool}")
    return bam
def _regtools_counts(bam,run_dir,cfg,family,sample,acc,tier):
    bed=run_dir/"junctions.bed"; _execute(["regtools","junctions","extract","-s",cfg["alignment"]["strand_code"],str(bam),"-o",str(bed)],run_dir/"commands.log")
    observed={}
    with bed.open() as h:
        for line in h:
            f=line.rstrip().split("\t"); observed[(f[0],f[1],f[2],f[5])]=int(f[4])
    records=[]
    for event in cfg["events"]:
        if event["family"]!=family: continue
        e=_definitions(cfg)[event["id"]]; key=lambda k:(e["chrom"],e[f"{k}_start"],e[f"{k}_end"],e["strand"])
        vals={k:observed.get(key(k),0) for k in ("inclusion_up","inclusion_down","skipping")}
        # A valid empty extraction is a verified biological zero only when the command succeeded and definitions were verified.
        records.append(_record(family,sample,acc,event["id"],vals,tier,"none"))
    return records
def _record(family,sample,acc,event,vals,tier,label,anomaly_threshold=0.01):
    inclusion=int(vals.get("inclusion_fragments",int(vals["inclusion_up"])+int(vals["inclusion_down"])))
    skipping=int(vals.get("skipping_fragments",vals["skipping"])); anomalies=int(vals.get("anomalous_fragments",0))
    total=inclusion+skipping; classified=total+anomalies; anomaly_rate=anomalies/classified if classified else 0.0
    status="missing_junction" if total==0 and tier=="targeted" else ("verified_zero" if total==0 else "observed")
    return {"family":family,"sample":sample,"run_accession":acc,"event":event,**vals,"inclusion_fragments":inclusion,"skipping_fragments":skipping,"anomalous_fragments":anomalies,"anomaly_rate":f"{anomaly_rate:.12g}","marker_specificity_warning":str(anomaly_rate>anomaly_threshold).lower(),"informative_reads":total,"count_status":status,"quantification_tier":tier,"approximation_label":label}
def rebuild_combined(cfg,manifest=None):
    count_dir=Path(cfg["paths"]["work_dir"])/"counts"; combined=count_dir/"all_counts.tsv"
    with file_lock(count_dir/".all_counts.lock"):
        rows=[]
        allowed=None if manifest is None else {r["run_accession"] for r in read_tsv(manifest) if r["endpoint_eligible"]=="true"}
        for p in sorted(count_dir.glob("*.tsv")):
            if allowed is not None and p.stem not in allowed: continue
            if p.name!="all_counts.tsv": rows.extend(read_tsv(p))
        if not rows: raise ValueError("no per-run counts exist; refusing empty combined counts")
        write_tsv(combined,rows,COUNT_FIELDS)
    return combined
def run(config_path,manifest=None):
    cfg=validate_config(config_path); paths=cfg["paths"]; tier=cfg["alignment"]["tier"]; manifest=manifest or str(Path(paths["work_dir"])/"manifests/all_runs.frozen.tsv"); cp=Checkpoints(paths["checkpoints"],"quantify",config_path,cfg)
    for row in read_tsv(manifest):
        if row["endpoint_eligible"]!="true": continue
        acc=row["run_accession"]
        with _accession_lock(paths["scratch_dir"],acc):
            _run_accession(cfg,paths,tier,manifest,cp,row,acc)
    return rebuild_combined(cfg,manifest)
def _run_accession(cfg,paths,tier,manifest,cp,row,acc):
        target=Path(paths["work_dir"])/"counts"/f"{acc}.tsv"; fastq_dir=Path(paths["scratch_dir"])/"fastq"/acc; run_dir=Path(paths["scratch_dir"])/"alignment"/acc
        if cp.valid(acc,[target],[manifest],row): shutil.rmtree(fastq_dir,ignore_errors=True); shutil.rmtree(run_dir,ignore_errors=True); return
        # Use only the frozen manifest's final archive names. Broad globs also match
        # resumable/corrupt debris left by an interrupted earlier fetch.
        fastqs=[fastq_dir/Path(url).name for url in row["fastq_ftp"].split(";")]
        missing=[str(path) for path in fastqs if not path.is_file()]
        if missing: raise FileNotFoundError(f"missing final FASTQ files for {acc}: {', '.join(missing)}")
        sample=row["sample_alias"] or row["sample_accession"] or acc; run_dir.mkdir(parents=True,exist_ok=True)
        if tier=="targeted":
            counts=[]
            for event in cfg["events"]:
                if event["family"]==row["family"]:
                    progress=run_dir/f"{event['id']}.targeted-progress.json"
                    counts.append(_record(row["family"],sample,acc,event["id"],_targeted_fastq(fastqs,_definitions(cfg)[event["id"]],progress,cfg.get("checkpoint",{}).get("heartbeat_seconds",600)),tier,"fragment-level inclusion/skipping junction-marker count",cfg.get("psi",{}).get("marker_specificity_anomaly_rate_threshold",0.01)))
        elif tier=="full": counts=_regtools_counts(_full(fastqs,run_dir,cfg,run_dir/"commands.log"),run_dir,cfg,row["family"],sample,acc,tier)
        else: raise RuntimeError("only targeted and full tiers may emit primary labels")
        expected=sum(e["family"]==row["family"] for e in cfg["events"])
        if not counts or len(counts)!=expected: raise ValueError(f"{acc}: empty/incomplete per-run counts")
        write_tsv(target,counts,COUNT_FIELDS); cp.complete(acc,[target],inputs=[manifest],resolved_metadata=row,tier=tier)
        shutil.rmtree(run_dir,ignore_errors=True)
        if cfg["fetch"].get("cleanup_fastq_after_quantification",True): shutil.rmtree(fastq_dir,ignore_errors=True)
