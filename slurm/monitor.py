#!/usr/bin/env python3
import argparse,fcntl,hashlib,json,os,shutil,subprocess,tempfile,time
from pathlib import Path
def atomic(path,data):
    path.parent.mkdir(parents=True,exist_ok=True); fd,tmp=tempfile.mkstemp(dir=path.parent,prefix=".monitor.")
    with os.fdopen(fd,"w") as h: json.dump(data,h,indent=2); h.flush(); os.fsync(h.fileno())
    os.replace(tmp,path)
def digest(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for block in iter(lambda:f.read(1048576),b""):h.update(block)
    return h.hexdigest()
def quarantine_invalid(root):
    bad=[]
    for cp in Path(root).glob("*/*.json"):
        try: state=json.loads(cp.read_text())
        except Exception: state=None
        if state is None: raise RuntimeError(f"corrupt checkpoint state: {cp}")
        valid=state.get("status")=="complete" and all(Path(n).is_file() and digest(n)==v for n,v in state.get("outputs",{}).items())
        if not valid:
            q=Path(root)/"quarantine"/f"{cp.parent.name}-{cp.name}-{int(time.time())}"; q.parent.mkdir(parents=True,exist_ok=True); shutil.move(cp,q); bad.append(str(q))
    return bad
def command(*args):return subprocess.run(args,text=True,capture_output=True,check=False)
def load_attempts(log):
    if not log.exists(): return []
    attempts=[]
    for number,line in enumerate(log.read_text().splitlines(),1):
        try: attempts.append(json.loads(line))
        except json.JSONDecodeError as exc: raise RuntimeError(f"corrupt durable retry log line {number}") from exc
    return attempts
def append_attempt(log,event):
    log.parent.mkdir(parents=True,exist_ok=True)
    with log.open("a") as h: h.write(json.dumps(event,sort_keys=True)+"\n"); h.flush(); os.fsync(h.fileno())
def evaluate(args):
    state_path=Path(args.state); attempts_path=state_path.with_suffix(state_path.suffix+".attempts.jsonl")
    if state_path.exists():
        try: state=json.loads(state_path.read_text())
        except Exception as exc: raise RuntimeError(f"corrupt monitor state: {state_path}; refusing retry reset") from exc
    else: state={"job_id":args.job_id,"history":[]}
    attempts=load_attempts(attempts_path); active=command("squeue","-h","-j",state["job_id"],"-o","%T").stdout.strip(); accounting=command("sacct","-n","-j",state["job_id"],"--format=State","-P").stdout.strip().splitlines(); status=active or (accounting[0].split("|")[0] if accounting else "UNKNOWN")
    heartbeat=Path("data/reanalysis/heartbeats")/f"{state['job_id']}.json"; stale=heartbeat.exists() and time.time()-heartbeat.stat().st_mtime>args.stale_seconds; bad=quarantine_invalid(args.checkpoints)
    running=status.upper() in ("RUNNING","PENDING","CONFIGURING","COMPLETING")
    retryable=any(x in status.upper() for x in ("PREEMPTED","NODE_FAIL","CANCELLED","TIMEOUT"))
    if stale and running and args.submit:
        # Cancellation is recorded; resubmission is forbidden until a later poll confirms terminal state.
        result=command("scancel",state["job_id"])
        if result.returncode: raise RuntimeError(result.stderr)
        append_attempt(attempts_path,{"time":int(time.time()),"action":"cancel_stale","job_id":state["job_id"]})
    elif retryable and len([a for a in attempts if a.get("action")=="submit"])<args.max_retries and args.submit:
        result=command("sbatch","--parsable",args.script)
        if result.returncode: raise RuntimeError(result.stderr)
        new=result.stdout.strip().split(";")[0]; append_attempt(attempts_path,{"time":int(time.time()),"action":"submit","from_job_id":state["job_id"],"job_id":new}); state["job_id"]=new
    state["retries"]=len([a for a in load_attempts(attempts_path) if a.get("action")=="submit"]); state["history"].append({"time":int(time.time()),"status":status,"stale":stale,"quarantined":bad}); state["history"]=state["history"][-1000:]; atomic(state_path,state); return state
def main():
    p=argparse.ArgumentParser(); p.add_argument("--job-id",required=True); p.add_argument("--script",default="slurm/run_preemptable.sbatch"); p.add_argument("--state",default="data/reanalysis/monitor.json"); p.add_argument("--checkpoints",default="data/reanalysis/checkpoints"); p.add_argument("--max-retries",type=int,default=3); p.add_argument("--stale-seconds",type=int,default=900); p.add_argument("--watch",action="store_true"); p.add_argument("--interval",type=int,default=60); p.add_argument("--submit",action="store_true"); args=p.parse_args(); lock=Path(args.state+".lock"); lock.parent.mkdir(parents=True,exist_ok=True)
    with lock.open("w") as h:
        try: fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError: raise SystemExit("another monitor owns this state")
        while True:
            print(json.dumps(evaluate(args),indent=2),flush=True)
            if not args.watch:break
            time.sleep(max(10,args.interval))
if __name__=="__main__":main()
