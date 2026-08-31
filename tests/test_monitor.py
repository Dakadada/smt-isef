import importlib.util,json,tempfile,time,unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
ROOT=Path(__file__).resolve().parents[1]; spec=importlib.util.spec_from_file_location("monitor",ROOT/"slurm/monitor.py"); monitor=importlib.util.module_from_spec(spec); spec.loader.exec_module(monitor)
def args(root,submit=True): return SimpleNamespace(state=str(root/"state.json"),job_id="1",script="job.sbatch",checkpoints=str(root/"cp"),max_retries=1,stale_seconds=1,submit=submit)
class MonitorTest(unittest.TestCase):
 def test_corrupt_state_fails_closed(self):
  with tempfile.TemporaryDirectory() as t:
   r=Path(t); (r/"state.json").write_text("bad")
   with self.assertRaisesRegex(RuntimeError,"corrupt monitor state"): monitor.evaluate(args(r))
 def test_stale_running_cancels_but_never_submits(self):
  with tempfile.TemporaryDirectory() as t:
   r=Path(t); (r/"cp").mkdir(); hb=Path("data/reanalysis/heartbeats/1.json"); hb.parent.mkdir(parents=True,exist_ok=True); hb.write_text("{}"); old=time.time()-10; import os; os.utime(hb,(old,old))
   calls=[]
   def cmd(*a): calls.append(a); return SimpleNamespace(stdout="RUNNING\n" if a[0]=="squeue" else "",returncode=0,stderr="")
   try:
    with patch.object(monitor,"command",side_effect=cmd): monitor.evaluate(args(r))
    self.assertIn(("scancel","1"),calls); self.assertFalse(any(c[0]=="sbatch" for c in calls))
   finally:
    if hb.exists(): hb.unlink()
 def test_durable_retry_cap_survives_state_loss(self):
  with tempfile.TemporaryDirectory() as t:
   r=Path(t); (r/"cp").mkdir(); log=r/"state.json.attempts.jsonl"; log.write_text(json.dumps({"action":"submit","job_id":"old"})+"\n")
   calls=[]
   def cmd(*a): calls.append(a); return SimpleNamespace(stdout="PREEMPTED\n" if a[0]=="sacct" else "",returncode=0,stderr="")
   with patch.object(monitor,"command",side_effect=cmd): state=monitor.evaluate(args(r))
   self.assertEqual(state["retries"],1); self.assertFalse(any(c[0]=="sbatch" for c in calls))
if __name__=="__main__":unittest.main()
