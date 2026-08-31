import hashlib,json,multiprocessing,os,queue,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
from rawread_reanalysis import cli,compute,fetch,harmonize,quantify,resolve
from rawread_reanalysis.common import Checkpoints,read_tsv,require_space,sha256,validate_config,write_tsv

def _hold_accession_lock(scratch,accession,events):
 with fetch._accession_lock(scratch,accession):
  events.put("acquired")
  events.get()

def _acquire_accession_lock(scratch,accession,events):
 events.put("attempting")
 with fetch._accession_lock(scratch,accession): events.put("acquired")

class PipelineTest(unittest.TestCase):
 def setUp(self):
  self.t=tempfile.TemporaryDirectory(); b=Path(self.t.name); self.b=b
  assets=b/"assets"; assets.mkdir(); (assets/"star").mkdir(); (assets/"hisat2").mkdir()
  for p,text in ((assets/"ref.fa",">chr1\nACGT\n"),(assets/"ref.gtf","chr1\tt\texon\t1\t4\t.\t+\t.\tgene_id \"x\";\n"),(assets/"container.sif","container"),(assets/"star/i","star"),(assets/"hisat2/i","hisat")): p.write_text(text)
  event=assets/"events.tsv"; event.write_text("event\tchrom\tstrand\tinclusion_up_start\tinclusion_up_end\tinclusion_down_start\tinclusion_down_end\tskipping_start\tskipping_end\tinclusion_up_marker\tinclusion_down_marker\tskipping_marker\nSMN\tchr1\t+\t10\t20\t30\t40\t10\t40\tAAAACCCC\tGGGGTTTT\tCCCCGGGG\nDMD\tchr2\t-\t50\t60\t70\t80\t50\t80\tTTTTAAAA\tAAAAGGGG\tGGGGAAAA\n")
  def bh(p):
   if p.is_file(): return sha256(p)
   d=hashlib.sha256()
   for f in sorted(x for x in p.rglob("*") if x.is_file()): d.update(str(f.relative_to(p)).encode()+b"\0"+bytes.fromhex(sha256(f)))
   return d.hexdigest()
  self.cfg={"paths":{"work_dir":str(b/"work"),"scratch_dir":str(b/"scratch/ns"),"scratch_namespace":str(b/"scratch"),"checkpoints":str(b/"work/cp")},"fetch":{"space_multiplier":1,"minimum_free_bytes":1,"aggregate_scratch_cap_bytes":1000000,"cleanup_fastq_after_quantification":True},"alignment":{"tier":"targeted","threads":1,"full_tool":"STAR","strand_code":"2"},"references":{"build":"GRCh38","annotation_release":"GENCODE v99 TEST","fasta":str(assets/"ref.fa"),"gtf":str(assets/"ref.gtf"),"container":str(assets/"container.sif"),"star_index":str(assets/"star"),"hisat2_index":str(assets/"hisat2"),"event_junctions_tsv":str(event)},"asset_sha256":{"reference_fasta":bh(assets/"ref.fa"),"annotation_gtf":bh(assets/"ref.gtf"),"container":bh(assets/"container.sif"),"star_index":bh(assets/"star"),"hisat2_index":bh(assets/"hisat2"),"event_junctions":bh(event)},"psi":{"minimum_informative_reads":2},"tier_concordance":{"passed":False},"studies":{"smn":{"project":"P","sample_map":{"S_T":{"intervention":"drug","condition":"treated","contrast":"main","pair_set_id":"donor1","batch_covariates":{"donor":"1","tissue":"cell"}},"S_C":{"intervention":"vehicle","condition":"control","contrast":"main","pair_set_id":"donor1","batch_covariates":{"donor":"1","tissue":"cell"}}},"matched_theta_proxy":False,"confound_flags":[]},"dmd":{"project":"D","sample_map":{"D_T":{"intervention":"U7","condition":"treated","contrast":"main","pair_set_id":"donor2-muscle","batch_covariates":{"donor":"2","tissue":"muscle"}},"D_C":{"intervention":"untreated","condition":"control","contrast":"main","pair_set_id":"donor2-muscle","batch_covariates":{"donor":"2","tissue":"muscle"}}},"matched_theta_proxy":False,"confound_flags":[]}},"events":[{"id":"SMN","family":"smn","primary":True,"therapeutic_polarity_multiplier":1},{"id":"DMD","family":"dmd","primary":True,"therapeutic_polarity_multiplier":-1}]}
  self.config=b/"config.json"; self.save()
 def tearDown(self): self.t.cleanup()
 def save(self): self.config.write_text(json.dumps(self.cfg))
 def counts(self,rows):
  p=self.b/"counts.tsv"; fields=["family","sample","run_accession","event","inclusion_up","inclusion_down","skipping","informative_reads","count_status","quantification_tier","approximation_label"]; write_tsv(p,rows,fields); return p
 def row(self,f,a,e,u,d,s,tier="targeted",status="observed"): return {"family":f,"sample":a,"run_accession":a,"event":e,"inclusion_up":u,"inclusion_down":d,"skipping":s,"informative_reads":u+d+s,"count_status":status,"quantification_tier":tier,"approximation_label":"test"}
 def test_asymmetric_smn_and_dmd_signs_polarity_and_pairing(self):
  p=self.counts([self.row("smn","S_C","SMN",2,2,8),self.row("smn","S_T","SMN",8,8,2),self.row("dmd","D_C","DMD",8,8,2),self.row("dmd","D_T","DMD",2,2,8)])
  out=read_tsv(compute.run(self.config,p)); smn=next(x for x in out if x["run_accession"]=="S_T"); dmd=next(x for x in out if x["run_accession"]=="D_T")
  self.assertGreater(float(smn["dpsi"]),0); self.assertLess(float(dmd["dpsi"]),0); self.assertEqual(dmd["therapeutic_polarity_multiplier"],"-1"); self.assertGreater(float(dmd["therapeutic_effect"]),0)
 def test_low_coverage_and_missing_junction_suppress_labels(self):
  p=self.counts([self.row("smn","S_C","SMN",0,0,0,status="verified_zero"),self.row("smn","S_T","SMN",1,0,0)])
  out=read_tsv(compute.run(self.config,p)); self.assertTrue(all(x["psi"]=="" and x["dpsi"]=="" and x["direction"]=="MISSING" for x in out))
 def test_ambiguous_regex_rejected(self):
  self.cfg["studies"]["smn"]["sample_map"]={}; self.cfg["studies"]["smn"]["sample_rules"]=[{"regex":"treated","intervention":"a","condition":"treated","contrast":"x","pair_set_id":"p"},{"regex":"treated","intervention":"b","condition":"treated","contrast":"x","pair_set_id":"p"}]; self.save()
  with self.assertRaisesRegex(ValueError,"matched 2"): compute.run(self.config,self.counts([self.row("smn","untreated_treated","SMN",2,2,2)]))
 def test_empty_event_and_placeholder_or_bad_hash_fail_closed(self):
  event=Path(self.cfg["references"]["event_junctions_tsv"]); event.write_text(event.read_text().splitlines()[0]+"\n"); self.cfg["asset_sha256"]["event_junctions"]=sha256(event); self.save()
  with self.assertRaisesRegex(ValueError,"empty"): validate_config(self.config)
  self.cfg["asset_sha256"]["event_junctions"]="REQUIRED-INPUT"; self.save()
  with self.assertRaisesRegex(ValueError,"SHA-256"): validate_config(self.config)
 def test_targeted_mode_needs_no_genome_or_index_assets(self):
  for key in ("reference_fasta","annotation_gtf","star_index","hisat2_index","container"): self.cfg["asset_sha256"][key]="REQUIRED-INPUT"
  for key in ("fasta","gtf","star_index","hisat2_index","container"): self.cfg["references"][key]="/path/to/REQUIRED-INPUT"
  self.save(); cfg=validate_config(self.config)
  self.assertEqual(set(cfg["_verified_asset_hashes"]),{"event_junctions"})
 def test_one_missing_junction_definition_fails_closed(self):
  event=Path(self.cfg["references"]["event_junctions_tsv"]); lines=event.read_text().splitlines(); event.write_text("\n".join(lines[:2])+"\n"); self.cfg["asset_sha256"]["event_junctions"]=sha256(event); self.save()
  with self.assertRaisesRegex(ValueError,"missing=.*DMD"): validate_config(self.config)
 def test_cross_tier_pooling_rejected(self):
  p=self.counts([self.row("smn","S_C","SMN",2,2,2,"targeted"),self.row("smn","S_T","SMN",2,2,2,"full")])
  with self.assertRaisesRegex(ValueError,"cross-tier"): compute.run(self.config,p)
 def test_targeted_fragment_and_tier_equivalent_definition(self):
  defs=quantify._definitions(validate_config(self.config)); fq=self.b/"x.fastq"; fq.write_text("@x\nAAAACCCC\n+\nIIIIIIII\n")
  self.assertEqual(quantify._targeted_fastq([fq],defs["SMN"])["inclusion_up"],1)
  fq.write_text("@u\nAAAACCCC\n+\nIIIIIIII\n@d\nGGGGTTTT\n+\nIIIIIIII\n@s\nCCCCGGGG\n+\nIIIIIIII\n")
  targeted=quantify._targeted_fastq([fq],defs["SMN"]); run_dir=self.b/"align"; run_dir.mkdir(); bam=run_dir/"x.bam"; bam.write_text("x")
  def fake_execute(command,log):
   Path(command[-1]).write_text("chr1\t10\t20\tJ\t1\t+\nchr1\t30\t40\tJ\t1\t+\nchr1\t10\t40\tJ\t1\t+\n")
  with patch.object(quantify,"_execute",side_effect=fake_execute): full=quantify._regtools_counts(bam,run_dir,validate_config(self.config),"smn","S_T","S_T","full")[0]
  self.assertEqual({k:targeted[k] for k in ("inclusion_up","inclusion_down","skipping")},{"inclusion_up":int(full["inclusion_up"]),"inclusion_down":int(full["inclusion_down"]),"skipping":int(full["skipping"])})
 def test_real_shaped_read_spanning_both_inclusion_junctions_counts_once(self):
  event=quantify._definitions(validate_config(self.config))["SMN"]; fq=self.b/"both.fastq"
  seq="TGCATGCA"+event["inclusion_up_marker"]+"ACGT"*14+event["inclusion_down_marker"]+"GATTACA"
  fq.write_text(f"@both\n{seq}\n+\n{'I'*len(seq)}\n")
  counts=quantify._targeted_fastq([fq],event)
  self.assertEqual((counts["inclusion_up"],counts["inclusion_down"],counts["inclusion_fragments"],counts["anomalous_fragments"]),(1,1,1,0))
 def test_inclusion_skipping_double_match_is_discarded_and_tallied(self):
  event=quantify._definitions(validate_config(self.config))["SMN"]; fq=self.b/"anomaly.fastq"
  seq="ACGT"*10+event["inclusion_up_marker"]+"TGCA"*8+event["skipping_marker"]+"ACGT"*5
  fq.write_text(f"@bad\n{seq}\n+\n{'I'*len(seq)}\n")
  counts=quantify._targeted_fastq([fq],event); record=quantify._record("smn","S","R","SMN",counts,"targeted","test",0.5)
  self.assertEqual((record["inclusion_fragments"],record["skipping_fragments"],record["anomalous_fragments"]),(0,0,1)); self.assertEqual(record["count_status"],"missing_junction"); self.assertEqual(record["marker_specificity_warning"],"true")
 def test_fragment_dedup_psi_math_with_double_junction_inclusion(self):
  event=quantify._definitions(validate_config(self.config))["SMN"]; fq=self.b/"mix.fastq"
  seqs=[event["inclusion_up_marker"]+"ACGT"*14+event["inclusion_down_marker"],event["inclusion_up_marker"]+"A"*92,event["skipping_marker"]+"C"*92]
  fq.write_text("".join(f"@r{i}\n{s}\n+\n{'I'*len(s)}\n" for i,s in enumerate(seqs)))
  record=quantify._record("smn","S","R","SMN",quantify._targeted_fastq([fq],event),"targeted","test")
  self.assertEqual((record["inclusion_fragments"],record["skipping_fragments"],record["informative_reads"]),(2,1,3))
  record.update({"sample":"S_T","run_accession":"S_T"}); count_file=self.b/"dedup-counts.tsv"; write_tsv(count_file,[record],quantify.COUNT_FIELDS)
  result=read_tsv(compute.run(self.config,count_file))[0]; self.assertAlmostEqual(float(result["psi"]),66.6667)
 def test_paired_mates_covering_different_inclusion_junctions_count_once(self):
  event=quantify._definitions(validate_config(self.config))["SMN"]; left=self.b/"pair_R1.fastq"; right=self.b/"pair_R2.fastq"
  for path,suffix,marker in ((left,"1",event["inclusion_up_marker"]),(right,"2",event["inclusion_down_marker"])):
   seq="ACGT"*12+marker+"TGCA"*12; path.write_text(f"@fragment/"+suffix+f"\n{seq}\n+\n{'I'*len(seq)}\n")
  counts=quantify._targeted_fastq([left,right],event)
  self.assertEqual((counts["inclusion_up"],counts["inclusion_down"],counts["inclusion_fragments"]),(1,1,1))
 def test_checkpoint_direct_input_staleness(self):
  cfg=validate_config(self.config); inp=self.b/"in"; out=self.b/"out"; inp.write_text("a"); out.write_text("x"); cp=Checkpoints(self.b/"cp","s",self.config,cfg); cp.complete("x",[out],[inp]); self.assertTrue(cp.valid("x",[out],[inp])); inp.write_text("b"); self.assertFalse(cp.valid("x",[out],[inp]))
 def test_cleanup_reconciliation_on_resume(self):
  cfg=validate_config(self.config); manifest=self.b/"m.tsv"; row={"family":"smn","run_accession":"S_T","sample_alias":"S_T","sample_accession":"","endpoint_eligible":"true"}; write_tsv(manifest,[row],row.keys()); target=self.b/"work/counts/S_T.tsv"; write_tsv(target,[self.row("smn","S_T","SMN",2,2,2)],quantify.COUNT_FIELDS); cp=Checkpoints(self.b/"work/cp","quantify",self.config,cfg); cp.complete("S_T",[target],[manifest],resolved_metadata=row); fq=self.b/"scratch/ns/fastq/S_T"; fq.mkdir(parents=True); (fq/"x.fastq").write_text("x"); quantify.run(self.config,manifest); self.assertFalse(fq.exists())
 def test_quantify_ignores_nonfinal_fastq_debris(self):
  fastq_dir=self.b/"scratch/ns/fastq/S_T"; fastq_dir.mkdir(parents=True)
  final=fastq_dir/"S_T.fastq"; final.write_text("@x\nAAAACCCC\n+\nIIIIIIII\n")
  (fastq_dir/"S_T.fastq.partial.corrupt").write_text("debris")
  row={"family":"smn","run_accession":"S_T","sample_alias":"S_T","sample_accession":"","endpoint_eligible":"true","fastq_ftp":"ftp/x/S_T.fastq"}
  manifest=self.b/"m.tsv"; write_tsv(manifest,[row],row.keys())
  out=read_tsv(quantify.run(self.config,manifest))
  self.assertEqual((len(out),out[0]["inclusion_up"]),(1,"1"))
 def test_aggregate_quota_cap(self):
  with self.assertRaisesRegex(RuntimeError,"aggregate cap"): require_space(self.b/"quota",10,5)
 def test_accession_lock_excludes_competing_process(self):
  ctx=multiprocessing.get_context("spawn"); holder_events=ctx.Queue(); contender_events=ctx.Queue()
  holder=ctx.Process(target=_hold_accession_lock,args=(self.b/"scratch/ns","S_T",holder_events)); holder.start()
  self.assertEqual(holder_events.get(timeout=5),"acquired")
  contender=ctx.Process(target=_acquire_accession_lock,args=(self.b/"scratch/ns","S_T",contender_events)); contender.start()
  self.assertEqual(contender_events.get(timeout=5),"attempting")
  with self.assertRaises(queue.Empty): contender_events.get(timeout=0.25)
  holder_events.put("release")
  self.assertEqual(contender_events.get(timeout=5),"acquired")
  holder.join(5); contender.join(5)
  self.assertEqual((holder.exitcode,contender.exitcode),(0,0))
 def test_process_isolates_accession_failure_and_marks_partial(self):
  manifest=self.b/"two.tsv"; rows=[{"family":"smn","run_accession":acc,"endpoint_eligible":"true"} for acc in ("BAD","GOOD")]; write_tsv(manifest,rows,rows[0].keys())
  processed=[]
  def fake_fetch(config,shard):
   acc=read_tsv(shard)[0]["run_accession"]
   if acc=="BAD": raise RuntimeError("empty reply")
  def fake_quantify(config,shard): processed.append(read_tsv(shard)[0]["run_accession"])
  with patch.object(cli.fetch,"run",side_effect=fake_fetch), patch.object(cli.quantify,"run",side_effect=fake_quantify), patch.object(cli.quantify,"rebuild_combined") as rebuild:
   cli.main(["process","--config",str(self.config),"--manifest",str(manifest)])
  self.assertEqual(processed,["GOOD"]); rebuild.assert_called_once()
  failures=read_tsv(self.b/"work/fetch_failures.tsv"); self.assertEqual([(x["run_accession"],x["stage"]) for x in failures],[("BAD","fetch")])
  marker=cli.mark_stage_a_completion(validate_config(self.config),self.b/"run")
  self.assertEqual(marker.name,"STAGE_A_PARTIAL"); self.assertTrue(marker.exists()); self.assertFalse((self.b/"run/STAGE_A_COMPLETE").exists())
 def test_offline_end_to_end_and_resume(self):
  # Reduce the synthetic project to one event, while retaining fully verified local assets.
  self.cfg["studies"]={"smn":self.cfg["studies"]["smn"]}; self.cfg["events"]=[self.cfg["events"][0]]
  event=Path(self.cfg["references"]["event_junctions_tsv"]); event.write_text("\n".join(event.read_text().splitlines()[:2])+"\n"); self.cfg["asset_sha256"]["event_junctions"]=sha256(event)
  archive=self.b/"archive"; archive.mkdir(); treated=archive/"t.fastq"; control=archive/"c.fastq"
  treated.write_text("@t1\nAAAACCCC\n+\nIIIIIIII\n@t2\nGGGGTTTT\n+\nIIIIIIII\n@t3\nAAAACCCC\n+\nIIIIIIII\n")
  control.write_text("@c1\nCCCCGGGG\n+\nIIIIIIII\n@c2\nCCCCGGGG\n+\nIIIIIIII\n@c3\nAAAACCCC\n+\nIIIIIIII\n")
  import hashlib as hl
  md=self.b/"metadata"; md.mkdir(); header="study_accession\trun_accession\tsample_accession\tsample_alias\texperiment_title\tfastq_ftp\tfastq_md5\tfastq_bytes\tlibrary_layout\tlibrary_strategy\n"
  body="".join(f"P\t{acc}\t{acc}\t{acc}\tx\t{p.as_uri()}\t{hl.md5(p.read_bytes()).hexdigest()}\t{p.stat().st_size}\tSINGLE\tRNA-Seq\n" for acc,p in (("S_T",treated),("S_C",control)))
  (md/"P.tsv").write_text(header+body); self.save()
  manifest=resolve.run(self.config,md)
  debris=self.b/"scratch/ns/fastq/S_T/t.fastq.partial.corrupt"; debris.parent.mkdir(parents=True); debris.write_text("old")
  fetch.run(self.config,manifest); self.assertFalse(debris.exists())
  debris.write_text("orphaned-after-checkpoint")
  fetch.run(self.config,manifest); self.assertFalse(debris.exists())
  counts=quantify.run(self.config,manifest); effects=compute.run(self.config,counts); tidy=harmonize.run(self.config,effects)
  row=next(x for x in read_tsv(tidy) if x["condition"]=="treated"); self.assertGreater(float(row["magnitude"]),0); before=os.stat(tidy).st_mtime_ns; harmonize.run(self.config,effects); self.assertEqual(before,os.stat(tidy).st_mtime_ns)

if __name__=="__main__":unittest.main()
