"""Read-only global J5 completion audit. No directory discovery or recomputation.

The pinned manifest and production shard partition determine every expected path.
A successful audit establishes artifact completeness, not a biological GO.
"""
import argparse
from collections import Counter
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import tempfile

import pyarrow as pa
import pyarrow.parquet as pq
from conference.model_v2 import j5_rbp_context as j5
from conference.model_v2.jobs.j5_benchmark import profile_worker


def audit(*, index, j6_genome_dir, manifest_sha256, full_out, shard_count=64,
          expected_windows=52424, reference=None):
    if shard_count<1: raise ValueError('positive shard count required')
    sources=j5.genome_queries(j6_genome_dir,0,1,expected_windows=expected_windows,
                              manifest_sha256=manifest_sha256)
    index_state=index.state()
    gaps=[]; peaks=Counter(); knockdown=Counter(); controls=Counter()
    valid=0; rows=0; output_bytes=0
    for shard in range(shard_count):
        for source in sources[len(sources)*shard//shard_count:len(sources)*(shard+1)//shard_count]:
            destination=j5.safe_path(full_out)/f'shard-{shard}'/'features'/j5.feature_name(source)
            try:
                before=source.stat()
                source_meta=pq.read_metadata(source)
                query=json.loads(source.with_suffix('.done').read_text())
                contract=j5.feature_contract(source,index,reference)
                if (not isinstance(query,dict) or type(query.get('rows')) is not int
                        or query['rows']!=source_meta.num_rows or query.get('sha256')!=contract['query_sha256']):
                    raise ValueError('J6 marker count/hash mismatch')
                done=j5._valid_done(destination,contract)
                if done is None: raise ValueError('missing or invalid output/marker hash/count')
                if done['rows']!=query['rows']: raise ValueError('J5/J6 row count mismatch')
                parquet=pq.ParquetFile(destination)
                if parquet.schema_arrow!=j5._schema(): raise ValueError('output schema mismatch')
                local_peaks=Counter(); local_kd=Counter(); local_controls=Counter()
                for batch in parquet.iter_batches(batch_size=128,columns=['rbp_context','matched_controls']):
                    for row in batch.to_pylist():
                        contexts=json.loads(row['rbp_context'])
                        if len(contexts)!=len(j5.RBPS)*len(j5.CELLS): raise ValueError('assay count mismatch')
                        for context in contexts:
                            local_peaks[context['status']]+=1
                            local_kd[context['knockdown']['status']]+=1
                        local_controls[json.loads(row['matched_controls'])['status']]+=1
                # Detect replacement during streaming, including changes to .done.
                if j5._valid_done(destination,contract)!=done: raise ValueError('output changed during gate')
                if source.stat()!=before: raise ValueError('J6 source changed during gate')
                peaks.update(local_peaks); knockdown.update(local_kd); controls.update(local_controls)
                valid+=1; rows+=done['rows']; output_bytes+=destination.stat().st_size
            except (OSError,ValueError,KeyError,TypeError,pa.ArrowException) as exc:
                gaps.append(dict(source=str(source),output=str(destination),error=str(exc)))
    if index.state()!=index_state: raise ValueError('index changed during gate')
    return dict(ok=not gaps,expected_windows=len(sources),valid_windows=valid,gaps=gaps,
        rows=rows,output_bytes=output_bytes,peak_status_counts=dict(peaks),
        knockdown_status_counts=dict(knockdown),control_status_counts=dict(controls),
        index_file_status_counts=dict(Counter(m['status'] for m in index.manifests())),
        manifest_sha256=manifest_sha256,index_sha256=index.contract(),reference=reference,
        shard_count=shard_count,utc_time=datetime.now(timezone.utc).isoformat(),
        gate_sha256=j5.digests(__file__)[0],
        limitation='Artifact audit only; submit with afterok on the full array to also require every Slurm task exit to succeed.')


def write_report(path, report):
    path=j5.safe_path(path); path.parent.mkdir(parents=True,exist_ok=True)
    sidecar=path.with_suffix(path.suffix+'.sha256')
    sidecar.unlink(missing_ok=True)
    j5.atomic_json(path,report)
    sha=j5.digests(path)[0]
    fd,temp=tempfile.mkstemp(dir=path.parent,suffix='.tmp')
    try:
        with os.fdopen(fd,'w') as f:
            f.write(sha+'\n'); f.flush(); os.fsync(f.fileno())
        os.replace(temp,sidecar)
    finally:
        if os.path.exists(temp): os.unlink(temp)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index-dir',type=Path,required=True)
    parser.add_argument('--assembly-lookup',type=Path,required=True)
    parser.add_argument('--lookup-sha256',required=True)
    parser.add_argument('--j6-genome-dir',type=Path,required=True)
    parser.add_argument('--j6-manifest-sha256',required=True)
    parser.add_argument('--full-out',type=Path,required=True)
    parser.add_argument('--report',type=Path,required=True)
    parser.add_argument('--shard-count',type=int,default=64)
    parser.add_argument('--fasta',type=Path)
    parser.add_argument('--genes',type=Path)
    args=parser.parse_args(argv)
    if args.shard_count<1: parser.error('--shard-count must be positive')
    try:
        with profile_worker(args) as (index,_,reference):
            report=audit(index=index,j6_genome_dir=args.j6_genome_dir,
                manifest_sha256=args.j6_manifest_sha256,full_out=args.full_out,
                shard_count=args.shard_count,reference=reference)
    except (OSError,ValueError,KeyError,TypeError,pa.ArrowException) as exc:
        report=dict(ok=False,error=str(exc),utc_time=datetime.now(timezone.utc).isoformat())
    write_report(args.report,report)
    print(j5.dumps({k:v for k,v in report.items() if k!='gaps'}),flush=True)
    return 0 if report['ok'] else 1


if __name__=='__main__': raise SystemExit(main())
