"""Bounded J5 throughput measurements; temporary outputs, never production artifacts.

Default: 128-site hot synthetic and legacy real sentinel timing. The legacy
adapter is explicitly timing-only: it reads old intervals without blessing old
provenance or changing the production reader's fail-closed integrity checks.
--index-dir instead profiles a genuinely verified current index on representative
J6 windows; run that mode in Slurm (j5_benchmark.sbatch).
"""
import argparse
import csv
from contextlib import nullcontext, contextmanager
import gzip
import json
import math
import os
from pathlib import Path
import resource
import re
import sqlite3
import subprocess
import sys
import tempfile
import time
from types import SimpleNamespace

import pyarrow as pa
import pyarrow.parquet as pq
from conference.model_v2 import j5_rbp_context as j5


def measure(index, source, destination, label, *, fetch=None, reference_contract=None,
            require_query_done=False):
    timing=destination.with_suffix('.timing.json')
    options=dict(fetch=fetch,reference_contract=reference_contract,require_query_done=require_query_done)
    profile=dict(reference=reference_contract,benchmark_sha256=j5.digests(__file__)[0])
    if timing.exists():
        # Verify the parquet/index contract before reusing checkpointed timing.
        previous=json.loads(timing.read_text())
        if any(previous.get(k)!=v for k,v in profile.items()): raise ValueError('benchmark contract changed')
        result=j5.write_features(source,destination,index,**options)
        if previous['output_sha256']!=result['sha256']: raise ValueError('benchmark checkpoint changed')
        return previous | dict(reused_timing=True)
    # Preemption may leave a completed parquet but no timing record. Recompute
    # it rather than report the near-zero duration of a checkpoint cache hit.
    destination.with_suffix('.done').unlink(missing_ok=True)
    statements=[0]
    def trace(_): statements[0]+=1
    index.db.set_trace_callback(trace)
    start=time.perf_counter()
    result=j5.write_features(source,destination,index,**options)
    seconds=time.perf_counter()-start
    index.db.set_trace_callback(None)
    measured=dict(label=label,sites=result['rows'],seconds=seconds,sites_per_second=result['rows']/seconds,
        sql_statements=statements[0],output_bytes=destination.stat().st_size,
        output_sha256=result['sha256'],maxrss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    artifacts=[destination,destination.with_suffix('.done')]
    measured.update(profile)
    measured.update(disk_bytes=sum(p.stat().st_size for p in artifacts),
                    disk_allocated_bytes=sum(p.stat().st_blocks*512 for p in artifacts),
                    index_bytes=index.path.stat().st_size)
    j5.atomic_json(timing,measured)
    return measured


def make_queries(path, count, start):
    rows=[dict(site_id=f'bench-{i}',chrom='chr5',start1=start+i*10,end1=start+i*10+9,
               strand='+',assembly='GRCh38',annotation='GENCODE v50') for i in range(count)]
    pq.write_table(pa.Table.from_pylist(rows),path)
    return rows


def synthetic(directory, count):
    root=directory/'encode'; rows=[]
    for rbp in j5.RBPS:
        for cell in j5.CELLS:
            path=root/'eclip_peaks'/rbp/cell/'ENCFF000AAA.bed.gz'
            # Distinct, stable fixture accessions, two peaks for each of 42 assays.
            path=path.with_name(f'ENCFF{len(rows):06d}.bed.gz'); path.parent.mkdir(parents=True)
            with gzip.open(path,'wt') as f: f.write('chr5\t0\t10\ta\t0\t+\nchr5\t100000\t100010\tb\t0\t+\n')
            rows.append(dict(kind='eclip_peaks',rbp=rbp,cell=cell,file=path.name,accession=path.name.split('.')[0],
                assembly='GRCh38',file_format='bed',format_type='narrowPeak',output_type='peaks',
                md5=j5.digests(path)[1],lookup='ok',dataset='/experiments/ENCSR000AAA/',
                target_label=rbp,biosample_term_name=cell,identity_status='verified'))
    lookup_path=directory/'lookup.tsv'
    with lookup_path.open('w') as f:
        writer=csv.DictWriter(f,rows[0].keys(),delimiter='\t'); writer.writeheader(); writer.writerows(rows)
    from conference.model_v2.jobs.j5_encode_lookup import write_provenance
    write_provenance(lookup_path,rows)
    lookup=j5.AssemblyLookup(lookup_path,expected_sha256=j5.digests(lookup_path)[0])
    out=directory/'index'; out.mkdir()
    with j5.ContextIndex(out/'context.sqlite',lookup=lookup) as index:
        j5._build_index(SimpleNamespace(encode_root=root,genes=None),index,out,{})
        source=directory/'queries.parquet'; queries=make_queries(source,count,1000)
        index.contexts_bulk(queries[:1])  # hot interval + availability cache, outside timer
        return measure(index,source,directory/'features.parquet','hot_synthetic_42_assays_no_KD_no_GC')


class LegacyTimingIndex(j5.ContextIndex):
    """Read-only timing adapter, intentionally unusable by the production CLI."""
    def __init__(self,path):
        self.path=j5.safe_path(path)
        self.readonly=True; self._sealed=None; self._bulk=False; self._memory=None
        self._cached_contract=None; self._availability_cache={}; self._cache_version=None
        self.db=sqlite3.connect(self.path.as_uri()+'?mode=ro',uri=True)
        self.db.row_factory=sqlite3.Row
        self.db.execute('PRAGMA query_only=ON')


def legacy_sentinel(directory,path,count):
    original_sha=j5.digests(path)[0]
    # Only the 11 MB known sentinel is allowed in this local timing path.
    if path.stat().st_size>64*1024*1024: raise ValueError('use Slurm for larger indexes')
    with LegacyTimingIndex(path) as index:
        source=directory/'queries.parquet'; queries=make_queries(source,count,70076584)
        index.contexts_bulk(queries[:1]); index.contract()
        result=measure(index,source,directory/'features.parquet','legacy_real_sentinel_timing_only_no_GC')
        result.update(sqlite_sha256=original_sha,intervals=index.db.execute('select count(*) from intervals').fetchone()[0])
    if j5.digests(path)[0]!=original_sha: raise ValueError('sentinel changed during benchmark')
    return result


def select_windows(sources, *, high=12, median=4):
    """Rank only completion manifests named by the pinned J6 inventory."""
    if high<1 or median<1: raise ValueError('positive high/median sample counts required')
    counts=[]
    for path in sources:
        marker=json.loads(path.with_suffix('.done').read_text())
        n=marker.get('rows')
        if type(n) is not int or n<0: raise ValueError('invalid J6 manifest row count')
        counts.append((path,n))
    ranked=sorted(counts,key=lambda pair:(-pair[1],pair[0].name))
    top=ranked[:high]
    ordered=sorted(counts,key=lambda pair:(pair[1],pair[0].name))
    # Consecutive ranks centered on the population median, excluding top-N.
    center=(len(ordered)-1)/2
    excluded={p for p,_ in top}
    middle=sorted([(i,p,n) for i,(p,n) in enumerate(ordered) if p not in excluded],
                  key=lambda item:(abs(item[0]-center),item[0]))[:median]
    middle=sorted(middle)
    return [(p,n,'high') for p,n in top]+[(p,n,'median') for _,p,n in middle]


@contextmanager
def profile_worker(args):
    lookup=j5.AssemblyLookup(args.assembly_lookup,expected_sha256=args.lookup_sha256)
    fetch,reference=None,None
    if args.fasta:
        from conference.model_v2.j6_site_universe import FaiFastaReader, reference_identity, file_identity
        genome=FaiFastaReader(j5.safe_path(args.fasta))
        reference=reference_identity(genome) | dict(fasta=file_identity(args.fasta))
        fetch=genome.fetch
    with j5.ContextIndex(args.index_dir/'context.sqlite',readonly=True,lookup=lookup) as index:
        complete=json.loads((args.index_dir/'index.done').read_text())
        if (complete.get('index_sha256')!=index.contract() or
                complete.get('provenance_sha256')!=j5.digests(args.index_dir/'provenance.parquet')[0]):
            raise ValueError('completed index contract mismatch')
        if args.genes:
            genes=index.db.execute('select sha256 from sources where name=?',('genes',)).fetchone()
            if not genes or genes[0]!=j5.digests(args.genes)[0]:
                raise ValueError('gene provenance does not match completed index')
        yield index,fetch,reference


def profile_window(args, source, destination):
    with profile_worker(args) as (index,fetch,reference):
        return measure(index,source,destination,source.name,fetch=fetch,
                       reference_contract=reference,require_query_done=True)


def profile_persistent(args, sources, root, *, total_windows, shard_windows):
    """One worker, contiguous production order; first cold window then warm caches.

    Always remeasure the whole sequence. Reusing window timing checkpoints would
    skip chromosome loads and invalidate the persistent-worker cache history.
    RSS is the cumulative process high-water mark, including all previous windows.
    """
    if not sources: raise ValueError('persistent profile requires windows')
    start=time.perf_counter(); windows=[]
    with profile_worker(args) as (index,fetch,reference):
        setup=time.perf_counter()-start
        for source in sources:
            destination=root/(source.stem+'.parquet')
            destination.with_suffix('.timing.json').unlink(missing_ok=True)
            measured=measure(index,source,destination,source.name,fetch=fetch,
                reference_contract=reference,require_query_done=True)
            measured.update(worker_pid=os.getpid(),
                rss_scope='persistent_worker_cumulative_peak',sequence=len(windows),
                cache_state='cold' if not windows else 'persistent')
            windows.append(measured)
            print(json.dumps(measured),flush=True)
        peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
        largest=math.ceil(total_windows/args.shard_count)
        # Worst measured window, including cold load; conservative for this sample.
        per_window=max(w['seconds'] for w in windows)
        shard_seconds=setup+largest*per_window
        projection=dict(array_tasks=args.shard_count,total_windows=total_windows,
            largest_shard_windows=largest,measured_shard_windows=shard_windows,
            seconds_per_window=per_window,setup_seconds=setup,shard_seconds=shard_seconds,
            array_worker_seconds=args.shard_count*shard_seconds,
            array_wall_seconds_at_concurrency_4=math.ceil(args.shard_count/4)*shard_seconds,
            output_bytes=total_windows*max(w['disk_bytes'] for w in windows),
            peak_rss_mib=peak,recommended_seconds=math.ceil(1.5*shard_seconds),
            recommended_mem_mib=math.ceil(1.5*peak),safety_factor=1.5,
            formula='shard = setup + ceil(total_windows/shards)*max(window_seconds); array_wall@4 = ceil(shards/4)*shard; disk = total_windows*max(window_disk_bytes); requests >= 1.5*shard/peak_RSS')
        result=dict(mode='persistent',shard_id=args.shard_id,windows=windows,
            reference=reference,genes_supplied=bool(args.genes),
            genes_sha256=j5.digests(args.genes)[0] if args.genes else None,
            gene_bounds_note='supplied and index-verified' if args.genes else 'NOT supplied; embedded J6/index bounds may still enable controls',
            index_sha256=index.contract(),lookup_sha256=args.lookup_sha256,
            benchmark_sha256=j5.digests(__file__)[0],projection=projection,
            limitations='Sample projection; chromosome density, later cache peaks and shared filesystem contention can exceed it; disk excludes temporary peak/quota competition.')
    j5.atomic_json(root/'persistent_profile.json',result)
    return result


def allocation_seconds(value):
    """scontrol TimeLimit uses [days-]HH:MM:SS; refuse unknown/unlimited units."""
    match=re.fullmatch(r'(?:(\d+)-)?(\d+):(\d{2}):(\d{2})',value)
    if not match: raise ValueError('invalid Slurm TimeLimit')
    days,hours,minutes,seconds=(int(v or 0) for v in match.groups())
    if minutes>=60 or seconds>=60: raise ValueError('invalid Slurm TimeLimit')
    return days*86400+hours*3600+minutes*60+seconds


def check_allocation(args, path, *, seconds, memory_mib):
    """Refuse full-array resources below the same-config persistent projection."""
    profile=json.loads(j5.safe_path(path).read_text())
    with profile_worker(args) as (index,_,reference):
        expected=dict(mode='persistent',index_sha256=index.contract(),reference=reference,
            genes_sha256=j5.digests(args.genes)[0] if args.genes else None,
            lookup_sha256=args.lookup_sha256,manifest_sha256=args.j6_manifest_sha256,
            benchmark_sha256=j5.digests(__file__)[0])
        if any(profile.get(k)!=v for k,v in expected.items()):
            raise ValueError('persistent profile configuration mismatch')
    if len(profile.get('windows',[]))<40:
        raise ValueError('full run requires at least 40 consecutive profile windows')
    projection=profile['projection']
    if projection['array_tasks']!=args.shard_count:
        raise ValueError('persistent profile shard configuration mismatch')
    if seconds<projection['recommended_seconds'] or memory_mib<projection['recommended_mem_mib']:
        raise ValueError('allocation below persistent projection 1.5x; raise SBATCH_TIMELIMIT/SBATCH_MEM_PER_NODE')
    return True


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sites',type=int,default=128)
    parser.add_argument('--sentinel',type=Path,default=Path.home()/'orcd/scratch/feng/model_v2/j5_sentinel/context.sqlite')
    parser.add_argument('--index-dir',type=Path)
    parser.add_argument('--out',type=Path,help='Persistent full-index benchmark checkpoints')
    parser.add_argument('--assembly-lookup',type=Path,default=j5.LOOKUP_PATH)
    parser.add_argument('--lookup-sha256',default=j5.LOOKUP_SHA256)
    parser.add_argument('--j6-genome-dir',type=Path,default=Path.home()/'orcd/scratch/feng/model_v2/j6_genome_v1')
    parser.add_argument('--fasta',type=Path,help='Same J5_FASTA as production; omit only for a no-GC production run')
    parser.add_argument('--genes',type=Path)
    parser.add_argument('--high-windows',type=int,default=12)
    parser.add_argument('--median-windows',type=int,default=4)
    parser.add_argument('--j6-manifest-sha256',default='b7b64df4609781af7dbb7aebe9f32091f7e22fca5480bd2371506585a246cdd5')
    parser.add_argument('--check-allocation',type=Path,help='Validate full-job resource request against persistent profile')
    parser.add_argument('--allocation-time')
    parser.add_argument('--allocation-mib',type=int)
    parser.add_argument('--persistent',action='store_true',help='One worker on consecutive windows of one production shard')
    parser.add_argument('--windows',type=int,default=40)
    parser.add_argument('--shard-id',type=int,default=0)
    parser.add_argument('--shard-count',type=int,default=64)
    parser.add_argument('--offset',type=int,default=0,help='First window within the selected shard')
    parser.add_argument('--window',type=Path,help=argparse.SUPPRESS)
    args=parser.parse_args(argv)
    if args.check_allocation:
        if not args.index_dir or not args.allocation_time or not args.allocation_mib:
            parser.error('--check-allocation requires --index-dir, --allocation-time and --allocation-mib')
        check_allocation(args,args.check_allocation,seconds=allocation_seconds(args.allocation_time),memory_mib=args.allocation_mib)
        return
    if args.persistent and not args.index_dir: parser.error('--persistent requires --index-dir')
    if not 1<=args.sites<=1024: parser.error('bounded local sample: 1..1024 sites')
    if args.out:
        if not args.index_dir: parser.error('--out is only for the full-index Slurm profile')
        args.out.mkdir(parents=True,exist_ok=True)
    with (nullcontext(args.out) if args.out else tempfile.TemporaryDirectory(prefix='j5-bench-')) as temporary:
        root=Path(temporary)
        if args.index_dir:
            if args.persistent:
                if args.window: parser.error('--persistent and --window are incompatible')
                if args.windows<2 or args.offset<0: parser.error('need >=2 windows and nonnegative offset')
                all_sources=j5.genome_queries(args.j6_genome_dir,0,1,manifest_sha256=args.j6_manifest_sha256)
                shard=j5.genome_queries(args.j6_genome_dir,args.shard_id,args.shard_count,manifest_sha256=args.j6_manifest_sha256)
                selected=shard[args.offset:args.offset+args.windows]
                if len(selected)!=args.windows: raise ValueError('not enough consecutive shard windows')
                result=profile_persistent(args,selected,root,total_windows=len(all_sources),shard_windows=len(shard))
                result.update(manifest_sha256=args.j6_manifest_sha256,offset=args.offset)
                j5.atomic_json(root/'persistent_profile.json',result)
                print(json.dumps({k:v for k,v in result.items() if k!='windows'}),flush=True)
                return
            if args.window:
                # One fresh interpreter per window: RSS is this worker's peak,
                # with no inherited chromosome cache from earlier windows.
                print(json.dumps(profile_window(args,args.window,root/'features.parquet')),flush=True)
                return
            sources=j5.genome_queries(args.j6_genome_dir,0,1,manifest_sha256=args.j6_manifest_sha256)
            selected=select_windows(sources,high=args.high_windows,median=args.median_windows)
            j5.atomic_json(root/'selection.json',dict(manifest_sha256=args.j6_manifest_sha256,
                fasta=str(args.fasta) if args.fasta else None,
                windows=[dict(source=str(p),rows=n,group=g) for p,n,g in selected]))
            for source,n,group in selected:
                destination=root/source.stem
                command=[sys.executable,'-m','conference.model_v2.jobs.j5_benchmark',
                    '--index-dir',str(args.index_dir),'--out',str(destination),
                    '--assembly-lookup',str(args.assembly_lookup),'--lookup-sha256',args.lookup_sha256,
                    '--window',str(source)]
                if args.fasta: command+=['--fasta',str(args.fasta)]
                if args.genes: command+=['--genes',str(args.genes)]
                start=time.perf_counter()
                completed=subprocess.run(command,check=True,capture_output=True,text=True)
                result=json.loads(completed.stdout)
                if result['sites']!=n: raise ValueError('J6 selection count changed')
                result.update(selection_group=group,manifest_rows=n,
                    worker_wall_seconds=time.perf_counter()-start,rss_scope='fresh_window_worker_peak')
                j5.atomic_json(destination/'profile.json',result)
                print(json.dumps(result),flush=True)
        else:
            synthetic_dir=root/'synthetic'; synthetic_dir.mkdir()
            print(json.dumps(synthetic(synthetic_dir,args.sites)),flush=True)
            real_dir=root/'real'; real_dir.mkdir()
            print(json.dumps(legacy_sentinel(real_dir,args.sentinel,args.sites)),flush=True)


if __name__=='__main__': main()
