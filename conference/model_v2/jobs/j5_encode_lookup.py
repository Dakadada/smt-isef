"""Refresh the *pinned inventory* with public ENCODE file + experiment biology.

No local inventory substitution: even missing source accessions remain expected.
A failed request retains its row with lookup_failed, never inferred biology.
Writes a new snapshot and SHA256 sidecar; pass both explicitly to J5 ingestion.
"""
import argparse
import csv
from datetime import datetime, timezone
import platform
import hashlib
import json
from pathlib import Path
import re
import sys
import urllib.request

DEFAULT = Path(__file__).resolve().parents[1]/'ops/ENCODE_ASSEMBLY_LOOKUP.tsv'
COLUMNS = ['kind','rbp','cell','file','accession','assembly','file_format','format_type',
           'output_type','status','dataset','md5','file_size','target_label',
           'biosample_term_name','identity_status','lookup','error']


def fetch_json(path):
    if not re.fullmatch(r'/(files|experiments|targets|biosample-types)/[A-Za-z0-9_.-]+/',path):
        raise ValueError('unexpected ENCODE object path')
    request=urllib.request.Request('https://www.encodeproject.org'+path+'?format=json&frame=object',
        headers={'Accept':'application/json','User-Agent':'model-v2-provenance-check'})
    with urllib.request.urlopen(request,timeout=20) as response:
        return json.load(response)


def refresh(records, *, fetch=fetch_json):
    cache={}
    def validate(obj, path, collection):
        expected_type={'files':'File','experiments':'Experiment','targets':'Target',
                       'biosample-types':'BiosampleType'}[collection]
        if not isinstance(path,str) or not re.fullmatch('/'+collection+r'/[A-Za-z0-9_.-]+/',path):
            raise ValueError('unexpected '+collection+' object path')
        if not isinstance(obj,dict) or obj.get('@id')!=path:
            raise ValueError('ENCODE object @id mismatch')
        types=obj.get('@type')
        if (not isinstance(types,list) or any(not isinstance(t,str) for t in types)
                or expected_type not in types
                or (set(types) & {'File','Experiment','Target','BiosampleType'})!={expected_type}):
            raise ValueError('ENCODE object @type mismatch')
        return obj
    def cached(path, collection):
        # Validate before fetching, including when fetch is an injected test client.
        if not isinstance(path,str) or not re.fullmatch('/'+collection+r'/[A-Za-z0-9_.-]+/',path):
            raise ValueError('unexpected '+collection+' object path')
        if path not in cache: cache[path]=fetch(path)
        return validate(cache[path],path,collection)
    def biology(value, collection, field):
        path=value.get('@id') if isinstance(value,dict) else value
        canonical=cached(path,collection)
        term=canonical.get(field)
        if not isinstance(term,str) or not term.strip():
            raise ValueError('ENCODE '+field+' unavailable')
        if isinstance(value,dict):
            validate(value,path,collection)
            if value.get(field)!=term: raise ValueError('embedded '+field+' contradiction')
        return term
    output=[]
    for source in records:
        rec={k:source[k] for k in ('kind','rbp','cell','file','accession')}
        try:
            obj=cached('/files/'+rec['accession']+'/', 'files')
            if obj.get('accession')!=rec['accession']: raise ValueError('file accession mismatch')
            dataset=obj.get('dataset')
            if not isinstance(dataset,str) or not re.fullmatch(r'/experiments/ENCSR[A-Z0-9]{6}/',dataset):
                raise ValueError('experiment dataset unavailable')
            experiment=cached(dataset,'experiments')
            # The fetched experiment's canonical id must be the file's dataset.
            if experiment['@id']!=obj['dataset']: raise ValueError('file dataset contradiction')
            label=biology(experiment.get('target'),'targets','label')
            biosample=biology(experiment.get('biosample_ontology'),'biosample-types','term_name')
            supplied=experiment.get('biosample_term_name')
            if supplied is not None and supplied!=biosample:
                raise ValueError('experiment biosample term contradiction')
            rec.update(assembly=obj.get('assembly'),file_format=obj.get('file_format'),
                format_type=str(obj.get('file_format_type')),output_type=obj.get('output_type'),
                status=obj.get('status'),dataset=dataset,md5=obj.get('md5sum'),file_size=obj.get('file_size'),
                target_label=label,biosample_term_name=biosample,lookup='ok',
                identity_status='verified' if (label,biosample)==(rec['rbp'],rec['cell']) else 'contradiction')
        except (OSError,ValueError,KeyError,TypeError) as exc:
            rec.update(lookup='lookup_failed',identity_status='unavailable',error=str(exc)[:240])
        output.append(rec)
    return output


def write_provenance(path, rows):
    """Bind a snapshot to this generator; safe for injected, offline test clients."""
    from conference.model_v2 import j5_rbp_context as j5
    record=dict(schema_version=1,script_sha256=j5.digests(__file__)[0],
        utc_time=datetime.now(timezone.utc).isoformat(),python_version=platform.python_version(),
        rows=len(rows),failures=sum(r.get('lookup')!='ok' for r in rows),
        contradictions=sum(r.get('identity_status')=='contradiction' for r in rows),
        tsv_sha256=j5.digests(path)[0])
    j5.atomic_json(Path(path).with_suffix('.tsv.json'),record)
    return record


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory',type=Path,default=DEFAULT)
    parser.add_argument('--out',type=Path,default=Path(__file__).with_name('j5_encode_lookup.tsv'))
    args=parser.parse_args(argv)
    with args.inventory.open() as f: records=list(csv.DictReader(f,delimiter='\t'))
    rows=refresh(records)
    temp=args.out.with_suffix('.tsv.tmp')
    with temp.open('w') as f:
        writer=csv.DictWriter(f,COLUMNS,delimiter='\t',lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    temp.replace(args.out)
    provenance=write_provenance(args.out,rows)
    sha=provenance['tsv_sha256']
    args.out.with_suffix('.tsv.sha256').write_text(sha+'\n')
    failures=sum(r['lookup']!='ok' for r in rows)
    contradictions=[r['accession'] for r in rows if r['identity_status']=='contradiction']
    print(json.dumps(dict(path=str(args.out),sha256=sha,rows=len(rows),failures=failures,contradictions=contradictions)))
    return 1 if failures or contradictions else 0


if __name__=='__main__': sys.exit(main())
