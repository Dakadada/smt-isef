import csv
import io
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

from .common import Checkpoints, read_tsv, validate_config, write_tsv

FIELDS = ["study_accession", "run_accession", "sample_accession", "sample_alias", "experiment_title",
          "fastq_ftp", "fastq_md5", "fastq_bytes", "library_layout", "library_strategy"]
MANIFEST_FIELDS = ["family", "project", "run_accession", "sample_accession", "sample_alias", "experiment_title",
                   "fastq_ftp", "fastq_md5", "fastq_bytes", "layout", "library_strategy", "endpoint_eligible",
                   "exclusion_reason", "retrieved_at", "metadata_source"]


def _ena(project, timeout=120):
    query = urllib.parse.urlencode({"accession": project, "result": "read_run", "fields": ",".join(FIELDS),
                                    "format": "tsv", "download": "true"})
    url = "https://www.ebi.ac.uk/ena/portal/api/filereport?" + query
    with urllib.request.urlopen(url, timeout=timeout) as response:  # nosec: fixed ENA HTTPS endpoint
        text = response.read().decode("utf-8")
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def run(config_path, metadata_dir=None):
    cfg = validate_config(config_path)
    out = Path(cfg["paths"]["work_dir"]); cp = Checkpoints(cfg["paths"]["checkpoints"], "resolve", config_path, cfg)
    all_rows = []
    for family, study in cfg["studies"].items():
        target = out / "manifests" / f"{family}.frozen.tsv"
        source = Path(metadata_dir) / f"{study['project']}.tsv" if metadata_dir else None
        inputs = [source] if source and source.exists() else []
        if cp.valid(family, [target], inputs):
            all_rows.extend(read_tsv(target)); continue
        rows = read_tsv(source) if source and source.exists() else _ena(study["project"])
        now = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        cooked = []
        for row in rows:
            strategy = row.get("library_strategy", "")
            searchable = " ".join(row.get(k, "") for k in ("experiment_title", "sample_alias", "library_strategy"))
            eligible = strategy.casefold() == "rna-seq" and not any(re.search(p, searchable, re.I) for p in study.get("exclude_regex", []))
            reason = "" if eligible else ("non-RNA-seq library_strategy" if strategy.casefold() != "rna-seq" else "excluded assay/sample rule")
            cooked.append({"family": family, "project": study["project"], "run_accession": row.get("run_accession", ""),
                "sample_accession": row.get("sample_accession", ""), "sample_alias": row.get("sample_alias", ""),
                "experiment_title": row.get("experiment_title", ""), "fastq_ftp": row.get("fastq_ftp", ""),
                "fastq_md5": row.get("fastq_md5", ""), "fastq_bytes": row.get("fastq_bytes", ""),
                "layout": row.get("library_layout", ""), "library_strategy": strategy,
                "endpoint_eligible": str(eligible).lower(), "exclusion_reason": reason, "retrieved_at": now,
                "metadata_source": str(source) if source and source.exists() else "ENA filereport result=read_run"})
        eligible_rows=[r for r in cooked if r["endpoint_eligible"]=="true"]
        sample_map=study.get("sample_map",{})
        missing=[r["run_accession"] for r in eligible_rows if r["run_accession"] not in sample_map]
        extra=sorted(set(sample_map)-{r["run_accession"] for r in eligible_rows})
        if missing or extra: raise ValueError(f"{family}: frozen accession sample_map mismatch; missing={missing}, extra={extra}")
        write_tsv(target,cooked,MANIFEST_FIELDS)
        cp.complete(family,[target],inputs=inputs,resolved_metadata=sample_map,project=study["project"],eligible_runs=len(eligible_rows))
        all_rows.extend(cooked)
    combined = out / "manifests" / "all_runs.frozen.tsv"
    write_tsv(combined, all_rows, MANIFEST_FIELDS)
    return combined
