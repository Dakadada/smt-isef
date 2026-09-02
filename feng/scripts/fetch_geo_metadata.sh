#!/usr/bin/env bash

# Fetch and validate the GSE307431 sample metadata needed by Phase 3.
# Run from the repository root on an internet-connected login node:
#
#   bash feng/scripts/fetch_geo_metadata.sh
#
# The output is published atomically only after every manifest row has a unique
# GEO sample, an exact recognized sample title, and a matching SRA RunInfo row.

set -euo pipefail

SUBMIT_DIR="${SLURM_SUBMIT_DIR:-$PWD}"
MANIFEST="${FENG_P3_MANIFEST:-${SUBMIT_DIR}/feng/artifacts/phase3/ena_manifest.tsv}"
OUTPUT="${FENG_P3_CONDITIONS:-${SUBMIT_DIR}/feng/artifacts/phase3/sample_conditions.tsv}"
GEO_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE307nnn/GSE307431/soft/GSE307431_family.soft.gz"
RUNINFO_URL="https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA1321660"

[[ -s "${MANIFEST}" ]] || { echo "ERROR: manifest missing: ${MANIFEST}" >&2; exit 1; }
command -v curl >/dev/null || { echo "ERROR: curl is required" >&2; exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 is required" >&2; exit 1; }

mkdir -p "$(dirname "${OUTPUT}")"
FETCH_DIR="$(mktemp -d "${TMPDIR:-/tmp}/feng-p3-metadata.XXXXXX")"
cleanup() { rm -rf -- "${FETCH_DIR}"; }
trap cleanup EXIT INT TERM

echo "Fetching GEO family SOFT: ${GEO_URL}"
curl -fL --retry 5 --connect-timeout 20 \
    "${GEO_URL}" -o "${FETCH_DIR}/GSE307431_family.soft.gz"
echo "Fetching SRA RunInfo: ${RUNINFO_URL}"
curl -fL --retry 5 --connect-timeout 20 \
    "${RUNINFO_URL}" -o "${FETCH_DIR}/runinfo.csv"

python3 - "${MANIFEST}" "${FETCH_DIR}/GSE307431_family.soft.gz" \
    "${FETCH_DIR}/runinfo.csv" "${OUTPUT}" "${GEO_URL}" "${RUNINFO_URL}" <<'PY'
import csv
import gzip
import os
import re
import sys
import tempfile
from collections import Counter

manifest_path, soft_path, runinfo_path, output_path, geo_url, runinfo_url = sys.argv[1:]

with open(manifest_path, newline="") as handle:
    manifest = list(csv.DictReader(handle, delimiter="\t"))
if len(manifest) != 20:
    raise SystemExit(f"ERROR: expected 20 manifest rows, found {len(manifest)}")
required_manifest = {"run_accession", "sample_accession"}
if not manifest or not required_manifest.issubset(manifest[0]):
    raise SystemExit("ERROR: manifest lacks run_accession/sample_accession")

# Parse sample blocks from the downloaded family SOFT. The BioSample relation
# is the join key to the repository manifest; titles supply condition metadata.
samples = {}
current = None
with gzip.open(soft_path, "rt", encoding="utf-8", errors="strict") as handle:
    for raw in handle:
        line = raw.rstrip("\n")
        if line.startswith("^SAMPLE = "):
            gsm = line.split("=", 1)[1].strip()
            current = {"gsm_accession": gsm, "title": "", "biosample": ""}
            samples[gsm] = current
        elif current is not None and line.startswith("!Sample_title = "):
            current["title"] = line.split("=", 1)[1].strip()
        elif current is not None and line.startswith("!Sample_relation = "):
            match = re.search(r"SAMN\d+", line)
            if match:
                current["biosample"] = match.group(0)

by_biosample = {}
for sample in samples.values():
    biosample = sample["biosample"]
    if not biosample:
        continue
    if biosample in by_biosample:
        raise SystemExit(f"ERROR: duplicate GEO BioSample relation: {biosample}")
    by_biosample[biosample] = sample

with open(runinfo_path, newline="") as handle:
    runinfo_rows = list(csv.DictReader(handle))
runinfo_by_run = {row.get("Run", ""): row for row in runinfo_rows if row.get("Run")}

title_patterns = (
    (re.compile(r"^Untreated sample, replicate (\d+)$"),
     ("control", "untreated", "untreated")),
    (re.compile(r"^Scrambled OMe-modified ASO, replicate (\d+)$"),
     ("control", "OMe", "scrambled")),
    (re.compile(r"^ISS-N1-targeting OMe-modified ASO, replicate (\d+)$"),
     ("F18OMe", "OMe", "targeted")),
    (re.compile(r"^Scrambled MOE-modified ASO, replicate (\d+)$"),
     ("control", "MOE", "scrambled")),
    (re.compile(r"^ISS-N1-targeting MOE-modified ASO, replicate (\d+)$"),
     ("F18MOE", "MOE", "targeted")),
    (re.compile(r"^Scrambled PMO ASO, replicate (\d+)$"),
     ("control", "PMO", "scrambled")),
    (re.compile(r"^ISS-N1-targeting PMO ASO, replicate (\d+)$"),
     ("F20PMO", "PMO", "targeted")),
)

def classify(title):
    matches = []
    for pattern, labels in title_patterns:
        match = pattern.fullmatch(title)
        if match:
            matches.append((*labels, int(match.group(1))))
    if len(matches) != 1:
        raise SystemExit(f"ERROR: unrecognized or ambiguous GEO title: {title!r}")
    return matches[0]

rows = []
seen_gsm = set()
for item in manifest:
    run = item["run_accession"]
    biosample = item["sample_accession"]
    sample = by_biosample.get(biosample)
    if sample is None:
        raise SystemExit(f"ERROR: no GEO sample maps to manifest BioSample {biosample}")
    if sample["gsm_accession"] in seen_gsm:
        raise SystemExit(f"ERROR: GEO sample reused: {sample['gsm_accession']}")
    seen_gsm.add(sample["gsm_accession"])

    info = runinfo_by_run.get(run)
    if info is None:
        raise SystemExit(f"ERROR: run absent from fetched RunInfo: {run}")
    if info.get("BioSample") != biosample:
        raise SystemExit(
            f"ERROR: RunInfo BioSample mismatch for {run}: "
            f"{info.get('BioSample')!r} != {biosample!r}"
        )

    condition, chemistry, control_type, replicate = classify(sample["title"])
    rows.append({
        "run_accession": run,
        "sample_accession": biosample,
        "gsm_accession": sample["gsm_accession"],
        "condition": condition,
        "chemistry": chemistry,
        "control_type": control_type,
        "replicate": replicate,
        "geo_title": sample["title"],
        "provenance": f"GSE307431:{geo_url};PRJNA1321660:{runinfo_url}",
    })

counts = Counter(row["condition"] for row in rows)
expected = Counter({"control": 11, "F18OMe": 3, "F18MOE": 3, "F20PMO": 3})
if counts != expected:
    raise SystemExit(f"ERROR: unexpected condition counts: {dict(counts)}; expected {dict(expected)}")

scrambled = Counter(
    row["chemistry"] for row in rows if row["control_type"] == "scrambled"
)
if scrambled != Counter({"OMe": 3, "MOE": 3, "PMO": 3}):
    raise SystemExit(f"ERROR: unexpected matched-control counts: {dict(scrambled)}")

fields = [
    "run_accession", "sample_accession", "gsm_accession", "condition",
    "chemistry", "control_type", "replicate", "geo_title", "provenance",
]
output_dir = os.path.dirname(os.path.abspath(output_path))
fd, temporary = tempfile.mkstemp(prefix=".sample_conditions.", suffix=".tmp", dir=output_dir)
try:
    with os.fdopen(fd, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(temporary, output_path)
except BaseException:
    try:
        os.unlink(temporary)
    except FileNotFoundError:
        pass
    raise

print(f"Wrote {len(rows)} validated rows to {output_path}")
print("Condition counts: " + ", ".join(f"{key}={counts[key]}" for key in sorted(counts)))
PY
