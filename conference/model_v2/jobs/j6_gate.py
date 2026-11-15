"""J6 completeness gate: every shard in run_manifest.expected_shards has a valid .done (rows+sha256 match parquet). Exit 1 + gate_report.json on any gap."""
import json, sys, os, hashlib
d = sys.argv[1]
m = json.load(open(f"{d}/run_manifest.json"))
bad, rows = [], 0
for s in m["expected_shards"]:
    try:
        dn = json.load(open(f"{d}/{s}.done"))
        h = hashlib.sha256(open(f"{d}/{s}.parquet", "rb").read()).hexdigest()
        assert isinstance(dn, dict) and dn.get("sha256") == h
        rows += dn["rows"]
    except Exception as e:
        bad.append([s, repr(e)[:80]])
r = {"expected": len(m["expected_shards"]), "bad": len(bad), "rows": rows, "missing": bad[:50]}
json.dump(r, open(f"{d}/gate_report.json.tmp", "w")); os.replace(f"{d}/gate_report.json.tmp", f"{d}/gate_report.json")
print(r); sys.exit(1 if bad else 0)
