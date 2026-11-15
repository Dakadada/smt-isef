"""J3 completeness gate: all N shard .done markers exist, sha256 matches parquet, rows sum == J6 gate rows, channel failures reported (not hidden)."""
import json, sys, os, hashlib
d, n, expect_rows = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
bad, rows, cc = [], 0, {}
for i in range(n):
    b = f"{d}/j3_shard{i:04d}_of{n:04d}"
    try:
        m = json.load(open(b + ".done"))
        h = hashlib.sha256(open(b + ".parquet", "rb").read()).hexdigest()
        assert isinstance(m, dict) and h in json.dumps(m), "sha mismatch"
        rows += int(m.get("rows", 0))
        for ch, c in m["channel_counts"].items():
            for k, v in c.items(): cc.setdefault(ch, {}).setdefault(k, 0); cc[ch][k] += v
    except Exception as e:
        bad.append([i, repr(e)[:80]])
r = {"shards": n, "bad": len(bad), "rows": rows, "expected_rows": expect_rows, "rows_match": rows == expect_rows, "channel_counts": cc, "missing": bad[:50]}
json.dump(r, open(f"{d}/gate_report.json.tmp", "w")); os.replace(f"{d}/gate_report.json.tmp", f"{d}/gate_report.json")
print(r); sys.exit(0 if not bad and rows == expect_rows else 1)
