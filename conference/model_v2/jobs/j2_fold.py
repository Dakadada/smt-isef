"""J2: local RNA accessibility (RNAplfold unpaired prob) for loci. proxy-class output.
Usage: j2_fold.py LOCI_TSV FASTA OUT_DIR SHARD NSHARDS
LOCI_TSV cols: locus_id chrom start1 end1 strand   (1-based inclusive)
Read-only tiny fasta slices via .fai. Idempotent: skips shard if .done exists."""
import sys, os, json, hashlib, itertools
import RNA, pyarrow as pa, pyarrow.parquet as pq
COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")
def fai(path):
    d = {}
    for l in open(path + ".fai"):
        n, ln, off, bl, bw = l.split()[:5]; d[n] = (int(ln), int(off), int(bl), int(bw))
    return d
def fetch(fa, idx, chrom, s, e):
    ln, off, bl, bw = idx[chrom]
    assert 1 <= s <= e <= ln
    a, b = s - 1, e
    start = off + a // bl * bw + a % bl
    end = off + (b - 1) // bl * bw + (b - 1) % bl + 1
    with open(fa, "rb") as f:
        f.seek(start); raw = f.read(end - start)
    return raw.replace(b"\n", b"").decode().upper()
def main(loci, fa, out, shard, n, windows=((80, 40), (120, 60), (200, 100)), ulens=(8, 12, 18)):
    os.makedirs(out, exist_ok=True)
    done = f"{out}/j2_shard{shard:04d}.done"
    if os.path.exists(done): print("skip", done); return
    idx = fai(fa); rows = []
    recs = [l.split("\t") for l in open(loci).read().splitlines() if l and not l.startswith("#")]
    for i, (lid, chrom, s, e, st) in enumerate(recs):
        if i % n != shard: continue
        seq = fetch(fa, idx, chrom, int(s), int(e))
        if st == "-": seq = seq.translate(COMP)[::-1]
        rna = seq.replace("T", "U")
        for (W, L), u in itertools.product(windows, ulens):
            if W > len(rna): continue
            # RNAplfold: unpaired prob of u-mers ending at each position, 5'->3' transcript orientation
            p = RNA.pfl_fold_up(rna, u, W, L)
            for pos in range(1, len(rna) + 1):
                v = p[pos][u] if pos < len(p) and u < len(p[pos]) else float("nan")
                rows.append((lid, st, pos, W, L, u, v))
    cols = list(zip(*rows)) if rows else [[]] * 7
    t = pa.table(dict(zip(["locus_id","strand","pos_end_1b_transcript","W","L","u","p_unpaired"], cols)))
    t = t.replace_schema_metadata({"evidence_class": "proxy", "tool": "ViennaRNA " + RNA.__version__, "method": "RNAplfold pfl_fold_up", "chemistry": "unmodified RNA (model_assumption)"})
    tmp = f"{out}/j2_shard{shard:04d}.parquet.tmp"; fin = tmp[:-4]
    pq.write_table(t, tmp); os.replace(tmp, fin)
    h = hashlib.sha256(open(fin, "rb").read()).hexdigest()
    json.dump({"rows": t.num_rows, "sha256": h, "loci": loci, "fasta": fa}, open(done + ".tmp", "w")); os.replace(done + ".tmp", done)
    print("done", fin, t.num_rows, h)
if __name__ == "__main__":
    a = sys.argv; main(a[1], a[2], a[3], int(a[4]), int(a[5]))
