"""Bounded synthetic scaling/RSS probe. No production references or Slurm jobs.

Run with python -m conference.model_v2.jobs.j6_profile and --genes 1000 or
--genes 8000 in separate processes for independent peak RSS.
"""
import argparse
import cProfile
import json
import pstats
import resource
import tempfile
import time

from conference.model_v2 import j6_site_universe as j6
from conference.model_v2.trans_intervention import ASO, InMemoryGenome


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", type=int, choices=(1000, 8000), default=1000)
    args = parser.parse_args()
    genes = [j6.Gene(f"g{i}", "chrT", 20*i + 1, 20*i + 12, "+") for i in range(args.genes)]
    genome = InMemoryGenome({"chrT": "A" * (20 * args.genes)})
    started = time.perf_counter()
    worker = j6.SiteWorker(genome, genes, [ASO("q", "TTTT")], synthetic_reference=True)
    setup_seconds = time.perf_counter() - started
    profile = cProfile.Profile()
    profile.enable()
    started = time.perf_counter()
    rows = 0
    for i in range(20):
        start = 20 * (args.genes // 2 + i) + 1
        rows += sum(1 for _ in worker.iter_sites(shard=j6.Shard("chrT", start, start), max_mismatch=0))
    query_seconds = time.perf_counter() - started
    profile.disable()
    stats = pstats.Stats(profile).stats
    interval_calls = sum(value[1] for key, value in stats.items() if key[2] == "_overlapping")
    assert rows == 20
    # Dense output probes the scanner-hit set and bounded Arrow writer together.
    dense = InMemoryGenome({"chrD": "A" * 20000})
    dense_worker = j6.SiteWorker(dense, [j6.Gene("dense", "chrD", 1, 20000, "+")],
                                 [ASO("dense", "TTTT")], synthetic_reference=True)
    with tempfile.TemporaryDirectory(prefix="j6-profile-") as output:
        path = j6.run_shard(worker=dense_worker, shard=j6.Shard("chrD", 1, 20000), output_dir=output,
                            identities={"fasta": {"synthetic": True}, "gtf": {"synthetic": True}},
                            max_mismatch=3)
        dense_rows = json.loads(path.with_suffix(".done").read_text())["rows"]
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    assert dense_rows == 19997
    assert peak < 2 * 1024 * 1024, f"peak RSS exceeds 2 GiB: {peak} KiB"
    print(json.dumps(dict(genes=args.genes, query_windows=20, query_rows=rows, dense_rows=dense_rows,
                          setup_seconds=round(setup_seconds, 4), query_seconds=round(query_seconds, 4),
                          interval_calls=interval_calls, peak_rss_kib=peak), sort_keys=True))


if __name__ == "__main__":
    main()
