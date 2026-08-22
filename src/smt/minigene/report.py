"""A6 CLI: run the compression algorithm on SMN2 (blind validation, item 2),
STMN2 and UNC13A (item 3), and print the numbers the HANDOFF-style report
needs. Also runs the log2fc_threshold sweep explicitly (item 2's "if you do
explore thresholds, report the full sweep" requirement) using the already-
computed per-gap measurement trace -- no refolding for the sweep.

Run with:  python -m smt.minigene.report [--out DIR]
"""

from __future__ import annotations

import argparse
from pathlib import Path

from smt.minigene.compress import CompressionParams, DEFAULT_PARAMS, GenePlan, build_construct, compress_gene, rescore_gap
from smt.minigene.spec import blocks, junctions, write_spec
from smt.seq import targets as T
from smt.seq.reference import SMN2_EXONS, SMN2_LENGTH, Interval, exon, load_smn2

THRESHOLD_SWEEP = (0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0)


def _smn2_plan(params: CompressionParams = DEFAULT_PARAMS) -> GenePlan:
    seq = load_smn2()
    exon_order = list(SMN2_EXONS.keys())
    exons = {k: exon(k) for k in exon_order}
    return compress_gene("SMN2", seq, exon_order, exons, params)


def _print_intron_table(plan: GenePlan) -> None:
    print(f"{'intron':<6}{'nat_len':>9}{'retained':>10}{'frac':>8}{'ndel':>6}  gaps(method)")
    for name, ip in plan.introns.items():
        gap_methods = ",".join(g.measurements[0].method if g.measurements else "n/a" for g in ip.gaps)
        print(
            f"{name:<6}{len(ip.interval):>9}{ip.retained_length:>10}"
            f"{ip.retained_fraction:>8.3f}{len(ip.deletions):>6}  {gap_methods}"
        )


def run_smn2_validation(out: Path) -> GenePlan:
    print("=" * 78)
    print("A6 VALIDATION: SMN2, blind (no gene-specific tuning, default threshold "
          f"log2fc_threshold={DEFAULT_PARAMS.log2fc_threshold})")
    print("=" * 78)
    plan = _smn2_plan()
    _print_intron_table(plan)

    c = build_construct(plan)
    print(f"\ncompressed length {len(c)} nt")
    print(f"natural transcribed unit (promoter..exon8) {plan.transcribed_unit_length()} nt")
    print(f"natural full pinned record (NG_008728) {plan.natural_length()} nt")
    print(f"compression ratio vs transcribed unit: {len(c)/plan.transcribed_unit_length():.3f}")

    print("\nTHRESHOLD SWEEP (re-decided from the SAME measurement trace, no refolding):")
    header = "log2fc_thr  " + "  ".join(f"i{n:<4}" for n in plan.introns)
    print(header)
    for thr in THRESHOLD_SWEEP:
        fracs = []
        for name, ip in plan.introns.items():
            retained = 0
            total = len(ip.interval)
            deleted = 0
            for g in ip.gaps:
                if g.gap is None:
                    continue
                _, deletion = rescore_gap(g, thr)
                if deletion is not None:
                    deleted += len(deletion)
            fracs.append((total - deleted) / total)
        row = f"{thr:<10.2f}  " + "  ".join(f"{f:.3f}" for f in fracs)
        print(row)

    paths = write_spec(c, out, plan)
    print(f"\nwrote {paths}")
    return plan


def run_stmn2(out: Path) -> GenePlan:
    print("=" * 78)
    print("A6 APPLICATION: STMN2 (full gene, cryptic exon 2a retained as a pinned "
          "regulatory anchor -- see smt.seq.targets)")
    print("=" * 78)
    seq = T.load_stmn2()
    exon_order = list(T.STMN2_EXONS.keys())
    exons = {k: Interval(*v) for k, v in T.STMN2_EXONS.items()}
    keep = Interval(min(T.STMN2_CRYPTIC_3SS.start, T.STMN2_CRYPTIC_EXON.start), T.STMN2_CRYPTIC_EXON.end)
    extra = {"1": [(keep, T.STMN2_CRYPTIC_3SS)]}  # no pinned exit feature -> auto-derived (see compress.py)
    plan = compress_gene("STMN2", seq, exon_order, exons, DEFAULT_PARAMS, extra_keep_by_intron=extra)
    _print_intron_table(plan)
    c = build_construct(plan, extra_features={"STMN2_CRYPTIC_EXON": T.STMN2_CRYPTIC_EXON, "STMN2_TDP43_GU_SITE": T.STMN2_TDP43_GU_SITE})
    print(f"\ncompressed length {len(c)} nt vs natural gene body {T.STMN2_LENGTH} nt "
          f"(ratio {len(c)/T.STMN2_LENGTH:.4f})")
    paths = write_spec(c, out, plan)
    print(f"wrote {paths}")
    return plan


def run_unc13a(out: Path) -> GenePlan:
    print("=" * 78)
    print("A6 APPLICATION: UNC13A -- PARTIAL. Only exons 20/21 and the intervening")
    print("cryptic-exon region are pinned in smt.seq.targets (not the full 44-exon")
    print("gene), so this compresses that segment only, not the 87019 nt locus.")
    print("=" * 78)
    seq = T.load_unc13a()
    exon_order = ["20", "21"]
    exons = {k: Interval(*T.UNC13A_EXONS[k]) for k in exon_order}
    keep = Interval(T.UNC13A_CRYPTIC_EXON_178.start, T.UNC13A_CRYPTIC_EXON_178.end)
    # entry (near keep.start) = the 178bp isoform's own acceptor; exit (near
    # keep.end) = the shared donor -- both pinned, unlike STMN2's exit side.
    extra = {"20": [(keep, T.UNC13A_CRYPTIC_ACCEPTOR_178, T.UNC13A_CRYPTIC_DONOR)]}
    # keep_upstream_from=exon20.start: compress_gene's default (1) retains
    # everything upstream of the first exon in exon_order as "promoter",
    # which is correct for a whole-gene run (SMN2, STMN2) but WRONG here --
    # UNC13A's reference sequence is the whole 87019 nt gene body, and exon
    # 20 is nowhere near its start. Only exon 20 onward is pinned/modelled;
    # the ~45 kb genuinely upstream of exon 20 (introns/exons 1-19) is out of
    # scope, not a promoter, and must not be retained wholesale.
    plan = compress_gene(
        "UNC13A", seq, exon_order, exons, DEFAULT_PARAMS,
        extra_keep_by_intron=extra, keep_upstream_from=exons["20"].start,
    )
    _print_intron_table(plan)
    c = build_construct(
        plan,
        extra_features={
            "UNC13A_CRYPTIC_EXON_128": T.UNC13A_CRYPTIC_EXON_128,
            "UNC13A_CRYPTIC_EXON_178": T.UNC13A_CRYPTIC_EXON_178,
            "UNC13A_CRYPTIC_DONOR": T.UNC13A_CRYPTIC_DONOR,
            "UNC13A_CRYPTIC_ACCEPTOR_128": T.UNC13A_CRYPTIC_ACCEPTOR_128,
            "UNC13A_CRYPTIC_ACCEPTOR_178": T.UNC13A_CRYPTIC_ACCEPTOR_178,
        },
    )
    segment_natural_len = exons["21"].end - exons["20"].start + 1
    print(f"\ncompressed length {len(c)} nt vs natural pinned segment (exon20..exon21) "
          f"{segment_natural_len} nt (ratio {len(c)/segment_natural_len:.4f})")
    print(f"for context, full UNC13A gene body is {T.UNC13A_LENGTH} nt (out of scope here)")
    paths = write_spec(c, out, plan)
    print(f"wrote {paths}")
    return plan


GENE_RUNNERS = {"SMN2": run_smn2_validation, "STMN2": run_stmn2, "UNC13A": run_unc13a}


def run_one_gene_idempotent(gene: str, out: Path) -> None:
    """For slurm/minigene.sbatch: run one gene's compression, skip if its
    marker file already exists, and write it atomically (.tmp then rename)
    so a preempted/killed run never leaves a half-written marker that a
    resubmission would trust. The actual outputs (FASTA/CSV, from
    write_spec) are already written under ``out`` by the runner itself;
    the marker just records that this gene's run completed."""
    marker = out / f"{gene}.done"
    if marker.exists():
        print(f"{gene}: marker exists, skipping -> {marker}")
        return
    GENE_RUNNERS[gene](out)
    tmp = marker.with_suffix(".done.tmp")
    tmp.write_text("ok\n")
    tmp.rename(marker)
    print(f"{gene}: wrote marker {marker}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=Path("/home/dashii/orcd/pool/smt-isef-data/minigene"))
    p.add_argument("--gene", choices=[*GENE_RUNNERS, "all"], default="all")
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.gene == "all":
        for gene in GENE_RUNNERS:
            run_one_gene_idempotent(gene, args.out)
    else:
        run_one_gene_idempotent(args.gene, args.out)


if __name__ == "__main__":
    main()
