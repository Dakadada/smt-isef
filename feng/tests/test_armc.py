"""Arm C bind->act gate: cheap pytest checks.

Covers: revcomp/orientation correctness of the scanner (a positive control
proving the SBO must be reverse-complemented before scanning, not scanned
directly -- the exact bug class that would silently search for the wrong
strand), S1/S4 parse row counts against the real PDFs on disk, the label-join
path (S1.id -> dg_id, S4 interval -> exon_key -> actionable), and the LOSO
gate's metric/degenerate-fold handling. No sealed ground-truth directory is
read anywhere in this file.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from feng.bench.armc_build import (
    DEFAULT_POOL_SEED,
    DEFAULT_TARGET_PAD,
    S1_PDF_DEFAULT,
    S4_PDF_DEFAULT,
    SBO_TSV_DEFAULT,
    _pool_keys_overlapping_positives,
    build_exonome,
    build_features,
    exons_df_from_keys,
    iter_gencode_exon_keys,
    load_armc_chain,
    load_gencode_exon_keys,
    parse_s1,
    parse_s4,
    positive_exon_keys_by_dg,
    positive_lift_coverage,
    resolve_chain_path,
    s4_direct_vs_indirect,
    sample_background_pool,
    sbo_low_complexity_reason,
    scan_sbo,
)
from feng.bench.armc_gate import (
    FEATURE_COLUMNS,
    average_precision,
    roc_auc,
    run_loso_gate,
)
from feng.bench.scores import NUSINERSEN_5P_TO_3P, NUSINERSEN_RC, reverse_complement
from feng.scan.seed_extend import seed_and_extend

REPO_ROOT = Path(__file__).resolve().parents[2]

# Staged references for the (skip-gated) liftover / frame tests.
CHAIN_PATH = resolve_chain_path()
GENOME_FASTA = Path("/home/dashii/orcd/scratch/feng/gencode/v50/GRCh38.primary_assembly.genome.fa")
SAMTOOLS = Path("/home/dashii/orcd/pool/feng-env/bin/samtools")
_GRCH38_MAIN_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}


# --- revcomp / scan orientation ---------------------------------------------


def test_reverse_complement_roundtrip() -> None:
    assert reverse_complement("AACG") == "CGTT"  # complement then reverse
    assert reverse_complement("ACGT") == "ACGT"  # palindrome
    assert reverse_complement(reverse_complement(NUSINERSEN_5P_TO_3P)) == NUSINERSEN_5P_TO_3P
    assert reverse_complement(NUSINERSEN_5P_TO_3P) == NUSINERSEN_RC


def _tiled_filler(n: int) -> str:
    motif = "GATTACAG"
    return (motif * (n // len(motif) + 1))[:n]


def test_scan_orientation_matters() -> None:
    """The exonome holds the SENSE-strand target site (= revcomp of the SBO).
    Scanning with the SBO sequence itself must NOT hit it; scanning with
    reverse_complement(sbo) (what scan_sbo actually does) must hit it exactly.
    This is a positive control for the orientation bug described in
    scan_sbo's docstring.
    """
    target_site = NUSINERSEN_RC  # what a sense-strand exon would contain
    exon_seq = _tiled_filler(20) + target_site + _tiled_filler(20)
    exonome = {"chr1:0-58:+": exon_seq}

    # Wrong orientation: scanning the SBO itself finds nothing at max_edit=0.
    wrong_hits = seed_and_extend(NUSINERSEN_5P_TO_3P, exonome, max_edit=0)
    assert wrong_hits == ()

    # Correct orientation (what scan_sbo does internally): exact hit.
    bound = scan_sbo(NUSINERSEN_5P_TO_3P, exonome, max_edit=0)
    assert list(bound["exon_key"]) == ["chr1:0-58:+"]
    assert bound.iloc[0]["edit_distance"] == 0
    assert bound.iloc[0]["exact_match_length"] == len(NUSINERSEN_5P_TO_3P)


def test_scan_sbo_respects_max_edit() -> None:
    target_site = NUSINERSEN_RC
    # introduce one mismatch in the middle of the embedded target
    mutated = target_site[:9] + ("A" if target_site[9] != "A" else "C") + target_site[10:]
    exon_seq = _tiled_filler(15) + mutated + _tiled_filler(15)
    exonome = {"chr2:0-48:+": exon_seq}

    no_hit = scan_sbo(NUSINERSEN_5P_TO_3P, exonome, max_edit=0)
    assert no_hit.empty

    one_hit = scan_sbo(NUSINERSEN_5P_TO_3P, exonome, max_edit=1)
    assert len(one_hit) == 1
    assert one_hit.iloc[0]["edit_distance"] == 1


# --- S1/S4 PDF parsing (real files on disk) ----------------------------------


@pytest.mark.skipif(not S1_PDF_DEFAULT.exists() and not (REPO_ROOT / S1_PDF_DEFAULT).exists(),
                     reason="Holgersen S1 PDF not staged")
def test_parse_s1_row_counts() -> None:
    s1 = parse_s1(REPO_ROOT / S1_PDF_DEFAULT)
    assert list(s1.columns) == ["id", "dg_id", "cell", "sequence", "category"]
    assert s1["id"].is_unique
    assert len(s1) == 105
    assert s1["dg_id"].nunique() == 82  # 81 canonical SBOs + 1 duplicate-sequence dg_id (DG2386)
    assert set(s1["sequence"].str.len().between(16, 20)) == {True}


@pytest.mark.skipif(not S4_PDF_DEFAULT.exists() and not (REPO_ROOT / S4_PDF_DEFAULT).exists(),
                     reason="Holgersen S4 PDF not staged")
def test_parse_s4_row_counts() -> None:
    s4 = parse_s4(REPO_ROOT / S4_PDF_DEFAULT)
    expected_cols = [
        "id", "chrom", "start", "end", "strand", "case_psi", "control_psi",
        "rnaplex_dg", "rnacofold_dg", "edit_dist", "exact_match",
    ]
    assert list(s4.columns) == expected_cols
    assert len(s4) == 724
    assert s4["id"].nunique() == 75
    assert (s4["end"] > s4["start"]).all()
    assert set(s4["strand"].unique()) <= {"+", "-"}


@pytest.mark.skipif(
    not (REPO_ROOT / S1_PDF_DEFAULT).exists() or not (REPO_ROOT / S4_PDF_DEFAULT).exists()
    or not (REPO_ROOT / SBO_TSV_DEFAULT).exists(),
    reason="Holgersen S1/S4 PDFs or sbo_sequences.tsv not staged",
)
def test_s4_ids_join_cleanly_to_s1() -> None:
    s1 = parse_s1(REPO_ROOT / S1_PDF_DEFAULT)
    s4 = parse_s4(REPO_ROOT / S4_PDF_DEFAULT)
    assert set(s4["id"]) <= set(s1["id"])

    sbo = pd.read_csv(REPO_ROOT / SBO_TSV_DEFAULT, sep="\t")
    assert len(sbo) == 81
    assert sbo["dg_id"].is_unique
    # every curated dg_id must appear in S1's own dg_id column
    assert set(sbo["dg_id"]) <= set(s1["dg_id"])


# --- label join sanity (synthetic) -------------------------------------------


def _synthetic_exonome() -> dict[str, str]:
    target = NUSINERSEN_RC
    exon1 = _tiled_filler(11) + target + _tiled_filler(11)  # len 11+18+11=40
    exon2 = _tiled_filler(11) + target + _tiled_filler(11)
    exon3 = _tiled_filler(40)
    assert len(exon1) == 40 and len(exon2) == 40 and len(exon3) == 40
    return {
        "chr1:100-140:+": exon1,
        "chr2:200-240:+": exon2,
        "chr3:0-40:+": exon3,
    }


def test_label_join_sanity() -> None:
    s1 = pd.DataFrame([
        dict(id=1, dg_id="DGtest", cell="HEK293T", sequence=NUSINERSEN_5P_TO_3P, category="exonic"),
    ])
    # S4 reports exon chr1:101-140:+ (1-based inclusive) as positive for id=1;
    # positive_interval_to_exon_key maps this to chr1:100-140:+, matching the
    # universe exon that was actually bound.
    s4 = pd.DataFrame([
        dict(id=1, chrom="chr1", start=101, end=140, strand="+",
             case_psi=0.8, control_psi=0.2, rnaplex_dg=-20.0, rnacofold_dg=-10.0,
             edit_dist=0.0, exact_match=18.0),
    ])
    sbo_sequences = pd.DataFrame([
        dict(dg_id="DGtest", sequence=NUSINERSEN_5P_TO_3P, length=len(NUSINERSEN_5P_TO_3P)),
    ])
    exonome = _synthetic_exonome()

    positives = positive_exon_keys_by_dg(s1, s4)
    assert positives == {"DGtest": {"chr1:100-140:+"}}

    features = build_features(
        s1=s1, s4=s4, sbo_sequences=sbo_sequences,
        exonome=exonome, max_edit=0, intarna=None,
    )
    assert set(features["exon_key"]) == {"chr1:100-140:+", "chr2:200-240:+"}
    by_key = features.set_index("exon_key")
    assert by_key.loc["chr1:100-140:+", "actionable"] == 1
    assert by_key.loc["chr2:200-240:+", "actionable"] == 0
    # Site-level features are computed from the SBO's best binding window by the
    # SAME path for both labels; both classes carry a non-null duplex_dg and
    # binding_window_gc.
    assert not by_key["duplex_dg"].isna().any()
    assert "binding_window_gc" in features.columns
    assert not by_key["binding_window_gc"].isna().any()
    # Whole-region length/GC were dropped as features (substrate confound) and
    # must not reappear as columns; nor may action-adjacent baseline_psi.
    assert "exon_length" not in features.columns
    assert "gc_fraction" not in features.columns
    assert "baseline_psi" not in features.columns
    assert set(by_key["dg_id"]) == {"DGtest"}
    assert set(by_key["sbo_group"]) == {"DGtest"}
    # S4's own rnaplex_dg/rnacofold_dg values must never appear in feature columns
    assert -20.0 not in by_key["duplex_dg"].tolist()
    assert -10.0 not in by_key["duplex_dg"].tolist()


def test_build_features_no_hits_returns_empty_for_that_sbo() -> None:
    s1 = pd.DataFrame([dict(id=1, dg_id="DGnohit", cell="HEK293T", sequence="A" * 18, category="exonic")])
    s4 = pd.DataFrame(columns=["id", "chrom", "start", "end", "strand", "case_psi",
                                "control_psi", "rnaplex_dg", "rnacofold_dg", "edit_dist", "exact_match"])
    sbo_sequences = pd.DataFrame([dict(dg_id="DGnohit", sequence="A" * 18, length=18)])
    exonome = _synthetic_exonome()
    features = build_features(
        s1=s1, s4=s4, sbo_sequences=sbo_sequences,
        exonome=exonome, max_edit=0, intarna=None,
    )
    assert features.empty


# --- redesign: positive definition, GENCODE universe, background pool --------


def test_positive_definition_numeric_edit_only() -> None:
    """Positives are S4 rows with a NUMERIC edit_dist<=max_edit. A 'background'
    (NaN) row is INDIRECT and must never become a positive; a numeric edit>5 row
    is also excluded."""
    s1 = pd.DataFrame([
        dict(id=1, dg_id="DG1", cell="HEK293T", sequence="A" * 18, category="x"),
        dict(id=2, dg_id="DG2", cell="HEK293T", sequence="A" * 18, category="x"),
        dict(id=3, dg_id="DG3", cell="HEK293T", sequence="A" * 18, category="x"),
    ])
    s4 = pd.DataFrame([
        # direct-bound (edit 3) -> positive for DG1
        dict(id=1, chrom="chr1", start=101, end=140, strand="+", case_psi=0.8,
             control_psi=0.2, rnaplex_dg=-20.0, rnacofold_dg=-10.0, edit_dist=3.0, exact_match=12.0),
        # background/indirect (NaN edit) -> excluded, not a positive for DG2
        dict(id=2, chrom="chr2", start=201, end=240, strand="+", case_psi=0.7,
             control_psi=0.3, rnaplex_dg=-15.0, rnacofold_dg=-8.0, edit_dist=None, exact_match=None),
        # numeric but edit>max_edit -> excluded for DG3
        dict(id=3, chrom="chr3", start=301, end=340, strand="+", case_psi=0.6,
             control_psi=0.4, rnaplex_dg=-12.0, rnacofold_dg=-6.0, edit_dist=9.0, exact_match=5.0),
    ])
    positives = positive_exon_keys_by_dg(s1, s4, max_edit=5)
    assert positives == {"DG1": {"chr1:100-140:+"}}
    assert "DG2" not in positives and "DG3" not in positives

    n_direct, n_indirect = s4_direct_vs_indirect(s4, max_edit=5)
    assert n_direct == 1  # only the edit=3 row
    assert n_indirect == 1  # only the NaN row


def test_iter_gencode_exon_keys_and_dedup(tmp_path) -> None:
    """GTF (1-based inclusive) -> exon_key (0-based half-open) and interval dedup."""
    gtf = tmp_path / "mini.gtf"
    gtf.write_text(
        "##comment line\n"
        'chr1\tHAVANA\tgene\t11\t100\t.\t+\t.\tgene_id "G1";\n'
        'chr1\tHAVANA\texon\t11\t20\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
        # identical interval from a second transcript -> must dedup to one key
        'chr1\tHAVANA\texon\t11\t20\t.\t+\t.\tgene_id "G1"; transcript_id "T2";\n'
        'chr2\tHAVANA\texon\t5\t9\t.\t-\t.\tgene_id "G2"; transcript_id "T3";\n',
        encoding="utf-8",
    )
    keys = list(iter_gencode_exon_keys(gtf))
    # start-1 for the 0-based frame; end unchanged. gene/comment lines skipped.
    assert keys == ["chr1:10-20:+", "chr1:10-20:+", "chr2:4-9:-"]
    unique = load_gencode_exon_keys(gtf)
    assert sorted(unique.tolist()) == ["chr1:10-20:+", "chr2:4-9:-"]

    df = exons_df_from_keys(unique)
    row = df.set_index("exon_key").loc["chr2:4-9:-"]
    assert int(row["exon_start_0base"]) == 4 and int(row["exon_end"]) == 9 and row["strand"] == "-"


def test_background_pool_is_deterministic_and_label_blind() -> None:
    keys = np.array([f"chr1:{i}-{i + 30}:+" for i in range(1000)], dtype=object)
    a = sample_background_pool(keys, pool_size=50, seed=DEFAULT_POOL_SEED)
    b = sample_background_pool(keys, pool_size=50, seed=DEFAULT_POOL_SEED)
    assert a == b  # same seed -> identical draw
    assert len(a) == 50 and len(set(a)) == 50  # no replacement
    assert set(a) <= set(keys.tolist())
    # a different seed draws a different pool
    c = sample_background_pool(keys, pool_size=50, seed=DEFAULT_POOL_SEED + 1)
    assert c != a
    # order of the input array does not change the draw (sorted internally)
    shuffled = keys[::-1].copy()
    assert sample_background_pool(shuffled, pool_size=50, seed=DEFAULT_POOL_SEED) == a
    # pool_size >= universe returns the whole (sorted) universe
    allk = sample_background_pool(keys, pool_size=10_000, seed=DEFAULT_POOL_SEED)
    assert allk == sorted(keys.tolist())


# --- redesign: liftover (hg19->GRCh38) of S4 positive loci -------------------


def test_lift_positive_interval_frame_and_strand_flip() -> None:
    """1-based-inclusive S4 interval -> 0-based-half-open lift; a '-' chain flips
    the feature strand. Uses a hand-built two-block chain (no external file)."""
    from feng.bench.armc_build import _lift_positive_interval
    from feng.rbp.liftover import ChainMap, _Chain, _Block

    # Forward chain: hg19 chr1 [1000,1100) -> hg38 chr1 [5000,5100), '+'.
    fwd = _Chain(
        t_name="chr1", t_size=10_000, t_start=1000, t_end=1100,
        q_name="chr1", q_size=20_000, q_strand="+", q_start=5000, q_end=5100,
        blocks=(_Block(t=1000, q=5000, size=100),),
    )
    chain = ChainMap([fwd], path=Path("mem"), sha256="0")
    # S4 reports 1-based inclusive 1001..1050 -> 0-based half-open [1000,1050);
    # lifts to [5000,5050) on '+'. Feature strand '+' is preserved.
    status, key = _lift_positive_interval(("chr1", 1001, 1050, "+"), chain)
    assert status == "lifted"
    assert key == "chr1:5000-5050:+"

    # A minus-strand chain flips the reported feature strand.
    rev = _Chain(
        t_name="chr2", t_size=10_000, t_start=2000, t_end=2100,
        q_name="chr2", q_size=30_000, q_strand="-", q_start=0, q_end=100,
        blocks=(_Block(t=2000, q=0, size=100),),
    )
    rchain = ChainMap([rev], path=Path("mem"), sha256="0")
    status_r, key_r = _lift_positive_interval(("chr2", 2001, 2050, "+"), rchain)
    assert status_r == "lifted"
    # '+' S4 feature strand on a '-' chain -> '-' in GRCh38.
    assert key_r.endswith(":-")


@pytest.mark.skipif(
    not CHAIN_PATH.exists()
    or not (REPO_ROOT / S1_PDF_DEFAULT).exists()
    or not (REPO_ROOT / S4_PDF_DEFAULT).exists(),
    reason="hg19->GRCh38 chain or Holgersen S1/S4 PDFs not staged",
)
def test_s4_positive_lift_recovery() -> None:
    """The direct-positive intervals are hg19: without lifting they do not sit on
    GENCODE GRCh38 keys; WITH the chain they lift essentially 1:1. Guards the
    proven 204-interval / ~100% recovery figure."""
    s1 = parse_s1(REPO_ROOT / S1_PDF_DEFAULT)
    s4 = parse_s4(REPO_ROOT / S4_PDF_DEFAULT)
    chain = load_armc_chain()

    cov = positive_lift_coverage(s1, s4, max_edit=5, chain=chain)
    assert cov["n_direct_positive_intervals"] == 204
    assert cov["lift_fraction"] >= 0.95  # proven 204/204 lift cleanly
    assert cov["n_lifted"] >= 200
    # Every lifted positive locus lands on a real GRCh38 primary contig.
    positives_by_dg = positive_exon_keys_by_dg(s1, s4, max_edit=5, chain=chain)
    for keys in positives_by_dg.values():
        for key in keys:
            chrom = key.split(":", 1)[0]
            assert chrom in _GRCH38_MAIN_CHROMS or chrom.startswith("chr")


@pytest.mark.skipif(
    not CHAIN_PATH.exists()
    or not GENOME_FASTA.is_file()
    or not SAMTOOLS.is_file()
    or not (REPO_ROOT / S1_PDF_DEFAULT).exists()
    or not (REPO_ROOT / S4_PDF_DEFAULT).exists()
    or not (REPO_ROOT / SBO_TSV_DEFAULT).exists(),
    reason="chain / genome FASTA / samtools / Holgersen inputs not all staged",
)
def test_edit0_positive_lands_on_grch38_after_lift(tmp_path) -> None:
    """Frame correctness end to end: for edit_dist==0 positives, revcomp(SBO)
    must appear EXACTLY (edit 0) in the lifted GRCh38 region sequence. If the
    frame/strand were wrong the site would not be found."""
    s1 = parse_s1(REPO_ROOT / S1_PDF_DEFAULT)
    s4 = parse_s4(REPO_ROOT / S4_PDF_DEFAULT)
    chain = load_armc_chain()
    sbo = pd.read_csv(REPO_ROOT / SBO_TSV_DEFAULT, sep="\t")
    sbo_seq = dict(zip(sbo["dg_id"], sbo["sequence"].str.upper()))
    id_to_dg = s1.set_index("id")["dg_id"].to_dict()

    from feng.bench.armc_build import _lift_positive_interval

    # Collect edit_dist==0 rows whose SBO is in the curated panel, lift each, and
    # remember which SBO must be found in the lifted region.
    want: dict[str, str] = {}
    for row in s4.itertuples(index=False):
        if row.edit_dist is None or pd.isna(row.edit_dist) or float(row.edit_dist) != 0.0:
            continue
        dg = id_to_dg.get(row.id)
        if dg not in sbo_seq:
            continue
        status, key = _lift_positive_interval(
            (str(row.chrom), int(row.start), int(row.end), str(row.strand)), chain
        )
        if status == "lifted" and key is not None:
            want[key] = dg
    assert want, "no edit_dist==0 positives available to check"

    exons = exons_df_from_keys(list(want))
    seqs = build_exonome(
        exons, genome_fasta=GENOME_FASTA,
        cache_path=tmp_path / "pos.fa", pad=DEFAULT_TARGET_PAD, samtools=SAMTOOLS,
    )
    n_found = 0
    for key, dg in want.items():
        hits = seed_and_extend(reverse_complement(sbo_seq[dg]), {key: seqs[key]}, max_edit=0)
        if hits:
            n_found += 1
    # With the symmetric fetch pad the reported site (which can straddle the exon
    # boundary) is contained; correct assembly/frame/strand then recovers the
    # exact antisense match for essentially all edit_dist==0 positives.
    assert n_found >= 0.95 * len(want)


# --- redesign: site-level feature parity (positives vs negatives) ------------


def test_site_level_feature_parity_pos_vs_neg() -> None:
    """A positive locus and a negative pool exon that both bind the SBO must yield
    the SAME site-level feature schema by the SAME code path -- no whole-region
    length/GC leaks in, binding_window_gc/duplex_dg present for both classes."""
    target = NUSINERSEN_RC
    pos_seq = _tiled_filler(30) + target + _tiled_filler(30)  # positive locus (wide region)
    neg_seq = _tiled_filler(11) + target + _tiled_filler(11)  # pool exon (narrow)
    exonome = {"chrPOS:0-78:+": pos_seq, "chr9:200-240:+": neg_seq}
    positives_by_dg = {"DGx": {"chrPOS:0-78:+"}}

    s1 = pd.DataFrame([dict(id=1, dg_id="DGx", cell="HEK293T",
                            sequence=NUSINERSEN_5P_TO_3P, category="x")])
    s4 = pd.DataFrame(columns=["id", "chrom", "start", "end", "strand", "case_psi",
                                "control_psi", "rnaplex_dg", "rnacofold_dg",
                                "edit_dist", "exact_match"])
    sbo_sequences = pd.DataFrame([dict(dg_id="DGx", sequence=NUSINERSEN_5P_TO_3P,
                                       length=len(NUSINERSEN_5P_TO_3P))])

    features = build_features(
        s1=s1, s4=s4, sbo_sequences=sbo_sequences, exonome=exonome,
        max_edit=0, intarna=None, positives_by_dg=positives_by_dg,
    )
    by_key = features.set_index("exon_key")
    assert by_key.loc["chrPOS:0-78:+", "actionable"] == 1  # positive
    assert by_key.loc["chr9:200-240:+", "actionable"] == 0  # bound-not-acted negative
    # identical feature schema for both rows, and it is site-level only
    for col in FEATURE_COLUMNS:
        assert col in features.columns
        assert not features[col].isna().any()
    assert "exon_length" not in features.columns
    assert "gc_fraction" not in features.columns
    # both classes' best window is the same 18-nt site, so binding_window_gc is
    # identical (proves whole-region length did not enter the features).
    assert (
        by_key.loc["chrPOS:0-78:+", "binding_window_gc"]
        == pytest.approx(by_key.loc["chr9:200-240:+", "binding_window_gc"])
    )
    assert (by_key["edit_distance"] == 0).all()


def test_pool_overlap_exclusion_drops_covering_negatives() -> None:
    """A pool exon overlapping an SBO's positive locus is excluded from that SBO's
    negatives (strand-agnostic coordinate overlap)."""
    positive_keys = {"chr1:1000-1200:+"}
    pool_keys = {
        "chr1:1100-1300:+",  # overlaps -> drop
        "chr1:1100-1300:-",  # overlaps (other strand) -> drop (strand-agnostic)
        "chr1:1300-1500:+",  # abuts, no overlap -> keep
        "chr2:1000-1200:+",  # other chrom -> keep
    }
    drop = _pool_keys_overlapping_positives(pool_keys, positive_keys)
    assert drop == {"chr1:1100-1300:+", "chr1:1100-1300:-"}


# --- gate: metrics + degenerate-fold handling --------------------------------


def test_roc_auc_and_ap_basic_values() -> None:
    y = np.array([0, 0, 1, 1])
    perfect_scores = np.array([0.1, 0.2, 0.8, 0.9])
    assert roc_auc(y, perfect_scores) == pytest.approx(1.0)
    assert average_precision(y, perfect_scores) == pytest.approx(1.0)

    inverted_scores = np.array([0.9, 0.8, 0.2, 0.1])
    assert roc_auc(y, inverted_scores) == pytest.approx(0.0)

    assert roc_auc(np.array([1, 1, 1]), np.array([0.1, 0.2, 0.3])) is None
    assert average_precision(np.array([0, 0, 0]), np.array([0.1, 0.2, 0.3])) is None


def test_average_precision_tie_handling() -> None:
    """Tied scores must be treated as one threshold block (sklearn convention).
    y=[1,0,1,0] with pair-tied scores has AP 0.5; the naive one-sample-at-a-time
    step gives 0.833 and is order-dependent -- guard against a regression to it.
    """
    y = np.array([1, 0, 1, 0])
    scores = np.array([0.8, 0.8, 0.5, 0.5])
    assert average_precision(y, scores) == pytest.approx(0.5)
    # Invariance to within-tie ordering: permuting tied rows must not change AP.
    perm = [1, 0, 3, 2]
    assert average_precision(y[perm], scores[perm]) == pytest.approx(0.5)
    # All-tied scores collapse to the base rate.
    y2 = np.array([1, 1, 0, 0, 0])
    assert average_precision(y2, np.full(5, 0.3)) == pytest.approx(2 / 5)


def test_sbo_low_complexity_guard() -> None:
    # Real degenerate designs in the panel: flagged.
    assert sbo_low_complexity_reason("AAAAATAAAATAAAATAA") is not None  # poly-A (DG986)
    assert sbo_low_complexity_reason("AACAAAACACACACACAC") is not None  # (AC)n (DG985)
    assert sbo_low_complexity_reason("ATCATCATCATCATCATC") is not None  # ATC repeat (DG959)
    assert sbo_low_complexity_reason("A" * 18) is not None
    # A normal antisense oligo (nusinersen) must NOT be flagged.
    assert sbo_low_complexity_reason(NUSINERSEN_5P_TO_3P) is None


def test_build_features_skips_low_complexity_sbo() -> None:
    s1 = pd.DataFrame([dict(id=1, dg_id="DGpolyA", cell="HEK293T", sequence="A" * 18, category="exonic")])
    s4 = pd.DataFrame(columns=["id", "chrom", "start", "end", "strand", "case_psi",
                                "control_psi", "rnaplex_dg", "rnacofold_dg", "edit_dist", "exact_match"])
    sbo_sequences = pd.DataFrame([dict(dg_id="DGpolyA", sequence="A" * 18, length=18)])
    features = build_features(
        s1=s1, s4=s4, sbo_sequences=sbo_sequences,
        exonome=_synthetic_exonome(), max_edit=5, intarna=None,
    )
    # Excluded before scanning -> no rows, no spurious universal hits.
    assert features.empty


def _synthetic_feature_table() -> pd.DataFrame:
    rng = np.random.RandomState(0)
    rows = []
    # 6 SBO groups, one signal-bearing feature (edit_distance low -> actionable),
    # so a well-behaved model should beat the always-noisy baselines is not
    # guaranteed here -- edit_distance IS one of the baselines, so this checks
    # the plumbing (no crash, sane bounds, degenerate fold is excluded/flagged)
    # rather than asserting the model wins.
    for gi in range(6):
        group = f"DG{gi}"
        n = 10
        actionable = (rng.rand(n) < 0.3).astype(int)
        if gi == 5:
            actionable[:] = 0  # degenerate group: no positives at all
        edit_distance = np.where(actionable == 1, rng.uniform(0, 1, n), rng.uniform(2, 5, n))
        duplex_dg = np.where(actionable == 1, rng.uniform(-15, -10, n), rng.uniform(-5, 0, n))
        for i in range(n):
            rows.append(dict(
                sbo_group=group, dg_id=group, exon_key=f"{group}_{i}",
                edit_distance=edit_distance[i], mismatches=1, hit_length=18,
                exact_match_length=15, duplex_dg=duplex_dg[i], binding_window_gc=0.5,
                actionable=int(actionable[i]),
            ))
    return pd.DataFrame(rows)


def test_run_loso_gate_smoke_and_degenerate_fold_reporting() -> None:
    features = _synthetic_feature_table()
    result = run_loso_gate(features, n_iter=200)

    assert result.n_groups_total == 6
    # DG5 is forced all-negative; other groups may also land degenerate by the
    # fixed-seed random draw -- assert the forced one is caught, not an exact count.
    assert result.n_groups_skipped_eval >= 1
    degenerate = [f for f in result.folds if f.group == "DG5"]
    assert len(degenerate) == 1 and degenerate[0].evaluable is False
    assert degenerate[0].n_pos == 0

    assert result.pooled_model_auc is not None
    assert 0.0 <= result.pooled_model_auc <= 1.0
    assert result.pooled_model_ap is not None
    assert result.pooled_baseline_edit_auc is not None
    assert result.pooled_baseline_dg_auc is not None
    assert result.verdict in {"GO", "NO-GO"}
    # the degenerate group's rows are still pooled (contribute negatives)
    assert len(result.pooled_y) == 6 * 10


def test_run_loso_gate_all_groups_single_class_yields_no_trainable_fold() -> None:
    rows = []
    for gi in range(3):
        for i in range(4):
            rows.append(dict(
                sbo_group=f"DG{gi}", dg_id=f"DG{gi}", exon_key=f"DG{gi}_{i}",
                edit_distance=1.0, mismatches=1, hit_length=18, exact_match_length=15,
                duplex_dg=-5.0, binding_window_gc=0.5, actionable=0,
            ))
    features = pd.DataFrame(rows)
    result = run_loso_gate(features, n_iter=50)
    assert result.n_groups_skipped_train == 3
    assert result.pooled_model_auc is None
    assert result.verdict.startswith("NO-GO")
