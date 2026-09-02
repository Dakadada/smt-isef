"""Freeze-guard + parsing tests for the Arm D detectable-exon universe."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
import pytest

from feng.bench.strata import (
    CONTRAST_FIELDS,
    DetectabilityThresholds,
    FreezeGuardError,
    GenomeFastaError,
    _resolve_fasta_contig,
    assert_no_contrast_fields,
    attach_gc_content,
    classify_detectable,
    collapse_to_unique_exons,
    load_se_events,
    positive_interval_to_exon_key,
    select_control_runs,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
SBATCH = REPO_ROOT / "feng" / "scripts" / "bench_detectable_universe.sbatch"
POSITIVES = REPO_ROOT / "feng" / "artifacts" / "bench" / "positives_v50.tsv"


def _conditions() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"run_accession": "SRR_ctrl_1", "condition": "control", "chemistry": "OMe", "control_type": "scrambled"},
            {"run_accession": "SRR_ctrl_2", "condition": "control", "chemistry": "untreated", "control_type": "untreated"},
            {"run_accession": "SRR_trt_1", "condition": "F18OMe", "chemistry": "OMe", "control_type": "targeted"},
            {"run_accession": "SRR_trt_2", "condition": "F20PMO", "chemistry": "PMO", "control_type": "targeted"},
        ]
    )


def test_select_control_runs_keeps_only_controls() -> None:
    selected = select_control_runs(_conditions())
    assert set(selected["run_accession"]) == {"SRR_ctrl_1", "SRR_ctrl_2"}
    assert not selected["control_type"].eq("targeted").any()
    assert not selected["condition"].isin({"F18OMe", "F18MOE", "F20PMO"}).any()


def test_select_control_runs_rejects_run_double_listed_as_treatment() -> None:
    # The same run accession appearing in both a control row and a targeted row
    # must abort the build rather than be trusted from its control row alone.
    bad = pd.concat(
        [
            _conditions(),
            pd.DataFrame(
                [{
                    "run_accession": "SRR_ctrl_1",
                    "condition": "F18OMe",
                    "chemistry": "OMe",
                    "control_type": "targeted",
                }]
            ),
        ],
        ignore_index=True,
    )
    with pytest.raises(FreezeGuardError):
        select_control_runs(bad)


def test_select_control_runs_empty_selection_aborts() -> None:
    treatment_only = _conditions()
    treatment_only = treatment_only[treatment_only["control_type"] == "targeted"]
    with pytest.raises(FreezeGuardError):
        select_control_runs(treatment_only)


def test_select_control_runs_rejects_duplicate_controls() -> None:
    duplicate = pd.concat([_conditions(), _conditions().iloc[[0]]], ignore_index=True)
    with pytest.raises(FreezeGuardError, match="duplicate control"):
        select_control_runs(duplicate)


def test_load_se_events_pools_both_control_groups(tmp_path) -> None:
    se = tmp_path / "SE.MATS.JC.txt"
    header = [
        "ID", "GeneID", "geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd",
        "upstreamES", "upstreamEE", "downstreamES", "downstreamEE",
        "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
        "IncFormLen", "SkipFormLen", "PValue", "FDR",
        "IncLevel1", "IncLevel2", "IncLevelDifference",
    ]
    row = [
        "1", '"ENSG0001.3"', '"GENEA"', "chr1", "+", "1000", "1100",
        "800", "900", "1200", "1300",
        "20,30", "1,1", "40,50", "0,2",
        "198", "99", "NA", "NA",
        "0.90,0.92", "0.95,0.80", "NA",
    ]
    se.write_text("\t".join(header) + "\n" + "\t".join(row) + "\n", encoding="utf-8")
    events = load_se_events(se, expected_sample_count=4)
    assert len(events) == 1
    e = events.iloc[0]
    assert e["gene_id_base"] == "ENSG0001"
    assert e["n_samples_total"] == 4  # both rMATS groups pooled
    assert e["n_samples_with_psi"] == 4
    assert e["exon_length"] == 100
    # mean total junction reads = mean([21, 31, 40, 52]) = 36.0
    assert e["mean_total_junction_reads"] == pytest.approx(36.0)
    assert e["baseline_psi"] == pytest.approx((0.90 + 0.92 + 0.95 + 0.80) / 4)
    assert e["exon_key"] == "chr1:1000-1100:+"
    assert CONTRAST_FIELDS.isdisjoint(events.columns)


def test_classify_detectable_thresholds() -> None:
    events = pd.DataFrame(
        [
            {"mean_total_junction_reads": 50.0, "n_samples_with_psi": 8, "baseline_psi": 0.9},
            {"mean_total_junction_reads": 5.0, "n_samples_with_psi": 8, "baseline_psi": 0.9},
            {"mean_total_junction_reads": 50.0, "n_samples_with_psi": 2, "baseline_psi": 0.9},
            {"mean_total_junction_reads": 50.0, "n_samples_with_psi": 8, "baseline_psi": float("nan")},
        ]
    )
    out = classify_detectable(events, DetectabilityThresholds())
    assert list(out["detectable"]) == [True, False, False, False]


def test_collapse_to_unique_exons_prefers_detectable_best_covered_event() -> None:
    events = pd.DataFrame(
        [
            {
                "event_id": "1", "chrom": "chr1", "strand": "+",
                "exon_start_0base": 100, "exon_end": 200,
                "detectable": False, "n_samples_with_psi": 11,
                "mean_total_junction_reads": 100.0,
                "gene_id": "ENSG1.1", "gene_id_base": "ENSG1",
                "gene_symbol": "GENEA", "control_mean_cpm": 10.0,
            },
            {
                "event_id": "2", "chrom": "chr1", "strand": "+",
                "exon_start_0base": 100, "exon_end": 200,
                "detectable": True, "n_samples_with_psi": 8,
                "mean_total_junction_reads": 30.0,
                "gene_id": "ENSG2.1", "gene_id_base": "ENSG2",
                "gene_symbol": "GENEB", "control_mean_cpm": 20.0,
            },
            {
                "event_id": "3", "chrom": "chr1", "strand": "+",
                "exon_start_0base": 300, "exon_end": 350,
                "detectable": True, "n_samples_with_psi": 9,
                "mean_total_junction_reads": 40.0,
                "gene_id": "ENSG3.1", "gene_id_base": "ENSG3",
                "gene_symbol": "GENEC", "control_mean_cpm": 30.0,
            },
        ]
    )
    collapsed = collapse_to_unique_exons(events)
    assert len(collapsed) == 2
    first = collapsed[collapsed["exon_start_0base"] == 100].iloc[0]
    assert first["event_id"] == "2"
    assert first["n_rmats_events"] == 2
    assert first["rmats_event_ids"] == "1,2"
    assert first["all_gene_ids"] == "ENSG1.1,ENSG2.1"
    assert first["all_gene_id_bases"] == "ENSG1,ENSG2"
    assert first["all_gene_symbols"] == "GENEA,GENEB"
    assert first["n_gene_ids"] == 2
    assert bool(first["multi_gene"])
    assert first["control_mean_cpm_by_gene"] == '{"ENSG1":10.0,"ENSG2":20.0}'
    assert not collapsed["exon_key"].duplicated().any()


def test_contrast_fields_are_rejected_before_publication() -> None:
    with pytest.raises(FreezeGuardError, match="contrast fields"):
        assert_no_contrast_fields(pd.DataFrame({"IncLevelDifference": [0.2]}))


def test_gc_content_fails_closed_without_fasta(tmp_path) -> None:
    events = pd.DataFrame(
        [{
            "chrom": "chr1", "exon_start_0base": 0, "exon_end": 4,
            "exon_length": 4, "exon_key": "chr1:0-4:+",
        }]
    )
    with pytest.raises(GenomeFastaError, match="missing or empty"):
        attach_gc_content(events, tmp_path / "missing.fa")


def test_gc_content_resolves_ucsc_prefix_for_unplaced_contigs() -> None:
    assert _resolve_fasta_contig("chrGL000008.2", {"GL000008.2"}) == "GL000008.2"
    with pytest.raises(GenomeFastaError, match="lacks contig"):
        _resolve_fasta_contig("chrMissing", {"chr1"})


def test_gc_content_never_returns_null(monkeypatch, tmp_path) -> None:
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nGCGT\n", encoding="utf-8")
    events = pd.DataFrame(
        [{
            "chrom": "chr1", "exon_start_0base": 0, "exon_end": 4,
            "exon_length": 4, "exon_key": "chr1:0-4:+",
        }]
    )
    monkeypatch.setattr(
        "feng.bench.strata._attach_gc_with_pysam", lambda events, fasta: [0.75]
    )
    out = attach_gc_content(events, fasta)
    assert out["gc_fraction"].tolist() == [0.75]
    assert out["gc_fraction"].isna().sum() == 0


def test_sbatch_pins_and_hashes_allowlist_before_rmats() -> None:
    script = SBATCH.read_text(encoding="utf-8")
    controls = [
        "SRR35298717", "SRR35298701", "SRR35298704", "SRR35298706",
        "SRR35298709", "SRR35298711", "SRR35298718", "SRR35298699",
        "SRR35298708", "SRR35298713", "SRR35298715",
    ]
    digest = hashlib.sha256(("\n".join(controls) + "\n").encode()).hexdigest()
    assert f'PINNED_CONTROL_SHA256="{digest}"' in script
    assert script.index("actual_allowlist_hash=") < script.index('"${PYTHON_BIN}" "${RMATS}"')
    assert "duplicate accession in pinned control allowlist" in script
    assert "duplicate control accession in conditions table" in script
    assert script.index("validate_published_artifact()") < script.index("DONE_FILE=")
    assert "duplicate exon_key values" in script
    assert "null gc_fraction values" in script
    assert "provenance SHA-256 mismatch" in script
    # The .done fast-path must re-check the freeze/labelling guarantees, not
    # just schema, so a stale sentinel cannot skip the applied fixes.
    assert "forbidden contrast fields" in script
    assert "does not match parquet" in script


def test_load_se_events_rejects_mismatched_or_wrong_sample_counts(tmp_path) -> None:
    se = tmp_path / "SE.MATS.JC.txt"
    header = [
        "ID", "GeneID", "geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd",
        "upstreamES", "upstreamEE", "downstreamES", "downstreamEE",
        "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
        "IncFormLen", "SkipFormLen", "PValue", "FDR",
        "IncLevel1", "IncLevel2", "IncLevelDifference",
    ]
    row = [
        "7", '"ENSG0001.3"', '"GENEA"', "chr1", "+", "1000", "1100",
        "800", "900", "1200", "1300", "20,30", "1", "40,50", "0,2",
        "198", "99", "NA", "NA", "0.9,0.8", "0.7,0.6", "NA",
    ]
    se.write_text("\t".join(header) + "\n" + "\t".join(row) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="IJC/SJC sample-count mismatch"):
        load_se_events(se, expected_sample_count=4)

    row[12] = "1,2"
    se.write_text("\t".join(header) + "\n" + "\t".join(row) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="expected 11 pooled control samples"):
        load_se_events(se)


def test_positive_intervals_match_figure_3b_exon_lengths() -> None:
    positives = pd.read_csv(POSITIVES, sep="\t")
    expected_lengths = {
        "POLR2H": 84, "PITHD1": 44, "SERPINB7": 117, "RTTN": 90,
        "REV3L": 85, "PRKRA": 170, "GOLGA4": 66, "PAK1": 211,
    }
    observed = (positives["end"] - positives["start"] + 1).groupby(
        positives["gene_symbol"]
    ).first().to_dict()
    assert observed == expected_lengths
    assert positives.loc[positives["gene_symbol"] == "POLR2H", "exon_number"].item() == 3
    assert positives.loc[positives["gene_symbol"] == "REV3L", "exon_number"].item() == 8
    golga4 = positives.loc[positives["gene_symbol"] == "GOLGA4"].iloc[0]
    assert golga4["gencode_v50_transcript_id"] == "ENST00001121839.1"
    assert golga4["gencode_v50_exon_id"] == "ENSE00003978486.1"
    assert (golga4["start"], golga4["end"]) == (37273536, 37273601)


def test_positive_interval_maps_to_universe_zero_based_key() -> None:
    # positives_v50.tsv is 1-based inclusive; the universe exon_key is 0-based
    # half-open. The converter must drop exactly one base at the start and leave
    # the end unchanged, preserving interval length.
    key = positive_interval_to_exon_key("chr3", 37273536, 37273601, "+")
    assert key == "chr3:37273535-37273601:+"
    start, end = key.split(":")[1].split("-")
    assert int(end) - int(start) == 37273601 - 37273536 + 1
    with pytest.raises(ValueError, match="non-positive interval"):
        positive_interval_to_exon_key("chr1", 100, 99, "+")


def test_positive_converter_matches_every_recorded_positive() -> None:
    positives = pd.read_csv(POSITIVES, sep="\t")
    for _, r in positives.iterrows():
        key = positive_interval_to_exon_key(
            r["chrom"], r["start"], r["end"], r["strand"]
        )
        assert key == f"{r['chrom']}:{int(r['start']) - 1}-{int(r['end'])}:{r['strand']}"
