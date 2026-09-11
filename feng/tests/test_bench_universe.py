"""Freeze-guard + parsing tests for the Arm D detectable-exon universe."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pandas as pd
import pytest
import feng.bench.strata as strata_module

from feng.bench.strata import (
    ALLOWED_GENCODE_TRANSCRIPT_TYPES,
    CONTRAST_FIELDS,
    DetectabilityThresholds,
    FreezeGuardError,
    GenomeFastaError,
    DETECTABILITY_POLICY,
    DETECTABILITY_POLICY_ID,
    FROZEN_MIN_MEAN_TOTAL_JUNCTION_READS,
    SCRMOE_MIN_SAMPLES_WITH_PSI,
    _resolve_fasta_contig,
    assert_no_contrast_fields,
    attach_gc_content,
    attach_annotated_junction_metrics,
    attach_input_provenance,
    classify_detectable,
    collapse_to_unique_exons,
    load_control_expression,
    load_annotated_internal_exons,
    load_se_events,
    positive_interval_to_exon_key,
    render_readme,
    require_frozen_thresholds,
    select_control_runs,
    validate_input_provenance,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
SBATCH = REPO_ROOT / "feng" / "scripts" / "bench_detectable_universe.sbatch"
POSITIVES = REPO_ROOT / "feng" / "artifacts" / "bench" / "positives_v50.tsv"


def _conditions() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {"run_accession": "SRR35298704", "condition": "control", "chemistry": "MOE", "control_type": "scrambled"},
            {"run_accession": "SRR35298699", "condition": "control", "chemistry": "MOE", "control_type": "scrambled"},
            {"run_accession": "SRR35298715", "condition": "control", "chemistry": "MOE", "control_type": "scrambled"},
            {"run_accession": "SRR_ctrl_2", "condition": "control", "chemistry": "untreated", "control_type": "untreated"},
            {"run_accession": "SRR_scr_ome", "condition": "control", "chemistry": "OMe", "control_type": "scrambled"},
            {"run_accession": "SRR_trt_1", "condition": "F18OMe", "chemistry": "OMe", "control_type": "targeted"},
            {"run_accession": "SRR_trt_2", "condition": "F20PMO", "chemistry": "PMO", "control_type": "targeted"},
        ]
    )


def test_select_control_runs_keeps_only_controls() -> None:
    selected = select_control_runs(_conditions())
    assert set(selected["run_accession"]) == {"SRR35298704", "SRR35298699", "SRR35298715"}
    assert selected["chemistry"].eq("MOE").all()
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
                    "run_accession": "SRR35298704",
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
    non_scrmoe = _conditions().query("not (control_type == 'scrambled' and chemistry == 'MOE')")
    with pytest.raises(FreezeGuardError):
        select_control_runs(non_scrmoe)


def test_select_control_runs_rejects_duplicate_controls() -> None:
    duplicate = pd.concat([_conditions(), _conditions().iloc[[0]]], ignore_index=True)
    with pytest.raises(FreezeGuardError, match="duplicate control"):
        select_control_runs(duplicate)


def test_select_control_runs_requires_exactly_three_scrmoe() -> None:
    with pytest.raises(FreezeGuardError, match="exactly 3 ScrMOE"):
        select_control_runs(_conditions().query("run_accession != 'SRR35298715'"))


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


def test_load_se_events_tolerates_missing_inclevel_cells(tmp_path) -> None:
    # pd.read_csv(dtype=str) yields float NaN for genuinely-empty cells; a
    # missing IncLevel must be treated as "no PSI" (like an NA cell), not crash.
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
        "0.90,0.92", "", "NA",  # IncLevel2 empty -> NaN
    ]
    se.write_text("\t".join(header) + "\n" + "\t".join(row) + "\n", encoding="utf-8")
    events = load_se_events(se, expected_sample_count=4)
    assert len(events) == 1
    e = events.iloc[0]
    assert e["n_samples_total"] == 4  # junction counts still pool all 4
    assert e["n_samples_with_psi"] == 2  # only IncLevel1's two values survive
    assert e["baseline_psi"] == pytest.approx((0.90 + 0.92) / 2)


def test_missing_star_count_fails_loud(tmp_path) -> None:
    paths = []
    for run in ("SRR35298704", "SRR35298699", "SRR35298715"):
        path = tmp_path / run / f"{run}.ReadsPerGene.out.tab"
        paths.append(path)
    with pytest.raises(FileNotFoundError, match="expected STAR count input"):
        load_control_expression(paths)


def _write_star_counts(tmp_path, run, rows):
    path = tmp_path / run / f"{run}.ReadsPerGene.out.tab"
    path.parent.mkdir(parents=True)
    path.write_text("\n".join("\t".join(map(str, row)) for row in rows) + "\n")
    return path


@pytest.mark.parametrize(
    "bad_row",
    [
        ("ENSG1", "nan", 1, 1),
        ("ENSG1", "inf", 1, 1),
        ("ENSG1", -1, 1, 1),
        ("ENSG1", 1.5, 1, 1),
        ("", 1, 1, 1),
        ("ENSG1", 1, 1, 1, 1),
    ],
)
def test_star_counts_require_finite_nonnegative_integer_schema(tmp_path, bad_row) -> None:
    paths = []
    for run in ("SRR35298704", "SRR35298699", "SRR35298715"):
        paths.append(_write_star_counts(tmp_path, run, [bad_row, ("ENSG2", 2, 2, 2)]))
    with pytest.raises(ValueError):
        load_control_expression(paths)


def test_star_counts_require_unique_gene_ids_and_exact_scrmoe_paths(tmp_path) -> None:
    paths = []
    duplicate_rows = [("ENSG1", 1, 1, 1), ("ENSG1", 2, 2, 2)]
    for run in ("SRR35298704", "SRR35298699", "SRR35298715"):
        paths.append(_write_star_counts(tmp_path, run, duplicate_rows))
    with pytest.raises(ValueError, match="duplicate gene IDs"):
        load_control_expression(paths)
    with pytest.raises(FreezeGuardError, match="exactly the three ScrMOE"):
        load_control_expression(paths[:2])


@pytest.mark.parametrize(
    "rows",
    [
        [("   ", 1, 1, 1), ("ENSG2", 2, 2, 2)],
        [("ENSG1", 1, 1, 1), (" ENSG1 ", 2, 2, 2)],
    ],
)
def test_star_gene_ids_validate_after_whitespace_normalization(tmp_path, rows) -> None:
    paths = [
        _write_star_counts(tmp_path, run, rows)
        for run in ("SRR35298704", "SRR35298699", "SRR35298715")
    ]
    with pytest.raises(ValueError, match="duplicate gene IDs"):
        load_control_expression(paths)


@pytest.mark.parametrize(
    "rows",
    [
        [("\u200b", 1, 1, 1), ("ENSG2", 2, 2, 2)],
        [("\ufeff", 1, 1, 1), ("ENSG2", 2, 2, 2)],
        [("ENSG1", 1, 1, 1), ("ENSG1\u200b", 2, 2, 2)],
        [("ENSG1", 1, 1, 1), ("ＥＮＳＧ１", 2, 2, 2)],
        [("ENSG1/invalid", 1, 1, 1), ("ENSG2", 2, 2, 2)],
    ],
)
def test_star_gene_ids_reject_invisible_unicode_and_non_gencode_grammar(
    tmp_path, rows
) -> None:
    paths = [
        _write_star_counts(tmp_path, run, rows)
        for run in ("SRR35298704", "SRR35298699", "SRR35298715")
    ]
    with pytest.raises(ValueError, match="after normalization"):
        load_control_expression(paths)


def test_classify_detectable_thresholds() -> None:
    events = pd.DataFrame(
        [
            {"mean_total_junction_reads": 20.0, "n_samples_with_psi": 2, "baseline_psi": 0.9},
            {"mean_total_junction_reads": 5.0, "n_samples_with_psi": 8, "baseline_psi": 0.9},
            {"mean_total_junction_reads": 50.0, "n_samples_with_psi": 1, "baseline_psi": 0.9},
            {"mean_total_junction_reads": 50.0, "n_samples_with_psi": 8, "baseline_psi": float("nan")},
        ]
    )
    out = classify_detectable(events, DetectabilityThresholds())
    assert list(out["detectable"]) == [True, False, False, False]


def test_frozen_threshold_constants_and_real_extraction_path_invariance(tmp_path) -> None:
    assert FROZEN_MIN_MEAN_TOTAL_JUNCTION_READS == 20.0
    assert SCRMOE_MIN_SAMPLES_WITH_PSI == 2
    assert DetectabilityThresholds() == DetectabilityThresholds(
        min_mean_total_junction_reads=FROZEN_MIN_MEAN_TOTAL_JUNCTION_READS,
        min_samples_with_psi=SCRMOE_MIN_SAMPLES_WITH_PSI,
    )
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text("".join(
        f'chr1\tT\texon\t{start}\t{start + 49}\t.\t+\t.\tgene_id "ENSG1.1"; '
        f'transcript_id "ENST1.1"; transcript_type "protein_coding"; '
        f'gene_name "GENE1"; exon_id "E{i}";\n'
        for i, start in enumerate((1, 101, 201, 301, 401), 1)
    ), encoding="utf-8")
    junctions = []
    boundary_counts = ((20, 10, 0, 10), (2, 2, 0, 2), (30, 30, 0, 0))
    for b1, b2, b3, b4 in boundary_counts:
        junctions.append({
            ("chr1", 51, 100, "+"): b1, ("chr1", 151, 200, "+"): b2,
            ("chr1", 251, 300, "+"): b3, ("chr1", 351, 400, "+"): b4,
        })
    annotated = load_annotated_internal_exons(gtf)
    annotated = classify_detectable(
        attach_annotated_junction_metrics(annotated, junctions)
    ).set_index("exon_key")

    se = tmp_path / "SE.MATS.JC.txt"
    header = [
        "ID", "GeneID", "geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd",
        "upstreamES", "upstreamEE", "downstreamES", "downstreamEE",
        "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
        "IncFormLen", "SkipFormLen", "PValue", "FDR", "IncLevel1", "IncLevel2",
        "IncLevelDifference",
    ]
    rows = []
    for event_id, start0, upstream_ee, downstream_es, inc in (
        ("E2", 100, 50, 200, (30, 4, 60)),
        ("E3", 200, 150, 300, (10, 2, 30)),
    ):
        psi = ["1" if value else "NA" for value in inc]
        rows.append([
            event_id, '"ENSG1.1"', '"GENE1"', "chr1", "+", str(start0),
            str(start0 + 50), str(upstream_ee - 50), str(upstream_ee),
            str(downstream_es), str(downstream_es + 50), f"{inc[0]},{inc[1]}",
            "0,0", str(inc[2]), "0", "149", "99", "NA", "NA",
            f"{psi[0]},{psi[1]}", psi[2], "NA",
        ])
    se.write_text("\t".join(header) + "\n" + "\n".join(
        "\t".join(row) for row in rows
    ) + "\n", encoding="utf-8")
    rmats = classify_detectable(load_se_events(se)).set_index("exon_key")
    shared_keys = rmats.index.intersection(annotated.index)
    assert len(shared_keys) == 2
    assert annotated.loc[shared_keys, "detectable"].equals(
        rmats.loc[shared_keys, "detectable"]
    )
    assert set(annotated.index).difference(rmats.index) == {"chr1:300-350:+"}


def test_publishable_thresholds_are_frozen() -> None:
    require_frozen_thresholds(DetectabilityThresholds())
    with pytest.raises(FreezeGuardError, match="requires frozen"):
        require_frozen_thresholds(DetectabilityThresholds(19.0, 2))


def test_authoritative_input_provenance_hashes_every_covariate_input(tmp_path) -> None:
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text("annotation-v1\n", encoding="utf-8")
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    sj_paths = []
    count_paths = []
    for run in ("SRR35298704", "SRR35298699", "SRR35298715"):
        path = tmp_path / run / f"{run}.SJ.out.tab"
        path.parent.mkdir()
        path.write_text(f"{run}\n", encoding="utf-8")
        sj_paths.append(path)
        count_path = path.parent / f"{run}.ReadsPerGene.out.tab"
        count_path.write_text(f"{run}-counts\n", encoding="utf-8")
        count_paths.append(count_path)
    out = attach_input_provenance(
        pd.DataFrame({"exon_key": ["x"]}),
        gtf,
        sj_paths,
        fasta,
        count_paths,
        101,
    )
    assert out["gencode_gtf_sha256"].item() == hashlib.sha256(gtf.read_bytes()).hexdigest()
    assert set(__import__("json").loads(out["scrmoe_sj_sha256"].item())) == {
        "SRR35298704", "SRR35298699", "SRR35298715"
    }
    assert out["genome_fasta_sha256"].item() == hashlib.sha256(
        fasta.read_bytes()
    ).hexdigest()
    assert __import__("json").loads(out["scrmoe_reads_per_gene_sha256"].item()) == {
        path.parent.name: hashlib.sha256(path.read_bytes()).hexdigest()
        for path in count_paths
    }
    assert out["junction_read_length"].item() == 101
    assert out["junction_read_length_sha256"].item() == hashlib.sha256(b"101\n").hexdigest()


def test_changed_covariate_input_fails_authoritative_provenance_validation(
    tmp_path,
) -> None:
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text("annotation-v1\n", encoding="utf-8")
    fasta = tmp_path / "genome.fa"
    original_fasta = ">chr1\nACGT\n"
    fasta.write_text(original_fasta, encoding="utf-8")
    sj_paths = []
    count_paths = []
    for run in ("SRR35298704", "SRR35298699", "SRR35298715"):
        run_dir = tmp_path / run
        run_dir.mkdir()
        sj_path = run_dir / f"{run}.SJ.out.tab"
        sj_path.write_text(f"{run}-junctions\n", encoding="utf-8")
        sj_paths.append(sj_path)
        count_path = run_dir / f"{run}.ReadsPerGene.out.tab"
        count_path.write_text(f"{run}-counts\n", encoding="utf-8")
        count_paths.append(count_path)

    artifact = attach_input_provenance(
        pd.DataFrame({"exon_key": ["x"]}),
        gtf,
        sj_paths,
        fasta,
        count_paths,
        100,
    )
    validate_input_provenance(
        artifact, gtf, sj_paths, fasta, count_paths, 100
    )

    fasta.write_text(">chr1\nTGCA\n", encoding="utf-8")
    with pytest.raises(FreezeGuardError, match="genome_fasta_sha256"):
        validate_input_provenance(
            artifact, gtf, sj_paths, fasta, count_paths, 100
        )

    fasta.write_text(original_fasta, encoding="utf-8")
    count_paths[1].write_text("changed-counts\n", encoding="utf-8")
    with pytest.raises(FreezeGuardError, match="scrmoe_reads_per_gene_sha256"):
        validate_input_provenance(
            artifact, gtf, sj_paths, fasta, count_paths, 100
        )


def _write_test_bundle(bundle: Path, parquet: bytes, readme: str) -> None:
    bundle.mkdir(parents=True)
    files = {
        "detectable_universe_v50.parquet": parquet,
        "detectable_universe_v50_README.md": readme.encode("utf-8"),
    }
    for name, content in files.items():
        (bundle / name).write_bytes(content)
    manifest = {
        "schema_version": 1,
        "files": {
            name: {
                "sha256": hashlib.sha256(content).hexdigest(),
                "size": len(content),
            }
            for name, content in files.items()
        },
    }
    (bundle / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def test_pointer_publish_crash_window_never_exposes_valid_mixed_bundle(
    tmp_path, monkeypatch
) -> None:
    bundles = tmp_path / "detectable_universe_v50.bundles"
    old_bundle = bundles / "old"
    staged_bundle = bundles / ".new.tmp"
    new_bundle = bundles / "new"
    _write_test_bundle(old_bundle, b"old-parquet", "old-readme\n")
    _write_test_bundle(staged_bundle, b"new-parquet", "new-readme\n")
    (old_bundle / ".done").write_text("old-manifest\n", encoding="utf-8")

    current = tmp_path / "detectable_universe_v50.current"
    current.symlink_to(old_bundle.relative_to(tmp_path), target_is_directory=True)
    compatibility_parquet = tmp_path / "detectable_universe_v50.parquet"
    compatibility_readme = tmp_path / "detectable_universe_v50_README.md"
    compatibility_parquet.write_bytes(b"old-parquet")
    compatibility_readme.write_text("old-readme\n", encoding="utf-8")
    legacy_done = tmp_path / "legacy.done"
    legacy_done.write_text("old\n", encoding="utf-8")
    real_replace = os.replace

    def crash_between_compatibility_updates(source, destination):
        if Path(destination) == compatibility_readme:
            visible = current.resolve(strict=True)
            assert visible == old_bundle
            assert (visible / "detectable_universe_v50.parquet").read_bytes() == b"old-parquet"
            assert (visible / "detectable_universe_v50_README.md").read_text() == "old-readme\n"
            assert compatibility_parquet.read_bytes() == b"old-parquet"
            assert compatibility_readme.read_text() == "old-readme\n"
            assert not (visible / ".done").exists()
            assert not legacy_done.exists()
            raise KeyboardInterrupt("simulated node loss between compatibility updates")
        return real_replace(source, destination)

    monkeypatch.setattr(
        "feng.bench.strata.os.replace", crash_between_compatibility_updates
    )
    with pytest.raises(KeyboardInterrupt, match="simulated node loss"):
        strata_module.publish_staged_bundle(
            staged_bundle, new_bundle, current, legacy_done
        )

    assert current.resolve(strict=True) == old_bundle
    assert (current / "detectable_universe_v50.parquet").read_bytes() == b"old-parquet"
    assert (current / "detectable_universe_v50_README.md").read_text() == "old-readme\n"

    monkeypatch.setattr("feng.bench.strata.os.replace", real_replace)
    strata_module.publish_staged_bundle(
        staged_bundle, new_bundle, current, legacy_done
    )
    assert current.resolve(strict=True) == new_bundle
    assert (current / "detectable_universe_v50.parquet").read_bytes() == b"new-parquet"
    assert (current / "detectable_universe_v50_README.md").read_text() == "new-readme\n"
    assert (current / "manifest.json").is_file()
    assert (current / ".done").is_file()


def test_pointer_publish_validation_failure_keeps_prior_bundle_valid(
    tmp_path,
) -> None:
    bundles = tmp_path / "detectable_universe_v50.bundles"
    old_bundle = bundles / "old"
    staged_bundle = bundles / ".new.tmp"
    new_bundle = bundles / "new"
    _write_test_bundle(old_bundle, b"old-parquet", "old-readme\n")
    _write_test_bundle(staged_bundle, b"new-parquet", "new-readme\n")
    (old_bundle / ".done").write_text("old-manifest\n", encoding="utf-8")
    (staged_bundle / "detectable_universe_v50.parquet").write_bytes(b"tampered")
    current = tmp_path / "detectable_universe_v50.current"
    current.symlink_to(old_bundle.relative_to(tmp_path), target_is_directory=True)
    legacy_done = tmp_path / "legacy.done"
    legacy_done.write_text("old\n", encoding="utf-8")

    with pytest.raises(ValueError, match="manifest .* mismatch"):
        strata_module.publish_staged_bundle(
            staged_bundle, new_bundle, current, legacy_done
        )

    assert current.resolve(strict=True) == old_bundle
    assert (current / ".done").is_file()
    assert legacy_done.is_file()
    assert staged_bundle.is_dir()
    assert not new_bundle.exists()


def test_pointer_publish_rename_then_interrupt_keeps_current_target_live(
    tmp_path, monkeypatch
) -> None:
    bundles = tmp_path / "detectable_universe_v50.bundles"
    old_bundle = bundles / "old"
    staged_bundle = bundles / ".new.tmp"
    new_bundle = bundles / "new"
    _write_test_bundle(old_bundle, b"old-parquet", "old-readme\n")
    _write_test_bundle(staged_bundle, b"new-parquet", "new-readme\n")
    (old_bundle / ".done").write_text("old-manifest\n", encoding="utf-8")

    current = tmp_path / "detectable_universe_v50.current"
    current.symlink_to(old_bundle.relative_to(tmp_path), target_is_directory=True)
    legacy_done = tmp_path / "legacy.done"
    legacy_done.write_text("old\n", encoding="utf-8")
    real_replace = os.replace

    def interrupt_after_pointer_rename(source, destination):
        result = real_replace(source, destination)
        if Path(destination) == current:
            raise KeyboardInterrupt("simulated interrupt after pointer rename")
        return result

    monkeypatch.setattr("feng.bench.strata.os.replace", interrupt_after_pointer_rename)
    with pytest.raises(KeyboardInterrupt, match="after pointer rename"):
        strata_module.publish_staged_bundle(
            staged_bundle, new_bundle, current, legacy_done
        )

    assert current.resolve(strict=True) == new_bundle
    assert new_bundle.is_dir()
    assert (current / ".done").is_file()
    assert not staged_bundle.exists()
    assert not legacy_done.exists()

    monkeypatch.setattr("feng.bench.strata.os.replace", real_replace)
    strata_module.publish_staged_bundle(
        staged_bundle, new_bundle, current, legacy_done
    )
    assert current.resolve(strict=True) == new_bundle


def test_pointer_publish_stage_rename_interrupt_restores_state_and_retries(
    tmp_path, monkeypatch
) -> None:
    bundles = tmp_path / "detectable_universe_v50.bundles"
    old_bundle = bundles / "old"
    staged_bundle = bundles / ".new.tmp"
    new_bundle = bundles / "new"
    _write_test_bundle(old_bundle, b"old-parquet", "old-readme\n")
    _write_test_bundle(staged_bundle, b"new-parquet", "new-readme\n")
    (old_bundle / ".done").write_text("old-manifest\n", encoding="utf-8")

    current = tmp_path / "detectable_universe_v50.current"
    current.symlink_to(old_bundle.relative_to(tmp_path), target_is_directory=True)
    legacy_done = tmp_path / "legacy.done"
    legacy_done.write_text("old\n", encoding="utf-8")
    real_replace = os.replace

    def interrupt_after_stage_rename(source, destination):
        result = real_replace(source, destination)
        if Path(source) == staged_bundle and Path(destination) == new_bundle:
            assert not staged_bundle.exists()
            assert new_bundle.is_dir()
            assert current.resolve(strict=True) == old_bundle
            raise KeyboardInterrupt("simulated interrupt after stage rename")
        return result

    monkeypatch.setattr("feng.bench.strata.os.replace", interrupt_after_stage_rename)
    with pytest.raises(KeyboardInterrupt, match="after stage rename"):
        strata_module.publish_staged_bundle(
            staged_bundle, new_bundle, current, legacy_done
        )

    assert current.resolve(strict=True) == old_bundle
    assert staged_bundle.is_dir()
    assert not new_bundle.exists()
    assert legacy_done.read_text(encoding="utf-8") == "old\n"

    monkeypatch.setattr("feng.bench.strata.os.replace", real_replace)
    strata_module.publish_staged_bundle(staged_bundle, new_bundle, current, legacy_done)
    assert current.resolve(strict=True) == new_bundle
    assert not staged_bundle.exists()
    assert (new_bundle / ".done").is_file()


def test_first_pointer_publish_mid_compat_interrupt_restores_legacy_and_retries(
    tmp_path, monkeypatch
) -> None:
    bundles = tmp_path / "detectable_universe_v50.bundles"
    staged_bundle = bundles / ".new.tmp"
    new_bundle = bundles / "new"
    _write_test_bundle(staged_bundle, b"new-parquet", "new-readme\n")

    current = tmp_path / "detectable_universe_v50.current"
    compatibility_parquet = tmp_path / "detectable_universe_v50.parquet"
    compatibility_readme = tmp_path / "detectable_universe_v50_README.md"
    compatibility_parquet.write_bytes(b"legacy-parquet")
    compatibility_readme.write_text("legacy-readme\n", encoding="utf-8")
    legacy_done = tmp_path / "legacy.done"
    legacy_done.write_text("legacy-valid\n", encoding="utf-8")
    real_replace = os.replace
    interrupted = False

    def interrupt_mid_compatibility(source, destination):
        nonlocal interrupted
        if Path(destination) == compatibility_readme and not interrupted:
            interrupted = True
            assert compatibility_parquet.is_symlink()
            assert not compatibility_parquet.exists()
            assert not compatibility_readme.is_symlink()
            assert not legacy_done.exists()
            assert not current.exists()
            raise KeyboardInterrupt("simulated first-migration compatibility interrupt")
        return real_replace(source, destination)

    monkeypatch.setattr("feng.bench.strata.os.replace", interrupt_mid_compatibility)
    with pytest.raises(KeyboardInterrupt, match="first-migration"):
        strata_module.publish_staged_bundle(
            staged_bundle, new_bundle, current, legacy_done
        )

    assert not current.exists()
    assert not current.is_symlink()
    assert staged_bundle.is_dir()
    assert not new_bundle.exists()
    assert not compatibility_parquet.is_symlink()
    assert compatibility_parquet.read_bytes() == b"legacy-parquet"
    assert not compatibility_readme.is_symlink()
    assert compatibility_readme.read_text(encoding="utf-8") == "legacy-readme\n"
    assert legacy_done.read_text(encoding="utf-8") == "legacy-valid\n"

    monkeypatch.setattr("feng.bench.strata.os.replace", real_replace)
    strata_module.publish_staged_bundle(staged_bundle, new_bundle, current, legacy_done)
    assert current.resolve(strict=True) == new_bundle
    assert compatibility_parquet.is_symlink()
    assert compatibility_parquet.read_bytes() == b"new-parquet"
    assert compatibility_readme.is_symlink()
    assert compatibility_readme.read_text(encoding="utf-8") == "new-readme\n"
    assert not legacy_done.exists()


def test_universe_stage_failed_validation_preserves_published_bundle(
    tmp_path,
) -> None:
    submit_dir = tmp_path / "submit"
    artifact_dir = submit_dir / "feng" / "artifacts" / "bench"
    conditions = submit_dir / "feng" / "artifacts" / "phase3" / "sample_conditions.tsv"
    artifact_dir.mkdir(parents=True)
    conditions.parent.mkdir(parents=True)
    runs = ("SRR35298704", "SRR35298699", "SRR35298715")
    conditions.write_text(
        "run_accession\tcondition\tchemistry\tcontrol_type\n"
        + "".join(f"{run}\tcontrol\tMOE\tscrambled\n" for run in runs),
        encoding="utf-8",
    )

    star_dir = tmp_path / "star"
    for run in runs:
        run_dir = star_dir / run
        run_dir.mkdir(parents=True)
        (run_dir / ".done").touch()
        (run_dir / f"{run}.Aligned.sortedByCoord.out.bam").write_bytes(b"bam")
        (run_dir / f"{run}.SJ.out.tab").write_text("junctions\n", encoding="utf-8")
        (run_dir / f"{run}.ReadsPerGene.out.tab").write_text(
            "counts\n", encoding="utf-8"
        )

    gtf = tmp_path / "annotation.gtf"
    fasta = tmp_path / "genome.fa"
    gtf.write_text("annotation\n", encoding="utf-8")
    fasta.write_text(">chr1\nACGT\n", encoding="utf-8")
    samtools = tmp_path / "samtools"
    samtools.write_text("#!/usr/bin/env bash\nexit 0\n", encoding="utf-8")
    samtools.chmod(0o755)

    out_parquet = artifact_dir / "detectable_universe_v50.parquet"
    out_readme = artifact_dir / "detectable_universe_v50_README.md"
    out_dir = tmp_path / "universe-work"
    out_dir.mkdir()
    legacy_done = out_dir / ".done"
    old_parquet = b"prior-published-parquet"
    old_readme = "prior-published-readme\n"
    out_parquet.write_bytes(old_parquet)
    out_readme.write_text(old_readme, encoding="utf-8")
    legacy_done.write_text("prior-valid-bundle\n", encoding="utf-8")

    rmats_stub = tmp_path / "rmats_stub.py"
    rmats_stub.write_text("# intercepted by python wrapper\n", encoding="utf-8")
    python_wrapper = tmp_path / "python-wrapper"
    python_wrapper.write_text(
        f"#!{sys.executable}\n"
        "import os\n"
        "from pathlib import Path\n"
        "import sys\n"
        f"real_python = {sys.executable!r}\n"
        f"rmats_stub = {str(rmats_stub)!r}\n"
        "args = sys.argv[1:]\n"
        "if args and args[0] == '-':\n"
        "    os.execv(real_python, [real_python, *args])\n"
        "if args and args[0] == rmats_stub:\n"
        "    out_dir = Path(args[args.index('--od') + 1])\n"
        "    out_dir.mkdir(parents=True, exist_ok=True)\n"
        "    (out_dir / 'SE.MATS.JC.txt').write_text('control-only SE\\n')\n"
        "    raise SystemExit(0)\n"
        "if args[:2] == ['-m', 'feng.bench.strata']:\n"
        "    Path(args[args.index('--out-parquet') + 1]).write_bytes(b'invalid parquet')\n"
        "    Path(args[args.index('--out-readme') + 1]).write_text('invalid readme\\n')\n"
        "    raise SystemExit(0)\n"
        "raise SystemExit('unexpected python-wrapper invocation: ' + repr(args))\n",
        encoding="utf-8",
    )
    python_wrapper.chmod(0o755)

    result = subprocess.run(
        ["bash", str(SBATCH)],
        cwd=REPO_ROOT,
        env={
            **os.environ,
            "SLURM_SUBMIT_DIR": str(submit_dir),
            "FENG_BENCH_CONDITIONS": str(conditions),
            "FENG_BENCH_STAR_DIR": str(star_dir),
            "FENG_BENCH_OUT_DIR": str(out_dir),
            "FENG_BENCH_GTF": str(gtf),
            "FENG_BENCH_FASTA": str(fasta),
            "FENG_BENCH_SAMTOOLS": str(samtools),
            "FENG_BENCH_PYTHON": str(python_wrapper),
            "FENG_BENCH_RMATS": str(rmats_stub),
            "PYTHONPATH": str(REPO_ROOT / "feng" / "src"),
        },
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "staged universe failed validation" in result.stderr
    assert out_parquet.read_bytes() == old_parquet
    assert out_readme.read_text(encoding="utf-8") == old_readme
    assert legacy_done.read_text(encoding="utf-8") == "prior-valid-bundle\n"
    assert not (artifact_dir / "detectable_universe_v50.current").exists()
    assert not list(artifact_dir.rglob("*.tmp"))


def test_annotation_not_rmats_seeds_constitutive_detectable_candidate(tmp_path) -> None:
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text(
        "chr1\tT\texon\t1\t50\t.\t+\t.\tgene_id \"ENSG1.1\"; transcript_id \"ENST1.1\"; transcript_type \"protein_coding\"; gene_name \"GENE1\"; exon_id \"E1\";\n"
        "chr1\tT\texon\t101\t150\t.\t+\t.\tgene_id \"ENSG1.1\"; transcript_id \"ENST1.1\"; transcript_type \"protein_coding\"; gene_name \"GENE1\"; exon_id \"E2\";\n"
        "chr1\tT\texon\t201\t250\t.\t+\t.\tgene_id \"ENSG1.1\"; transcript_id \"ENST1.1\"; transcript_type \"protein_coding\"; gene_name \"GENE1\"; exon_id \"E3\";\n",
        encoding="utf-8",
    )
    annotated = load_annotated_internal_exons(gtf)
    assert annotated["exon_key"].tolist() == ["chr1:100-150:+"]
    samples = [
        {("chr1", 51, 100, "+"): left, ("chr1", 151, 200, "+"): right}
        for left, right in [(12, 11), (15, 10), (20, 10)]
    ]
    measured = attach_annotated_junction_metrics(annotated, samples)
    detected = classify_detectable(measured)
    assert detected["mean_total_junction_reads"].item() == pytest.approx(26.0)
    assert detected["n_samples_with_psi"].item() == 3
    assert detected["baseline_psi"].item() == pytest.approx(1.0)
    assert bool(detected["detectable"].item())


def test_star_junction_estimator_matches_rmats_with_skipping_and_asymmetric_flanks(
    tmp_path,
) -> None:
    exon = pd.DataFrame([{
        "event_id": "annotation:ENST1:E2", "gene_id": "ENSG1.1",
        "gene_id_base": "ENSG1", "gene_symbol": "GENE1", "chrom": "chr1",
        "strand": "+", "exon_start_0base": 100, "exon_end": 150,
        "upstream_ee": 50, "downstream_es": 200, "exon_length": 50,
        "exon_key": "chr1:100-150:+",
    }])
    flank_counts = [(20, 10, 5), (30, 6, 11), (0, 0, 0)]
    junctions = [
        {("chr1", 51, 100, "+"): left,
         ("chr1", 151, 200, "+"): right,
         ("chr1", 51, 200, "+"): skipped}
        for left, right, skipped in flank_counts
    ]
    annotated = attach_annotated_junction_metrics(exon, junctions, read_length=100).iloc[0]
    inc = [left + right for left, right, _ in flank_counts]
    skip = [value for _, _, value in flank_counts]
    inc_len, skip_len = 149, 99
    psi = [
        (i / inc_len) / ((i / inc_len) + (s / skip_len))
        for i, s in zip(inc[:2], skip[:2])
    ]
    se = tmp_path / "SE.MATS.JC.txt"
    header = [
        "ID", "GeneID", "geneSymbol", "chr", "strand", "exonStart_0base", "exonEnd",
        "upstreamES", "upstreamEE", "downstreamES", "downstreamEE",
        "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
        "IncFormLen", "SkipFormLen", "PValue", "FDR",
        "IncLevel1", "IncLevel2", "IncLevelDifference",
    ]
    row = [
        "1", '"ENSG1.1"', '"GENE1"', "chr1", "+", "100", "150",
        "0", "50", "200", "250", ",".join(map(str, inc[:2])),
        ",".join(map(str, skip[:2])), str(inc[2]), str(skip[2]), str(inc_len),
        str(skip_len), "NA", "NA", ",".join(map(str, psi)), "NA", "NA",
    ]
    se.write_text("\t".join(header) + "\n" + "\t".join(row) + "\n", encoding="utf-8")
    rmats = load_se_events(se).iloc[0]
    assert annotated["mean_inclusion_junction_reads"] == pytest.approx(
        rmats["mean_inclusion_junction_reads"]
    )
    assert annotated["mean_skipping_junction_reads"] == pytest.approx(
        rmats["mean_skipping_junction_reads"]
    )
    assert annotated["mean_total_junction_reads"] == pytest.approx(
        rmats["mean_total_junction_reads"]
    )
    assert annotated["n_samples_with_psi"] == rmats["n_samples_with_psi"] == 2
    assert annotated["baseline_psi"] == pytest.approx(rmats["baseline_psi"])


def test_annotated_candidates_exclude_retained_intron_transcripts(tmp_path) -> None:
    assert "protein_coding" in ALLOWED_GENCODE_TRANSCRIPT_TYPES
    gtf = tmp_path / "classes.gtf"
    lines = []
    for transcript, transcript_type, tag, offset in [
        ("ENST_CODING", "protein_coding", "basic", 0),
        ("ENST_RI_TYPE", "retained_intron", "basic", 1000),
        ("ENST_RI_TAG", "protein_coding", "retained_intron_CDS", 2000),
    ]:
        for number, (start, end) in enumerate([(1, 50), (101, 150), (201, 250)], 1):
            attrs = (
                f'gene_id "ENSG{offset}"; transcript_id "{transcript}"; '
                f'transcript_type "{transcript_type}"; gene_name "G{offset}"; '
                f'exon_id "{transcript}_E{number}"; tag "{tag}";'
            )
            lines.append(f"chr1\tT\texon\t{start + offset}\t{end + offset}\t.\t+\t.\t{attrs}")
    gtf.write_text("\n".join(lines) + "\n", encoding="utf-8")
    candidates = load_annotated_internal_exons(gtf)
    assert candidates["event_id"].tolist() == ["annotation:ENST_CODING:ENST_CODING_E2"]


def test_generated_documentation_describes_transcript_position_filter(tmp_path) -> None:
    se = tmp_path / "SE.MATS.JC.txt"
    fasta = tmp_path / "genome.fa"
    se.write_text("source\n", encoding="utf-8")
    fasta.write_text(">chr1\nA\n", encoding="utf-8")
    events = pd.DataFrame([{
        "detectable": True, "chrom": "chr1", "baseline_psi": 0.5,
        "gene_id_base": "ENSG1", "junction_read_length": 100,
    }])
    readme = render_readme(
        events, select_control_runs(_conditions()), se,
        DetectabilityThresholds(), fasta,
    )
    assert "Sorting each transcript" in readme
    assert "retaining only positions between its first and" in readme
    assert "an SE event always has a defined flanking exon" not in readme
    assert (
        f"Reference FASTA SHA-256: `{hashlib.sha256(fasta.read_bytes()).hexdigest()}`"
        in readme
    )


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
    assert first["n_annotation_events"] == 2
    assert first["annotation_event_ids"] == "1,2"
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
    controls = ["SRR35298704", "SRR35298699", "SRR35298715"]
    digest = hashlib.sha256(("\n".join(controls) + "\n").encode()).hexdigest()
    assert f'PINNED_CONTROL_SHA256="{digest}"' in script
    assert script.index("actual_allowlist_hash=") < script.index('"${PYTHON_BIN}" "${RMATS}"')
    assert "duplicate accession in pinned control allowlist" in script
    assert "duplicate control accession in conditions table" in script
    assert script.index("validate_published_artifact()") < script.index("DONE_FILE=")
    assert "duplicate exon_key values" in script
    assert "null gc_fraction values" in script
    assert "provenance SHA-256 mismatch" in script
    assert "validate_input_provenance" in script
    assert '"genome_fasta_sha256"' in script
    assert '"scrmoe_reads_per_gene_sha256"' in script
    assert 'f"{run}.ReadsPerGene.out.tab"' in script
    # The .done fast-path must re-check the freeze/labelling guarantees, not
    # just schema, so a stale sentinel cannot skip the applied fixes.
    assert "forbidden contrast fields" in script
    assert "does not match parquet" in script
    assert "rMATS b1/b2 accessions differ from pinned ScrMOE allowlist" in script
    assert "published parquet is still gated by rMATS SE enumeration" in script
    assert '--gtf "${GTF}"' in script
    assert '--read-length "${READ_LENGTH}"' in script
    assert 'CURRENT_POINTER="${ARTIFACT_DIR}/detectable_universe_v50.current"' in script
    assert 'STAGED_BUNDLE="${BUNDLE_ROOT}/.${BUILD_ID}.tmp"' in script
    assert 'VERSIONED_BUNDLE="${BUNDLE_ROOT}/${BUILD_ID}"' in script
    assert "write_bundle_manifest" in script
    assert "validate_bundle_manifest" in script
    assert 'validate_published_artifact "${STAGED_PARQUET}" "${STAGED_README}"' in script
    assert "publish_staged_bundle" in script


def test_detectability_policy_records_migration_rationale() -> None:
    policy = __import__("json").loads(DETECTABILITY_POLICY)
    assert policy["policy_id"] == DETECTABILITY_POLICY_ID
    assert policy["min_samples_with_psi"] == 2
    assert policy["migration_from"] == "pooled-11-control minimum 6 of 11"
    assert "strict majority" in policy["rationale"]


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
    with pytest.raises(ValueError, match="expected 3 pooled control samples"):
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
