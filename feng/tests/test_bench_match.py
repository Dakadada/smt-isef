"""Matching + kill-gate + control-only tests for the Arm D benchmark strata."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from feng.bench.match import (
    COVARIATES,
    COVARIATE_TRANSFORMS,
    LINEAR,
    MatchConfig,
    build_strata,
    check_kill_gates,
    map_positives_to_universe,
)
from feng.bench.strata import (
    FreezeGuardError,
    INPUT_PROVENANCE_COLUMNS,
    UNIVERSE_PROVENANCE,
    build_universe,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
UNIVERSE = REPO_ROOT / "feng" / "artifacts" / "bench" / "detectable_universe_v50.parquet"
POSITIVES = REPO_ROOT / "feng" / "artifacts" / "bench" / "positives_v50.tsv"
MATCH_SCRIPT = REPO_ROOT / "feng" / "scripts" / "bench_match.sbatch"
TEST_INPUT_PROVENANCE = {
    "gencode_gtf_sha256": "1" * 64,
    "scrmoe_sj_sha256": json.dumps({"run": "2" * 64}),
    "genome_fasta_sha256": "3" * 64,
    "scrmoe_reads_per_gene_sha256": json.dumps({"run": "4" * 64}),
    "junction_read_length": 100,
    "junction_read_length_sha256": hashlib.sha256(b"100\n").hexdigest(),
}


def _universe_row(exon_key: str, chrom: str, start0: int, end: int, strand: str, **cov) -> dict:
    row = {
        "exon_key": exon_key,
        "chrom": chrom,
        "exon_start_0base": start0,
        "exon_end": end,
        "strand": strand,
        "detectable": cov.pop("detectable", True),
        **UNIVERSE_PROVENANCE,
        **TEST_INPUT_PROVENANCE,
    }
    # default covariate values, overridable
    defaults = {
        "baseline_psi": 0.5,
        "mean_total_junction_reads": 50.0,
        "control_mean_cpm": 10.0,
        "exon_length": end - start0,
        "gc_fraction": 0.5,
    }
    defaults.update(cov)
    row.update(defaults)
    return row


def _synthetic_universe(n_close: int, n_far: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Universe with one positive, ``n_close`` matchable + ``n_far`` unmatchable backgrounds."""
    rows = []
    # positive at chr1:100-200:+ -> exon_key chr1:99-200:+ (1-based 100..200)
    pos_key = "chr1:99-200:+"
    rows.append(_universe_row(pos_key, "chr1", 99, 200, "+"))
    # close backgrounds: identical covariates to the positive -> standardised dist 0
    for i in range(n_close):
        s0 = 1000 + i * 500
        rows.append(_universe_row(f"chr1:{s0}-{s0 + 101}:+", "chr1", s0, s0 + 101, "+"))
    # far backgrounds: extreme covariates -> outside the caliper
    for i in range(n_far):
        s0 = 500000 + i * 500
        rows.append(
            _universe_row(
                f"chr1:{s0}-{s0 + 300}:+", "chr1", s0, s0 + 300, "+",
                baseline_psi=0.99, mean_total_junction_reads=5000.0,
                control_mean_cpm=5000.0, exon_length=300, gc_fraction=0.9,
            )
        )
    universe = pd.DataFrame(rows)
    positives = pd.DataFrame(
        [{"event_id": "POS1", "chrom": "chr1", "start": 100, "end": 200, "strand": "+"}]
    )
    return universe, positives


def test_map_positives_to_universe_flags_presence() -> None:
    universe, positives = _synthetic_universe(5, 0)
    mapping = map_positives_to_universe(universe, positives)
    assert mapping.loc[0, "exon_key"] == "chr1:99-200:+"
    assert bool(mapping.loc[0, "in_universe"])
    assert bool(mapping.loc[0, "detectable"])


def test_map_positive_absent_from_universe() -> None:
    universe, positives = _synthetic_universe(5, 0)
    positives.loc[0, "end"] = 250  # shift so the exon_key no longer exists
    mapping = map_positives_to_universe(universe, positives)
    assert not bool(mapping.loc[0, "in_universe"])
    assert not bool(mapping.loc[0, "detectable"])


def test_duplicate_positive_event_id_is_rejected() -> None:
    universe, positives = _synthetic_universe(5, 0)
    positives = pd.concat([positives, positives], ignore_index=True)
    with pytest.raises(FreezeGuardError, match="event_id values must be unique"):
        map_positives_to_universe(universe, positives)


def test_build_strata_excludes_positive_and_selects_close_backgrounds() -> None:
    universe, positives = _synthetic_universe(n_close=60, n_far=40)
    strata, summary = build_strata(universe, positives, MatchConfig())
    # positive appears once as a flagged anchor, never as a background row
    backgrounds = strata.loc[~strata["is_positive"]]
    assert "chr1:99-200:+" not in set(backgrounds["exon_key"])
    anchors = strata.loc[strata["is_positive"]]
    assert set(anchors["exon_key"]) == {"chr1:99-200:+"}
    # only the close group is eligible; selection capped at n_background
    row = summary.iloc[0]
    assert int(row["eligible_background"]) == 60
    assert int(row["n_selected"]) == 60  # fewer than n_background=100 available
    assert bool(row["kill_gate_pass"])
    # every selected background is from the close group (exon_length 101), never
    # the far group (exon_length 300) -- eligible_background == 60 already proves
    # only the close group cleared the caliper.
    assert (strata["exon_length"] == 101).all()


def test_kill_gate_fails_below_floor() -> None:
    universe, positives = _synthetic_universe(n_close=10, n_far=40)
    strata, summary = build_strata(universe, positives, MatchConfig())
    row = summary.iloc[0]
    assert int(row["eligible_background"]) == 10
    assert not bool(row["kill_gate_pass"])
    assert not check_kill_gates(summary)
    assert "10 eligible" in row["kill_reason"]


def test_undetectable_positive_is_killed_with_no_stratum() -> None:
    # n_far>0 so the background pool has real covariate spread; the positive
    # itself is flipped undetectable below.
    universe, positives = _synthetic_universe(n_close=60, n_far=10)
    universe.loc[universe["exon_key"] == "chr1:99-200:+", "detectable"] = False
    strata, summary = build_strata(universe, positives, MatchConfig())
    row = summary.iloc[0]
    assert not bool(row["kill_gate_pass"])
    assert int(row["n_selected"]) == 0
    assert strata.empty


def test_contrast_field_in_universe_is_refused() -> None:
    universe, positives = _synthetic_universe(n_close=60, n_far=0)
    universe["IncLevelDifference"] = 0.0
    with pytest.raises(FreezeGuardError):
        build_strata(universe, positives, MatchConfig())


def test_matching_rejects_unbound_or_wrong_control_provenance() -> None:
    universe, positives = _synthetic_universe(n_close=5, n_far=1)
    with pytest.raises(FreezeGuardError, match="ScrMOE provenance"):
        build_strata(universe.drop(columns=["control_allowlist_sha256"]), positives)
    universe.loc[0, "control_run_ids"] = "pooled-controls"
    with pytest.raises(FreezeGuardError, match="invalid ScrMOE provenance"):
        build_strata(universe, positives)


def test_selection_is_deterministic() -> None:
    universe, positives = _synthetic_universe(n_close=60, n_far=40)
    a, _ = build_strata(universe, positives, MatchConfig(seed=1))
    b, _ = build_strata(universe, positives, MatchConfig(seed=999))
    # selection does not depend on the seed (deterministic nearest + tiebreak)
    assert list(a["exon_key"]) == list(b["exon_key"])


def test_build_strata_carries_exact_universe_input_provenance() -> None:
    universe, positives = _synthetic_universe(60, 10)
    universe["gencode_gtf_sha256"] = "a" * 64
    universe["genome_fasta_sha256"] = "b" * 64
    strata, _ = build_strata(universe, positives, MatchConfig())
    for column in INPUT_PROVENANCE_COLUMNS:
        assert column in strata.columns
        assert strata[column].eq(universe[column].iloc[0]).all()


def _input_provenance(token: str) -> dict:
    def digest(label: str) -> str:
        return hashlib.sha256(f"{token}:{label}".encode()).hexdigest()

    return {
        "gencode_gtf_sha256": digest("gtf"),
        "scrmoe_sj_sha256": json.dumps({"run": digest("sj")}),
        "genome_fasta_sha256": digest("fasta"),
        "scrmoe_reads_per_gene_sha256": json.dumps({"run": digest("counts")}),
        "junction_read_length": 100,
        "junction_read_length_sha256": hashlib.sha256(b"100\n").hexdigest(),
    }


def _write_controlled_match_universe(
    path: Path, positives: pd.DataFrame, token: str
) -> dict:
    provenance = _input_provenance(token)
    rows = []
    for _, positive in positives.reset_index(drop=True).iterrows():
        start0 = int(positive["start"]) - 1
        rows.append(_universe_row(
            f'{positive["chrom"]}:{start0}-{positive["end"]}:{positive["strand"]}',
            str(positive["chrom"]), start0, int(positive["end"]),
            str(positive["strand"]), **provenance,
        ))
    for index in range(800):
        start0 = 10_000_000 + index * 200
        rows.append(_universe_row(
            f"chrB:{start0}-{start0 + 100}:+", "chrB", start0, start0 + 100, "+",
            baseline_psi=0.2 + (index % 60) / 100,
            mean_total_junction_reads=50.0 + index,
            control_mean_cpm=10.0 + index / 10,
            exon_length=100,
            gc_fraction=0.3 + (index % 40) / 100,
            **provenance,
        ))
    pd.DataFrame(rows).to_parquet(path, index=False)
    return provenance


def _write_universe_bundle_metadata(bundle: Path, readme_text: str) -> None:
    readme = bundle / "detectable_universe_v50_README.md"
    readme.write_text(readme_text, encoding="utf-8")
    universe = bundle / "detectable_universe_v50.parquet"
    files = {universe.name: universe.read_bytes(), readme.name: readme.read_bytes()}
    manifest = bundle / "manifest.json"
    manifest.write_text(json.dumps({
        "schema_version": 1,
        "files": {
            name: {"sha256": hashlib.sha256(content).hexdigest(), "size": len(content)}
            for name, content in files.items()
        },
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (bundle / ".done").write_text(
        hashlib.sha256(manifest.read_bytes()).hexdigest() + "\n", encoding="ascii"
    )


def _write_match_python_stub(path: Path) -> None:
    path.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import os\n"
        "import sys\n"
        "if sys.argv[1:3] != ['-m', 'feng.bench.match']:\n"
        "    os.execv(sys.executable, [sys.executable, *sys.argv[1:]])\n"
        "import pandas as pd\n"
        "from feng.bench.match import COORD_COLUMNS, COVARIATES\n"
        "from feng.bench.strata import INPUT_PROVENANCE_COLUMNS, UNIVERSE_PROVENANCE_COLUMNS\n"
        "args = sys.argv[3:]\n"
        "def value(flag):\n"
        "    return Path(args[args.index(flag) + 1])\n"
        "switch_to = os.environ.get('MATCH_SWITCH_POINTER_TO')\n"
        "if switch_to:\n"
        "    pointer = Path(os.environ['MATCH_CURRENT_POINTER'])\n"
        "    temporary = pointer.with_name(pointer.name + '.switch')\n"
        "    temporary.unlink(missing_ok=True)\n"
        "    temporary.symlink_to(Path(switch_to), target_is_directory=True)\n"
        "    os.replace(temporary, pointer)\n"
        "positives = pd.read_csv(value('--positives'), sep='\\t')\n"
        "universe_path = value('--universe')\n"
        "universe = pd.read_parquet(universe_path)\n"
        "with Path(os.environ['MATCH_STUB_LOG']).open('a', encoding='utf-8') as handle:\n"
        "    handle.write(__import__('hashlib').sha256(universe_path.read_bytes()).hexdigest() + '\\n')\n"
        "positive_keys = [f'{row.chrom}:{int(row.start) - 1}-{int(row.end)}:{row.strand}' for row in positives.itertuples()]\n"
        "backgrounds = universe.loc[~universe['exon_key'].isin(positive_keys)].reset_index(drop=True)\n"
        "carry = ['exon_key', *COORD_COLUMNS, *COVARIATES, *UNIVERSE_PROVENANCE_COLUMNS, *INPUT_PROVENANCE_COLUMNS]\n"
        "rows = []\n"
        "summary = []\n"
        "mode = os.environ['MATCH_STUB_MODE']\n"
        "for event_index, positive in enumerate(positives.itertuples(index=False)):\n"
        "    event_id = str(positive.event_id)\n"
        "    positive_key = positive_keys[event_index]\n"
        "    anchor = universe.loc[universe['exon_key'].eq(positive_key)].iloc[0]\n"
        "    selected = backgrounds.iloc[event_index * 100:(event_index + 1) * 100]\n"
        "    if mode == 'short' and event_index == 0:\n"
        "        selected = selected.iloc[:-1]\n"
        "    anchor_row = {column: anchor[column] for column in carry}\n"
        "    anchor_row.update({'positive_event_id': event_id, 'positive_exon_key': positive_key, 'is_positive': True, 'match_distance': 0.0})\n"
        "    rows.append(anchor_row)\n"
        "    for distance, (_, background) in enumerate(selected.iterrows(), 1):\n"
        "        row = {column: background[column] for column in carry}\n"
        "        row.update({'positive_event_id': event_id, 'positive_exon_key': positive_key, 'is_positive': False, 'match_distance': distance / 100.0})\n"
        "        rows.append(row)\n"
        "    summary.append({\n"
        "        'event_id': event_id,\n"
        "        'positive_exon_key': positive_key,\n"
        "        'in_universe': True,\n"
        "        'detectable': True,\n"
        "        'eligible_background': 100.5 if (mode == 'fractional' and event_index == 0) else 50 if (mode == 'impossible' and event_index == 0) else 100,\n"
        "        'n_selected': 100,\n"
        "        'kill_gate_pass': True,\n"
        "        'kill_reason': '',\n"
        "    })\n"
        "strata = pd.DataFrame(rows)\n"
        "if mode == 'nullable':\n"
        "    extra = strata.iloc[1].copy()\n"
        "    extra['is_positive'] = pd.NA\n"
        "    strata = pd.concat([strata, extra.to_frame().T], ignore_index=True)\n"
        "    strata['is_positive'] = strata['is_positive'].astype('boolean')\n"
        "strata.to_parquet(value('--out-strata'), index=False)\n"
        "pd.DataFrame(summary).to_csv(value('--out-summary'), sep='\\t', index=False)\n",
        encoding="utf-8",
    )
    path.chmod(0o755)


def _match_stage_fixture(tmp_path: Path) -> dict:
    submit_dir = tmp_path / "submit"
    artifact_dir = submit_dir / "feng" / "artifacts" / "bench"
    artifact_dir.mkdir(parents=True)
    bundle = artifact_dir / "detectable_universe_v50.bundles" / "fixture"
    bundle.mkdir(parents=True)
    universe = bundle / "detectable_universe_v50.parquet"
    positives_path = artifact_dir / "positives_v50.tsv"
    positives_path.write_bytes(POSITIVES.read_bytes())
    positives = pd.read_csv(positives_path, sep="\t")
    provenance = _write_controlled_match_universe(universe, positives, "initial")
    _write_universe_bundle_metadata(bundle, "fixture universe\n")
    current = artifact_dir / "detectable_universe_v50.current"
    current.symlink_to(bundle.relative_to(artifact_dir), target_is_directory=True)
    universe_dir = tmp_path / "scratch" / "feng" / "bench" / "detectable_universe"
    universe_dir.mkdir(parents=True)
    match_dir = tmp_path / "scratch" / "feng" / "bench" / "matched_strata"
    match_dir.mkdir(parents=True)
    python_stub = tmp_path / "python-stub"
    _write_match_python_stub(python_stub)
    stub_log = tmp_path / "matcher.log"
    out_strata = artifact_dir / "strata_v50.parquet"
    out_summary = artifact_dir / "strata_summary_v50.tsv"
    env = {
        **os.environ,
        "SLURM_SUBMIT_DIR": str(submit_dir),
        "SCRATCH": str(tmp_path / "scratch"),
        "FENG_BENCH_PYTHON": str(python_stub),
        "FENG_BENCH_MATCH_SOURCE": str(REPO_ROOT / "feng" / "src" / "feng" / "bench" / "match.py"),
        "FENG_BENCH_ARTIFACT_DIR": str(artifact_dir),
        "FENG_BENCH_POSITIVES": str(positives_path),
        "FENG_BENCH_UNIVERSE_DIR": str(universe_dir),
        "FENG_BENCH_MATCH_DIR": str(match_dir),
        "FENG_BENCH_STRATA": str(out_strata),
        "FENG_BENCH_SUMMARY": str(out_summary),
        "MATCH_STUB_LOG": str(stub_log),
    }
    return {
        "universe": universe,
        "bundle": bundle,
        "current": current,
        "positives_path": positives_path,
        "positives": positives,
        "provenance": provenance,
        "done": match_dir / ".done",
        "out_strata": out_strata,
        "out_summary": out_summary,
        "stub_log": stub_log,
        "env": env,
    }


def _run_match_stage(fixture: dict, mode: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["bash", str(MATCH_SCRIPT)], cwd=REPO_ROOT,
        env={**fixture["env"], "MATCH_STUB_MODE": mode},
        capture_output=True, text=True,
    )


def test_match_stage_validates_before_atomic_publish_and_done(tmp_path: Path) -> None:
    fixture = _match_stage_fixture(tmp_path)
    out_strata = fixture["out_strata"]
    out_summary = fixture["out_summary"]
    done = fixture["done"]
    old_strata = b"prior-published-strata"
    old_summary = "prior-published-summary\n"
    out_strata.write_bytes(old_strata)
    out_summary.write_text(old_summary, encoding="utf-8")
    done.touch()
    failed = _run_match_stage(fixture, "short")
    assert failed.returncode != 0
    assert "staged matched strata failed validation" in failed.stderr
    assert out_strata.read_bytes() == old_strata
    assert out_summary.read_text(encoding="utf-8") == old_summary
    assert not done.exists()

    passed = _run_match_stage(fixture, "valid")
    assert passed.returncode == 0, passed.stderr
    assert done.is_file()
    published_strata = pd.read_parquet(out_strata)
    published_summary = pd.read_csv(out_summary, sep="\t")
    assert len(published_strata) == 8 * 101
    assert published_summary["event_id"].tolist() == fixture["positives"]["event_id"].tolist()
    for column in INPUT_PROVENANCE_COLUMNS:
        assert published_strata[column].eq(fixture["provenance"][column]).all()


def test_match_stage_default_universe_reads_atomic_current_pointer(
    tmp_path: Path,
) -> None:
    fixture = _match_stage_fixture(tmp_path)
    result = _run_match_stage(fixture, "valid")
    assert result.returncode == 0, result.stderr
    assert fixture["done"].is_file()


def test_match_stage_rejects_legacy_universe_override_bypassing_current(
    tmp_path: Path,
) -> None:
    fixture = _match_stage_fixture(tmp_path)
    legacy = tmp_path / "legacy-universe.parquet"
    legacy.write_bytes(fixture["universe"].read_bytes())
    legacy_done = tmp_path / "legacy.done"
    legacy_done.write_text("legacy-valid\n", encoding="utf-8")
    fixture["env"]["FENG_BENCH_UNIVERSE"] = str(legacy)
    fixture["env"]["FENG_BENCH_UNIVERSE_DONE"] = str(legacy_done)

    result = _run_match_stage(fixture, "valid")

    assert result.returncode != 0
    assert "FENG_BENCH_UNIVERSE" in result.stderr
    assert not fixture["stub_log"].exists()


def test_match_stage_pins_one_bundle_across_concurrent_pointer_switch(
    tmp_path: Path,
) -> None:
    fixture = _match_stage_fixture(tmp_path)
    old_provenance = fixture["provenance"]
    new_bundle = fixture["bundle"].parent / "replacement"
    new_bundle.mkdir()
    new_universe = new_bundle / fixture["universe"].name
    _write_controlled_match_universe(new_universe, fixture["positives"], "replacement")
    _write_universe_bundle_metadata(new_bundle, "replacement universe\n")
    fixture["env"]["MATCH_SWITCH_POINTER_TO"] = str(new_bundle)
    fixture["env"]["MATCH_CURRENT_POINTER"] = str(fixture["current"])

    result = _run_match_stage(fixture, "valid")

    assert result.returncode == 0, result.stderr
    assert fixture["current"].resolve(strict=True) == new_bundle
    published = pd.read_parquet(fixture["out_strata"])
    for column in INPUT_PROVENANCE_COLUMNS:
        assert published[column].eq(old_provenance[column]).all()


@pytest.mark.parametrize("mode", ["nullable", "fractional", "impossible"])
def test_match_stage_rejects_malformed_roles_and_counts(tmp_path: Path, mode: str) -> None:
    fixture = _match_stage_fixture(tmp_path)
    result = _run_match_stage(fixture, mode)
    assert result.returncode != 0
    assert "staged matched strata failed validation" in result.stderr
    assert not fixture["out_strata"].exists()
    assert not fixture["out_summary"].exists()
    assert not fixture["done"].exists()


@pytest.mark.parametrize("changed_input", ["universe", "positives"])
def test_match_stage_rebuilds_valid_done_when_inputs_change(
    tmp_path: Path, changed_input: str
) -> None:
    fixture = _match_stage_fixture(tmp_path)
    first = _run_match_stage(fixture, "valid")
    assert first.returncode == 0, first.stderr
    first_done = fixture["done"].read_bytes()
    if changed_input == "universe":
        replacement = fixture["bundle"].parent / "rebuilt"
        replacement.mkdir()
        rebuilt_universe = replacement / fixture["universe"].name
        expected_provenance = _write_controlled_match_universe(
            rebuilt_universe, fixture["positives"], "rebuilt"
        )
        _write_universe_bundle_metadata(replacement, "rebuilt universe\n")
        replacement_pointer = fixture["current"].with_name(
            fixture["current"].name + ".replacement"
        )
        replacement_pointer.symlink_to(replacement, target_is_directory=True)
        os.replace(replacement_pointer, fixture["current"])
    else:
        changed_positives = fixture["positives"].iloc[::-1].reset_index(drop=True)
        changed_positives.to_csv(fixture["positives_path"], sep="\t", index=False)
        expected_provenance = fixture["provenance"]

    second = _run_match_stage(fixture, "valid")
    assert second.returncode == 0, second.stderr
    assert "SKIP:" not in second.stdout
    assert len(fixture["stub_log"].read_text(encoding="utf-8").splitlines()) == 2
    assert fixture["done"].read_bytes() != first_done
    published = pd.read_parquet(fixture["out_strata"])
    assert published["gencode_gtf_sha256"].eq(expected_provenance["gencode_gtf_sha256"]).all()
    if changed_input == "positives":
        summary = pd.read_csv(fixture["out_summary"], sep="\t")
        assert summary["event_id"].tolist() == fixture["positives"]["event_id"].iloc[::-1].tolist()


def test_widened_heterogeneous_universe_preserves_six_positive_match_quality() -> None:
    positives = pd.read_csv(POSITIVES, sep="\t")
    profiles = {
        "OT01_POLR2H_e2": (0.31, 29.0, 3.2, 84, 0.41),
        "OT02_PITHD1_e2": (0.54, 43.0, 8.5, 44, 0.48),
        "OT03_SERPINB7_e4": (0.96, 61.0, 54.0, 117, 0.63),
        "OT04_RTTN_e4": (0.93, 62.0, 126.0, 90, 0.57),
        "OT05_REV3L_e9": (0.67, 88.0, 19.0, 85, 0.45),
        "OT06_PRKRA_e2": (0.42, 137.0, 37.0, 170, 0.52),
        "OT07_GOLGA4_e4": (0.78, 72.0, 14.0, 66, 0.36),
        "OT08_PAK1_e2": (0.59, 211.0, 82.0, 211, 0.60),
    }
    positive_rows = []
    for positive in positives.itertuples(index=False):
        start0 = int(positive.start) - 1
        positive_rows.append(_universe_row(
            f"{positive.chrom}:{start0}-{positive.end}:{positive.strand}",
            positive.chrom, start0, int(positive.end), positive.strand,
            **dict(zip(COVARIATES, profiles[positive.event_id])),
        ))

    background_rows = []
    for profile_index, (event_id, profile) in enumerate(profiles.items()):
        anchor = dict(zip(COVARIATES, profile))
        for index in range(80):
            centred = (index - 39.5) / 39.5
            covariates = {
                "baseline_psi": float(np.clip(
                    anchor["baseline_psi"] + 0.055 * centred
                    + 0.004 * ((index % 5) - 2), 0.01, 0.99,
                )),
                "mean_total_junction_reads": anchor["mean_total_junction_reads"]
                * float(np.exp(0.22 * centred + 0.018 * ((index % 7) - 3))),
                "control_mean_cpm": anchor["control_mean_cpm"]
                * float(np.exp(0.26 * centred + 0.025 * ((index % 6) - 2.5))),
                "exon_length": max(20, round(
                    anchor["exon_length"]
                    * float(np.exp(0.13 * centred + 0.012 * ((index % 4) - 1.5)))
                )),
                "gc_fraction": float(np.clip(
                    anchor["gc_fraction"] + 0.045 * centred
                    + 0.004 * ((index % 3) - 1), 0.05, 0.95,
                )),
            }
            start0 = 10_000_000 + profile_index * 1_000_000 + index * 1_000
            background_rows.append(_universe_row(
                f"chrB:{start0}-{start0 + covariates['exon_length']}:+",
                "chrB", start0, start0 + covariates["exon_length"], "+",
                **covariates,
            ))
    # Added annotated-exon candidates span the full covariate range and alter
    # every centre/scale used by matching; they are not mirror copies.
    for index in range(120):
        start0 = 30_000_000 + index * 2_000
        background_rows.append(_universe_row(
            f"chrW:{start0}-{start0 + 35 + index}:+",
            "chrW", start0, start0 + 35 + index, "+",
            baseline_psi=0.08 + 0.84 * ((index * 17) % 119) / 118,
            mean_total_junction_reads=float(np.exp(2.9 + index / 42)),
            control_mean_cpm=float(np.exp(0.4 + index / 28)),
            exon_length=35 + index,
            gc_fraction=0.23 + 0.54 * ((index * 29) % 119) / 118,
        ))
    widened = pd.DataFrame([*positive_rows, *background_rows])
    config = MatchConfig(n_background=60, min_background=50)
    strata, summary = build_strata(widened, positives, config)

    prior_ids = set(profiles).difference({"OT03_SERPINB7_e4", "OT04_RTTN_e4"})
    prior_summary = summary.loc[summary["event_id"].isin(prior_ids)]
    assert len(prior_summary) == 6
    assert prior_summary["kill_gate_pass"].all()
    positive_keys = set(summary["positive_exon_key"])
    backgrounds = strata.loc[~strata["is_positive"]]
    assert positive_keys.isdisjoint(backgrounds["exon_key"])

    pool = widened.loc[
        widened["detectable"] & ~widened["exon_key"].isin(positive_keys)
    ].copy()
    for event_id in prior_ids:
        stratum = strata.loc[strata["positive_event_id"].eq(event_id)]
        anchor = stratum.loc[stratum["is_positive"]].iloc[0]
        selected = stratum.loc[~stratum["is_positive"]]
        for name in COVARIATES:
            if COVARIATE_TRANSFORMS[name] == LINEAR:
                pool_values = pool[name].astype(float)
                anchor_value = float(anchor[name])
                selected_values = selected[name].astype(float)
            else:
                pool_values = np.log1p(pool[name].astype(float))
                anchor_value = float(np.log1p(anchor[name]))
                selected_values = np.log1p(selected[name].astype(float))
            scale = float(pool_values.std(ddof=0))
            selected_deviations = (selected_values - anchor_value) / scale
            assert selected_deviations.abs().max() <= config.caliper_sd + 1e-12
            assert abs(float(selected_deviations.mean())) < 0.25


def test_widened_universe_selected_backgrounds_respect_transformed_calipers() -> None:
    anchor = _universe_row("chr1:99-200:+", "chr1", 99, 200, "+")
    # Each transformed covariate independently places different candidates
    # near, exactly on, and beyond the inclusive 1-SD boundary.
    boundary = 1.0
    filler = ((8.0 - 2.0 * (1.1**2 + boundary**2 + 0.9**2)) / 2.0) ** 0.5
    offsets = [-1.1, -boundary, -0.9, -filler, filler, 0.9, boundary, 1.1]
    assert sum(value * value for value in offsets) / len(offsets) == pytest.approx(1.0)
    rows = [anchor]
    shifts = {name: index % 3 for index, name in enumerate(COVARIATES)}
    per_covariate_offsets = {
        name: offsets[shift:3] + offsets[:shift] + offsets[3:5]
        + offsets[5 + shift:8] + offsets[5:5 + shift]
        for name, shift in shifts.items()
    }
    for index in range(len(offsets)):
        covariates = {}
        for name in COVARIATES:
            offset = per_covariate_offsets[name][index]
            centre = float(anchor[name])
            if COVARIATE_TRANSFORMS[name] == LINEAR:
                covariates[name] = centre + 0.1 * offset
            else:
                covariates[name] = float(
                    np.expm1(np.log1p(centre) + 0.1 * offset)
                )
        start = 1000 + index * 500
        rows.append(_universe_row(
            f"chr2:{start}-{start + 101}:+", "chr2", start, start + 101, "+",
            **covariates,
        ))
    universe = pd.DataFrame(rows)
    positives = pd.DataFrame([{
        "event_id": "POS1", "chrom": "chr1", "start": 100, "end": 200, "strand": "+",
    }])
    config = MatchConfig(n_background=8, caliper_sd=1.0, min_background=1)
    strata, summary = build_strata(universe, positives, config)
    pool = universe.loc[universe["exon_key"] != anchor["exon_key"]]
    selected = strata.loc[~strata["is_positive"]]
    for name in COVARIATES:
        transform = (
            (lambda values: values.astype(float))
            if COVARIATE_TRANSFORMS[name] == LINEAR
            else (lambda values: np.log1p(values.astype(float)))
        )
        transformed = transform(pool[name])
        scale = transformed.std(ddof=0)
        anchor_value = transform(pd.Series([anchor[name]])).iloc[0]
        candidate_deviations = ((transformed - anchor_value) / scale).abs()
        assert (candidate_deviations < config.caliper_sd).any()
        assert any(
            value == pytest.approx(config.caliper_sd, abs=1e-12)
            for value in candidate_deviations
        )
        assert (candidate_deviations > config.caliper_sd).any()
        selected_deviations = (
            (transform(selected[name]) - anchor_value) / scale
        ).abs()
        assert (selected_deviations <= config.caliper_sd + 1e-12).all()
    assert int(summary.loc[0, "eligible_background"]) == 2


def _write_eight_positive_build_inputs(tmp_path, positives):
    gtf = tmp_path / "annotation.gtf"
    gtf_lines = []
    junction_rows_by_run = [[], [], []]
    rmats_rows = []
    rmats_header = [
        "ID", "GeneID", "geneSymbol", "chr", "strand", "exonStart_0base",
        "exonEnd", "upstreamES", "upstreamEE", "downstreamES", "downstreamEE",
        "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2",
        "IncFormLen", "SkipFormLen", "PValue", "FDR", "IncLevel1",
        "IncLevel2", "IncLevelDifference",
    ]
    annotation_only = {"OT03_SERPINB7_e4", "OT04_RTTN_e4"}
    for index, positive in positives.reset_index(drop=True).iterrows():
        start = int(positive["start"])
        end = int(positive["end"])
        upstream_start, upstream_end = start - 200, start - 101
        downstream_start, downstream_end = end + 101, end + 200
        transcript = str(positive["gencode_v50_transcript_id"])
        gene = str(positive["gencode_v50_gene_id"])
        symbol = str(positive["gene_symbol"])
        for exon_number, (exon_start, exon_end) in enumerate(
            ((upstream_start, upstream_end), (start, end), (downstream_start, downstream_end)),
            1,
        ):
            exon_id = (
                str(positive["gencode_v50_exon_id"])
                if exon_number == 2 else f"fixture-{index}-{exon_number}"
            )
            gtf_lines.append(
                f'{positive["chrom"]}\tfixture\texon\t{exon_start}\t{exon_end}'
                f'\t.\t{positive["strand"]}\t.\tgene_id "{gene}"; '
                f'transcript_id "{transcript}"; transcript_type "protein_coding"; '
                f'gene_name "{symbol}"; exon_id "{exon_id}";\n'
            )
        strand_code = 1 if positive["strand"] == "+" else 2
        per_sample_totals = (
            (66, 51, 66) if symbol == "SERPINB7"
            else (78, 33, 74) if symbol == "RTTN"
            else (30 + index, 32 + index, 34 + index)
        )
        for run_index, total in enumerate(per_sample_totals):
            left = total // 2
            right = total - left
            junction_rows_by_run[run_index].extend([
                [positive["chrom"], upstream_end + 1, start - 1, strand_code, 0, 1, left, 0, 20],
                [positive["chrom"], end + 1, downstream_start - 1, strand_code, 0, 1, right, 0, 20],
            ])
        if positive["event_id"] not in annotation_only:
            exon_length = end - start + 1
            inc_form_len = min(198, 99 + exon_length)
            rmats_rows.append([
                positive["event_id"], f'"{gene}"', f'"{symbol}"', positive["chrom"],
                positive["strand"], start - 1, end, upstream_start - 1, upstream_end,
                downstream_start - 1, downstream_end, "30,32", "0,0", "34", "0",
                inc_form_len, 99, "NA", "NA", "1,1", "1", "NA",
            ])
    gtf.write_text("".join(gtf_lines), encoding="utf-8")

    se_table = tmp_path / "SE.MATS.JC.txt"
    se_table.write_text(
        "\t".join(rmats_header) + "\n"
        + "\n".join("\t".join(map(str, row)) for row in rmats_rows) + "\n",
        encoding="utf-8",
    )
    star_dir = tmp_path / "star"
    sj_paths = []
    count_paths = []
    runs = ("SRR35298704", "SRR35298699", "SRR35298715")
    for run_index, run in enumerate(runs):
        run_dir = star_dir / run
        run_dir.mkdir(parents=True)
        sj_path = run_dir / f"{run}.SJ.out.tab"
        sample_rows = junction_rows_by_run[run_index]
        sj_path.write_text(
            "\n".join("\t".join(map(str, row)) for row in sample_rows) + "\n",
            encoding="utf-8",
        )
        sj_paths.append(sj_path)
        count_path = run_dir / f"{run}.ReadsPerGene.out.tab"
        count_path.write_text(
            "\n".join(
                [
                    "N_unmapped\t0\t0\t0",
                    "N_multimapping\t0\t0\t0",
                    "N_noFeature\t0\t0\t0",
                    "N_ambiguous\t0\t0\t0",
                    *[
                        f'{row["gencode_v50_gene_id"]}\t0\t0\t{100 + gene_index + run_index}'
                        for gene_index, (_, row) in enumerate(positives.iterrows())
                    ],
                ]
            ) + "\n",
            encoding="utf-8",
        )
        count_paths.append(count_path)

    fasta = tmp_path / "genome.fa"
    contigs = sorted(set(positives["chrom"].astype(str)))
    fasta.write_text("".join(f">{chrom}\nA\n" for chrom in contigs), encoding="utf-8")
    samtools = tmp_path / "samtools"
    fai_text = "".join(
        f"{chrom}\t1000000000\t0\t1\t2\n" for chrom in contigs
    )
    samtools.write_text(
        f"#!{sys.executable}\n"
        "from pathlib import Path\n"
        "import sys\n"
        "args = sys.argv[1:]\n"
        "fai = Path(args[args.index('--fai-idx') + 1])\n"
        "if '-r' not in args:\n"
        f"    fai.write_text({fai_text!r})\n"
        "else:\n"
        "    regions = Path(args[args.index('-r') + 1])\n"
        "    for region in regions.read_text().splitlines():\n"
        "        span = region.rsplit(':', 1)[1]\n"
        "        start, end = map(int, span.split('-'))\n"
        "        length = end - start + 1\n"
        "        print('>' + region)\n"
        "        print(('GC' * ((length + 1) // 2))[:length])\n",
        encoding="utf-8",
    )
    samtools.chmod(0o755)
    return se_table, count_paths, sj_paths, gtf, fasta, samtools


def test_build_universe_detects_all_eight_recorded_positives_end_to_end(tmp_path) -> None:
    positives = pd.read_csv(POSITIVES, sep="\t")
    se_table, count_paths, sj_paths, gtf, fasta, samtools = (
        _write_eight_positive_build_inputs(tmp_path, positives)
    )
    universe = build_universe(
        se_table,
        count_paths,
        sj_paths,
        gtf,
        fasta,
        samtools=samtools,
    )
    mapping = map_positives_to_universe(universe, positives)
    assert mapping["event_id"].tolist() == positives["event_id"].tolist()
    assert (mapping["in_universe"] & mapping["detectable"]).all()
    by_key = universe.set_index("exon_key")
    annotation_only_keys = {
        "chr18:63793160-63793277:+",
        "chr18:70201893-70201983:-",
    }
    assert set(by_key.index[by_key["n_rmats_events"].eq(0)]) == annotation_only_keys
    assert by_key.loc[list(annotation_only_keys), "annotation_event_ids"].str.startswith(
        "annotation:"
    ).all()
    assert universe["gc_fraction"].eq(1.0).all()
    assert universe["control_mean_cpm"].gt(0).all()
    assert universe["genome_fasta_sha256"].eq(
        hashlib.sha256(fasta.read_bytes()).hexdigest()
    ).all()


@pytest.mark.skipif(not UNIVERSE.is_file(), reason="committed universe parquet required")
def test_committed_universe_carries_scrmoe_provenance() -> None:
    universe = pd.read_parquet(UNIVERSE)
    for column, expected in UNIVERSE_PROVENANCE.items():
        assert universe[column].eq(expected).all()
