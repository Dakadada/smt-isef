"""J9 guard tests: scope, labelling, exclusions, refused metrics, rank stability.

Everything here runs without ViennaRNA and without the production reference:
the J3 backend is stubbed and the genome is a tiny in-memory contig. Real
ViennaRNA calls happen only through conference/model_v2/jobs/j9_smn2.sbatch.
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from conference.model_v2 import j9_smn2_prediction as j9  # noqa: E402
from conference.model_v2 import trans_intervention as ti  # noqa: E402

MODULE_PATH = Path(j9.__file__).resolve()

# Outcome / sealed artefacts this module must never name or open. Written out here,
# in the checker, on purpose: the assertion is that the MODULE text lacks them.
FORBIDDEN_FILENAMES = (
    "SMN2_Exon7_7_27+28_Complete_Data.xlsx",
    "SMN2_Exon7_7_27+28_Complete_Data",
    "Complete_Data.xlsx",
    "psi25",
    "7_27+28",
    "7-27",
    "7-28",
    "FENG_PHASE9_UNLOCK",
    "feng-freeze",
    "terminal",
)
FORBIDDEN_PATTERNS = (
    r"\bgel\w*\.csv\b",
    r"\boutcomes?\.(csv|tsv|parquet|xlsx)\b",
    r"\bdpsi\b",
    r"\.xlsx\b",
)


# --------------------------------------------------------------------------- fixtures
SMN2 = ti.LOCI["SMN2_exon7"]
SMN1 = ti.LOCI["SMN1_exon7"]


class SparseGenome:
    """chr5 stand-in: filler bases plus explicit patches, built per slice.

    A dense 70-Mb string would be wasteful on a shared node, and the real FASTA
    is deliberately not read by the unit tests.
    """

    def __init__(self, length: int, patches, filler: str = "T"):
        self.contig_length = length
        self.filler = filler
        self.patches = tuple(sorted(patches))  # (start_1based, sequence)

    def length(self, chrom: str) -> str:
        if chrom != "chr5":
            raise ti.GenomeAccessError(f"unknown contig {chrom!r}")
        return self.contig_length

    def fetch(self, chrom: str, start: int, end: int) -> str:
        n = self.length(chrom)
        if not (1 <= start <= end <= n):
            raise ti.GenomeAccessError(f"slice {chrom}:{start}-{end} outside 1..{n}")
        out = bytearray(self.filler.encode() * (end - start + 1))
        for patch_start, sequence in self.patches:
            lo = max(start, patch_start)
            hi = min(end, patch_start + len(sequence) - 1)
            if lo <= hi:
                out[lo - start:hi - start + 1] = sequence[lo - patch_start:hi - patch_start + 1].encode()
        return out.decode()


def _patches():
    """Real F18 target at its real offset; F14/F10 targets are prefixes of it.

    3UP8's target is placed so that exactly one candidate exists in each locus
    window, matching the F1 enumeration (SMN2 70076581-70076588).
    """
    f18 = ti.ASOS["F18"].target_sequence          # 70076584..70076601
    up8 = ti.ASOS["3UP8"].target_sequence         # 70076581..70076588
    smn1_offset = SMN1.exon_start - SMN2.exon_start
    return [
        (SMN2.exon_start, "C" * 54),
        (70076581, up8),
        (70076584, f18),
        (SMN1.exon_start, "C" * 54),
        (70076581 + smn1_offset, up8),
        (70076584 + smn1_offset, f18),
    ]


@pytest.fixture(scope="module")
def genome():
    return SparseGenome(SMN1.search_end + 1000, _patches())


@pytest.fixture(scope="module")
def registry():
    return j9.aso_registry()


@pytest.fixture(scope="module")
def built(genome, registry, tmp_path_factory):
    out = tmp_path_factory.mktemp("j9") / "j9_smn2_v1"
    backend = j9.StubBackend()
    try:
        features, provenance = j9.build_features(
            genome=genome, backend=backend,
            atlas_dir=j9.ATLAS_SENTINEL_DIR, cache_path=out, registry=registry)
    finally:
        backend.close()
    ranks = j9.build_ranks(features)
    return {"features": features, "ranks": ranks, "provenance": provenance, "out": out}


# --------------------------------------------------------------------------- scope guards
def test_module_never_names_an_outcome_or_sealed_artifact():
    text = MODULE_PATH.read_text(encoding="utf-8")
    lowered = text.lower()
    offenders = [name for name in FORBIDDEN_FILENAMES if name.lower() in lowered]
    assert offenders == [], f"module names forbidden artefact(s): {offenders}"
    for pattern in FORBIDDEN_PATTERNS:
        assert re.search(pattern, lowered) is None, f"module matches {pattern!r}"


def test_module_opens_no_file_outside_the_allowlist(tmp_path):
    for bad in (tmp_path / "outcomes.parquet", Path("/etc/passwd"),
                Path.home() / "orcd" / "scratch" / "feng" / "terminal_sealed"):
        with pytest.raises(j9.ForbiddenSourceError):
            j9.check_source_allowed(bad)
    allowed = j9.check_source_allowed(REPO_ROOT / "conference" / "aso_pipeline" / "sequences.py")
    assert allowed.name == "sequences.py"
    assert j9.check_source_allowed(j9.ATLAS_SENTINEL_DIR).name == "j4_sentinel_v2"


def test_sha256_refuses_a_path_outside_the_allowlist(tmp_path):
    victim = tmp_path / "x.txt"
    victim.write_text("x", encoding="utf-8")
    with pytest.raises(j9.ForbiddenSourceError):
        j9._sha256(victim)


def test_scientific_lock_prohibitions_are_restated_in_the_status_block():
    status = j9.STATISTICAL_STATUS
    assert status["verified_negatives"] is False
    assert j9.VERIFIED_NEGATIVES_AVAILABLE is False
    assert status["n_asos_ranked"] == 3
    assert "previously seen" in status["july_gels"].lower()
    assert "not prospective" in status["july_gels"].lower()
    assert status["supervision"].startswith("none")
    joined = " ".join(status["prohibited_outputs"]).lower()
    for banned in ("p-value", "auroc", "auprc", "fdr", "accuracy"):
        assert banned in joined


# --------------------------------------------------------------------------- forbidden metrics
@pytest.mark.parametrize("name", [
    "AUPRC", "auprc", "average_precision", "FDR", "fdr", "q_value", "AUC", "auroc",
    "roc_auc", "accuracy", "precision", "recall", "f1_score", "sklearn.metrics.auprc",
    "balanced-accuracy", "brier",
])
def test_request_metric_refuses_every_supervised_metric(name):
    with pytest.raises(j9.ForbiddenMetricError) as excinfo:
        j9.request_metric(name)
    assert "verified negatives" in str(excinfo.value)


def test_request_metric_still_refuses_when_negatives_are_merely_claimed():
    with pytest.raises(NotImplementedError):
        j9.request_metric("auprc", negatives={"verified": True, "n": 7})
    with pytest.raises(j9.ForbiddenMetricError):
        j9.request_metric("fdr", negatives={"verified": False, "n": 7})


def test_request_metric_rejects_unknown_names():
    with pytest.raises(ValueError):
        j9.request_metric("spearman_rho")


# --------------------------------------------------------------------------- 3UP8
def test_3up8_is_excluded_and_never_gets_a_unique_placement(genome, registry):
    assert "3UP8" in j9.EXCLUDED_ASOS
    assert "3UP8" not in j9.RANKED_ASO_IDS
    record = j9.place_asos(genome, registry)["3UP8"]
    assert record["status"] == "needs_human_confirmation"
    assert record["placement"] is None
    assert len(record["candidates"]) >= 1


def test_3up8_rows_are_not_available_or_hypothesis_only(built):
    rows = [r for r in built["features"] if r["aso_id"] == "3UP8"]
    assert rows, "3UP8 must still appear as explicit rows"
    assert all(r["placement_status"] != "placed" for r in rows)
    assert {r["baseline"] for r in rows} >= set(j9.BASELINES)
    excluded = [r for r in rows if r["feature_name"].endswith("_excluded")]
    assert len(excluded) == len(j9.BASELINES)
    assert all(r["evidence_class"] == "not_available" for r in excluded)
    assert all(r["status"] == "needs_human_confirmation" for r in excluded)
    hypotheses = [r for r in rows if r["is_hypothesis"]]
    assert hypotheses, "candidate placements must be carried as labelled hypotheses"
    assert all(r["hypothesis_label"] for r in hypotheses)


def test_3up8_and_hypotheses_enter_no_rank_table(built):
    assert all(r["aso_id"] != "3UP8" for r in built["ranks"])
    assert set(r["aso_id"] for r in built["ranks"]) <= set(j9.RANKED_ASO_IDS)


def test_place_asos_raises_if_an_excluded_aso_were_ever_placed(genome, monkeypatch):
    unflagged = ti.ASO("3UP8", ti.ASOS["3UP8"].sequence_5to3,
                       requires_human_confirmation=False)
    monkeypatch.setitem(ti.ASOS, "3UP8", unflagged)
    with pytest.raises(AssertionError):
        j9.place_asos(genome, {"3UP8": j9.aso_registry()["3UP8"]})


# --------------------------------------------------------------------------- labelling
def test_every_value_carries_an_evidence_label(built):
    for row in built["features"]:
        assert row["evidence_class"] in ti.EVIDENCE_CLASSES, row
        assert row["status"], row
        assert row["reason"], row
        assert row["source"], row
        has_value = any(row[field] is not None
                        for field in ("value_double", "value_int", "value_string"))
        if has_value and row["status"] == "ok":
            assert row["evidence_class"] != "not_available", row
        if not has_value:
            assert row["evidence_class"] == "not_available", row
    for row in built["ranks"]:
        assert row["evidence_class"] in ti.EVIDENCE_CLASSES, row


def test_row_builder_rejects_an_unlabeled_or_mislabeled_value(registry):
    aso = registry["F18"]
    with pytest.raises(ValueError):
        j9._row(aso=aso, baseline="M1", feature_name="x", evidence_class="bogus",
                status="ok", reason="r", source="s", value_int=1)
    with pytest.raises(ValueError):
        j9._row(aso=aso, baseline="M1", feature_name="x", evidence_class="derived",
                status="ok", reason="r", source="s")  # ok without a value
    with pytest.raises(ValueError):
        j9._row(aso=aso, baseline="M1", feature_name="x", evidence_class="not_available",
                status="ok", reason="r", source="s", value_double=1.0)
    with pytest.raises(ValueError):
        j9._row(aso=aso, baseline="M9", feature_name="x", evidence_class="derived",
                status="ok", reason="r", source="s", value_int=1)


def test_duplex_is_model_assumption_and_accessibility_is_proxy():
    assert j9._J3_EVIDENCE["duplex"] == "model_assumption"
    assert j9._J3_EVIDENCE["accessibility"] == "proxy"
    assert "chemistry NOT modelled" in j3_assumption()


def j3_assumption() -> str:
    from conference.model_v2 import j3_scoring
    return j3_scoring.ASSUMPTION


# --------------------------------------------------------------------------- registry / scrambles
def test_registry_matches_the_locked_sequence_module(registry):
    from conference.aso_pipeline import sequences

    assert set(registry) == set(sequences.ASOS)
    for aso_id, entry in registry.items():
        assert entry["sequence_5to3"] == sequences.ASOS[aso_id].dna_sequence.upper()
        assert entry["target_sequence"] == ti.reverse_complement(entry["sequence_5to3"])
        assert entry["chemistry_family"] == "unknown"
        assert entry["chemistry_evidence"] == "not_available"


def test_registry_rejects_a_divergent_sequence(monkeypatch):
    from conference.aso_pipeline import sequences

    monkeypatch.setitem(sequences.ASOS, "F18",
                        sequences.ASO("F18", "AAAAAAAAAAAAAAAAAA", "tampered"))
    with pytest.raises(ValueError, match="disagree"):
        j9.aso_registry()


def test_scrambles_are_absent_and_recorded_as_not_available(built, registry):
    controls = j9.scramble_controls(registry)
    assert controls["available"] is False
    assert controls["aso_ids"] == ()
    rows = [r for r in built["features"] if r["feature_name"] == "matched_scramble_controls"]
    assert len(rows) == 1
    assert rows[0]["evidence_class"] == "not_available"
    assert rows[0]["status"] == "not_available"
    assert "sequences.py" in rows[0]["reason"]


def test_scrambles_would_be_detected_if_they_existed(registry):
    extended = dict(registry)
    extended["F18_SCR"] = dict(registry["F18"], aso_id="F18_SCR")
    controls = j9.scramble_controls(extended)
    assert controls["available"] is True
    assert controls["aso_ids"] == ("F18_SCR",)


# --------------------------------------------------------------------------- baselines
def test_all_six_baselines_are_present_for_every_ranked_aso(built):
    for aso_id in j9.RANKED_ASO_IDS:
        seen = {r["baseline"] for r in built["features"] if r["aso_id"] == aso_id}
        assert seen == set(j9.BASELINES), (aso_id, sorted(seen))


def test_m0_is_a_seeded_null_and_is_reproducible(registry):
    first = j9.m0_rows(registry["F18"], order_rank=1, n_ordered=3)
    second = j9.m0_rows(registry["F18"], order_rank=1, n_ordered=3)
    assert first == second
    values = {r["feature_name"]: r for r in first}
    assert values["m0_n_permutations"]["value_int"] == j9.M0_PERMUTATIONS
    assert values["m0_observed_mismatch_count"]["value_int"] == 0
    assert values["m0_permuted_mismatch_mean"]["value_double"] > 0
    assert all(r["evidence_label"] in ("seeded_permutation_null", "seeded_random_order")
               for r in first)


def test_m1_complementarity_is_perfect_for_every_ranked_aso(built):
    for aso_id in j9.RANKED_ASO_IDS:
        rows = {r["feature_name"]: r for r in built["features"]
                if r["aso_id"] == aso_id and r["baseline"] == "M1"
                and r["variant_id"] == "default"}
        assert rows["m1_mismatch_count"]["value_int"] == 0
        assert rows["m1_gu_pairs"]["value_int"] == 0
        assert rows["m1_wc_pairs"]["value_int"] == rows["m1_aso_length"]["value_int"]
        assert rows["m1_wc_pairs"]["evidence_class"] == "derived"


def test_m4_geometry_is_derived_and_strand_aware(built):
    rows = {r["feature_name"]: r for r in built["features"]
            if r["aso_id"] == "F18" and r["baseline"] == "M4" and not r["is_hypothesis"]
            and r["feature_name"].startswith("m4_")}
    assert rows["m4_donor_last_exonic_base"]["value_int"] == 70076574
    assert rows["m4_acceptor_first_exonic_base"]["value_int"] == 70076521
    # F18 target 70076584-70076601 => +10..+27 from the last exonic base (F1 result).
    assert rows["m4_donor_dist_5p_from_last_exonic"]["value_int"] == 10
    assert rows["m4_donor_dist_3p_from_last_exonic"]["value_int"] == 27
    assert rows["m4_donor_dist_5p_from_first_intronic"]["value_int"] == 9
    assert rows["m4_region"]["value_string"] == "intronic"
    assert all(r["evidence_class"] == "derived" for r in rows.values())


def test_m5_reads_the_atlas_only_where_it_covers_the_footprint(built):
    atlas_status = built["provenance"]["atlas_status"]
    rows = [r for r in built["features"] if r["baseline"] == "M5"
            and r["aso_id"] in j9.RANKED_ASO_IDS]
    assert rows
    covered = [r for r in rows if r["status"] == "ok"]
    if atlas_status != "committed":
        assert covered == [], "no atlas => nothing may be labelled ok"
        return
    assert all(r["evidence_label"] == j9.ATLAS_EVIDENCE_LABEL for r in rows)
    assert all(r["evidence_class"] in ("proxy", "model_assumption", "not_available")
               for r in rows)
    # The j4_sentinel_v2 grid is length 18 only -> F18 covered, F14/F10 not.
    by_aso = {aso_id: [r for r in covered if r["aso_id"] == aso_id]
              for aso_id in j9.RANKED_ASO_IDS}
    assert by_aso["F18"], "F18 (18 nt) lies inside the sentinel footprint grid"
    for aso_id in ("F14", "F10"):
        assert by_aso[aso_id] == []
        missing = [r for r in rows if r["aso_id"] == aso_id]
        assert all(r["evidence_class"] == "not_available" for r in missing)
        assert any("does not cover" in r["reason"] for r in missing)


def test_m5_delta_is_labelled_raw_model_delta_not_dpsi(built):
    rows = [r for r in built["features"]
            if r["baseline"] == "M5" and r["status"] == "ok"
            and r["feature_name"].startswith("m5_raw_model_delta")]
    if not rows:
        pytest.skip("J4 atlas unavailable in this environment")
    assert all("not delta-PSI" in r["reason"] for r in rows)
    assert all(r["evidence_label"] == j9.ATLAS_EVIDENCE_LABEL for r in rows)


def test_atlas_inputs_reject_an_uncommitted_directory(tmp_path):
    with pytest.raises(j9.ForbiddenSourceError):
        j9.atlas_inputs(tmp_path)


# --------------------------------------------------------------------------- stub discipline
def test_stub_backend_values_are_never_a_result_and_are_never_ranked(built):
    rows = [r for r in built["features"] if r["baseline"] in ("M2", "M3")
            and r["aso_id"] in j9.RANKED_ASO_IDS
            and r["feature_name"] != "m2_model_assumption"]
    assert rows
    assert all(r["backend_kind"] == "stub" for r in rows)
    assert all(r["status"] == "stub_not_a_result" for r in rows)
    assert all(r["evidence_class"] == "not_available" for r in rows)
    assert all(r["value_double"] is None and r["value_string"] is None for r in rows)
    assert all("j9_smn2.sbatch" in r["reason"] for r in rows)
    assert [r for r in built["ranks"] if r["baseline"] in ("M2", "M3")] == []


def test_m2_assumption_row_is_carried_verbatim(built):
    rows = [r for r in built["features"] if r["feature_name"] == "m2_model_assumption"]
    assert rows
    assert all(r["evidence_class"] == "model_assumption" for r in rows)
    assert all(r["value_string"] == j3_assumption() for r in rows)


def test_null_backend_marks_m2_m3_not_available_and_points_at_the_sbatch(registry, tmp_path, genome):
    backend = j9.NullBackend()
    features, provenance = j9.build_features(
        genome=genome, backend=backend, atlas_dir=j9.ATLAS_SENTINEL_DIR,
        cache_path=tmp_path, registry=registry)
    assert provenance["backend_kind"] == "unavailable"
    rows = [r for r in features if r["baseline"] in ("M2", "M3")
            and r["aso_id"] in j9.RANKED_ASO_IDS
            and r["feature_name"] != "m2_model_assumption"]
    assert rows
    assert all(r["evidence_class"] == "not_available" for r in rows)
    assert all("j9_smn2.sbatch" in r["reason"] for r in rows)
    ranks = j9.build_ranks(features)
    assert [r for r in ranks if r["baseline"] in ("M2", "M3")] == []
    assert [r for r in ranks if r["baseline"] in ("M0", "M1", "M4")]


def test_missing_reference_degrades_instead_of_guessing(registry, tmp_path):
    backend = j9.StubBackend()
    try:
        features, provenance = j9.build_features(
            genome=None, backend=backend, atlas_dir=j9.ATLAS_SENTINEL_DIR,
            cache_path=tmp_path, registry=registry)
    finally:
        backend.close()
    assert provenance["genome_available"] is False
    assert all(r["placement_status"] != "placed" for r in features)
    geometry = [r for r in features if r["feature_name"] == "m4_splice_geometry"]
    assert geometry and all(r["evidence_class"] == "not_available" for r in geometry)
    assert j9.build_ranks(features)  # M0/M1 still rank; no placement needed


# --------------------------------------------------------------------------- ranks
def test_ranks_are_ordinal_only_and_declare_a_direction(built):
    rank_rows = [r for r in built["ranks"] if r["section"] == "rank"]
    assert rank_rows
    for row in rank_rows:
        assert row["rank_method"] == j9.RANK_METHOD
        assert row["direction"] in ("asc", "desc")
        assert 1 <= row["rank"] <= row["n_ranked"] <= len(j9.RANKED_ASO_IDS)
        assert "no significance" in row["reason"]
    assert all(not any(bad in r["feature_name"] for bad in ("auc", "fdr", "pvalue"))
               for r in rank_rows)


def test_rank_direction_is_applied(built):
    rows = {r["aso_id"]: r for r in built["ranks"]
            if r["section"] == "rank" and r["feature_name"] == "m1_wc_pairs"}
    assert rows["F18"]["rank"] == 1          # 18 nt, desc
    assert rows["F10"]["rank"] == 3          # 10 nt
    lengths = {r["aso_id"]: r for r in built["ranks"]
               if r["section"] == "rank" and r["feature_name"] == "m1_n_sites_in_locus_window"
               and r["variant_id"] == "mm0"}
    assert all(row["direction"] == "asc" for row in lengths.values())


def test_rank_stability_section_spans_the_prespecified_alternatives(built):
    stability = [r for r in built["ranks"] if r["section"] == "stability"]
    assert stability
    for row in stability:
        assert row["variant_id"] == "all"
        assert row["rank_range"] == row["max_rank"] - row["min_rank"]
        assert row["stable"] == (row["rank_range"] == 0)
        assert row["n_variants"] >= 1
    mismatch = [r for r in stability if r["feature_name"] == "m1_n_sites_in_locus_window"]
    assert mismatch and all(r["n_variants"] == len(j9.MAX_MISMATCH_ALTERNATIVES)
                            for r in mismatch)


def test_prespecified_alternatives_match_j2_and_j6(built):
    assert [(w, s) for _v, _f, w, s in j9.VIENNA_ALTERNATIVES] == [(80, 40), (120, 60), (200, 100)]
    assert j9.MAX_MISMATCH_ALTERNATIVES == (0, 1, 2, 3)
    source = (REPO_ROOT / "conference" / "model_v2" / "jobs" / "j2_fold.py").read_text()
    assert "((80, 40), (120, 60), (200, 100))" in source


def test_mismatch_sensitivity_is_monotone_nondecreasing(built):
    for aso_id in j9.RANKED_ASO_IDS:
        counts = {}
        for row in built["features"]:
            if (row["aso_id"] == aso_id
                    and row["feature_name"] == "m1_n_sites_in_locus_window"
                    and row["status"] == "ok"):
                counts[json.loads(row["variant_json"])["max_mismatch"]] = row["value_int"]
        assert set(counts) == set(j9.MAX_MISMATCH_ALTERNATIVES), (aso_id, counts)
        ordered = [counts[m] for m in sorted(counts)]
        assert ordered == sorted(ordered), (aso_id, ordered)
        assert counts[0] >= 1, "the confirmed exact site must be found at max_mismatch=0"


def test_paralog_is_reported_as_a_separate_locus_never_joined(built):
    rows = [r for r in built["features"] if r["feature_name"] == "m4_paralog_placement_status"]
    assert rows
    assert all(r["locus_id"] == "SMN1_exon7" for r in rows)
    assert all("never joined" in r["reason"] for r in rows)


# --------------------------------------------------------------------------- output artefacts
@pytest.fixture(scope="module")
def written(built, monkeypatch_module_allowlist):
    out = built["out"]
    manifest = j9.write_dataset(out, built["features"], built["ranks"],
                                built["provenance"], command="pytest")
    return {"dir": out, "manifest": manifest}


@pytest.fixture(scope="module")
def monkeypatch_module_allowlist(tmp_path_factory, built):
    """Allow the temporary output directory for the duration of the write tests."""
    original = j9._allowed_roots
    extra = built["out"].resolve()
    j9._allowed_roots = lambda: original() + (extra,)
    yield
    j9._allowed_roots = original


def test_write_dataset_emits_all_four_artifacts(written):
    directory = written["dir"]
    for name in ("features.parquet", "ranks.parquet", "manifest.json", "README.md"):
        assert (directory / name).is_file(), name


def test_manifest_records_hashes_sources_and_git_commit(written):
    manifest = json.loads((written["dir"] / "manifest.json").read_text())
    assert manifest["dataset_id"] == "j9_smn2_v1"
    assert manifest["statistical_status"] == j9.STATISTICAL_STATUS
    assert set(manifest["excluded_asos"]) == {"3UP8"}
    assert manifest["ranked_aso_ids"] == list(j9.RANKED_ASO_IDS)
    names = {o["name"] for o in manifest["outputs"]}
    assert names == {"features.parquet", "ranks.parquet"}
    for output in manifest["outputs"]:
        assert re.fullmatch(r"[0-9a-f]{64}", output["sha256"])
        assert output["bytes"] > 0
    sources = {s["path"]: s["sha256"] for s in manifest["sources"]}
    assert "conference/model_v2/j9_smn2_prediction.py" in sources
    assert "conference/aso_pipeline/sequences.py" in sources
    assert "conference/model_v2/SCIENTIFIC_LOCK.md" in sources
    assert all(re.fullmatch(r"[0-9a-f]{64}", v) for v in sources.values())
    assert manifest["git"]["commit"] == "unavailable" or re.fullmatch(
        r"[0-9a-f]{40}", manifest["git"]["commit"])
    assert "no timestamp" in manifest["note"]


def test_manifest_and_parquet_rebuild_identically(written, built):
    before = {name: (written["dir"] / name).read_bytes()
              for name in ("features.parquet", "ranks.parquet", "manifest.json")}
    j9.write_dataset(written["dir"], built["features"], built["ranks"],
                     built["provenance"], command="pytest")
    for name, blob in before.items():
        assert (written["dir"] / name).read_bytes() == blob, name


def test_parquet_carries_the_statistical_status_in_its_metadata(written):
    import pyarrow.parquet as pq

    for name in ("features.parquet", "ranks.parquet"):
        schema = pq.read_schema(written["dir"] / name)
        metadata = schema.metadata or {}
        assert metadata[b"dataset_id"] == b"j9_smn2_v1"
        assert b"no verified negatives" in metadata[b"supervision"]
        status = json.loads(metadata[b"statistical_status"])
        assert status["verified_negatives"] is False
        assert status["n_asos_ranked"] == 3


def test_features_parquet_round_trips_with_every_label_present(written):
    import pyarrow.parquet as pq

    rows = pq.read_table(written["dir"] / "features.parquet").to_pylist()
    assert rows
    assert all(r["evidence_class"] in ti.EVIDENCE_CLASSES for r in rows)
    assert all(r["dataset_version"] == j9.DATASET_VERSION for r in rows)
    assert {r["aso_id"] for r in rows} >= set(j9.RANKED_ASO_IDS) | {"3UP8"}


def test_readme_states_the_statistical_limits(written):
    text = (written["dir"] / "README.md").read_text()
    assert "Not a trained" in text
    assert "3UP8 is excluded" in text
    assert "ForbiddenMetricError" in text
    for banned in ("AUPRC", "FDR", "accuracy"):
        assert banned in text
    assert j9.ATLAS_EVIDENCE_LABEL in text


# --------------------------------------------------------------------------- CLI / sbatch
def test_cli_status_and_plan_are_machine_readable(capsys):
    assert j9.main(["status"]) == 0
    status = json.loads(capsys.readouterr().out)
    assert status["verified_negatives"] is False
    assert j9.main(["plan"]) == 0
    plan = json.loads(capsys.readouterr().out)
    assert plan["baselines"] == list(j9.BASELINES)
    assert plan["excluded_asos"] == dict(j9.EXCLUDED_ASOS)
    assert any(root.endswith("j4_sentinel_v2") for root in plan["allowed_source_roots"])


def test_sbatch_script_is_valid_quicktest_and_does_not_submit_anything():
    script = REPO_ROOT / "conference" / "model_v2" / "jobs" / "j9_smn2.sbatch"
    text = script.read_text()
    assert "--partition=mit_quicktest" in text
    assert "--gres" not in text and "gpu" not in text.lower()
    assert "sbatch" not in text.replace("j9_smn2.sbatch", "")
    assert "j9_smn2_prediction" in text
    assert "/home/dashii/orcd/pool/feng-env/bin/python" in text
    assert subprocess.run(["bash", "-n", str(script)], capture_output=True).returncode == 0
