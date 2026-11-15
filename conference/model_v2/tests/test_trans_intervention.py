"""Tests for conference.model_v2.trans_intervention.

Synthetic-genome tests always run. Golden tests read tiny slices of the GRCh38
fasta via its .fai and skip (with reason) when the fasta is absent.
"""
import dataclasses
import json
import random
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from conference.model_v2 import trans_intervention as ti  # noqa: E402

# ------------------------------------------------------------------ synthetic fixtures
_rng = random.Random(7)


def _rand(n):
    return "".join(_rng.choice("ACGT") for _ in range(n))


TX = _rand(90)  # transcript sense: 30 nt intron | 20 nt exon | 40 nt intron
PAD_L, PAD_R = _rand(10), _rand(10)
GENOME = ti.InMemoryGenome({
    "chrP": PAD_L + TX + PAD_R,                           # + strand gene, tx at 11..100
    "chrM": PAD_L + ti.reverse_complement(TX) + PAD_R,    # - strand gene, tx at 100..11
    "chrQ": PAD_R + TX + PAD_L,                           # sequence-identical paralog
})
TARGET = TX[59:71]                      # intron +10..+21 downstream of the donor
SYN_ASO = ti.ASO("SYN12", ti.reverse_complement(TARGET))


def _locus(locus_id, chrom, strand, exon, role="reference"):
    return ti.Locus(locus_id=locus_id, gene="SYN", paralog_group="SYN", role=role, chrom=chrom,
                    strand=strand, search_start=11, search_end=100, exon_label="e1",
                    exon_start=exon[0], exon_end=exon[1], coords_evidence="model_assumption")


PLUS = _locus("SYN_plus", "chrP", "+", (41, 60))
MINUS = _locus("SYN_minus", "chrM", "-", (51, 70))
PARALOG = _locus("SYN_paralog", "chrQ", "+", (41, 60), role="alternate")


def _fasta_or_skip():
    try:
        return ti.FaiFastaReader()
    except ti.GenomeAccessError as exc:
        pytest.skip(f"GRCh38 fasta/.fai not available: {exc}")


@pytest.fixture(scope="module")
def genome():
    return _fasta_or_skip()


SMN2, SMN1 = ti.LOCI["SMN2_exon7"], ti.LOCI["SMN1_exon7"]


# ------------------------------------------------------------------ sequence / strand
def test_reverse_complement_and_normalize():
    assert ti.reverse_complement("TCACTTTCATAATGCTGG") == "CCAGCATTATGAAAGTGA"
    assert ti.normalize_sequence("ucacuu") == "TCACTT"
    with pytest.raises(ti.SchemaError):
        ti.normalize_sequence("ACGX")
    with pytest.raises(ti.SchemaError):
        ti.normalize_sequence("")


def test_five_prime_three_prime_flip_is_not_placed():
    flipped = ti.ASO("SYN12_flipped", SYN_ASO.sequence_5to3[::-1])  # 3'->5' written as 5'->3'
    assert ti.place(SYN_ASO, PLUS, GENOME).status == "placed"
    assert ti.place(flipped, PLUS, GENOME).status == "no_match"
    sense = ti.ASO("SYN12_sense", TARGET)  # sense copy of the target, not antisense
    assert ti.place(sense, PLUS, GENOME).status == "no_match"


def test_plus_strand_placement_and_geometry():
    res = ti.place(SYN_ASO, PLUS, GENOME)
    p = res.placement
    assert (res.status, p.start, p.end, p.strand) == ("placed", 70, 81, "+")
    assert p.target_sequence == TARGET
    g = p.geometry
    assert (g.target_5p_pos, g.target_3p_pos) == (70, 81)
    assert (g.donor_last_exonic_base, g.donor_first_intronic_base) == (60, 61)
    assert (g.acceptor_first_exonic_base, g.acceptor_last_intronic_base) == (41, 40)
    assert (g.donor_dist_5p_from_last_exonic, g.donor_dist_3p_from_last_exonic) == (10, 21)
    assert (g.donor_dist_5p_from_first_intronic, g.donor_dist_3p_from_first_intronic) == (9, 20)
    assert (g.acceptor_dist_5p_from_first_exonic, g.acceptor_dist_3p_from_first_exonic) == (29, 40)
    assert (g.region_5p, g.region_3p, g.region) == ("intronic_downstream", "intronic_downstream", "intronic")


def test_minus_strand_gene_mirrors_plus_strand():
    plus = ti.place(SYN_ASO, PLUS, GENOME).placement
    res = ti.place(SYN_ASO, MINUS, GENOME)
    m = res.placement
    assert (res.status, m.start, m.end, m.strand) == ("placed", 30, 41, "-")
    assert m.target_sequence == TARGET  # transcript sense on both strands
    assert ti.reverse_complement(GENOME.fetch("chrM", 30, 41)) == TARGET
    gm, gp = m.geometry, plus.geometry
    assert (gm.target_5p_pos, gm.target_3p_pos) == (41, 30)  # 5' end is the higher coordinate
    assert (gm.donor_last_exonic_base, gm.donor_first_intronic_base) == (51, 50)
    assert (gm.acceptor_first_exonic_base, gm.acceptor_last_intronic_base) == (70, 71)
    for name in ("donor_dist_5p_from_last_exonic", "donor_dist_3p_from_last_exonic",
                 "donor_dist_5p_from_first_intronic", "donor_dist_3p_from_first_intronic",
                 "acceptor_dist_5p_from_first_exonic", "acceptor_dist_3p_from_first_exonic",
                 "acceptor_dist_5p_from_last_intronic", "acceptor_dist_3p_from_last_intronic",
                 "region_5p", "region_3p", "region"):
        assert getattr(gm, name) == getattr(gp, name), name


@pytest.mark.parametrize("start,end,region,r5,r3", [
    (45, 56, "exonic", "exonic", "exonic"),
    (55, 66, "exon_intron_spanning", "exonic", "intronic_downstream"),
    (35, 46, "exon_intron_spanning", "intronic_upstream", "exonic"),
    (20, 31, "intronic", "intronic_upstream", "intronic_upstream"),
    (35, 66, "exon_intron_spanning", "intronic_upstream", "intronic_downstream"),
])
def test_region_labels(start, end, region, r5, r3):
    g = ti.splice_geometry(PLUS, start, end)
    assert (g.region, g.region_5p, g.region_3p) == (region, r5, r3)


def test_donor_boundary_offsets_have_no_off_by_one():
    last_exonic = ti.splice_geometry(PLUS, 60, 60)
    first_intronic = ti.splice_geometry(PLUS, 61, 61)
    assert (last_exonic.donor_dist_5p_from_last_exonic, last_exonic.donor_dist_5p_from_first_intronic) == (0, -1)
    assert (first_intronic.donor_dist_5p_from_last_exonic, first_intronic.donor_dist_5p_from_first_intronic) == (1, 0)
    assert (last_exonic.region, first_intronic.region) == ("exonic", "intronic")
    acc_exonic = ti.splice_geometry(PLUS, 41, 41)
    acc_intronic = ti.splice_geometry(PLUS, 40, 40)
    assert (acc_exonic.acceptor_dist_5p_from_first_exonic, acc_exonic.acceptor_dist_5p_from_last_intronic) == (0, 1)
    assert (acc_intronic.acceptor_dist_5p_from_first_exonic, acc_intronic.acceptor_dist_5p_from_last_intronic) == (-1, 0)


# ------------------------------------------------------------------ window / match handling
def test_aso_longer_than_window():
    res = ti.place(SYN_ASO, PLUS, GENOME, window=(70, 80))  # 11-nt window, 12-nt ASO
    assert (res.status, res.reason, res.placement, res.candidates) == ("no_match", "aso_longer_than_window", None, ())


def test_aso_shorter_than_window_and_exact_window():
    assert ti.place(SYN_ASO, PLUS, GENOME, window=(60, 95)).placement.start == 70
    exact = ti.place(SYN_ASO, PLUS, GENOME, window=(70, 81))
    assert (exact.status, exact.placement.start, exact.placement.end) == ("placed", 70, 81)


def test_window_outside_locus_rejected():
    with pytest.raises(ValueError):
        ti.place(SYN_ASO, PLUS, GENOME, window=(5, 50))


def test_single_mismatch_is_no_match_not_fuzzy():
    seq = SYN_ASO.sequence_5to3
    mutated = ti.ASO("SYN12_mm", ("A" if seq[5] != "A" else "C").join([seq[:5], seq[6:]]))
    res = ti.place(mutated, PLUS, GENOME)
    assert (res.status, res.reason) == ("no_match", "no_exact_match")
    at = ti.place_at(mutated, PLUS, 70, 81, GENOME)
    assert (at.status, at.reason, at.placement) == ("mismatch", "sequence_mismatch:n=1", None)


def test_place_at_off_by_one_and_length():
    assert ti.place_at(SYN_ASO, PLUS, 70, 81, GENOME).status == "placed"
    for s, e in ((69, 80), (71, 82)):
        assert ti.place_at(SYN_ASO, PLUS, s, e, GENOME).status == "mismatch"
    assert ti.place_at(SYN_ASO, PLUS, 70, 80, GENOME).reason == "length_mismatch"
    assert ti.place_at(SYN_ASO, MINUS, 30, 41, GENOME).status == "placed"
    assert ti.place_at(SYN_ASO, MINUS, 31, 42, GENOME).status == "mismatch"


def test_multiple_matches_never_silently_picked():
    g = ti.InMemoryGenome({"chrR": "ACGT" * 3 + "TTGACCGGA" + "ACGT" * 2 + "TTGACCGGA" + "ACGT" * 3})
    loc = ti.Locus("REP", "REP", "REP", "reference", "chrR", "+", 1, 50, "e1", 20, 30, "model_assumption")
    res = ti.place(ti.ASO("REP9", ti.reverse_complement("TTGACCGGA")), loc, g)
    assert (res.status, res.reason, res.placement) == ("needs_human_confirmation", "multiple_exact_matches", None)
    assert [(c.start, c.end) for c in res.candidates] == [(13, 21), (30, 38)]
    assert len({c.site_id for c in res.candidates}) == 2


def test_flagged_aso_returns_candidates_even_when_unique():
    flagged = dataclasses.replace(SYN_ASO, requires_human_confirmation=True)
    res = ti.place(flagged, PLUS, GENOME)
    assert (res.status, res.placement) == ("needs_human_confirmation", None)
    assert [(c.start, c.end) for c in res.candidates] == [(70, 81)]
    assert ti.place_at(flagged, PLUS, 70, 81, GENOME).status == "needs_human_confirmation"


# ------------------------------------------------------------------ paralogs / IDs
def test_site_id_deterministic_and_locus_scoped():
    a = ti.place(SYN_ASO, PLUS, GENOME).placement
    b = ti.place(SYN_ASO, PLUS, GENOME).placement
    assert a.site_id == b.site_id == ti.target_site_id("SYN_plus", "GRCh38", "chrP", 70, 81, "+", TARGET)
    assert a.site_id.startswith("ts_") and len(a.site_id) == 23
    assert ti.target_site_id("SYN_plus", "GRCh38", "chrP", 71, 82, "+", TARGET) != a.site_id
    assert ti.target_site_id("SYN_plus", "GRCh38", "chrP", 70, 81, "-", TARGET) != a.site_id
    with pytest.raises(ti.SchemaError):
        dataclasses.replace(a, start=71, end=82)  # stale site_id rejected


def test_paralog_join_prevented_despite_identical_sequence():
    ref = ti.place(SYN_ASO, PLUS, GENOME).placement
    alt = ti.place(SYN_ASO, PARALOG, GENOME).placement
    assert ref.target_sequence == alt.target_sequence
    assert ref.site_id != alt.site_id and ref.locus_id != alt.locus_id
    assert ti.join_placements(ref, ref) is True
    with pytest.raises(ti.ParalogJoinError):
        ti.join_placements(ref, alt)
    with pytest.raises(ti.SchemaError):  # a result cannot carry another locus' placement
        ti.PlacementResult("SYN_plus", "reference", "placed", "x", 11, 100, placement=alt)


def test_intervention_reference_and_alternate_roles():
    iv = ti.build_intervention(SYN_ASO, PLUS, GENOME, alternate_loci=[PARALOG])
    assert (iv.reference.locus_id, iv.reference.role) == ("SYN_plus", "reference")
    assert [(a.locus_id, a.role, a.status) for a in iv.alternates] == [("SYN_paralog", "alternate", "placed")]
    assert ti.validate(iv) == []
    with pytest.raises(ti.SchemaError):
        ti.build_intervention(SYN_ASO, PARALOG, GENOME)  # alternate cannot be the reference
    with pytest.raises(ti.SchemaError):
        ti.build_intervention(SYN_ASO, PLUS, GENOME, alternate_loci=[PLUS])
    with pytest.raises(ti.SchemaError):
        dataclasses.replace(iv, alternates=(dataclasses.replace(iv.reference, role="alternate"),))


# ------------------------------------------------------------------ evidence / missingness
def test_unidentifiable_predictions_never_assigned():
    iv = ti.build_intervention(SYN_ASO, PLUS, GENOME)
    preds = {p.name: p for p in iv.predictions}
    assert set(preds) == {"molecular_occupancy", "observed_psi"}
    for p in preds.values():
        assert (p.value, p.evidence, p.identifiability) == (None, "not_available", "unidentifiable")
    with pytest.raises(ti.SchemaError):
        ti.Prediction("molecular_occupancy", 0.9, "derived", "identifiable")
    with pytest.raises(ti.SchemaError):
        ti.Prediction("observed_psi", 0.5, "model_assumption", "identifiable")
    with pytest.raises(ti.SchemaError):
        ti.Prediction("observed_psi", 0.5, "not_available", "unidentifiable")
    with pytest.raises(ti.SchemaError):
        ti.Prediction("anything", None, "derived", "identifiable")  # missing value must say so
    assert ti.Prediction("observed_psi", 0.5, "observed", "identifiable").value == 0.5
    with pytest.raises(ti.SchemaError):
        dataclasses.replace(iv, predictions=())  # explicit entries are mandatory


def test_evidence_classes_and_chemistry_rules():
    assert ti.EVIDENCE_CLASSES == ("observed", "derived", "proxy", "model_assumption", "unlabeled", "not_available")
    unknown = ti.Chemistry()
    assert (unknown.family, unknown.evidence, unknown.modifications) == ("unknown", "not_available", ())
    with pytest.raises(ti.SchemaError):
        ti.Chemistry("unknown", (), "observed")
    with pytest.raises(ti.SchemaError):
        ti.Chemistry("MOE-PS", (), "certain")
    mods = tuple(ti.Modification(i, "2'-MOE", "PS", "unmodified") for i in (1, 2, 12))
    aso = ti.ASO("SYN12_moe", SYN_ASO.sequence_5to3, ti.Chemistry("MOE-PS", mods, "model_assumption"))
    assert [m.position for m in aso.chemistry.modifications] == [1, 2, 12]
    with pytest.raises(ti.SchemaError):
        ti.ASO("bad", "ACGT", ti.Chemistry("MOE-PS", (ti.Modification(5),), "model_assumption"))
    with pytest.raises(ti.SchemaError):
        ti.Chemistry("MOE-PS", (ti.Modification(2), ti.Modification(2)), "model_assumption")
    for res in (ti.place(SYN_ASO, PLUS, GENOME),):
        assert res.placement.evidence == "derived" and res.placement.geometry.evidence == "derived"


def test_records_are_immutable():
    iv = ti.build_intervention(SYN_ASO, PLUS, GENOME)
    for obj, attr in ((iv, "schema_version"), (iv.aso, "sequence_5to3"), (iv.reference, "status"),
                      (iv.reference.placement, "start"), (iv.reference.placement.geometry, "region"),
                      (iv.predictions[0], "value"), (PLUS, "exon_end")):
        with pytest.raises(dataclasses.FrozenInstanceError):
            setattr(obj, attr, "x")
    assert isinstance(iv.alternates, tuple) and isinstance(iv.predictions, tuple)


# ------------------------------------------------------------------ serialization
def _examples():
    flagged = dataclasses.replace(SYN_ASO, requires_human_confirmation=True)
    return [
        ti.build_intervention(SYN_ASO, PLUS, GENOME, alternate_loci=[PARALOG]),
        ti.build_intervention(SYN_ASO, MINUS, GENOME),
        ti.build_intervention(flagged, PLUS, GENOME),
        ti.build_intervention(ti.ASO("NOHIT", "ACGTACGTACGTACGTACGT"), PLUS, GENOME),
    ]


def test_json_round_trip_is_deterministic():
    for iv in _examples():
        text = ti.to_json(iv)
        back = ti.from_json(text)
        assert back == iv
        assert ti.to_json(back) == text
        assert text.endswith("\n") and json.loads(text)["schema_version"] == ti.SCHEMA_VERSION
        assert list(ti.to_dict(iv)) == ["schema_version", "aso", "reference", "alternates", "predictions"]
        shuffled = json.dumps(dict(reversed(list(json.loads(text).items()))))
        assert ti.to_json(ti.from_json(shuffled)) == text  # key order of input is irrelevant


def test_from_dict_is_strict():
    good = ti.to_dict(_examples()[0])
    for mutate in (
        lambda d: d.pop("predictions"),
        lambda d: d.update(extra=1),
        lambda d: d.update(schema_version="0.0.1"),
        lambda d: d["reference"]["placement"].update(start="70"),
        lambda d: d["reference"]["placement"].update(start=True),
        lambda d: d["reference"]["placement"].update(start=71),
        lambda d: d["reference"].update(status="guessed"),
        lambda d: d["predictions"][0].update(value=0.7),
        lambda d: d["aso"].update(sequence_5to3="ACGTACGTACGT"),
        lambda d: d["alternates"][0].update(locus_id="SYN_plus"),
    ):
        bad = json.loads(json.dumps(good))
        mutate(bad)
        with pytest.raises(ti.SchemaError):
            ti.from_dict(bad)
    assert ti.from_dict(json.loads(json.dumps(good))) == _examples()[0]


def test_parquet_round_trip(tmp_path):
    pytest.importorskip("pyarrow", reason="pyarrow not installed; parquet round trip skipped")
    import pyarrow.parquet as pq

    ivs = _examples()
    p1, p2 = tmp_path / "a.parquet", tmp_path / "b.parquet"
    ti.write_parquet(ivs, p1)
    assert ti.read_parquet(p1) == ivs
    ti.write_parquet(ti.read_parquet(p1), p2)
    t1, t2 = pq.read_table(str(p1)), pq.read_table(str(p2))
    assert t1.equals(t2) and t1.column_names == list(ti._PARQUET_COLUMNS)
    assert t1.column("start").to_pylist() == [70, 30, None, None]  # missingness is null, not a sentinel
    assert t1.column("status").to_pylist() == ["placed", "placed", "needs_human_confirmation", "no_match"]


def test_schema_file_matches_generator():
    schema = ti.json_schema()
    assert schema["x-schema-version"] == ti.SCHEMA_VERSION
    defs = schema["$defs"]
    assert set(defs) == {"Modification", "Chemistry", "ASO", "SpliceGeometry", "Placement",
                         "PlacementResult", "Prediction", "Intervention"}
    for name, d in defs.items():
        assert d["additionalProperties"] is False and d["required"] == list(d["properties"]), name
    assert defs["Prediction"]["properties"]["evidence"]["enum"] == list(ti.EVIDENCE_CLASSES)
    assert defs["PlacementResult"]["properties"]["status"]["enum"] == list(ti.STATUSES)
    on_disk = ti.SCHEMA_DIR / ti.SCHEMA_FILENAME
    assert on_disk.is_file(), "run: python -m conference.model_v2.trans_intervention schema"
    assert on_disk.read_text(encoding="ascii") == ti.schema_text()
    for iv in _examples():  # every emitted key is declared in the schema
        d = ti.to_dict(iv)
        assert set(d) == set(defs["Intervention"]["properties"])
        assert set(d["aso"]) == set(defs["ASO"]["properties"])
        assert set(d["reference"]) == set(defs["PlacementResult"]["properties"])


# ------------------------------------------------------------------ CLI
def test_cli_schema_and_validate(tmp_path, capsys):
    assert ti.main(["schema", "--out-dir", str(tmp_path)]) == 0
    assert (tmp_path / ti.SCHEMA_FILENAME).read_text(encoding="ascii") == ti.schema_text()
    assert ti.main(["schema", "--stdout"]) == 0
    assert ti.schema_text() in capsys.readouterr().out
    good = tmp_path / "good.json"
    good.write_text(ti.to_json(_examples()[0]), encoding="ascii")
    assert ti.main(["validate", str(good)]) == 0
    assert capsys.readouterr().out.startswith("OK ")
    bad = tmp_path / "bad.json"
    bad.write_text(ti.to_json(_examples()[0]).replace('"placed"', '"guessed"'), encoding="ascii")
    assert ti.main(["validate", str(bad)]) == 1
    assert capsys.readouterr().out.startswith("INVALID")
    assert ti.main(["validate", str(tmp_path / "absent.json")]) == 1


def test_cli_place_error_paths(tmp_path, capsys):
    assert ti.main(["place", "--aso", "NOPE"]) == 2
    assert ti.main(["place", "--aso", "F18", "--locus", "NOPE"]) == 2
    assert ti.main(["place", "--aso", "F18", "--fasta", str(tmp_path / "absent.fa")]) == 2
    capsys.readouterr()


def test_fai_reader_on_tiny_fasta(tmp_path):
    seq = "".join(random.Random(3).choice("ACGT") for _ in range(137))
    body = "\n".join(seq[i:i + 10] for i in range(0, len(seq), 10))
    fa = tmp_path / "t.fa"
    fa.write_text(">c1\nACGTACGTAC\nGT\n>c2\n" + body + "\n", encoding="ascii")
    (tmp_path / "t.fa.fai").write_text("c1\t12\t4\t10\t11\nc2\t137\t22\t10\t11\n", encoding="ascii")
    r = ti.FaiFastaReader(str(fa))
    assert r.fetch("c1", 1, 12) == "ACGTACGTACGT" and r.length("c2") == 137
    for s, e in ((1, 1), (10, 11), (9, 31), (128, 137), (1, 137), (60, 60)):
        assert r.fetch("c2", s, e) == seq[s - 1:e], (s, e)
    for s, e in ((0, 5), (5, 4), (130, 138)):
        with pytest.raises(ti.GenomeAccessError):
            r.fetch("c2", s, e)
    with pytest.raises(ti.GenomeAccessError):
        r.fetch("c3", 1, 2)
    with pytest.raises(ti.GenomeAccessError):
        ti.FaiFastaReader(str(tmp_path / "absent.fa"))


# ------------------------------------------------------------------ golden (real GRCh38)
def test_golden_anchor_is_revcomp_f18(genome):
    assert genome.fetch("chr5", *ti.ANCHOR_WINDOW) == "CCAGCATTATGAAAGTGA"
    assert ti.reverse_complement(ti.ASOS["F18"].sequence_5to3) == "CCAGCATTATGAAAGTGA"
    res = ti.place(ti.ASOS["F18"], SMN2, genome, window=ti.ANCHOR_WINDOW)
    p = res.placement
    assert (res.status, p.chrom, p.start, p.end, p.strand) == ("placed", "chr5", 70076584, 70076601, "+")
    assert (p.assembly, p.annotation, p.locus_id) == ("GRCh38", "GENCODE v50", "SMN2_exon7")
    wide = ti.place(ti.ASOS["F18"], SMN2, genome)  # also unique in the whole locus window
    assert wide.status == "placed" and wide.placement.site_id == p.site_id


@pytest.mark.parametrize("aso_id,start,end", [("F14", 70076584, 70076597), ("F10", 70076584, 70076593)])
def test_golden_truncations_unique_in_anchor_window(genome, aso_id, start, end):
    aso = ti.ASOS[aso_id]
    assert ti.ASOS["F18"].sequence_5to3.endswith(aso.sequence_5to3)  # 3'-anchored truncation
    res = ti.place(aso, SMN2, genome, window=ti.ANCHOR_WINDOW)
    assert res.status == "placed" and res.candidates == ()  # exactly one match
    assert (res.placement.start, res.placement.end) == (start, end)
    assert res.placement.start == ti.ANCHOR_WINDOW[0]  # ASO 3' end == target 5' end stays fixed
    assert ti.place_at(aso, SMN2, start, end, genome).status == "placed"
    assert ti.place_at(aso, SMN2, start + 1, end + 1, genome).status == "mismatch"
    assert ti.place_at(aso, SMN2, start - 1, end - 1, genome).status == "mismatch"


def test_golden_3up8_needs_human_confirmation(genome):
    res = ti.place(ti.ASOS["3UP8"], SMN2, genome)
    assert (res.status, res.placement) == ("needs_human_confirmation", None)
    cands = [(c.start, c.end) for c in res.candidates]
    assert (70076581, 70076588) in cands and len(cands) >= 1
    for c in res.candidates:
        assert genome.fetch("chr5", c.start, c.end) == "CTGCCAGC"
    iv = ti.build_intervention(ti.ASOS["3UP8"], SMN2, genome, alternate_loci=[SMN1])
    assert iv.reference.status == "needs_human_confirmation" and iv.reference.placement is None
    assert all(a.status != "placed" for a in iv.alternates)
    assert ti.from_json(ti.to_json(iv)) == iv


def test_golden_donor_both_conventions(genome):
    # Repo convention: exon-7 last base 70076574, first intron-7 base 70076575.
    # The "donor at 70076574" note == the last-exonic-base convention.
    assert (SMN2.donor_last_exonic_base, SMN2.donor_first_intronic_base) == (70076574, 70076575)
    assert genome.fetch("chr5", 70076575, 70076576) == "GT"   # intron starts at ...575
    assert genome.fetch("chr5", 70076574, 70076575) != "GT"   # off-by-one guard
    assert genome.fetch("chr5", 70076576, 70076577) != "GT"
    g = ti.place(ti.ASOS["F18"], SMN2, genome, window=ti.ANCHOR_WINDOW).placement.geometry
    assert (g.donor_dist_5p_from_last_exonic, g.donor_dist_3p_from_last_exonic) == (10, 27)      # intron +10..+27
    assert (g.donor_dist_5p_from_first_intronic, g.donor_dist_3p_from_first_intronic) == (9, 26)
    assert g.donor_dist_5p_from_last_exonic - g.donor_dist_5p_from_first_intronic == 1
    assert g.region == "intronic" and g.region_5p == "intronic_downstream"


def test_golden_acceptor_and_exon_length(genome):
    assert (SMN2.acceptor_first_exonic_base, SMN2.acceptor_last_intronic_base) == (70076521, 70076520)
    assert SMN2.exon_end - SMN2.exon_start + 1 == 54
    assert genome.fetch("chr5", 70076519, 70076520) == "AG"
    g = ti.place(ti.ASOS["F18"], SMN2, genome, window=ti.ANCHOR_WINDOW).placement.geometry
    assert (g.acceptor_dist_5p_from_first_exonic, g.acceptor_dist_5p_from_last_intronic) == (63, 64)


def test_golden_smn1_is_separate_locus(genome):
    iv = ti.build_intervention(ti.ASOS["F18"], SMN2, genome, alternate_loci=[SMN1], window=ti.ANCHOR_WINDOW)
    ref, alt = iv.reference.placement, iv.alternates[0].placement
    assert (alt.locus_id, alt.start, alt.end, alt.strand) == ("SMN1_exon7", 70952004, 70952021, "+")
    assert ref.target_sequence == alt.target_sequence and ref.site_id != alt.site_id
    with pytest.raises(ti.ParalogJoinError):
        ti.join_placements(ref, alt)
    assert alt.geometry.donor_dist_5p_from_last_exonic == ref.geometry.donor_dist_5p_from_last_exonic == 10
    assert genome.fetch("chr5", 70951995, 70951996) == "GT" and genome.fetch("chr5", 70951939, 70951940) == "AG"
    # exon-7 position 6 distinguishes the paralogs: loci are real, distinct sequence contexts
    assert genome.fetch("chr5", SMN2.exon_start + 5, SMN2.exon_start + 5) == "T"
    assert genome.fetch("chr5", SMN1.exon_start + 5, SMN1.exon_start + 5) == "C"


def test_golden_cli_place(genome, tmp_path):
    out = tmp_path / "f18.json"
    rc = ti.main(["place", "--aso", "F18", "--window", "70076584-70076601", "--alt-locus", "SMN1_exon7", "--out", str(out)])
    assert rc == 0 and ti.main(["validate", str(out)]) == 0
    assert ti.from_json(out.read_text(encoding="ascii")).reference.placement.start == 70076584
    assert ti.main(["place", "--aso", "3UP8", "--out", str(tmp_path / "u.json")]) == 3
    assert ti.main(["place", "--seq", "ACGTACGTACGTACGTAC", "--out", str(tmp_path / "n.json")]) == 4
