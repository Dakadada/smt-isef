"""J3 contracts: independent scientific channels and crash-safe artifacts."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from conference.model_v2 import trans_intervention as ti


def module():
    from conference.model_v2 import j3_scoring
    return j3_scoring


class Backend:
    identity = {"tool": "test", "version": "1"}

    def __init__(self, failures=()):
        self.failures = failures
        self.calls = []

    def call(self, method, **kwargs):
        self.calls.append((method, kwargs))
        if method in self.failures:
            raise RuntimeError(method + " deliberately failed")
        if method == "duplex":
            return {"energy_kcal_mol": -5.8, "structure": "((((&))))"}
        return {"p_unpaired": 0.25}


def row(site_id="s", **kwargs):
    return dict(site_id=site_id, aso_id="q", chrom="chrT", start1=5, end1=8,
                strand="+", target_sequence="CCCC", assembly="GRCh38",
                annotation="GENCODE v50", **kwargs)


def scorer(tmp_path, backend=None, genome=None, **kwargs):
    return module().Scorer(backend=backend or Backend(),
                           genome=genome or ti.InMemoryGenome({"chrT": "AAAACCCCUUUU"}),
                           asos={"q": "GGGG"}, cache_path=tmp_path / "cache.sqlite",
                           flank=2, **kwargs)


def input_shard(tmp_path, rows):
    path = tmp_path / "j6.parquet"
    pq.write_table(pa.Table.from_pylist(rows), path)
    path.with_suffix(".done").write_text(json.dumps({
        "rows": len(rows), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}))
    return path


def test_antiparallel_complementarity_and_wobble():
    j = module()
    assert j.complementarity("GGGG", "CCCC") == dict(
        wc_pairs=4, gu_pairs=0, mismatch_count=0, mismatch_positions_1b=[])
    # target is 5'->3'; mismatch positions refer to ASO 5'->3'.
    assert j.complementarity("GGAU", "ACUC") == dict(
        wc_pairs=2, gu_pairs=1, mismatch_count=2, mismatch_positions_1b=[2, 3])


@pytest.mark.parametrize("failed", ["duplex", "accessibility"])
def test_failure_rows_and_channel_independence(tmp_path, failed):
    backend = Backend([failed])
    with scorer(tmp_path, backend) as s:
        r = s.score(row())
    other = "accessibility" if failed == "duplex" else "duplex"
    assert r[failed + "_status"] == "failed"
    assert "deliberately failed" in r[failed + "_reason"]
    assert r[other + "_status"] == r["complementarity_status"] == "ok"
    assert r["duplex_energy_kcal_mol" if failed == "duplex" else "accessibility_p_unpaired"] is None
    assert r["duplex_model_assumption"] == "unmodified RNA; chemistry NOT modelled"
    assert r["accessibility_evidence_class"] == "proxy"


def test_normalized_cache_key_and_flank_boundaries():
    j = module()
    assert j.cache_key(" ggtt ", "aacc", ("tt", " a ")) == j.cache_key("GGUU", "AACC", ("UU", "A"))
    assert j.cache_key("GG", "CC", ("A", "UU")) != j.cache_key("GG", "CC", ("AU", "U"))


def test_cache_reuses_success_and_retries_failure(tmp_path):
    b = Backend(["duplex"])
    with scorer(tmp_path, b) as s:
        assert s.score(row())["duplex_status"] == "failed"
    b = Backend()
    with scorer(tmp_path, b) as s:
        assert s.score(row())["duplex_status"] == "ok"
        assert s.score(row("different-id"))["duplex_status"] == "ok"
    assert [m for m, _ in b.calls] == ["duplex"]


def test_cache_configuration_isolation(tmp_path):
    b = Backend()
    with scorer(tmp_path, b, window=80) as s:
        s.score(row())
    with scorer(tmp_path, b, window=120) as s:
        s.score(row())
    assert len(b.calls) == 4


def test_missing_aso_does_not_block_accessibility(tmp_path):
    r = row()
    r["aso_id"] = "unknown"
    with scorer(tmp_path) as s:
        result = s.score(r)
    assert result["complementarity_status"] == result["duplex_status"] == "missing"
    assert result["accessibility_status"] == "ok"


def test_bad_aso_does_not_block_accessibility(tmp_path):
    with scorer(tmp_path) as s:
        s.asos["q"] = "NNNN"
        r = s.score(row())
    assert r["complementarity_status"] == r["duplex_status"] == "failed"
    assert r["accessibility_status"] == "ok"


def test_minus_strand_context_and_exact_footprint(tmp_path):
    b = Backend()
    r = row()
    r["strand"] = "-"
    with scorer(tmp_path, b, ti.InMemoryGenome({"chrT": "AAAAGGGGTTTT"})) as s:
        result = s.score(r)
    assert result["target_with_flank"] == "AACCCCUU"
    access = dict(b.calls)["accessibility"]
    assert (access["start1"], access["length"]) == (3, 4)


def test_reference_mismatch_preserves_complementarity(tmp_path):
    with scorer(tmp_path, genome=ti.InMemoryGenome({"chrT": "AAAAAAAAAAAA"})) as s:
        r = s.score(row())
    assert r["complementarity_status"] == "ok"
    assert r["duplex_status"] == r["accessibility_status"] == "failed"
    assert "reference" in r["accessibility_reason"]


def test_assembly_fence(tmp_path):
    r = row()
    r["assembly"] = "GRCh37"
    with scorer(tmp_path) as s, pytest.raises(ValueError, match="GRCh38"):
        s.score(r)


@pytest.mark.parametrize("method,field,value", [
    ("duplex", "energy_kcal_mol", float("nan")),
    ("duplex", "energy_kcal_mol", 100000.0),
    ("accessibility", "p_unpaired", float("nan")),
    ("accessibility", "p_unpaired", -0.1),
    ("accessibility", "p_unpaired", 1.2),
])
def test_invalid_backend_values_are_explicit_failures(tmp_path, method, field, value):
    class Invalid(Backend):
        def call(self, name, **kwargs):
            if name == method:
                return {field: value}
            return super().call(name, **kwargs)
    with scorer(tmp_path, Invalid()) as s:
        r = s.score(row())
    assert r[method + "_status"] == "failed"
    assert r[method + "_reason"]


def test_shard_union_restart_and_checksum(tmp_path):
    rows = [row(f"s{i}") for i in range(17)]
    path = input_shard(tmp_path, rows)
    j = module()
    with scorer(tmp_path) as s:
        whole = j.run_shard([path], tmp_path / "whole", s)
        first = pq.read_table(whole).to_pylist()
        mtime = whole.stat().st_mtime_ns
        assert j.run_shard([path], tmp_path / "whole", s) == whole
        assert whole.stat().st_mtime_ns == mtime
        parts = [j.run_shard([path], tmp_path / "parts", s, shard_index=i, shard_count=4)
                 for i in range(4)]
        combined = [r for p in parts for r in pq.read_table(p).to_pylist()]
        assert sorted(first, key=lambda r: r["site_id"]) == sorted(combined, key=lambda r: r["site_id"])
        manifest = json.loads(whole.with_suffix(".done").read_text())
        assert manifest["sha256"] == hashlib.sha256(whole.read_bytes()).hexdigest()
        assert manifest["rows"] == 17
        whole.write_bytes(b"corrupt")
        j.run_shard([path], tmp_path / "whole", s)
        assert pq.read_table(whole).num_rows == 17
        whole.with_suffix(".done").unlink()
        j.run_shard([path], tmp_path / "whole", s)
        assert whole.with_suffix(".done").is_file()


def test_failed_channels_survive_parquet_as_null_with_status(tmp_path):
    path = input_shard(tmp_path, [row()])
    with scorer(tmp_path, Backend(["duplex"])) as s:
        output = module().run_shard([path], tmp_path / "out", s)
    r = pq.read_table(output).to_pylist()[0]
    assert r["duplex_energy_kcal_mol"] is None
    assert r["duplex_status"] == "failed" and r["duplex_reason"]
    assert r["accessibility_status"] == "ok"


def test_uncommitted_and_tampered_inputs_rejected(tmp_path):
    path = input_shard(tmp_path, [row()])
    path.with_suffix(".done").unlink()
    with scorer(tmp_path) as s, pytest.raises(ValueError, match="committed"):
        module().run_shard([path], tmp_path / "out", s)
    path = input_shard(tmp_path, [row()])
    path.write_bytes(b"changed")
    with scorer(tmp_path) as s, pytest.raises(ValueError, match="checksum"):
        module().run_shard([path], tmp_path / "out", s)


def test_malformed_completion_marker_is_rebuilt(tmp_path):
    path = input_shard(tmp_path, [row()])
    with scorer(tmp_path) as s:
        output = module().run_shard([path], tmp_path / "out", s)
        output.with_suffix(".done").write_text("null")
        module().run_shard([path], tmp_path / "out", s)
    assert json.loads(output.with_suffix(".done").read_text())["rows"] == 1


def test_stale_output_rebuilt_after_input_change(tmp_path):
    path = input_shard(tmp_path, [row()])
    with scorer(tmp_path) as s:
        output = module().run_shard([path], tmp_path / "out", s)
        input_shard(tmp_path, [row(), row("new")])
        module().run_shard([path], tmp_path / "out", s)
    assert pq.read_table(output).num_rows == 2


def test_empty_shard_retains_explicit_types(tmp_path):
    path = input_shard(tmp_path, [row()])
    with scorer(tmp_path) as s:
        output = module().run_shard([path], tmp_path / "out", s, shard_index=20, shard_count=21)
    table = pq.read_table(output)
    assert table.num_rows == 0
    assert table.schema.field("duplex_energy_kcal_mol").type == pa.float64()


def test_real_vienna_hand_checkable_duplex_and_footprint():
    j = module()
    with j.ViennaBackend() as b:
        duplex = b.call("duplex", aso="GGGG", context="CCCC")
        assert duplex["structure"] == "((((&))))"
        assert duplex["energy_kcal_mol"] == pytest.approx(-5.8, abs=0.15)
        # With no possible base pairing, a complete 4-nt segment is unpaired.
        # start=1 tests the endpoint index (entry 1,4 would incorrectly be zero).
        access = b.call("accessibility", context="AAAAAA", start1=1, length=4, window=6, span=6)
        assert access["p_unpaired"] == pytest.approx(1.0)
        with pytest.raises(j.NotApplicable, match="duplex"):
            b.call("duplex", aso="AAAA", context="AAAA")


def test_missing_vienna_is_reported_per_channel(tmp_path):
    j = module()
    with j.ViennaBackend(python="/no/such/python") as b, scorer(tmp_path, b) as s:
        r = s.score(row())
    assert r["complementarity_status"] == "ok"
    assert r["duplex_status"] == r["accessibility_status"] == "missing"


def test_crash_before_marker_resumes_from_committed_channel_cache(tmp_path, monkeypatch):
    j = module()
    path = input_shard(tmp_path, [row()])
    publish = j._atomic_json

    def crash(*args):
        raise RuntimeError("simulated eviction before marker")

    with scorer(tmp_path) as s:
        monkeypatch.setattr(j, "_atomic_json", crash)
        with pytest.raises(RuntimeError, match="eviction"):
            j.run_shard([path], tmp_path / "out", s)
    assert len(list((tmp_path / "out").glob("*.parquet"))) == 1
    assert not list((tmp_path / "out").glob("*.done"))
    b = Backend(["duplex", "accessibility"])
    monkeypatch.setattr(j, "_atomic_json", publish)
    with scorer(tmp_path, b) as s:
        output = j.run_shard([path], tmp_path / "out", s)
    assert b.calls == []
    assert pq.read_table(output).to_pylist()[0]["duplex_status"] == "ok"
    assert not list((tmp_path / "out").glob("*.tmp"))


def test_cli_real_vienna_atomic_output_and_idempotence(tmp_path):
    target = ti.ASOS["F10"].target_sequence
    sequence = "A" * 12 + target + "A" * 12
    fasta = tmp_path / "tiny.fa"
    fasta.write_text(">chrT\n" + sequence + "\n")
    Path(str(fasta) + ".fai").write_text(f"chrT\t{len(sequence)}\t6\t{len(sequence)}\t{len(sequence) + 1}\n")
    inputs = tmp_path / "input"
    inputs.mkdir()
    r = row()
    r.update(aso_id="F10", start1=13, end1=12 + len(target), target_sequence=target)
    input_shard(inputs, [r])
    command = [sys.executable, "-m", "conference.model_v2.j3_scoring", "--input-dir", str(inputs),
               "--output-dir", str(tmp_path / "out"), "--fasta", str(fasta),
               "--flank", "12", "--requeue-safe"]
    env = {**os.environ, "PYTHONPATH": str(ROOT)}
    subprocess.run(command, check=True, capture_output=True, text=True, env=env, cwd=ROOT)
    output = next((tmp_path / "out").glob("*.parquet"))
    result = pq.read_table(output).to_pylist()[0]
    assert all(result[channel + "_status"] == "ok" for channel in module().CHANNELS)
    assert result["complementarity_mismatch_count"] == 0
    assert 0 <= result["accessibility_p_unpaired"] <= 1
    assert result["duplex_energy_kcal_mol"] < 0
    mtime = output.stat().st_mtime_ns
    subprocess.run(command, check=True, capture_output=True, text=True, env=env, cwd=ROOT)
    assert output.stat().st_mtime_ns == mtime


@pytest.mark.parametrize("payload", [dict(status="ok", reason="computed", values={}),
                                      dict(status="not_applicable", reason="no duplex")])
def test_bad_cached_payload_is_recomputed(tmp_path, payload):
    b = Backend()
    with scorer(tmp_path, b) as s:
        r = s.score(row())
        with s.cache.db:
            s.cache.db.execute("UPDATE channels SET value=? WHERE channel='duplex'",
                               (json.dumps(payload),))
        assert s.score(row())["duplex_status"] == "ok"
    assert [m for m, _ in b.calls].count("duplex") == 2


def test_cache_does_not_require_host_local_wal_memory(tmp_path):
    # Shards/requeues may use a shared HPC filesystem on different nodes.
    with scorer(tmp_path) as s:
        assert s.cache.db.execute("PRAGMA journal_mode").fetchone()[0] == "delete"


@pytest.mark.parametrize("method", ["duplex", "accessibility"])
def test_malformed_backend_output_is_failed_not_missing(tmp_path, method):
    class Malformed(Backend):
        def call(self, name, **kwargs):
            return {} if name == method else super().call(name, **kwargs)
    with scorer(tmp_path, Malformed()) as s:
        assert s.score(row())[method + "_status"] == "failed"


def test_no_genome_preserves_rows_and_descriptors(tmp_path):
    j = module()
    with j.Scorer(backend=Backend(), genome=None, cache_path=tmp_path / "cache.sqlite", asos={"q": "GGGG"}) as s:
        r = s.score(row())
    assert r["complementarity_status"] == "ok"
    assert r["duplex_status"] == r["accessibility_status"] == "missing"
    assert r["duplex_energy_kcal_mol"] is r["accessibility_p_unpaired"] is None


def test_length_mismatch_descriptor_not_applicable_but_other_channels_run(tmp_path):
    with scorer(tmp_path) as s:
        s.asos["q"] = "GGG"
        r = s.score(row())
    assert r["complementarity_status"] == "not_applicable"
    assert r["complementarity_wc_pairs"] is None
    assert r["duplex_status"] == r["accessibility_status"] == "ok"
