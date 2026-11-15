"""Independent J3 channels over committed J6 site shards.

Complementarity is an ungapped antiparallel descriptor; GU pairs count as
Watson-Crick mismatches and are also reported separately. Duplexfold minimizes
RNA:RNA energy anywhere in the strand-oriented target plus flank (no footprint
anchor). Accessibility is the RNAplfold probability of the *whole footprint*
being unpaired, not a product/mean of marginal probabilities. Neither energy nor
accessibility models modified ASO chemistry, occupancy, or potency.

Cache commits every completed channel. Requeue may replay the output stream but
reuses cached computations. A verified .done file alone commits a parquet shard.
"""
from __future__ import annotations

import argparse
from contextlib import AbstractContextManager
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import select
import signal
import sqlite3
import subprocess
import tempfile

from conference.model_v2.trans_intervention import (
    ANNOTATION, ASSEMBLY, ASOS, DEFAULT_FASTA, FaiFastaReader,
)

VERSION = "1.0.0"
VIENNA_PYTHON = "/home/dashii/orcd/pool/feng-env/bin/python"
WORKER = Path(__file__).with_name("jobs") / "j3_vienna.py"
ASSUMPTION = "unmodified RNA; chemistry NOT modelled"
CHANNELS = ("complementarity", "duplex", "accessibility")


class NotApplicable(ValueError):
    """A defined computation has no applicable result, e.g. no duplex."""


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _digest(value):
    return hashlib.sha256(_json(value).encode()).hexdigest()


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def normalize(sequence, *, empty=False):
    if not isinstance(sequence, str):
        raise ValueError("sequence must be text")
    sequence = sequence.strip().upper().replace("T", "U")
    if (not sequence and not empty) or set(sequence) - set("ACGU"):
        raise ValueError("empty or non-ACGU sequence")
    return sequence


def cache_key(aso_seq, target_seq, flank):
    """Normalized RNA triple; keep left/right flank boundary in the key."""
    left, right = flank
    return _digest([normalize(aso_seq), normalize(target_seq),
                    [normalize(left, empty=True), normalize(right, empty=True)]])


def complementarity(aso, target):
    aso, target = normalize(aso), normalize(target)
    if len(aso) != len(target):
        raise NotApplicable("ungapped descriptors require equal ASO/target lengths")
    pairs = list(zip(aso, reversed(target)))
    wc = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G")}
    gu = {("G", "U"), ("U", "G")}
    positions = [i for i, pair in enumerate(pairs, 1) if pair not in wc]
    return dict(wc_pairs=len(pairs) - len(positions),
                gu_pairs=sum(p in gu for p in pairs), mismatch_count=len(positions),
                mismatch_positions_1b=positions)


class ViennaBackend(AbstractContextManager):
    """One persistent isolated RNA worker, bounded per-request timeout."""

    def __init__(self, python=VIENNA_PYTHON, timeout=30):
        self.process = None
        self.stderr = None
        self.timeout = timeout
        self.error = None
        self.identity = {"python": str(Path(python).expanduser().absolute()),
                         "worker_sha256": _sha256(WORKER)}
        try:
            self.stderr = tempfile.TemporaryFile(mode="w+b")
            self.process = subprocess.Popen(
                [str(Path(python).expanduser()), "-u", str(WORKER)],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                stderr=self.stderr, text=True,
                env={**os.environ, "OMP_NUM_THREADS": "1"})
            self.identity.update(self._read())
        except Exception as exc:
            self.error = f"ViennaRNA worker unavailable: {type(exc).__name__}: {exc}"
            self.identity["unavailable"] = self.error
            self.close()

    def _read(self):
        if not select.select([self.process.stdout], [], [], self.timeout)[0]:
            raise TimeoutError("ViennaRNA worker timeout")
        line = self.process.stdout.readline()
        if not line:
            raise RuntimeError("ViennaRNA worker exited without a result")
        return json.loads(line)

    def call(self, method, **kwargs):
        if self.error:
            raise FileNotFoundError(self.error)
        try:
            self.process.stdin.write(_json(dict(method=method, **kwargs)) + "\n")
            self.process.stdin.flush()
            result = self._read()
        except Exception:
            self.close()
            # A timed-out worker must never answer the next request out of order.
            self.error = "ViennaRNA worker interrupted; restart to retry"
            raise
        if result["status"] == "not_applicable":
            raise NotApplicable(result["reason"])
        if result["status"] != "ok":
            raise RuntimeError(result["reason"])
        return result["result"]

    def close(self):
        if self.process is not None:
            if self.process.poll() is None:
                self.process.terminate()
                try:
                    self.process.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    self.process.kill()
                    self.process.wait()
            self.process.stdin.close()
            self.process.stdout.close()
        if self.stderr is not None:
            self.stderr.close()

    def __exit__(self, *args):
        self.close()


class ChannelCache:
    """Single-file transactional cache, including shared array-worker writers."""

    def __init__(self, path, namespace):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        self.db = sqlite3.connect(path, timeout=60)
        # Rollback journal avoids WAL's host-local shared-memory requirement.
        self.db.execute("PRAGMA journal_mode=DELETE")
        self.db.execute("PRAGMA synchronous=FULL")
        self.db.execute("CREATE TABLE IF NOT EXISTS channels "
                        "(namespace TEXT, key TEXT, channel TEXT, value TEXT, "
                        "PRIMARY KEY(namespace,key,channel))")
        self.db.commit()
        self.namespace = namespace

    def get(self, key, channel):
        if key is None:
            return None
        row = self.db.execute("SELECT value FROM channels WHERE namespace=? AND key=? AND channel=?",
                              (self.namespace, key, channel)).fetchone()
        if row is None:
            return None
        try:
            result = json.loads(row[0])
            if (result["status"] in ("ok", "not_applicable") and
                    isinstance(result["reason"], str) and result["reason"] and
                    isinstance(result["values"], dict) and
                    (result["status"] == "ok" or not result["values"])):
                return result
        except (ValueError, KeyError, TypeError):
            pass
        return None

    def put(self, key, channel, result):
        if key is not None and result["status"] in ("ok", "not_applicable"):
            with self.db:
                self.db.execute("INSERT OR REPLACE INTO channels VALUES (?,?,?,?)",
                                (self.namespace, key, channel, _json(result)))

    def close(self):
        self.db.close()


def _attempt(function):
    """Preserve exception text for every non-ok channel without swallowing signals."""
    try:
        return dict(status="ok", reason="computed", values=function())
    except Exception as exc:
        status = ("not_applicable" if isinstance(exc, NotApplicable) else
                  "missing" if isinstance(exc, (FileNotFoundError, KeyError)) else "failed")
        return dict(status=status, reason=f"{type(exc).__name__}: {exc}", values={})


def _validate_values(channel, values):
    field = "energy_kcal_mol" if channel == "duplex" else "p_unpaired"
    if not isinstance(values, dict) or field not in values:
        raise ValueError(f"malformed {channel} result: {field} missing")
    if channel == "duplex":
        value = values["energy_kcal_mol"]
        if not isinstance(value, (float, int)) or not math.isfinite(value) or abs(value) >= 99999:
            raise ValueError("invalid duplex energy")
        if not isinstance(values.get("structure"), str) or "(" not in values["structure"]:
            raise ValueError("duplex structure missing")
    else:
        value = values["p_unpaired"]
        if not isinstance(value, (float, int)) or not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError("invalid unpaired probability")
    return values


class Scorer(AbstractContextManager):
    def __init__(self, *, backend, genome, cache_path, asos=None, flank=80, window=120, span=60):
        if type(flank) is not int or not 0 <= flank <= 1000:
            raise ValueError("flank must be 0..1000 nt")
        if any(type(x) is not int or x <= 0 for x in (window, span)) or span > window:
            raise ValueError("require 1 <= span <= window")
        self.backend, self.genome = backend, genome
        self.asos = dict(asos if asos is not None else {k: v.sequence_5to3 for k, v in ASOS.items()})
        self.flank, self.window, self.span = flank, window, span
        self.params = dict(version=VERSION, backend=backend.identity, flank=flank, window=window,
                           span=span, asos=self.asos.copy(), implementation_sha256=_sha256(__file__),
                           evidence=dict(complementarity="derived", duplex="model_assumption", accessibility="proxy"))
        self.cache = ChannelCache(cache_path, _digest(self.params))

    def _context(self, row):
        if self.genome is None:
            raise FileNotFoundError("indexed GRCh38 FASTA unavailable; flanks not supplied by J6")
        start, end, strand = row["start1"], row["end1"], row["strand"]
        if (type(start) is not int or type(end) is not int or not 1 <= start <= end or
                strand not in ("+", "-") or end - start + 1 != len(normalize(row["target_sequence"]))):
            raise ValueError("invalid site coordinates/strand/footprint")
        length = self.genome.length(row["chrom"])
        lo, hi = max(1, start - self.flank), min(length, end + self.flank)
        context = normalize(self.genome.fetch(row["chrom"], lo, hi))
        if strand == "-":
            context = context.translate(str.maketrans("ACGU", "UGCA"))[::-1]
        offset = start - lo if strand == "+" else hi - end
        target = normalize(row["target_sequence"])
        if context[offset:offset + len(target)] != target:
            raise ValueError("target does not match strand-oriented GRCh38 reference")
        return context, offset, context[:offset], context[offset + len(target):]

    def score(self, row):
        if row.get("assembly") != ASSEMBLY or row.get("annotation") != ANNOTATION:
            raise ValueError("J3 requires GRCh38 / GENCODE v50; liftover before scoring")
        aso = _attempt(lambda: dict(sequence=normalize(self.asos[row["aso_id"]])))
        target = _attempt(lambda: dict(sequence=normalize(row["target_sequence"])))
        comp = (_attempt(lambda: complementarity(aso["values"]["sequence"], target["values"]["sequence"]))
                if aso["status"] == target["status"] == "ok" else
                aso if aso["status"] != "ok" else target)
        context = _attempt(lambda: dict(zip(("sequence", "offset", "left", "right"), self._context(row))))
        ctx = context["values"]
        key = (cache_key(aso["values"]["sequence"], target["values"]["sequence"], (ctx["left"], ctx["right"]))
               if aso["status"] == target["status"] == context["status"] == "ok" else None)
        results = {"complementarity": comp}
        for channel in ("duplex", "accessibility"):
            dependency = (aso if channel == "duplex" and aso["status"] != "ok" else context)
            if dependency["status"] != "ok":
                results[channel] = dependency
                continue
            cached = self.cache.get(key, channel)
            if cached and cached["status"] == "ok":
                check = _attempt(lambda: _validate_values(channel, cached["values"]))
                if check["status"] != "ok":
                    cached = None
            if cached is not None:
                results[channel] = cached
                continue
            kwargs = (dict(aso=aso["values"]["sequence"], context=ctx["sequence"]) if channel == "duplex" else
                      dict(context=ctx["sequence"], start1=ctx["offset"] + 1,
                           length=len(target["values"]["sequence"]), window=self.window, span=self.span))
            result = _attempt(lambda: _validate_values(channel, self.backend.call(channel, **kwargs)))
            results[channel] = result
            self.cache.put(key, channel, result)
        output = dict(score_cache_key=key, aso_sequence=aso["values"].get("sequence"),
                      target_with_flank=ctx.get("sequence"), footprint_start1_in_context=(ctx["offset"] + 1 if ctx else None),
                      flank_left=ctx.get("left"), flank_right=ctx.get("right"),
                      duplex_model_assumption=ASSUMPTION,
                      duplex_scope="minimum duplex anywhere in target+flank; not footprint constrained",
                      accessibility_method="RNAplfold whole-footprint unpaired probability")
        for channel, result in results.items():
            output[channel + "_status"] = result["status"]
            output[channel + "_reason"] = result["reason"]
            output[channel + "_evidence_class"] = self.params["evidence"][channel]
        for name in ("wc_pairs", "gu_pairs", "mismatch_count", "mismatch_positions_1b"):
            output["complementarity_" + name] = comp["values"].get(name)
        for channel, field, name in (("duplex", "energy_kcal_mol", "duplex_energy_kcal_mol"),
                                      ("duplex", "structure", "duplex_structure"),
                                      ("accessibility", "p_unpaired", "accessibility_p_unpaired")):
            output[name] = results[channel]["values"].get(field)
        return output

    def __exit__(self, *args):
        self.cache.close()


def _score_schema():
    import pyarrow as pa
    strings = ["score_cache_key", "aso_sequence", "target_with_flank", "flank_left", "flank_right",
               "duplex_model_assumption", "duplex_scope", "accessibility_method", "duplex_structure"]
    strings += [c + suffix for c in CHANNELS for suffix in ("_status", "_reason", "_evidence_class")]
    return pa.schema([(name, pa.string()) for name in strings] + [
        ("footprint_start1_in_context", pa.int32()),
        ("complementarity_wc_pairs", pa.int32()), ("complementarity_gu_pairs", pa.int32()),
        ("complementarity_mismatch_count", pa.int32()),
        ("complementarity_mismatch_positions_1b", pa.list_(pa.int32())),
        ("duplex_energy_kcal_mol", pa.float64()), ("accessibility_p_unpaired", pa.float64())])


def _committed_inputs(paths):
    import pyarrow.parquet as pq
    identities = []
    for path in paths:
        try:
            manifest = json.loads(path.with_suffix(".done").read_text())
            if not isinstance(manifest, dict):
                raise ValueError("completion marker is not an object")
            sha = _sha256(path)
        except (OSError, ValueError) as exc:
            raise ValueError(f"J6 input must be committed: {path}: {exc}") from exc
        if manifest.get("sha256") != sha:
            raise ValueError(f"J6 input checksum mismatch: {path}")
        if manifest.get("rows") != pq.read_metadata(path).num_rows:
            raise ValueError(f"J6 input row count mismatch: {path}")
        identities.append(dict(path=str(path.resolve()), sha256=sha, rows=manifest["rows"]))
    return identities


def _fsync_dir(path):
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _atomic_json(path, value):
    fd, temporary = tempfile.mkstemp(prefix=path.name + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(_json(value) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        _fsync_dir(path.parent)
    finally:
        Path(temporary).unlink(missing_ok=True)


def _reference_identity(genome):
    if genome is None:
        return {"status": "missing"}
    if hasattr(genome, "path"):
        path = Path(genome.path)
        stat = path.stat()
        # Never hash/read an entire reference genome on the login node.
        return dict(path=str(path.resolve()), size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                    fai_sha256=_sha256(genome.fai_path))
    if hasattr(genome, "_seqs"):
        return {"synthetic_sha256": _digest(genome._seqs)}
    raise ValueError("reference must expose indexed FASTA identity")


def run_shard(paths, output_dir, scorer, *, shard_index=0, shard_count=1):
    """Bounded streaming writer; interrupted runs retain committed cache entries."""
    import pyarrow as pa
    import pyarrow.parquet as pq
    if (type(shard_index) is not int or type(shard_count) is not int or
            not 0 <= shard_index < shard_count):
        raise ValueError("require 0 <= shard_index < shard_count")
    paths = sorted({Path(p).resolve() for p in paths})
    if not paths:
        raise ValueError("no J6 input parquet shards")
    inputs = _committed_inputs(paths)
    reference = _reference_identity(scorer.genome)
    expected = dict(params=scorer.params, inputs=inputs, reference=reference,
                    shard_index=shard_index, shard_count=shard_count)
    schema = pq.read_schema(paths[0])
    required = {"site_id", "aso_id", "chrom", "start1", "end1", "strand", "target_sequence", "assembly", "annotation"}
    if not required <= set(schema.names) or set(schema.names) & set(_score_schema().names):
        raise ValueError("expected unscored J6 site schema")
    if any(pq.read_schema(p) != schema for p in paths):
        raise ValueError("J6 input schemas differ")
    schema = pa.schema(list(schema) + list(_score_schema()), metadata={
        b"j3_provenance": _json(expected).encode(), b"duplex_model_assumption": ASSUMPTION.encode()})
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"j3_shard{shard_index:04d}_of{shard_count:04d}.parquet"
    done = path.with_suffix(".done")
    with path.with_suffix(".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            manifest = json.loads(done.read_text())
            if not isinstance(manifest, dict):
                raise ValueError("completion marker is not an object")
            if (all(manifest.get(k) == v for k, v in expected.items()) and
                    manifest["sha256"] == _sha256(path) and manifest["rows"] == pq.read_metadata(path).num_rows):
                return path
        except (OSError, ValueError, KeyError, TypeError, pa.ArrowException):
            pass
        done.unlink(missing_ok=True)
        _fsync_dir(output_dir)
        fd, temp = tempfile.mkstemp(prefix=path.stem + ".", suffix=".tmp", dir=output_dir)
        os.close(fd)
        count = 0
        counts = {c: {s: 0 for s in ("ok", "failed", "missing", "not_applicable")} for c in CHANNELS}
        try:
            with pq.ParquetWriter(temp, schema, compression="zstd") as writer:
                for source in paths:
                    for batch in pq.ParquetFile(source).iter_batches(batch_size=256):
                        rows = []
                        for row in batch.to_pylist():
                            if not isinstance(row["site_id"], str) or not row["site_id"]:
                                raise ValueError("nonempty site_id required")
                            if int(hashlib.sha256(row["site_id"].encode()).hexdigest(), 16) % shard_count != shard_index:
                                continue
                            values = scorer.score(row)
                            rows.append(dict(row, **values))
                            for channel in CHANNELS:
                                counts[channel][values[channel + "_status"]] += 1
                        if rows:
                            writer.write_table(pa.Table.from_pylist(rows, schema=schema))
                            count += len(rows)
            if _committed_inputs(paths) != inputs or _reference_identity(scorer.genome) != reference:
                raise ValueError("input/reference changed during scoring")
            with open(temp, "rb") as handle:
                os.fsync(handle.fileno())
            sha = _sha256(temp)
            os.replace(temp, path)
            _fsync_dir(output_dir)
            _atomic_json(done, dict(**expected, rows=count, sha256=sha, channel_counts=counts))
        finally:
            Path(temp).unlink(missing_ok=True)
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--cache")
    parser.add_argument("--fasta", default=DEFAULT_FASTA)
    parser.add_argument("--vienna-python", default=VIENNA_PYTHON)
    parser.add_argument("--flank", type=int, default=80)
    parser.add_argument("--window", type=int, default=120)
    parser.add_argument("--span", type=int, default=60)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--requeue-safe", action="store_true", help="exit cleanly on TERM/USR1; retain per-channel checkpoints")
    args = parser.parse_args(argv)
    if args.requeue_safe:
        def interrupted(signum, frame):
            raise SystemExit(75)
        signal.signal(signal.SIGTERM, interrupted)
        signal.signal(signal.SIGUSR1, interrupted)
    out = Path(args.output_dir).expanduser()
    paths = sorted(Path(args.input_dir).expanduser().glob("*.parquet"))
    try:
        genome = FaiFastaReader(args.fasta)
    except RuntimeError as exc:
        # Preserve rows and complementarity even when the reference is unavailable.
        print(f"reference unavailable: {exc}", flush=True)
        genome = None
    with ViennaBackend(args.vienna_python) as backend, Scorer(
            backend=backend, genome=genome, cache_path=Path(args.cache).expanduser() if args.cache else out / "j3_cache.sqlite",
            flank=args.flank, window=args.window, span=args.span) as scorer:
        path = run_shard(paths, out, scorer, shard_index=args.shard_index, shard_count=args.shard_count)
        print(path.with_suffix(".done").read_text(), end="", flush=True)


if __name__ == "__main__":
    main()
