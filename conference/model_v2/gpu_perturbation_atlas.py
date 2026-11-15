"""Deterministic reference/Pangolin N-ablation atlas (GRCh38 / GENCODE v50).

N masking measures cis-sequence sensitivity, NOT physical ASO occupancy or PSI.
Each of four tissue splice-site probability channels is reported separately, averaged over its
three released models, matching feng.titration.gpu_scores' channel selection.
All public coordinates are 1-based inclusive. FaiFastaReader converts at I/O;
mask slices convert in footprint_mask. Minus-strand input is reverse-complemented
once; readout indices then follow transcript orientation.

CLI: --smoke --output-dir DIR, or --config JSON --output-dir DIR.
JSON: {"assembly":"GRCh38", "annotation":"GENCODE v50", "windows":[
 {"locus":"id", "chrom":"chr5", "strand":"+", "center":70076574,
  "distance":50, "sites":[{"position":70076574,"kind":"donor"}]}],
 "grid":{"lengths":[18],"start":70076584,"end":70076601,"stride":1}}
Sites and grid must be prespecified independently of measured outcomes. One
absolute genomic tiling interval per config; use separate configs for distant
loci. No inferred or outcome-selected splice sites are added by this engine.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import asdict, dataclass
import fcntl
import hashlib
import json
import os
from pathlib import Path
from typing import Iterator

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from conference.model_v2 import trans_intervention as ti
from feng.gpu import checkpoint as ck
from feng.titration import gpu_scores as gs

VERSION = 1
TISSUES = (0, 2, 4, 6)  # published weight-group indices, not inferred tissue labels
# GRCh38 primary chromosomes, including mitochondrial sequence. Extra scaffolds
# are allowed. Source: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
GRCH38_PRIMARY_LENGTHS = dict(zip(
    [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"],
    [248956422, 242193529, 198295559, 190214555, 181538259,
     170805979, 159345973, 145138636, 138394717, 133797422,
     135086622, 133275309, 114364328, 107043718, 101991189,
     90338345, 83257441, 80373285, 58617616, 64444167,
     46709983, 50818468, 156040895, 57227415, 16569],
))


def digest(value) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def file_sha(path) -> str:
    h = hashlib.sha256()
    with Path(path).expanduser().resolve().open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def grch38_reference(fasta):
    """Authenticate the index's assembly signature; no whole-genome scan.

    This verifies contig lengths, not every reference base. The full index hash
    and the separately hashed input slices bind subsequent outputs to inputs.
    """
    path = Path(fasta).expanduser()
    fai = Path(str(path) + ".fai").resolve()
    index_bytes = fai.read_bytes()
    lengths = {}
    try:
        for line in index_bytes.decode("ascii").splitlines():
            fields = line.split("\t")
            if len(fields) != 5 or not fields[0] or fields[0] in lengths:
                raise ValueError("malformed or duplicate contig")
            length, offset, bases, width = map(int, fields[1:])
            if length <= 0 or offset < 0 or not 0 < bases <= width:
                raise ValueError("invalid index geometry")
            lengths[fields[0]] = length
    except (ValueError, UnicodeError) as exc:
        raise ValueError(f"GRCh38 reference index invalid: {fai}: {exc}") from exc
    mismatches = [chrom for chrom, length in GRCH38_PRIMARY_LENGTHS.items() if lengths.get(chrom) != length]
    if mismatches:
        raise ValueError(f"GRCh38 primary contig-length signature mismatch: {', '.join(mismatches)}")
    genome = ti.FaiFastaReader(str(path))
    genome.path, genome.fai_path = genome.path.resolve(), genome.fai_path.resolve()
    index_sha = hashlib.sha256(index_bytes).hexdigest()
    if file_sha(genome.fai_path) != index_sha:
        raise ValueError("GRCh38 reference index changed during authentication")
    identity = {
        "assembly": "GRCh38", "verification": "FAI-primary-contig-length-signature",
        "primary_contig_lengths": dict(GRCH38_PRIMARY_LENGTHS),
        "primary_signature_sha256": digest(GRCH38_PRIMARY_LENGTHS),
        "fasta_path": str(genome.path), "fai_path": str(genome.fai_path), "fai_sha256": index_sha,
    }
    return genome, identity


def half_open(start: int, end: int) -> tuple[int, int]:
    if type(start) is not int or type(end) is not int or not 1 <= start <= end:
        raise ValueError("invalid 1-based inclusive interval")
    return start - 1, end


def inclusive(start: int, end: int) -> tuple[int, int]:
    if type(start) is not int or type(end) is not int or not 0 <= start < end:
        raise ValueError("invalid 0-based half-open interval")
    return start + 1, end


@dataclass(frozen=True)
class Site:
    position: int
    kind: str

    def __post_init__(self):
        half_open(self.position, self.position)
        if self.kind not in ("donor", "acceptor"):
            raise ValueError("site kind must be donor or acceptor; position is exonic base")


@dataclass(frozen=True)
class Window:
    locus: str
    chrom: str
    strand: str
    center: int
    distance: int
    sites: tuple[Site, ...]

    def __post_init__(self):
        if not self.locus or not self.chrom or self.strand not in ("+", "-"):
            raise ValueError("invalid locus/chrom/strand")
        if type(self.distance) is not int or not 0 <= self.distance <= 45000:
            raise ValueError("distance must be an integer in 0..45000")
        half_open(self.center - self.distance, self.center + self.distance)
        if not self.sites or len(set(self.sites)) != len(self.sites):
            raise ValueError("nonempty unique splice sites required")
        for site in self.sites:
            gs.readout_index(self.center, site.position, self.distance)


@dataclass(frozen=True)
class Grid:
    lengths: tuple[int, ...]
    start: int
    end: int
    stride: int

    def __post_init__(self):
        half_open(self.start, self.end)
        if type(self.stride) is not int or self.stride <= 0:
            raise ValueError("positive integer stride required")
        if (not self.lengths or len(set(self.lengths)) != len(self.lengths)
                or any(type(n) is not int or not 1 <= n <= self.end - self.start + 1 for n in self.lengths)):
            raise ValueError("unique positive lengths must fit tiling interval")


def output_index(window: Window, position: int) -> int:
    index = gs.readout_index(window.center, position, window.distance)
    return index if window.strand == "+" else 2 * window.distance - index


def oriented_sequence(genome, window: Window, context: int, footprint=None) -> str:
    start, end = gs.window_bounds_1b(window.center, window.distance, context)
    half_open(start, end)
    seq = genome.fetch(window.chrom, start, end).upper()
    if len(seq) != end - start + 1 or set(seq) - set("ACGTN"):
        raise ValueError("invalid reference sequence or context length")
    if footprint is not None:
        seq = gs.footprint_mask(seq, footprint[0], footprint[1], start)
    return seq if window.strand == "+" else ti.reverse_complement(seq)


def tasks(windows, grid: Grid, shard=0, shards=1) -> Iterator[tuple[int, int, int]]:
    if type(shards) is not int or type(shard) is not int or not 0 <= shard < shards:
        raise ValueError("require 0 <= shard < shards")
    for wi, window in enumerate(windows):
        for length in sorted(grid.lengths):
            for start in range(grid.start, grid.end - length + 2, grid.stride):
                identity = [asdict(window), start, length]
                if int(digest(identity), 16) % shards == shard:
                    yield wi, start, length


class PangolinScorer:
    """Real, deterministic batched inference. Construction requires Slurm.

    Models remain on CPU between calls; only one is resident on GPU at a time.
    This bounds accelerator memory independently of ensemble size. A failed
    model call invalidates the entire batch, never a partial ensemble average.
    """
    context = gs.PANGOLIN_CONTEXT

    def __init__(self, device="cuda"):
        if not os.environ.get("SLURM_JOB_ID"):
            raise RuntimeError("real Pangolin inference requires a Slurm allocation")
        os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
        from feng import source_auth
        identity = source_auth.verify_pangolin_package()
        import torch
        import pangolin.model as pm
        if device not in ("cpu", "cuda") or (device == "cuda" and not torch.cuda.is_available()):
            raise RuntimeError("requested CUDA device unavailable or invalid device")
        self.torch, self.device = torch, device
        torch.set_num_threads(max(1, int(os.environ.get("SLURM_CPUS_PER_TASK", "1"))))
        torch.manual_seed(0)
        torch.use_deterministic_algorithms(True)
        torch.backends.cudnn.benchmark = False
        torch.backends.cudnn.deterministic = True
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
        mdir = Path(pm.__file__).resolve().parent / "models"
        hashes, self.models = {}, []
        for filename in gs.PANGOLIN_WEIGHT_FILENAMES:
            path = mdir / filename
            before = file_sha(path)
            model = pm.Pangolin(pm.L, pm.W, pm.AR)
            model.load_state_dict(torch.load(path, map_location="cpu", weights_only=True))
            if before != file_sha(path):
                raise ValueError("weights changed during load")
            model.eval()
            self.models.append(model)
            hashes[filename] = before
        root = Path(__file__).resolve().parents[2]
        source_paths = [Path(__file__), Path(ti.__file__), Path(gs.__file__), Path(ck.__file__),
                        Path(source_auth.__file__), root / "feng/src/feng/pangolin_package_seed_provenance.json"]
        sources = {str(p.resolve().relative_to(root)): file_sha(p.resolve()) for p in source_paths}
        self.provenance = {
            "model": "Pangolin-1.0.2-splice-probability-four-tissue-three-replicate",
            "output_channels": [1, 4, 7, 10],
            "model_sha256": file_sha(pm.__file__), "weights_sha256": digest(hashes),
            "weight_files": hashes, "source_sha256": digest(sources), "source_files": sources,
            **identity, "model_evidence": True, "torch_version": torch.__version__,
            "device": device, "deterministic": True, "tf32": False,
        }

    def clear_oom(self):
        if self.device == "cuda":
            self.torch.cuda.empty_cache()

    def predict(self, sequences):
        torch = self.torch
        if not sequences or len({len(s) for s in sequences}) != 1:
            raise ValueError("batch requires equal, nonempty sequence lengths")
        if any(set(s) - set("ACGTN") for s in sequences):
            raise ValueError("non-DNA model input")
        size = len(sequences[0]) - 2 * self.context
        if size < 1:
            raise ValueError("sequence shorter than model context")
        # N is all-zero, matching the existing feng scorer exactly.
        encoded = np.array([[[float(b == base) for b in s] for base in "ACGT"] for s in sequences], dtype=np.float32)
        x = torch.from_numpy(encoded).to(self.device)
        scores = np.zeros((len(sequences), 4, size), dtype=np.float64)
        try:
            with torch.inference_mode():
                for i, model in enumerate(self.models):
                    try:
                        model.to(self.device)
                        prediction = model(x)[:, (1, 4, 7, 10)[i // 3], :].cpu().numpy()
                        if prediction.shape != (len(sequences), size):
                            raise ValueError("unexpected Pangolin output shape")
                        scores[:, i // 3, :] += prediction / 3.0
                    finally:
                        model.cpu()
        finally:
            del x
        return scores


def score_batch(scorer, sequences):
    """Return (scores-or-None, reason) per input, retrying smaller OOM batches."""
    if not sequences:
        return []
    try:
        values = np.asarray(scorer.predict(sequences), dtype=float)
        shape = (len(sequences), 4, len(sequences[0]) - 2 * scorer.context)
        if values.shape != shape:
            raise ValueError(f"invalid scorer shape {values.shape}, expected {shape}")
        return [(row, None) if np.isfinite(row).all() else (None, "nonfinite model output") for row in values]
    except Exception as exc:
        reason = f"{type(exc).__name__}: {exc}"
        oom = "out of memory" in str(exc).lower()
    # Outside except: release traceback references to accelerator tensors before retry.
    if oom:
        scorer.clear_oom()
    if len(sequences) > 1:
        middle = len(sequences) // 2
        return score_batch(scorer, sequences[:middle]) + score_batch(scorer, sequences[middle:])
    return [(None, reason)]


def row_schema(provenance):
    strings = ["row_id", "locus", "chrom", "strand", "site_kind", "status", "error",
               "evidence_class", "evidence_category"]
    integers = ["window_start", "window_end", "site_position", "tissue", "footprint_start", "footprint_end", "length"]
    fields = [(s, pa.string()) for s in strings] + [(s, pa.int64()) for s in integers]
    fields += [(s, pa.float64()) for s in ("reference_score", "masked_score", "raw_model_delta")]
    fields += [("model_evidence", pa.bool_())]
    return pa.schema(fields, metadata={b"provenance": json.dumps(provenance, sort_keys=True).encode(),
                                      b"coordinates": b"GRCh38; GENCODE v50; 1-based inclusive",
                                      b"interpretation": b"N-ablation cis-feature sensitivity; not occupancy or delta PSI"})


@contextmanager
def file_lock(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a") as stream:
        fcntl.flock(stream, fcntl.LOCK_EX)
        yield


def atomic_parquet(path, rows, provenance):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    try:
        pq.write_table(pa.Table.from_pylist(rows, schema=row_schema(provenance)), temp)
        with open(temp, "rb") as stream:
            os.fsync(stream.fileno())
        os.replace(temp, path)
    finally:
        temp.unlink(missing_ok=True)
    return artifact(path, len(rows))


def artifact(path, rows):
    return {"path": str(path.resolve()), "sha256": file_sha(path), "rows": rows}


def save_manifest(data, path):
    """Detect torn/manual edits as well as corruption of referenced artifacts."""
    ck.save_checkpoint(dict(data, integrity_sha256=digest(data)), path)


def load_manifest(path):
    data = ck.load_checkpoint(path)
    if data is not None:
        actual = data.pop("integrity_sha256", None)
        if actual != digest(data):
            raise ValueError(f"manifest integrity failed: {path}")
    return data


def verify_artifact(item, provenance):
    path = Path(item["path"]).resolve()
    if not path.is_file() or file_sha(path) != item["sha256"]:
        raise ValueError(f"artifact integrity failed: {path}")
    metadata = pq.read_metadata(path)
    if metadata.num_rows != item["rows"] or json.loads(metadata.metadata[b"provenance"]) != provenance:
        raise ValueError(f"artifact integrity/metadata failed: {path}")


def rows_for(window, provenance, reference, alt=None, footprint=None, error=None):
    rows = []
    start, end = window.center - window.distance, window.center + window.distance
    for site in sorted(window.sites, key=lambda s: (s.position, s.kind)):
        for tii, tissue in enumerate(TISSUES):
            idx = output_index(window, site.position)
            ref = None if reference is None else float(reference[tii, idx])
            masked = None if alt is None else float(alt[tii, idx])
            failed = error is not None or ref is None or (footprint is not None and masked is None)
            rows.append({
                "row_id": digest([asdict(window), site.position, site.kind, tissue, footprint]),
                "locus": window.locus, "chrom": window.chrom, "strand": window.strand,
                "window_start": start, "window_end": end, "site_position": site.position,
                "site_kind": site.kind, "tissue": tissue,
                "footprint_start": None if footprint is None else footprint[0],
                "footprint_end": None if footprint is None else footprint[1],
                "length": None if footprint is None else footprint[1] - footprint[0] + 1,
                "reference_score": ref, "masked_score": masked,
                "raw_model_delta": None if failed or footprint is None else masked - ref,
                "status": "failed" if failed else "ok", "error": error or ("missing score" if failed else None),
                "evidence_class": "derived" if footprint is None else "sequence_ablation_proxy",
                "evidence_category": "not_available" if failed or not provenance["model_evidence"] else ("derived" if footprint is None else "proxy"),
                "model_evidence": provenance["model_evidence"] and not failed,
            })
    return rows


def reference_scores(window, sequence, scorer, directory, provenance):
    key = digest({"window": asdict(window), "context": scorer.context,
                  "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(), "provenance": provenance})
    marker = directory / f"{key}.done"
    with file_lock(directory / f"{key}.lock"):
        if marker.exists():
            record = load_manifest(marker)
            if record["key"] != key or record["provenance"] != provenance:
                raise ValueError("reference cache integrity failed")
            verify_artifact(record["artifact"], provenance)
            rows = pq.read_table(record["artifact"]["path"]).to_pylist()
        else:
            values, reason = score_batch(scorer, [sequence])[0]
            rows = rows_for(window, provenance, values, error=reason)
            # Failed attempts cannot overwrite artifacts cited by another run.
            path = directory / f"{key}_{digest(rows)}.parquet"
            record = {"key": key, "provenance": provenance,
                      "artifact": atomic_parquet(path, rows, provenance)}
            # Failed references remain explicit artifacts, but are not reusable cache hits.
            if reason is None:
                save_manifest(record, marker)
    return reference_from_rows(window, rows, record["artifact"])


def reference_from_rows(window, rows, item):
    values = np.full((4, 2 * window.distance + 1), np.nan)
    for row in rows:
        if row["status"] != "ok":
            return None, row["error"], item
        values[TISSUES.index(row["tissue"]), output_index(window, row["site_position"])] = row["reference_score"]
    return values, None, item


def run_atlas(windows, grid, genome, scorer, output_dir, *, batch_size=4, shard=0, shards=1, reference_dir=None):
    if type(batch_size) is not int or not 1 <= batch_size <= 64:
        raise ValueError("batch_size must be 1..64")
    windows = sorted(windows, key=lambda w: (w.locus, w.chrom, w.center, w.strand))
    if not windows or len({digest(asdict(w)) for w in windows}) != len(windows):
        raise ValueError("unique nonempty windows required")
    # Cheap bounded FASTA slices only; never hash or scan a whole genome.
    sequences = [oriented_sequence(genome, w, scorer.context) for w in windows]
    for w in windows:
        low, high = gs.window_bounds_1b(w.center, w.distance, scorer.context)
        if not low <= grid.start <= grid.end <= high:
            raise ValueError("grid outside input context")
    work = list(tasks(windows, grid, shard, shards))
    provenance = dict(scorer.provenance)
    plan = {"version": VERSION, "windows": [asdict(w) for w in windows], "grid": asdict(grid),
            "batch_size": batch_size,
            "context": scorer.context, "shard": shard, "shards": shards, "provenance": provenance,
            "sequence_sha256": [hashlib.sha256(s.encode()).hexdigest() for s in sequences]}
    run_id = digest(plan)
    out = Path(output_dir).resolve()
    references = Path(reference_dir).resolve() if reference_dir else out / "reference_scores"
    marker, checkpoint = out / f"shard_{shard}.done", out / f"shard_{shard}.checkpoint.json"
    expected = sum(len(windows[wi].sites) * 4 for wi, _, _ in work)
    with file_lock(out / f"shard_{shard}.lock"):
        if marker.exists():
            done = load_manifest(marker)
            if done["run_id"] != run_id or done["expected_rows"] != expected or done["plan"] != json.loads(json.dumps(plan)):
                raise ValueError("done manifest integrity/config mismatch")
            verify_artifact(done["ablation"], provenance)
            if done["ablation"]["rows"] != expected:
                raise ValueError("done row count integrity failed")
            for item in done["references"]:
                verify_artifact(item, provenance)
            return done
        state = load_manifest(checkpoint) or {"run_id": run_id, "next": 0, "parts": [], "references": []}
        if state["run_id"] != run_id or not 0 <= state["next"] <= len(work):
            raise ValueError("checkpoint config/integrity mismatch")
        for part in state["parts"]:
            verify_artifact(part, provenance)
        if sum(p["rows"] for p in state["parts"]) != sum(len(windows[wi].sites) * 4 for wi, _, _ in work[:state["next"]]):
            raise ValueError("checkpoint row count integrity failed")
        save_manifest(state, checkpoint)
        with ck.checkpoint_on_signal(lambda _sig: save_manifest(state, checkpoint), requeue=False):
            if state["references"]:
                if len(state["references"]) != len(windows):
                    raise ValueError("checkpoint reference count integrity failed")
                refs = []
                for w, item in zip(windows, state["references"]):
                    verify_artifact(item, provenance)
                    refs.append(reference_from_rows(w, pq.read_table(item["path"]).to_pylist(), item))
            else:
                refs = [reference_scores(w, seq, scorer, references, provenance) for w, seq in zip(windows, sequences)]
                state = dict(state, references=[r[2] for r in refs])
                save_manifest(state, checkpoint)
            while state["next"] < len(work):
                first = state["next"]
                batch = work[first:first + batch_size]
                # Equal-length contexts are mandatory for batching. Split at window boundary.
                batch = batch[:next((i for i, t in enumerate(batch) if t[0] != batch[0][0]), len(batch))]
                wi = batch[0][0]
                w = windows[wi]
                low, high = gs.window_bounds_1b(w.center, w.distance, scorer.context)
                # Use the hashed reference snapshot, never fetch a potentially changed FASTA.
                masked = []
                for _, start, n in batch:
                    offset = start - low if w.strand == "+" else high - (start + n - 1)
                    masked.append(gs.footprint_mask(sequences[wi], offset + 1, offset + n, 1))
                scored = score_batch(scorer, masked) if refs[wi][0] is not None else [(None, refs[wi][1])] * len(batch)
                rows = []
                for (_, start, n), (values, reason) in zip(batch, scored):
                    rows.extend(rows_for(w, provenance, refs[wi][0], values, (start, start + n - 1), reason))
                part = atomic_parquet(out / "parts" / f"{run_id}_{first}.parquet", rows, provenance)
                # Single assignment keeps signal snapshots internally consistent.
                state = dict(state, next=first + len(batch), parts=state["parts"] + [part])
                save_manifest(state, checkpoint)
            target = out / "ablation_maps" / f"shard_{shard}.parquet"
            target.parent.mkdir(parents=True, exist_ok=True)
            temp = target.with_name(f"{target.name}.tmp.{os.getpid()}")
            failures, total = 0, 0
            try:
                with pq.ParquetWriter(temp, row_schema(provenance)) as writer:
                    for part in state["parts"]:
                        table = pq.read_table(part["path"])
                        writer.write_table(table)
                        total += len(table)
                        failures += sum(s == "failed" for s in table["status"].to_pylist())
                if total != expected:
                    raise ValueError("final row count mismatch")
                with open(temp, "rb") as stream:
                    os.fsync(stream.fileno())
                os.replace(temp, target)
            finally:
                temp.unlink(missing_ok=True)
            done = {"run_id": run_id, "plan": plan, "provenance": provenance, "expected_rows": expected,
                    "completed_masks": len(work), "failures": failures,
                    "status": "complete_with_failures" if failures else "complete",
                    "ablation": artifact(target, total), "references": [r[2] for r in refs]}
            save_manifest(done, marker)
            return json.loads(json.dumps(done))


def cpu_gpu_tolerance(sequences, *, atol=1e-5, rtol=1e-4):
    scorer = PangolinScorer("cuda")
    gpu = scorer.predict(sequences)
    scorer.device = "cpu"
    cpu = scorer.predict(sequences)
    return {"passed": bool(np.allclose(cpu, gpu, atol=atol, rtol=rtol)),
            "max_absolute_error": float(np.max(np.abs(cpu - gpu))), "atol": atol, "rtol": rtol,
            "provenance": scorer.provenance}


def smoke_config():
    # Exon-7 donor is LAST EXONIC base, not intron first base 70076575.
    return [Window("SMN2_exon7", "chr5", "+", gs.EXON7_LAST_BASE_1B, 50,
                   (Site(gs.EXON7_LAST_BASE_1B, "donor"),))], Grid((18,), gs.ISS_N1_START_1B - 2, gs.ISS_N1_END_1B + 2, 1)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--smoke", action="store_true")
    mode.add_argument("--config", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fasta", default=ti.DEFAULT_FASTA)
    parser.add_argument("--batch-size", type=int, default=4)
    parser.add_argument("--shard", type=int, default=0)
    parser.add_argument("--shards", type=int, default=1)
    parser.add_argument("--check-tolerance", action="store_true")
    args = parser.parse_args(argv)
    if args.smoke:
        windows, grid = smoke_config()
    else:
        config = json.loads(args.config.read_text())
        if config.get("assembly") != ti.ASSEMBLY or config.get("annotation") != ti.ANNOTATION:
            raise ValueError("GRCh38 / GENCODE v50 assembly fence required")
        windows = [Window(**dict(w, sites=tuple(Site(**s) for s in w["sites"]))) for w in config["windows"]]
        grid = Grid(**dict(config["grid"], lengths=tuple(config["grid"]["lengths"])))
    # Startup has no completed work to persist. Shell remains sole requeue owner.
    with ck.checkpoint_on_signal(lambda _sig: None, requeue=False):
        # Publish only AFTER signal handlers exist and BEFORE any output work.
        # Shell checks this marker after wait, avoiding a readiness/TERM race.
        if ready_file := os.environ.get("J4_READY_FILE"):
            Path(ready_file).resolve().touch(exist_ok=False)
        genome, reference_identity = grch38_reference(args.fasta)
        scorer = PangolinScorer()
        scorer.provenance = dict(scorer.provenance, reference=reference_identity)
        done = run_atlas(windows, grid, genome, scorer, args.output_dir,
                         batch_size=args.batch_size, shard=args.shard, shards=args.shards)
        if args.check_tolerance:
            del scorer
            report = cpu_gpu_tolerance([oriented_sequence(genome, windows[0], gs.PANGOLIN_CONTEXT)])
            ck.save_checkpoint(report, args.output_dir / "cpu_gpu_tolerance.json")
            if not report["passed"]:
                return 1
        print(json.dumps({k: done[k] for k in ("status", "expected_rows", "failures", "run_id")}, sort_keys=True))
        return int(done["failures"] > 0)


if __name__ == "__main__":
    raise SystemExit(main())
