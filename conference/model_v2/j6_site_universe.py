"""Derived, substitution-only ASO sites in GRCh38 / GENCODE v50 pre-mRNA.

Coverage: complete contiguous gene bodies, including exons AND introns, on the
annotated gene strand; no intergenic sequence, no exon-union index, and no
splice-concatenated RNA. A hit must fit entirely inside its gene. Ambiguous
reference bases count as substitutions, never wildcards. K defaults to 3 (use
0, 1, 2, or 3 for the usual nested universes); edit distance is not implemented.

The GTF loader reads gene records only. Consequently region_class and transcript
membership are not_available/null unless callers supply Transcript annotations.
With annotations, junction-spanning means a *genomic exon/intron boundary*, not
a spliced exon/exon junction. Alternative-transcript classes are preserved.
Same ASO/coordinates/strand within one gene collapse to one record, with sorted
transcript IDs; different genes (including paralogs) NEVER collapse. Gene-only
records have transcript_ids=None, not a fabricated empty transcript list.

Shards own genomic leftmost starts in closed, fixed chromosome windows; a right
halo retains boundary hits on both strands. Region bounds restrict starts, not
ends. IDs depend on coordinates/gene/ASO/sequence, not sharding or K. SMN ambiguity
means the identical transcript-sense sequence was verified in the other SMN gene
in the supplied gene universe, even if it lies outside the current shard.
If the other paralog is unavailable the flag is null, never an unchecked false.
CLI workers retain the full GTF universe regardless of the requested region.

Outputs are explicit-schema Parquet, then an atomic .done JSON commit marker.
Restart verifies parameters, input identities, row count, and Parquet SHA256.
Input identities use resolved path/size/mtime_ns/inode plus FAI SHA256; they are
filesystem identities, NOT full FASTA/GTF content hashes. Inputs must be immutable
during a run. Production requires the pinned 25-chromosome GRCh38 length
signature; this is not base-content authentication. --synthetic-reference is
explicit, stamps assembly=synthetic and different IDs, and cannot share a run
with production. No full reference hashing or unbounded fetch is performed.
An immutable atomic run_manifest.json lists every expected shard; incompatible
or legacy directory reuse is rejected. Use a fresh directory for a changed run.
Parquet writes use 512-row batches. iter_sites streams unsharded output;
enumerate_sites intentionally retains a sorted list for API compatibility.

CLI (from repo root, with PYTHONPATH=feng/src):
  python -m conference.model_v2.j6_site_universe --output-dir OUT [--region chr5:70000000-71100000]
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
from typing import Mapping, Sequence
from types import MappingProxyType
from itertools import islice
from contextlib import closing
import importlib

from conference.model_v2.trans_intervention import (
    ANNOTATION, ASSEMBLY, ASOS, ASO, DEFAULT_FASTA, FaiFastaReader,
    MAX_SLICE_BASES, reverse_complement,
)
from feng.scan.seed_extend import seed_and_extend

VERSION = "1.1.0"
DEFAULT_GTF = "~/orcd/scratch/feng/gencode/v50/gencode.v50.primary_assembly.annotation.gtf"
DEFAULT_WINDOW = 50_000


def _interval(start, end):
    if type(start) is not int or type(end) is not int or not 1 <= start <= end:
        raise ValueError("coordinates must be integers with 1 <= start <= end")


@dataclass(frozen=True)
class Transcript:
    transcript_id: str
    start1: int
    end1: int
    exons: tuple[tuple[int, int], ...]
    utrs: tuple[tuple[int, int], ...] = ()

    def __post_init__(self):
        _interval(self.start1, self.end1)
        if not self.transcript_id or not self.exons:
            raise ValueError("transcript ID and complete exon intervals required")
        for lo, hi in self.exons + self.utrs:
            _interval(lo, hi)
            if not self.start1 <= lo <= hi <= self.end1:
                raise ValueError("feature outside transcript")
        exons = sorted(self.exons)
        if any(a[1] >= b[0] for a, b in zip(exons, exons[1:])):
            raise ValueError("overlapping exon intervals within a transcript")
        for lo, hi in self.utrs:
            if not any(a <= lo <= hi <= b for a, b in self.exons):
                raise ValueError("UTR outside exons")


@dataclass(frozen=True)
class Gene:
    gene_id: str
    chrom: str
    start1: int
    end1: int
    strand: str
    gene_name: str = ""
    transcripts: tuple[Transcript, ...] | None = None

    def __post_init__(self):
        _interval(self.start1, self.end1)
        if not self.gene_id or not self.chrom or self.strand not in ("+", "-"):
            raise ValueError("gene ID, chromosome and +/- strand required")
        if self.transcripts is not None:
            ids = [t.transcript_id for t in self.transcripts]
            if not ids or len(ids) != len(set(ids)):
                raise ValueError("supply distinct transcripts, or None for unavailable")
            if any(not self.start1 <= t.start1 <= t.end1 <= self.end1 for t in self.transcripts):
                raise ValueError("transcript outside gene")


@dataclass(frozen=True, order=True)
class Shard:
    chrom: str
    start1: int
    end1: int

    def __post_init__(self):
        _interval(self.start1, self.end1)
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", self.chrom):
            raise ValueError("unsafe chromosome name")

    @property
    def shard_id(self):
        return f"{self.chrom}.{self.start1:012d}-{self.end1:012d}"


@dataclass(frozen=True)
class Site:
    site_id: str
    aso_id: str
    chrom: str
    start1: int
    end1: int
    strand: str
    gene_id: str
    gene_name: str
    region_class: str
    region_classes: tuple[str, ...]
    transcript_ids: tuple[str, ...] | None
    mismatch_count: int
    mismatch_positions: tuple[int, ...]
    target_sequence: str
    paralog_ambiguity: bool | None
    evidence_class: str = "derived"
    assembly: str = ASSEMBLY
    annotation: str = ANNOTATION


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def _digest(value):
    return hashlib.sha256(_json(value).encode()).hexdigest()


def _genes(genes):
    unique = {}
    for gene in genes:
        if gene.gene_id in unique and unique[gene.gene_id] != gene:
            raise ValueError(f"conflicting records for gene {gene.gene_id}")
        unique[gene.gene_id] = gene
    return tuple(sorted(unique.values(), key=lambda g: (g.chrom, g.start1, g.end1, g.gene_id)))


def _budget(value):
    if type(value) is not int:
        raise TypeError("max_mismatch must be an integer")
    if value < 0:
        raise ValueError("max_mismatch must be nonnegative")


_GTF_CACHE_CHUNK = 8 * 1024 * 1024


def _gtf_lines(path):
    """Stream GTF lines and release consumed file cache every 8 MiB.

    Slurm cgroup peaks may include file cache: the v50 GTF alone is 4.37 GiB.
    DONTNEED is advisory, read-only, and limited to bytes this reader consumed.
    Platforms/filesystems without this advisory keep the ordinary streaming path.
    """
    with Path(path).expanduser().open("rb") as handle:
        consumed = released = 0
        def release():
            nonlocal released
            if consumed > released and hasattr(os, "posix_fadvise"):
                try:
                    os.posix_fadvise(handle.fileno(), released, consumed - released, os.POSIX_FADV_DONTNEED)
                except OSError:
                    pass  # An optional cache hint must not alter annotation parsing.
            released = consumed
        try:
            for raw in handle:
                consumed += len(raw)
                if consumed - released >= _GTF_CACHE_CHUNK:
                    release()
                yield raw.decode("utf-8")
        finally:
            release()


def read_gene_bodies(path, *, region=None) -> tuple[Gene, ...]:
    """Stream one known GTF file; parse only gene lines, retain full gene bounds.

    Header must explicitly identify GRCh38 and version 50. No feature/transcript
    inference is attempted. Optional region bounds filter overlapping genes.
    """
    if region is not None:
        Shard(*region)
    genes = []
    valid_header = False
    with closing(_gtf_lines(path)) as lines:
        for line in lines:
            if line.startswith("#"):
                if line.startswith("##description:"):
                    if not (re.search(r"\bGRCh38\b", line) and re.search(r"\bversion 50\b", line)):
                        raise ValueError("GTF must declare GRCh38 / GENCODE version 50")
                    valid_header = True
                continue
            prefix = line.split("\t", 3)
            if len(prefix) < 4 or prefix[2] != "gene":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError("malformed GTF gene record")
            chrom, _, _, lo, hi, _, strand, _, attrs = fields
            lo, hi = int(lo), int(hi)
            if region and (chrom != region[0] or hi < region[1] or lo > region[2]):
                continue
            attributes = dict(re.findall(r'(\w+) "([^"\n]*)";', attrs))
            genes.append(Gene(attributes["gene_id"], chrom, lo, hi, strand,
                              attributes.get("gene_name", "")))
    if not valid_header:
        raise ValueError("missing GRCh38 / GENCODE v50 GTF description")
    return _genes(genes)


def make_shards(genes, *, window_size=DEFAULT_WINDOW, region=None) -> tuple[Shard, ...]:
    """Stable chromosome windows; region clipping affects the explicit shard ID."""
    if type(window_size) is not int or window_size < 1 or window_size > MAX_SLICE_BASES:
        raise ValueError("window_size must be 1..100000")
    if region is not None:
        Shard(*region)
    windows = set()
    for gene in _genes(genes):
        lo, hi = gene.start1, gene.end1
        if region:
            if gene.chrom != region[0]:
                continue
            lo, hi = max(lo, region[1]), min(hi, region[2])
        if lo > hi:
            continue
        for index in range((lo - 1) // window_size, (hi - 1) // window_size + 1):
            start, end = index * window_size + 1, (index + 1) * window_size
            if region:
                start, end = max(start, region[1]), min(end, region[2])
            windows.add(Shard(gene.chrom, start, end))
    return tuple(sorted(windows))


def site_sort_key(site):
    return site.chrom, site.start1, site.end1, site.strand, site.gene_id, site.aso_id


def _annotation(gene, start, end):
    if gene.transcripts is None:
        return "not_available", (), None
    classes, transcripts = set(), []
    for tx in gene.transcripts:
        if not tx.start1 <= start <= end <= tx.end1:
            continue
        transcripts.append(tx.transcript_id)
        exon_bases = sum(max(0, min(end, b) - max(start, a) + 1) for a, b in tx.exons)
        if exon_bases == 0:
            classes.add("intron")
        elif exon_bases < end - start + 1:
            classes.add("junction-spanning")
        elif any(a <= start <= end <= b for a, b in tx.utrs):
            classes.add("UTR")
        else:
            classes.add("exon")
    regions = tuple(sorted(classes))
    region = regions[0] if len(regions) == 1 else "transcript_dependent" if regions else "not_available"
    return region, regions, tuple(sorted(transcripts))


def _fetch(genome, chrom, start, end):
    """Bound every individual public FASTA fetch, including API-only large shards."""
    parts = []
    for lo in range(start, end + 1, MAX_SLICE_BASES):
        hi = min(end, lo + MAX_SLICE_BASES - 1)
        part = genome.fetch(chrom, lo, hi)
        if len(part) != hi - lo + 1:
            raise ValueError("short reference fetch")
        parts.append(part)
    return "".join(parts).upper()


# GRCh38 primary chromosome lengths (GRC accession GCA_000001405.15).
# This is an assembly signature, not a checksum of FASTA sequence content.
GRCH38_PRIMARY_LENGTHS = MappingProxyType(dict(zip(
    [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"],
    (248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
     159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
     114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
     58617616, 64444167, 46709983, 50818468, 156040895, 57227415, 16569))))


def reference_identity(genome, *, synthetic_reference=False):
    """Fail closed before fetching production sequence; synthetic is opt-in.

    All 25 primary chromosomes must match the pinned GRCh38 length signature.
    Unplaced/unlocalized scaffolds may also be present. This does not authenticate
    base content; the immutable FASTA/FAI filesystem identities fence the run.
    """
    if synthetic_reference:
        return dict(mode="synthetic")
    if not isinstance(genome, FaiFastaReader):
        raise ValueError("GRCh38 requires an indexed FASTA; use explicit synthetic_reference for tests")
    contigs = {}
    with genome.fai_path.open() as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) != 5 or fields[0] in contigs:
                raise ValueError("invalid GRCh38 FAI contig table")
            contigs[fields[0]] = int(fields[1])
    primary = {name: contigs.get(name) for name in GRCH38_PRIMARY_LENGTHS}
    if primary != dict(GRCH38_PRIMARY_LENGTHS):
        raise ValueError("reference does not match trusted GRCh38 primary chromosome length signature")
    # Alternate loci/patch chromosomes are outside the primary assembly universe.
    if any(name.endswith(("_alt", "_fix", "_novel")) for name in contigs):
        raise ValueError("GRCh38 primary assembly required, not alternate/patch loci")
    if any(genome.length(name) != length for name, length in contigs.items()):
        raise ValueError("GRCh38 reader/index identity mismatch")
    return dict(mode="GRCh38_primary", primary_lengths_sha256=_digest(primary),
                contig_lengths_sha256=_digest(contigs))


@dataclass(frozen=True, slots=True)
class _IntervalNode:
    gene: Gene
    max_end: int
    left: _IntervalNode | None = None
    right: _IntervalNode | None = None


def _interval_tree(genes):
    if not genes:
        return None
    mid = len(genes) // 2
    left, right = _interval_tree(genes[:mid]), _interval_tree(genes[mid + 1:])
    return _IntervalNode(genes[mid], max(genes[mid].end1,
                         left.max_end if left else 0, right.max_end if right else 0), left, right)


def _overlapping(node, start, end):
    if node is None or node.max_end < start:
        return
    yield from _overlapping(node.left, start, end)
    if node.gene.start1 > end:
        return
    if node.gene.end1 >= start:
        yield node.gene
    yield from _overlapping(node.right, start, end)


@dataclass(frozen=True)
class SiteWorker:
    """One immutable universe/index/fingerprint per worker, shared by all shards.

    The only mutable member is a private, bounded SMN sequence cache. Keep input
    files immutable for the worker lifetime. iter_sites streams in deterministic
    window/gene/ASO order; enumerate_sites remains the sorted-list compatibility API.
    """
    genome: object
    genes: Sequence[Gene]
    asos: Sequence[ASO]
    synthetic_reference: bool = False

    def __post_init__(self):
        reference = reference_identity(self.genome, synthetic_reference=self.synthetic_reference)
        genes = _genes(self.genes)
        asos = tuple(sorted(self.asos, key=lambda a: a.aso_id))
        if not asos or len({a.aso_id for a in asos}) != len(asos):
            raise ValueError("distinct ASO IDs required")
        grouped = {}
        digest = hashlib.sha256()
        digest.update(b"[")
        for i, gene in enumerate(genes):
            if i:
                digest.update(b",")
            digest.update(_json(asdict(gene)).encode())
            grouped.setdefault(gene.chrom, []).append(gene)
        digest.update(b"]")
        scanner = importlib.import_module("feng.scan.seed_extend")
        params = dict(version=VERSION, assembly="synthetic" if self.synthetic_reference else ASSEMBLY,
                      annotation=ANNOTATION, mode="substitutions", reference=reference,
                      asos=[dict(aso_id=a.aso_id, sequence=a.sequence_5to3) for a in asos],
                      gene_universe_sha256=digest.hexdigest(), implementation_sha256=_sha256(__file__),
                      scanner_sha256=_sha256(scanner.__file__))
        object.__setattr__(self, "genes", genes)
        object.__setattr__(self, "asos", asos)
        object.__setattr__(self, "params_json", _json(params))
        object.__setattr__(self, "intervals", MappingProxyType({c: _interval_tree(gs) for c, gs in grouped.items()}))
        object.__setattr__(self, "paralogs", MappingProxyType({name: tuple(g for g in genes if g.gene_name == name)
                                                           for name in ("SMN1", "SMN2")}))
        object.__setattr__(self, "_paralog_sequences", {})

    def iter_sites(self, *, max_mismatch=3, shard=None):
        _budget(max_mismatch)
        windows = make_shards(self.genes) if shard is None else (shard,)
        for window in windows:
            # Bound scanner hit/candidate sets even for a caller-supplied huge shard.
            for lo in range(window.start1, window.end1 + 1, DEFAULT_WINDOW):
                yield from self._iter_window(Shard(window.chrom, lo, min(window.end1, lo + DEFAULT_WINDOW - 1)),
                                             max_mismatch)

    def _iter_window(self, shard, max_mismatch):
        selected = tuple(_overlapping(self.intervals.get(shard.chrom), shard.start1, shard.end1))
        if not selected:
            return
        length = self.genome.length(shard.chrom)
        if any(g.end1 > length for g in selected):
            raise ValueError("gene outside reference chromosome")
        halo_end = min(length, shard.end1 + max(a.length for a in self.asos) - 1)
        sequence = _fetch(self.genome, shard.chrom, shard.start1, halo_end)
        for gene in selected:
            start = max(gene.start1, shard.start1)
            end = min(gene.end1, halo_end)
            reference = sequence[start - shard.start1:end - shard.start1 + 1]
            for aso in self.asos:
                query = aso.target_sequence if gene.strand == "+" else aso.sequence_5to3
                for hit in seed_and_extend(query, reference, max_mismatch=max_mismatch):
                    lo, hi = start + hit.start, start + hit.end - 1
                    if lo > shard.end1:
                        continue
                    target = reference[hit.start:hit.end]
                    if gene.strand == "-":
                        target = reverse_complement(target)
                    positions = tuple(sorted(aso.length - i for i, (a, b) in
                                             enumerate(zip(aso.target_sequence, target)) if a != b))
                    if len(positions) != hit.mismatches:
                        raise ValueError("scanner/mismatch-position disagreement")
                    ambiguous = False
                    if gene.gene_name in self.paralogs:
                        other_name = "SMN1" if gene.gene_name == "SMN2" else "SMN2"
                        others = self.paralogs[other_name]
                        ambiguous = False if others else None
                        for other in others:
                            if other.gene_id not in self._paralog_sequences:
                                seq = _fetch(self.genome, other.chrom, other.start1, other.end1)
                                self._paralog_sequences[other.gene_id] = seq if other.strand == "+" else reverse_complement(seq)
                            if target in self._paralog_sequences[other.gene_id]:
                                ambiguous = True
                                break
                    region, regions, transcripts = _annotation(gene, lo, hi)
                    assembly = "synthetic" if self.synthetic_reference else ASSEMBLY
                    identity = [assembly, ANNOTATION, aso.aso_id, aso.sequence_5to3,
                                gene.gene_id, gene.chrom, lo, hi, gene.strand, target]
                    yield Site("j6-" + _digest(identity), aso.aso_id, gene.chrom, lo, hi,
                               gene.strand, gene.gene_id, gene.gene_name, region, regions,
                               transcripts, hit.mismatches, positions, target, ambiguous, assembly=assembly)


def iter_sites(genome, genes, asos, *, max_mismatch=3, shard=None, synthetic_reference=False):
    """Stream an unsharded or sharded universe without retaining output records."""
    return SiteWorker(genome, genes, asos, synthetic_reference).iter_sites(max_mismatch=max_mismatch, shard=shard)


def enumerate_sites(genome, genes: Sequence[Gene], asos: Sequence[ASO], *,
                    max_mismatch=3, shard: Shard | None = None, synthetic_reference=False) -> list[Site]:
    """Sorted-list compatibility API; use iter_sites for memory-bounded output."""
    return sorted(iter_sites(genome, genes, asos, max_mismatch=max_mismatch, shard=shard,
                             synthetic_reference=synthetic_reference), key=site_sort_key)


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_identity(path):
    path = Path(path).expanduser().resolve(strict=True)
    info = path.stat()
    result = dict(path=str(path), size=info.st_size, mtime_ns=info.st_mtime_ns, inode=info.st_ino)
    fai = Path(str(path) + ".fai")
    if fai.is_file():
        result["fai_sha256"] = _sha256(fai)
    return result


def _check_identities(identities):
    for identity in identities.values():
        if "path" in identity and file_identity(identity["path"]) != identity:
            raise ValueError("reference input identity changed during run")


def _schema():
    import pyarrow as pa
    return pa.schema([
        ("site_id", pa.string()), ("aso_id", pa.string()), ("chrom", pa.string()),
        ("start1", pa.int64()), ("end1", pa.int64()), ("strand", pa.string()),
        ("gene_id", pa.string()), ("gene_name", pa.string()), ("region_class", pa.string()),
        ("region_classes", pa.list_(pa.string())), ("transcript_ids", pa.list_(pa.string())),
        ("mismatch_count", pa.int32()), ("mismatch_positions", pa.list_(pa.int32())),
        ("target_sequence", pa.string()), ("paralog_ambiguity", pa.bool_()),
        ("evidence_class", pa.string()), ("assembly", pa.string()), ("annotation", pa.string()),
    ])


def _fsync_dir(path):
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


@dataclass(frozen=True)
class RunConfig:
    output_dir: Path
    manifest_json: str
    worker_params_json: str
    identities_json: str
    expected_shards: frozenset[str]
    max_mismatch: int
    manifest_sha256: str
    manifest_stat: tuple[int, ...]


def _manifest_stat(path):
    info = path.stat()
    return info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns, info.st_ctime_ns


def _check_run(run):
    if _manifest_stat(run.output_dir / "run_manifest.json") != run.manifest_stat:
        raise ValueError("incompatible run: immutable manifest changed")


def _atomic_json(path, value):
    fd, temp = tempfile.mkstemp(prefix=path.stem + ".", suffix=".tmp", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(_json(value) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temp, path)
        _fsync_dir(path.parent)
    finally:
        Path(temp).unlink(missing_ok=True)


def prepare_run(output_dir, worker, identities, shards, *, max_mismatch=3,
                window_size=None, region=None, shard_count=1):
    """Atomically seal the exact expected shard set; incompatible reuse is rejected.

    Legacy directories without a run manifest are rejected, never silently adopted.
    Choose a fresh output directory rather than deleting/mixing prior generations.
    Worker index is deliberately excluded: all array workers share one generation.
    """
    _budget(max_mismatch)
    if not {"fasta", "gtf"} <= identities.keys():
        raise ValueError("both FASTA and GTF identities are required")
    _check_identities(identities)
    shard_ids = sorted(s.shard_id for s in shards)
    if len(shard_ids) != len(set(shard_ids)):
        raise ValueError("duplicate expected shards")
    config = dict(worker=json.loads(worker.params_json), identities=identities,
                  max_mismatch=max_mismatch, window_size=window_size, region=region,
                  shard_count=shard_count)
    manifest = dict(schema_version=1, config=config, expected_shards=shard_ids)
    manifest_json = _json(manifest) + "\n"
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / "run_manifest.json"
    with (output_dir / "run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if path.exists():
            if path.read_text() != manifest_json:
                raise ValueError("incompatible run configuration; use a fresh output directory")
        else:
            if any(output_dir.glob("*.done")) or any(output_dir.glob("*.parquet")):
                raise ValueError("incompatible run: legacy shards lack an immutable run manifest")
            _atomic_json(path, manifest)
        allowed = set(shard_ids)
        for suffix in ("*.done", "*.parquet"):
            if any(p.stem not in allowed for p in output_dir.glob(suffix)):
                raise ValueError("incompatible run: unexpected shard files")
    return RunConfig(output_dir, manifest_json, worker.params_json, _json(identities),
                     frozenset(shard_ids), max_mismatch, hashlib.sha256(manifest_json.encode()).hexdigest(),
                     _manifest_stat(path))


def _valid_done(manifest, *, expected=None):
    """Validate JSON types before any field access (bool is not an integer count)."""
    if not isinstance(manifest, dict):
        return False
    if type(manifest.get("rows")) is not int or manifest["rows"] < 0:
        return False
    sha = manifest.get("sha256")
    if not isinstance(sha, str) or re.fullmatch(r"[0-9a-f]{64}", sha) is None:
        return False
    return expected is None or all(manifest.get(k) == v for k, v in expected.items())


def run_shard(*, shard, output_dir, identities: Mapping, genome=None, genes=None, asos=None,
              max_mismatch=3, synthetic_reference=False, worker=None, run=None) -> Path:
    """Commit/reuse one listed shard, writing bounded Arrow batches atomically.

    Standalone callers get a single-shard run. Multi-shard callers must prepare_run
    once with the complete shard list and reuse the same worker and RunConfig.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq
    _budget(max_mismatch)
    if worker is None:
        worker = SiteWorker(genome, genes, asos, synthetic_reference)
    elif any(value is not None for value in (genome, genes, asos)):
        raise ValueError("pass a worker or genome/genes/asos, not both")
    _check_identities(identities)
    if run is None:
        run = prepare_run(output_dir, worker, identities, [shard], max_mismatch=max_mismatch)
    if (run.output_dir != Path(output_dir).resolve() or run.worker_params_json != worker.params_json or
            run.identities_json != _json(identities) or run.max_mismatch != max_mismatch or
            shard.shard_id not in run.expected_shards):
        raise ValueError("incompatible run/shard configuration")
    _check_run(run)
    output_dir = run.output_dir
    path = output_dir / (shard.shard_id + ".parquet")
    done = path.with_suffix(".done")
    params = dict(json.loads(worker.params_json), shard=asdict(shard), max_mismatch=max_mismatch)
    expected = dict(params=params, identities=json.loads(run.identities_json),
                    run_sha256=run.manifest_sha256)
    with path.with_suffix(".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        _check_identities(identities)
        try:
            manifest = json.loads(done.read_text())
            if (_valid_done(manifest, expected=expected) and manifest["sha256"] == _sha256(path) and
                    manifest["rows"] == pq.read_metadata(path).num_rows and
                    pq.read_schema(path).equals(_schema())):
                return path
        except (OSError, ValueError, KeyError, TypeError, pa.ArrowException):
            pass
        done.unlink(missing_ok=True)
        _fsync_dir(output_dir)
        fd, temp = tempfile.mkstemp(prefix=shard.shard_id + ".", suffix=".tmp", dir=output_dir)
        os.close(fd)
        try:
            stream = worker.iter_sites(max_mismatch=max_mismatch, shard=shard)
            rows = 0
            with pq.ParquetWriter(temp, _schema(), compression="zstd") as writer:
                while batch := list(islice(stream, 512)):
                    table = pa.Table.from_pylist([asdict(s) for s in batch], schema=_schema())
                    writer.write_table(table)
                    rows += table.num_rows
            _check_identities(identities)
            _check_run(run)
            with open(temp, "rb") as handle:
                os.fsync(handle.fileno())
            sha = _sha256(temp)
            os.replace(temp, path)
            _fsync_dir(output_dir)
            _check_identities(identities)
            _check_run(run)
            _atomic_json(done, dict(**expected, rows=rows, sha256=sha))
        finally:
            Path(temp).unlink(missing_ok=True)
    return path


def _region(text):
    match = re.fullmatch(r"([A-Za-z0-9_.-]+):(\d+)-(\d+)", text)
    if not match:
        raise argparse.ArgumentTypeError("region must be chrom:start1-end1")
    chrom, start, end = match.groups()
    try:
        return_tuple = chrom, int(start), int(end)
        Shard(*return_tuple)
        return return_tuple
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc


def check_smn_anchors(output_dir, *, paths=None) -> int:
    """Sentinel acceptance: six exact F18/F14/F10 records, two distinct genes.

    Inspect only committed shards in this specific output directory; verify each
    checksum before reading. No new genomic search or reference traversal.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq
    expected = {(gene, aso): (start, start + length - 1, target)
                for gene, start in (("SMN2", 70076584), ("SMN1", 70952004))
                for aso, length, target in (("F18", 18, "CCAGCATTATGAAAGTGA"),
                                            ("F14", 14, "CCAGCATTATGAAA"),
                                            ("F10", 10, "CCAGCATTAT"))}
    found = {}
    markers = (Path(output_dir).glob("*.done") if paths is None else
               (Path(path).with_suffix(".done") for path in paths))
    for done in sorted(markers):
        path = done.with_suffix(".parquet")
        try:
            manifest = json.loads(done.read_text())
            if not _valid_done(manifest):
                raise ValueError(f"invalid sentinel commit schema: {path}")
            if manifest["sha256"] != _sha256(path):
                raise ValueError(f"sentinel checksum mismatch: {path}")
            if manifest["rows"] != pq.read_metadata(path).num_rows:
                raise ValueError(f"sentinel row-count mismatch: {path}")
            rows = pq.read_table(path, filters=[("start1", "in", [70076584, 70952004])]).to_pylist()
        except (OSError, KeyError, TypeError, pa.ArrowException) as exc:
            raise ValueError(f"invalid committed sentinel shard: {path}") from exc
        for row in rows:
            key = row["gene_name"], row["aso_id"]
            if key not in expected or row["start1"] != expected[key][0]:
                continue
            if (row["chrom"] != "chr5" or row["strand"] != "+" or
                    (row["start1"], row["end1"], row["target_sequence"]) != expected[key] or
                    row["mismatch_count"] != 0 or row["mismatch_positions"] != [] or
                    row["paralog_ambiguity"] is not True or row["evidence_class"] != "derived"):
                raise ValueError(f"invalid sentinel anchor: {key}")
            if key in found:
                raise ValueError(f"duplicate sentinel anchor: {key}")
            found[key] = row
    if set(found) != set(expected):
        raise ValueError(f"missing sentinel anchors: {sorted(set(expected) - set(found))}")
    if len({r["site_id"] for r in found.values()}) != 6:
        raise ValueError("sentinel site IDs collapse distinct ASOs/paralogs")
    if {found[("SMN1", a)]["gene_id"] for a in ("F18", "F14", "F10")} & {
            found[("SMN2", a)]["gene_id"] for a in ("F18", "F14", "F10")}:
        raise ValueError("sentinel gene IDs collapse paralogs")
    return 6


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", default=DEFAULT_FASTA)
    parser.add_argument("--gtf", default=DEFAULT_GTF)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--region", type=_region)
    parser.add_argument("--asos", nargs="+", choices=sorted(ASOS), default=["F18", "F14", "F10"])
    parser.add_argument("--max-mismatch", type=int, default=3)
    parser.add_argument("--window-size", type=int, default=DEFAULT_WINDOW)
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    parser.add_argument("--check-smn-anchors", action="store_true")
    parser.add_argument("--synthetic-reference", action="store_true",
                        help="test-only reference: synthetic assembly/IDs and isolated run configuration")
    parser.add_argument("--profile-memory", action="store_true", help="emit peak process RSS in KiB at stages")
    args = parser.parse_args(argv)
    _budget(args.max_mismatch)
    if not 0 <= args.shard_index < args.shard_count:
        parser.error("require 0 <= shard-index < shard-count")
    if args.check_smn_anchors and args.shard_count != 1:
        parser.error("sentinel acceptance requires a single complete shard worker")
    if args.synthetic_reference and args.check_smn_anchors:
        parser.error("synthetic references cannot pass production sentinel acceptance")
    def memory(stage):
        if args.profile_memory:
            import resource
            import sys
            print(f"peak_rss_kib stage={stage} value={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}",
                  file=sys.stderr, flush=True)
    memory("startup")
    genome = FaiFastaReader(args.fasta)
    reference_identity(genome, synthetic_reference=args.synthetic_reference)
    identities = dict(fasta=file_identity(genome.path), gtf=file_identity(args.gtf))
    genes = read_gene_bodies(args.gtf)  # Retain both SMN paralogs regardless of output region.
    if file_identity(args.gtf) != identities["gtf"]:
        raise ValueError("GTF changed while reading")
    if not genes:
        parser.error("no gene bodies overlap the requested universe")
    asos = [ASOS[name] for name in args.asos]
    shards = make_shards(genes, window_size=args.window_size, region=args.region)
    if not shards:
        parser.error("no gene bodies overlap the requested universe")
    worker = SiteWorker(genome, genes, asos, args.synthetic_reference)
    run = prepare_run(args.output_dir, worker, identities, shards, max_mismatch=args.max_mismatch,
                      window_size=args.window_size, region=args.region, shard_count=args.shard_count)
    memory("prepared")
    completed = []
    for index, shard in enumerate(shards):
        if index % args.shard_count != args.shard_index:
            continue
        if (file_identity(genome.path) != identities["fasta"] or
                file_identity(args.gtf) != identities["gtf"]):
            raise ValueError("reference input identity changed during run")
        path = run_shard(worker=worker, run=run, shard=shard,
                         output_dir=args.output_dir, identities=identities,
                         max_mismatch=args.max_mismatch)
        completed.append(path)
        print(path, flush=True)
        memory(shard.shard_id)
    if args.check_smn_anchors:
        print(f"verified_smn_anchors={check_smn_anchors(args.output_dir, paths=completed)}", flush=True)
    memory("finished")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
