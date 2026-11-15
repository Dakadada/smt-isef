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

Outputs are explicit-schema Parquet, then an atomic .done JSON commit marker.
Restart verifies parameters, input identities, row count, and Parquet SHA256.
Input identities use resolved path/size/mtime_ns/inode plus FAI SHA256; they are
filesystem identities, NOT full FASTA/GTF content hashes. Inputs must be immutable
during a run. No full reference hashing or unbounded reference fetch is performed.

CLI (from repo root, with PYTHONPATH=feng/src):
  python -m conference.model_v2.j6_site_universe --output-dir OUT [--region chr5:70000000-71100000]
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import asdict, dataclass
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
from typing import Mapping, Sequence

from conference.model_v2.trans_intervention import (
    ANNOTATION, ASSEMBLY, ASOS, ASO, DEFAULT_FASTA, FaiFastaReader,
    MAX_SLICE_BASES, reverse_complement,
)
from feng.scan.seed_extend import seed_and_extend

VERSION = "1.0.0"
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
    paralog_ambiguity: bool
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


def read_gene_bodies(path, *, region=None) -> tuple[Gene, ...]:
    """Stream one known GTF file; parse only gene lines, retain full gene bounds.

    Header must explicitly identify GRCh38 and version 50. No feature/transcript
    inference is attempted. Optional region bounds filter overlapping genes.
    """
    if region is not None:
        Shard(*region)
    genes = []
    valid_header = False
    with Path(path).expanduser().open() as handle:
        for line in handle:
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


def enumerate_sites(genome, genes: Sequence[Gene], asos: Sequence[ASO], *,
                    max_mismatch=3, shard: Shard | None = None) -> list[Site]:
    """Complete Hamming enumeration, using Feng seed-and-extend read-only.

    Memory is bounded by a chromosome window plus output rows for that window.
    Unsharded calls internally use the same windows, retaining identical site IDs.
    """
    _budget(max_mismatch)
    genes = _genes(genes)
    asos = tuple(sorted(asos, key=lambda a: a.aso_id))
    if len({a.aso_id for a in asos}) != len(asos):
        raise ValueError("duplicate ASO IDs")
    if not asos:
        raise ValueError("at least one ASO is required")
    if shard is None:
        return sorted([site for window in make_shards(genes)
                       for site in enumerate_sites(genome, genes, asos,
                                                   max_mismatch=max_mismatch, shard=window)],
                      key=site_sort_key)
    selected = [g for g in genes if g.chrom == shard.chrom and
                g.start1 <= shard.end1 and g.end1 >= shard.start1]
    if not selected:
        return []
    length = genome.length(shard.chrom)
    if any(g.end1 > length for g in selected):
        raise ValueError("gene outside reference chromosome")
    halo_end = min(length, shard.end1 + max(a.length for a in asos) - 1)
    if shard.start1 > length:
        return []
    sequence = _fetch(genome, shard.chrom, shard.start1, halo_end)
    sites = []
    paralog_sequences = {}
    for gene in selected:
        start = max(gene.start1, shard.start1)
        end = min(gene.end1, halo_end)
        reference = sequence[start - shard.start1:end - shard.start1 + 1]
        for aso in asos:
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
                if gene.gene_name in ("SMN1", "SMN2"):
                    other_name = "SMN1" if gene.gene_name == "SMN2" else "SMN2"
                    for other in genes:
                        if other.gene_name != other_name:
                            continue
                        if other.gene_id not in paralog_sequences:
                            seq = _fetch(genome, other.chrom, other.start1, other.end1)
                            paralog_sequences[other.gene_id] = seq if other.strand == "+" else reverse_complement(seq)
                        if target in paralog_sequences[other.gene_id]:
                            ambiguous = True
                            break
                region, regions, transcripts = _annotation(gene, lo, hi)
                identity = [ASSEMBLY, ANNOTATION, aso.aso_id, aso.sequence_5to3,
                            gene.gene_id, gene.chrom, lo, hi, gene.strand, target]
                sites.append(Site("j6-" + _digest(identity), aso.aso_id, gene.chrom, lo, hi,
                                  gene.strand, gene.gene_id, gene.gene_name, region, regions,
                                  transcripts, hit.mismatches, positions, target, ambiguous))
    return sorted(sites, key=site_sort_key)


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


def run_shard(*, genome, genes, asos, shard, output_dir, identities: Mapping,
              max_mismatch=3) -> Path:
    """Commit/reuse one shard; .done is the sole completeness marker.

    Per-shard advisory locks serialize concurrent retries. Orphan temporary files
    are ignored; absent/invalid markers or changed parameters force a rebuild.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq
    _budget(max_mismatch)
    if not {"fasta", "gtf"} <= identities.keys():
        raise ValueError("both FASTA and GTF identities are required")
    genes = _genes(genes)
    asos = tuple(sorted(asos, key=lambda a: a.aso_id))
    if not asos or len({a.aso_id for a in asos}) != len(asos):
        raise ValueError("distinct ASO IDs required")
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / (shard.shard_id + ".parquet")
    done = path.with_suffix(".done")
    import feng.scan.seed_extend as scanner
    params = dict(version=VERSION, assembly=ASSEMBLY, annotation=ANNOTATION,
                  shard=asdict(shard), max_mismatch=max_mismatch, mode="substitutions",
                  asos=[dict(aso_id=a.aso_id, sequence=a.sequence_5to3) for a in asos],
                  gene_universe_sha256=_digest([asdict(g) for g in genes]),
                  implementation_sha256=_sha256(__file__), scanner_sha256=_sha256(scanner.__file__))
    expected = json.loads(_json(dict(params=params, identities=identities)))
    with path.with_suffix(".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            manifest = json.loads(done.read_text())
            if (all(manifest.get(k) == v for k, v in expected.items()) and
                    manifest["sha256"] == _sha256(path) and
                    manifest["rows"] == pq.read_metadata(path).num_rows):
                return path
        except (OSError, ValueError, KeyError, TypeError, pa.ArrowException):
            pass
        # Invalidate the old commit marker BEFORE replacing its data file.
        done.unlink(missing_ok=True)
        _fsync_dir(output_dir)
        sites = enumerate_sites(genome, genes, asos, max_mismatch=max_mismatch, shard=shard)
        table = pa.Table.from_pylist([asdict(s) for s in sites], schema=_schema())
        temps = []
        try:
            fd, temp = tempfile.mkstemp(prefix=shard.shard_id + ".", suffix=".tmp", dir=output_dir)
            os.close(fd)
            temps.append(Path(temp))
            pq.write_table(table, temp, compression="zstd")
            with open(temp, "rb") as handle:
                os.fsync(handle.fileno())
            sha = _sha256(temp)
            os.replace(temp, path)
            _fsync_dir(output_dir)
            manifest = dict(**expected, rows=table.num_rows, sha256=sha)
            fd, temp = tempfile.mkstemp(prefix=shard.shard_id + ".", suffix=".tmp", dir=output_dir)
            temps.append(Path(temp))
            with os.fdopen(fd, "w") as handle:
                handle.write(_json(manifest) + "\n")
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temp, done)
            _fsync_dir(output_dir)
        finally:
            for temp in temps:
                temp.unlink(missing_ok=True)
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
    args = parser.parse_args(argv)
    _budget(args.max_mismatch)
    if not 0 <= args.shard_index < args.shard_count:
        parser.error("require 0 <= shard-index < shard-count")
    genome = FaiFastaReader(args.fasta)
    identities = dict(fasta=file_identity(genome.path), gtf=file_identity(args.gtf))
    genes = read_gene_bodies(args.gtf, region=args.region)
    if file_identity(args.gtf) != identities["gtf"]:
        raise ValueError("GTF changed while reading")
    if not genes:
        parser.error("no gene bodies overlap the requested universe")
    asos = [ASOS[name] for name in args.asos]
    shards = make_shards(genes, window_size=args.window_size, region=args.region)
    for index, shard in enumerate(shards):
        if index % args.shard_count != args.shard_index:
            continue
        if (file_identity(genome.path) != identities["fasta"] or
                file_identity(args.gtf) != identities["gtf"]):
            raise ValueError("reference input identity changed during run")
        path = run_shard(genome=genome, genes=genes, asos=asos, shard=shard,
                         output_dir=args.output_dir, identities=identities,
                         max_mismatch=args.max_mismatch)
        print(path, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
