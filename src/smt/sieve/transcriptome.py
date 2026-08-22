"""A5 step 1: transcriptome acquisition and exon/intron structural index.

Data source and pin
--------------------
GENCODE human release 44 (genome build GRCh38.p14), "basic" annotation set,
fetched from the EBI mirror. Two files, both pinned by URL + MD5 (GENCODE
publishes MD5, not SHA256, in its release MD5SUMS manifest, so we record both
theirs and our own SHA256 computed on download -- the same "pin and checksum"
discipline as :mod:`smt.seq.reference`, adapted to the checksum algorithm the
upstream actually provides):

  gencode.v44.basic.annotation.gtf.gz
      https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/
      MD5 (GENCODE) 7450ef42cf9cb3d29625320b22d4bb45
      29,570,410 bytes compressed.

  GRCh38.primary_assembly.genome.fa.gz
      same directory
      MD5 (GENCODE) 9c3fc2ca260a767530dddb0f26721a6b
      844,691,642 bytes compressed (~806 MiB); 3,151,417,447 bytes (~2.94 GiB)
      decompressed. Documented here BEFORE fetching per project policy: this
      is a bounded, one-time bulk download, not an unscoped scan of anything.

Why GENCODE over RefSeq: GENCODE's GTF gives exon/UTR/CDS coordinates in one
self-consistent, transcript-hierarchical file (gene -> transcript -> exon in
file order) that is trivial to stream in one pass, and "basic" is itself a
pre-filtered representative subset (not the full ~250k-transcript
comprehensive set), which keeps the index a manageable size without an extra
filtering pass on our side. RefSeq would work equally well scientifically;
GENCODE was picked for the file-format convenience, not a biological reason,
and that is flagged here explicitly so it isn't mistaken for one.

Why the *primary_assembly* genome FASTA and not just GENCODE's pre-built
transcript FASTA: the transcript FASTA GENCODE ships contains only *mature*
(spliced) mRNA -- no introns. This project's whole point is that the same
motif is a silencer in an intron and an enhancer in an exon, so intronic
sequence is not optional. Only the genome FASTA has it.

Scope decisions (assumptions, not measurements) -- read before using the index
--------------------------------------------------------------------------
1. ONE TRANSCRIPT PER GENE. We keep only transcripts tagged
   ``Ensembl_canonical`` (GENCODE's one-per-gene representative isoform,
   MANE Select where one exists) among ``gene_type protein_coding`` genes.
   Non-coding biotypes and alternative isoforms are OUT of this index. This is
   a tractability choice: it turns "every isoform of every gene" (order
   100,000+ transcripts in the comprehensive set) into ~20,000 canonical
   transcripts, which is what makes k-mer seeding and per-candidate
   thermodynamics affordable on a shared cluster. A real off-target could in
   principle sit on a non-canonical isoform only; this index will miss it.
2. INTRONS ARE REPRESENTED ONLY NEAR THEIR SPLICE SITES, not in full. For each
   intron we keep up to ``flank_len`` nt (default 200, see ``build_index``)
   from the donor (5') end and up to ``flank_len`` nt from the acceptor (3')
   end; if the intron is shorter than ``2 * flank_len`` we keep the whole
   thing as one unit instead of two overlapping ones. Human introns range from
   tens of nt to over a megabase (some genes span >1 Mb); storing full introns
   transcriptome-wide is not tractable here, and every documented intronic
   splicing regulatory element this project is aware of (ISS-N1 itself sits at
   intron 7 +10..+24 of SMN2, i.e. within 24 nt of the 5'ss) is near a splice
   site. This is a MODELING ASSUMPTION: an intronic off-target site deep in an
   intron's interior (hundreds of nt from either splice site) is invisible to
   this index by construction, not ruled out by any evidence.
3. Only the 25 canonical chromosomes/contigs (chr1-22, chrX, chrY, chrM) are
   indexed; alternate haplotypes and unplaced scaffolds in the primary
   assembly FASTA are skipped. GENCODE "basic" annotation is overwhelmingly on
   these already.

Coordinate conventions
-----------------------
Genomic coordinates are 1-based inclusive on the pinned GRCh38.p14 FASTA
records, exactly like the GTF's own convention and like
:mod:`smt.seq.reference`. Transcript-space (mRNA) coordinates are 1-based
inclusive on the spliced mRNA sequence, 5'->3'. Flank sequences are stored
5'->3' in transcript orientation (i.e. already reverse-complemented for
minus-strand genes), with their genomic span recorded separately so a flank
hit can be mapped back to a genomic position.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Optional

DATA = Path("/home/dashii/orcd/pool/smt-isef-data/transcriptome")

GENCODE_RELEASE = "44"
GENCODE_BASE_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44"
)

GTF_NAME = "gencode.v44.basic.annotation.gtf.gz"
GTF_URL = f"{GENCODE_BASE_URL}/{GTF_NAME}"
GTF_PATH = DATA / GTF_NAME
GTF_MD5 = "7450ef42cf9cb3d29625320b22d4bb45"  # from GENCODE's own MD5SUMS manifest
GTF_SIZE_BYTES = 29_570_410

FASTA_NAME = "GRCh38.primary_assembly.genome.fa.gz"
FASTA_URL = f"{GENCODE_BASE_URL}/{FASTA_NAME}"
FASTA_GZ_PATH = DATA / FASTA_NAME
FASTA_PATH = DATA / "GRCh38.primary_assembly.genome.fa"  # decompressed, for faidx
FASTA_MD5 = "9c3fc2ca260a767530dddb0f26721a6b"
FASTA_GZ_SIZE_BYTES = 844_691_642
FASTA_SIZE_BYTES = 3_151_417_447

# Canonical chromosomes/contigs to index. Everything else in the primary
# assembly (unplaced scaffolds, alt haplotypes) is skipped -- see module
# docstring scope decision 3.
CANONICAL_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}

INDEX_PATH = DATA / "transcripts.jsonl"
FAI_PATH = DATA / "GRCh38.primary_assembly.genome.fa.fai"

#: Default intronic flank length, nt, from each splice site. See module
#: docstring scope decision 2. This is a tractability/coverage assumption,
#: not a literature-measured value.
DEFAULT_FLANK_LEN = 200

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    """Reverse complement, DNA alphabet, case-preserving."""
    return seq.translate(_COMPLEMENT)[::-1]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def md5_file(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# fetch (idempotent, checksum-verified)
# ---------------------------------------------------------------------------


def _fetch_one(url: str, dest: Path, expected_md5: str, expected_size: int) -> None:
    if dest.exists() and dest.stat().st_size == expected_size:
        got = md5_file(dest)
        if got != expected_md5:
            raise ValueError(
                f"{dest} exists with the pinned size but MD5 {got} != expected "
                f"{expected_md5}; corrupt or a different file, re-fetch by hand."
            )
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".tmp")
    import urllib.request

    with urllib.request.urlopen(url, timeout=600) as resp, tmp.open("wb") as out:
        while True:
            chunk = resp.read(1 << 20)
            if not chunk:
                break
            out.write(chunk)
    got_size = tmp.stat().st_size
    got_md5 = md5_file(tmp)
    if got_size != expected_size or got_md5 != expected_md5:
        tmp.unlink(missing_ok=True)
        raise ValueError(
            f"download of {url} failed verification: size {got_size} "
            f"(expected {expected_size}), md5 {got_md5} (expected {expected_md5})"
        )
    os.replace(tmp, dest)  # atomic


def fetch_gtf() -> Path:
    """Download and MD5-verify the pinned GTF, skipping if already present."""
    _fetch_one(GTF_URL, GTF_PATH, GTF_MD5, GTF_SIZE_BYTES)
    return GTF_PATH


def fetch_genome() -> Path:
    """Download+verify the pinned genome FASTA, then decompress for faidx.

    Decompression is required (not just convenience): plain gzip does not
    support random-access seeking, and the whole point of this file is
    O(1) slicing of arbitrary genomic intervals via :class:`GenomeIndex`.
    """
    _fetch_one(FASTA_URL, FASTA_GZ_PATH, FASTA_MD5, FASTA_GZ_SIZE_BYTES)
    if not FASTA_PATH.exists() or FASTA_PATH.stat().st_size != FASTA_SIZE_BYTES:
        tmp = FASTA_PATH.with_suffix(".fa.tmp")
        with gzip.open(FASTA_GZ_PATH, "rb") as src, tmp.open("wb") as out:
            while True:
                chunk = src.read(1 << 20)
                if not chunk:
                    break
                out.write(chunk)
        os.replace(tmp, FASTA_PATH)
    return FASTA_PATH


# ---------------------------------------------------------------------------
# minimal faidx: O(1) random access into a plain (non-bgzipped) FASTA
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class _FaiRecord:
    length: int
    offset: int  # byte offset of first sequence base
    linebases: int  # bases per line (all but possibly the last line)
    linewidth: int  # bytes per line including the newline


def build_fai(fasta_path: Path) -> dict[str, _FaiRecord]:
    """Build a samtools-.fai-equivalent index by one linear scan.

    Assumes (as samtools faidx does) that every line of a record has the same
    length except the last. GENCODE/Ensembl genome FASTAs satisfy this.
    """
    fai: dict[str, _FaiRecord] = {}
    name: Optional[str] = None
    seqlen = 0
    offset = 0
    linebases: Optional[int] = None
    linewidth: Optional[int] = None
    with fasta_path.open("rb") as fh:
        while True:
            line = fh.readline()
            if not line:
                if name is not None:
                    fai[name] = _FaiRecord(seqlen, offset, linebases or 0, linewidth or 0)
                break
            if line.startswith(b">"):
                if name is not None:
                    fai[name] = _FaiRecord(seqlen, offset, linebases or 0, linewidth or 0)
                name = line[1:].split()[0].decode()
                seqlen = 0
                linebases = None
                linewidth = None
                offset = fh.tell()
            else:
                l = len(line.rstrip(b"\r\n"))
                if linebases is None:
                    linebases = l
                    linewidth = len(line)
                seqlen += l
    return fai


def load_or_build_fai(fasta_path: Path = FASTA_PATH, fai_path: Path = FAI_PATH) -> dict[str, _FaiRecord]:
    if fai_path.exists():
        fai: dict[str, _FaiRecord] = {}
        with fai_path.open() as fh:
            for line in fh:
                name, length, offset, linebases, linewidth = line.rstrip("\n").split("\t")
                fai[name] = _FaiRecord(int(length), int(offset), int(linebases), int(linewidth))
        return fai
    fai = build_fai(fasta_path)
    tmp = fai_path.with_suffix(".fai.tmp")
    with tmp.open("w") as fh:
        for name, r in fai.items():
            fh.write(f"{name}\t{r.length}\t{r.offset}\t{r.linebases}\t{r.linewidth}\n")
    os.replace(tmp, fai_path)
    return fai


class GenomeIndex:
    """Random-access reader over the decompressed pinned genome FASTA."""

    def __init__(self, fasta_path: Path = FASTA_PATH, fai_path: Path = FAI_PATH):
        self.fasta_path = fasta_path
        self.fai = load_or_build_fai(fasta_path, fai_path)
        self._fh = fasta_path.open("rb")

    def close(self) -> None:
        self._fh.close()

    def __enter__(self) -> "GenomeIndex":
        return self

    def __exit__(self, *exc) -> None:
        self.close()

    def fetch(self, chrom: str, start: int, end: int) -> str:
        """1-based inclusive slice of ``chrom``, uppercase DNA."""
        rec = self.fai[chrom]
        if start < 1 or end > rec.length or start > end:
            raise ValueError(f"{chrom}:{start}-{end} outside record of length {rec.length}")
        start0 = start - 1
        n = end - start0
        line_gap = rec.linewidth - rec.linebases
        start_line = start0 // rec.linebases
        start_col = start0 % rec.linebases
        file_pos = rec.offset + start_line * rec.linewidth + start_col
        # Read generously past `n` bytes to account for embedded newlines,
        # then strip them. +2 lines of slack covers boundary rounding.
        extra_lines = n // rec.linebases + 2
        to_read = n + extra_lines * line_gap
        self._fh.seek(file_pos)
        raw = self._fh.read(to_read)
        seq = raw.replace(b"\n", b"").replace(b"\r", b"")[:n]
        if len(seq) != n:
            raise ValueError(f"short read for {chrom}:{start}-{end}: got {len(seq)} of {n}")
        return seq.decode("ascii").upper()


# ---------------------------------------------------------------------------
# GTF parsing
# ---------------------------------------------------------------------------


def _parse_attrs(field: str) -> dict[str, object]:
    attrs: dict[str, object] = {}
    for part in field.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        key, _, val = part.partition(" ")
        val = val.strip().strip('"')
        if key == "tag":
            attrs.setdefault("tag", []).append(val)
        else:
            attrs[key] = val
    return attrs


@dataclass
class _RawTranscript:
    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    biotype: str
    exon_spans: list[tuple[int, int]] = field(default_factory=list)  # genomic, 1-based incl
    utr_spans: list[tuple[int, int]] = field(default_factory=list)
    cds_spans: list[tuple[int, int]] = field(default_factory=list)


def iter_canonical_protein_coding_transcripts(gtf_path: Path = GTF_PATH) -> Iterator[_RawTranscript]:
    """Stream the GTF once, yielding one ``_RawTranscript`` per qualifying gene.

    Qualifying = ``gene_type protein_coding`` and transcript tagged
    ``Ensembl_canonical`` (GENCODE's one-per-gene representative isoform) on
    a canonical chromosome. See module docstring scope decision 1.

    GENCODE GTF is hierarchically ordered (gene, then its transcript, then
    that transcript's exon/CDS/UTR lines, in file order), so a single forward
    pass with a "currently tracked transcript_id" is sufficient -- no need to
    buffer the whole file.
    """
    current: Optional[_RawTranscript] = None
    tracked_id: Optional[str] = None
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, _src, feature, start, end, _score, strand, _frame, attr_field = f
            if chrom not in CANONICAL_CHROMS:
                continue
            if feature == "transcript":
                if current is not None:
                    yield current
                    current = None
                attrs = _parse_attrs(attr_field)
                tags = attrs.get("tag", [])
                if (
                    attrs.get("gene_type") == "protein_coding"
                    and "Ensembl_canonical" in tags
                ):
                    current = _RawTranscript(
                        transcript_id=attrs["transcript_id"],
                        gene_id=attrs["gene_id"],
                        gene_name=attrs.get("gene_name", attrs["gene_id"]),
                        chrom=chrom,
                        strand=strand,
                        biotype=attrs.get("gene_type", "unknown"),
                    )
                    tracked_id = current.transcript_id
                else:
                    tracked_id = None
            elif feature in ("exon", "UTR", "CDS") and current is not None:
                attrs = _parse_attrs(attr_field)
                if attrs.get("transcript_id") != tracked_id:
                    continue
                span = (int(start), int(end))
                if feature == "exon":
                    current.exon_spans.append(span)
                elif feature == "UTR":
                    current.utr_spans.append(span)
                else:
                    current.cds_spans.append(span)
        if current is not None:
            yield current


# ---------------------------------------------------------------------------
# structural index: exon/intron in both transcript AND genomic coordinates
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Exon:
    order: int  # 1-based, transcript 5'->3' order
    g_start: int  # genomic, 1-based inclusive, g_start <= g_end always
    g_end: int
    t_start: int  # transcript-space (mRNA), 1-based inclusive
    t_end: int


@dataclass(frozen=True)
class IntronFlank:
    """One scannable window near a splice site (or a whole short intron).

    ``side`` is ``"donor"`` (5' end of the intron, immediately after the
    upstream exon), ``"acceptor"`` (3' end, immediately before the downstream
    exon), or ``"whole"`` when the intron is shorter than ``2 * flank_len``
    and is kept in full instead of as two overlapping pieces.
    """

    intron_index: int  # 1-based, between exon[i] and exon[i+1] in transcript order
    side: str
    g_start: int  # genomic span of THIS flank, 1-based inclusive
    g_end: int
    seq: str  # transcript-oriented 5'->3' (revcomp'd already if strand '-')
    intron_len: int  # length of the FULL intron this flank belongs to

    def dist_to_splice_site(self, local_pos: int, span: int = 1) -> int:
        """Nearest-splice-site distance (nt) for a ``span``-nt window starting at
        ``local_pos`` (1-based, within :attr:`seq`, transcript orientation).

        ``side="donor"``: the donor (5') splice site sits immediately before
        this flank (at intron position 0); distance grows with ``local_pos``.
        ``side="acceptor"``: the acceptor (3') splice site sits immediately
        after this flank's last base; distance grows as ``local_pos`` moves
        away from the end.
        ``side="whole"``: the intron has splice sites at BOTH ends (it is
        shorter than ``2 * flank_len``); distance is to whichever is nearer.
        """
        n = len(self.seq)
        if self.side == "donor":
            return local_pos  # 1 == the base right after the 5'ss
        if self.side == "acceptor":
            return n - (local_pos + span - 1) + 1  # 1 == the base right before the 3'ss
        # "whole": distance to the nearer end.
        d_donor = local_pos
        d_acceptor = n - (local_pos + span - 1) + 1
        return min(d_donor, d_acceptor)


@dataclass(frozen=True)
class Transcript:
    transcript_id: str
    gene_id: str
    gene_name: str
    chrom: str
    strand: str
    biotype: str
    exons: tuple[Exon, ...]  # transcript order
    five_utr: tuple[tuple[int, int], ...]  # genomic spans
    three_utr: tuple[tuple[int, int], ...]  # genomic spans
    mrna_seq: str  # spliced, 5'->3'
    flanks: tuple[IntronFlank, ...]  # transcript order, donor then acceptor per intron

    def mrna_pos_to_genomic(self, t_pos: int) -> int:
        """Transcript-space 1-based position -> genomic 1-based position."""
        for ex in self.exons:
            if ex.t_start <= t_pos <= ex.t_end:
                offset = t_pos - ex.t_start
                return ex.g_start + offset if self.strand == "+" else ex.g_end - offset
        raise ValueError(f"{t_pos} outside {self.transcript_id} mRNA (len {len(self.mrna_seq)})")

    def flank_pos_to_genomic(self, flank: IntronFlank, local_pos: int) -> int:
        """Flank-local 1-based (transcript-oriented) position -> genomic 1-based position."""
        if self.strand == "+":
            return flank.g_start + (local_pos - 1)
        return flank.g_end - (local_pos - 1)

    def genomic_region_kind(self, g_pos: int) -> str:
        """"five_utr" / "three_utr" / "exonic" / "outside_exon" for a genomic pos in an exon."""
        for s, e in self.five_utr:
            if s <= g_pos <= e:
                return "five_utr"
        for s, e in self.three_utr:
            if s <= g_pos <= e:
                return "three_utr"
        for ex in self.exons:
            if ex.g_start <= g_pos <= ex.g_end:
                return "exonic"
        return "outside_exon"


def _build_transcript(raw: _RawTranscript, genome: GenomeIndex, flank_len: int) -> Transcript:
    spans = sorted(raw.exon_spans)  # ascending genomic order
    if raw.strand == "-":
        ordered = list(reversed(spans))  # 5'->3' transcript order
    else:
        ordered = spans

    exons: list[Exon] = []
    t_cursor = 1
    pieces: list[str] = []
    for i, (g_start, g_end) in enumerate(ordered, start=1):
        raw_seq = genome.fetch(raw.chrom, g_start, g_end)
        seq = raw_seq if raw.strand == "+" else revcomp(raw_seq)
        pieces.append(seq)
        t_start = t_cursor
        t_end = t_cursor + (g_end - g_start + 1) - 1
        exons.append(Exon(i, g_start, g_end, t_start, t_end))
        t_cursor = t_end + 1
    mrna_seq = "".join(pieces)

    flanks: list[IntronFlank] = []
    for i in range(len(ordered) - 1):
        exon_a = ordered[i]
        exon_b = ordered[i + 1]
        if raw.strand == "+":
            intron_g_start = exon_a[1] + 1
            intron_g_end = exon_b[0] - 1
        else:
            intron_g_start = exon_b[1] + 1
            intron_g_end = exon_a[0] - 1
        intron_len = intron_g_end - intron_g_start + 1
        if intron_len <= 0:
            continue  # adjacent/overlapping exon annotations; skip defensively
        idx = i + 1
        if intron_len <= 2 * flank_len:
            g_s, g_e = intron_g_start, intron_g_end
            raw_seq = genome.fetch(raw.chrom, g_s, g_e)
            seq = raw_seq if raw.strand == "+" else revcomp(raw_seq)
            flanks.append(IntronFlank(idx, "whole", g_s, g_e, seq, intron_len))
        else:
            if raw.strand == "+":
                donor_g = (intron_g_start, intron_g_start + flank_len - 1)
                acceptor_g = (intron_g_end - flank_len + 1, intron_g_end)
            else:
                # Transcript 5' end of the intron is at the HIGH genomic coordinate.
                donor_g = (intron_g_end - flank_len + 1, intron_g_end)
                acceptor_g = (intron_g_start, intron_g_start + flank_len - 1)
            for side, (g_s, g_e) in (("donor", donor_g), ("acceptor", acceptor_g)):
                raw_seq = genome.fetch(raw.chrom, g_s, g_e)
                seq = raw_seq if raw.strand == "+" else revcomp(raw_seq)
                flanks.append(IntronFlank(idx, side, g_s, g_e, seq, intron_len))

    return Transcript(
        transcript_id=raw.transcript_id,
        gene_id=raw.gene_id,
        gene_name=raw.gene_name,
        chrom=raw.chrom,
        strand=raw.strand,
        biotype=raw.biotype,
        exons=tuple(exons),
        five_utr=tuple(_utr_spans(raw, five=True)),
        three_utr=tuple(_utr_spans(raw, five=False)),
        mrna_seq=mrna_seq,
        flanks=tuple(flanks),
    )


def _utr_spans(raw: _RawTranscript, five: bool) -> list[tuple[int, int]]:
    """Split GENCODE's undifferentiated ``UTR`` feature into 5'/3' by CDS position.

    GENCODE GTF marks both UTRs with the single feature type ``UTR``; genomic
    position relative to the CDS span disambiguates which side each is on.
    A UTR interval strictly upstream of the CDS (in transcript direction) is
    5'; strictly downstream is 3'.
    """
    if not raw.cds_spans:
        return []
    cds_g_start = min(s for s, _ in raw.cds_spans)
    cds_g_end = max(e for _, e in raw.cds_spans)
    out = []
    for s, e in raw.utr_spans:
        if raw.strand == "+":
            is_five = e < cds_g_start
        else:
            is_five = s > cds_g_end
        if is_five == five:
            out.append((s, e))
    return out


# ---------------------------------------------------------------------------
# build / load the on-disk index (idempotent, atomic write)
# ---------------------------------------------------------------------------


def _to_jsonable(t: Transcript) -> dict:
    return {
        "transcript_id": t.transcript_id,
        "gene_id": t.gene_id,
        "gene_name": t.gene_name,
        "chrom": t.chrom,
        "strand": t.strand,
        "biotype": t.biotype,
        "exons": [[e.order, e.g_start, e.g_end, e.t_start, e.t_end] for e in t.exons],
        "five_utr": [list(s) for s in t.five_utr],
        "three_utr": [list(s) for s in t.three_utr],
        "mrna_seq": t.mrna_seq,
        "flanks": [
            [f.intron_index, f.side, f.g_start, f.g_end, f.seq, f.intron_len]
            for f in t.flanks
        ],
    }


def _from_jsonable(d: dict) -> Transcript:
    return Transcript(
        transcript_id=d["transcript_id"],
        gene_id=d["gene_id"],
        gene_name=d["gene_name"],
        chrom=d["chrom"],
        strand=d["strand"],
        biotype=d["biotype"],
        exons=tuple(Exon(*e) for e in d["exons"]),
        five_utr=tuple(tuple(s) for s in d["five_utr"]),
        three_utr=tuple(tuple(s) for s in d["three_utr"]),
        mrna_seq=d["mrna_seq"],
        flanks=tuple(IntronFlank(*f) for f in d["flanks"]),
    )


def build_index(
    out_path: Path = INDEX_PATH,
    gtf_path: Path = GTF_PATH,
    fasta_path: Path = FASTA_PATH,
    flank_len: int = DEFAULT_FLANK_LEN,
) -> Path:
    """Build the transcript structural index as newline-delimited JSON.

    Idempotent: skips if ``out_path`` already exists. Atomic: writes to a
    ``.tmp`` file and ``os.replace``s it into place, so a killed build never
    leaves a half-written file that a later run would mistake for done --
    same discipline as :mod:`smt.validate.length_sweep`.
    """
    if out_path.exists():
        print(f"index exists, skipping -> {out_path}")
        return out_path
    genome = GenomeIndex(fasta_path)
    tmp = out_path.with_suffix(".jsonl.tmp")
    n = 0
    try:
        with tmp.open("w") as fh:
            for raw in iter_canonical_protein_coding_transcripts(gtf_path):
                if not raw.exon_spans:
                    continue
                t = _build_transcript(raw, genome, flank_len)
                fh.write(json.dumps(_to_jsonable(t)) + "\n")
                n += 1
                if n % 2000 == 0:
                    print(f"  indexed {n} transcripts...", flush=True)
    finally:
        genome.close()
    os.replace(tmp, out_path)
    print(f"wrote {n} transcripts -> {out_path}")
    return out_path


def load_index(path: Path = INDEX_PATH) -> list[Transcript]:
    out = []
    with path.open() as fh:
        for line in fh:
            out.append(_from_jsonable(json.loads(line)))
    return out


def main() -> None:
    import argparse

    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--flank-len", type=int, default=DEFAULT_FLANK_LEN)
    args = p.parse_args()
    fetch_gtf()
    fetch_genome()
    build_index(flank_len=args.flank_len)


if __name__ == "__main__":
    main()
