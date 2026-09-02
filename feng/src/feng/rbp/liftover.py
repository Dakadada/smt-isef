"""Arm B hg19->GRCh38 (GENCODE v50 frame) liftover bridge for the RBP layer.

Claim boundary: this module only *re-coordinates* processed ENCODE data from the
``hg19`` frame the files are distributed on to the ``GRCh38`` frame the rest of
the arm works in. It computes no concordance statistic and asserts no
interaction; those belong to a separate, post-freeze step. What it does provide
is the single missing bridge that lets an eCLIP peak frame (read by
:mod:`feng.rbp.encode`, which stamps every row ``assembly='hg19'`` and refuses
cross-assembly joins) legitimately pass :func:`feng.rbp.encode.assert_same_assembly`
against the GRCh38 exon universe: the coordinates are actually lifted and the
frame is re-stamped ``GRCh38`` with chain-file provenance. The hard fence in
``encode.py`` is never weakened; a still-hg19 frame stays rejected.

Liftover engine: a self-contained parser + interval mapper for the UCSC ``chain``
format (the same format ``liftOver``/``pyliftover``/``CrossMap`` consume). No
external tool is required. Supply the standard UCSC chain
``hg19ToHg38.over.chain(.gz)`` (see :data:`UCSC_HG19_TO_HG38_CHAIN_URL`); its
path and SHA-256 are recorded on every lifted frame for provenance.

Coordinate conventions
----------------------
* Input peak coordinates are 0-based half-open BED intervals (as read by
  :func:`feng.rbp.encode.read_eclip_peaks`).
* In a UCSC ``chain``, the *target* (``t*``) side is the OLD genome (hg19) and
  the *query* (``q*``) side is the NEW genome (hg38). ``tStrand`` is always
  ``+``; ``qStrand`` may be ``+`` or ``-`` (a lifted interval on a ``-`` chain
  has its feature strand flipped and its endpoints reordered).
* An interval lifts only if BOTH endpoints map, land on the same query contig,
  and via a single unambiguous chain result. Anything else is dropped and
  counted (``unmapped`` or ``multimapped``) -- never silently kept.
"""

from __future__ import annotations

from dataclasses import dataclass
import gzip
import hashlib
from pathlib import Path

import pandas as pd

from . import encode

# The production chain is the standard UCSC over.chain. It is a small (~230 KB
# gzipped) public download; it is NOT vendored into the repo. Stage it on scratch
# (or any readable path) and pass the path to :func:`load_chain_file`.
UCSC_HG19_TO_HG38_CHAIN_URL = (
    "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz"
)

# Provenance/liftover columns added to a lifted frame (on top of encode.py's
# source_* columns). ``assembly`` is re-stamped GRCh38 so encode.py's fence
# accepts a downstream join; ``needs_liftover`` is cleared.
LIFTED_PROVENANCE_COLUMNS = (
    "target_assembly",
    "liftover_chain_path",
    "liftover_chain_sha256",
)


class LiftoverError(RuntimeError):
    """Raised on malformed chain input or an invalid liftover request."""


@dataclass(frozen=True, slots=True)
class _Block:
    """One ungapped alignment block: target [t, t+size), query [q, q+size)."""

    t: int
    q: int
    size: int


@dataclass(frozen=True, slots=True)
class _Chain:
    """One UCSC chain record (target = hg19, query = hg38)."""

    t_name: str
    t_size: int
    t_start: int
    t_end: int
    q_name: str
    q_size: int
    q_strand: str
    q_start: int
    q_end: int
    blocks: tuple[_Block, ...]

    def lift_point(self, pos: int) -> tuple[str, int, str] | None:
        """Lift a single 0-based target position to (q_name, q_pos_plus, strand).

        ``q_pos_plus`` is always reported on the forward (``+``) genome strand,
        regardless of the chain's query strand; ``strand`` is the resulting
        feature strand (flipped for a ``-`` chain).
        """
        if not (self.t_start <= pos < self.t_end):
            return None
        for block in self.blocks:
            if block.t <= pos < block.t + block.size:
                q_in_chain = block.q + (pos - block.t)
                if self.q_strand == "+":
                    q_forward = q_in_chain
                else:  # minus-strand chain: convert to the forward coordinate
                    q_forward = self.q_size - 1 - q_in_chain
                return (self.q_name, q_forward, self.q_strand)
        # Position falls in a chain gap (deleted in the new assembly).
        return None


class ChainMap:
    """Parsed UCSC chain file indexed by target (hg19) contig."""

    def __init__(self, chains: list[_Chain], *, path: Path, sha256: str) -> None:
        self._by_contig: dict[str, list[_Chain]] = {}
        for chain in chains:
            self._by_contig.setdefault(chain.t_name, []).append(chain)
        for contig_chains in self._by_contig.values():
            contig_chains.sort(key=lambda c: c.t_start)
        self.path = path
        self.sha256 = sha256

    @property
    def n_chains(self) -> int:
        return sum(len(v) for v in self._by_contig.values())

    def _candidate_chains(self, contig: str, pos: int) -> list[_Chain]:
        # Chains may overlap on the target side (the source of multi-mapping),
        # so every chain whose target span covers ``pos`` is a candidate.
        return [
            c
            for c in self._by_contig.get(contig, ())
            if c.t_start <= pos < c.t_end
        ]

    def lift_point(self, contig: str, pos: int) -> list[tuple[str, int, str]]:
        """All distinct forward-strand lifts of one 0-based target position."""
        results: list[tuple[str, int, str]] = []
        for chain in self._candidate_chains(contig, pos):
            lifted = chain.lift_point(pos)
            if lifted is not None and lifted not in results:
                results.append(lifted)
        return results

    def lift_interval(
        self, contig: str, start: int, end: int
    ) -> list[tuple[str, int, int, str]]:
        """Lift a 0-based half-open interval; return distinct (chrom,start,end,strand).

        Both the first (``start``) and last (``end-1``) covered bases must lift
        through the SAME chain to the same query contig and strand; the result is
        rebuilt as a half-open interval (endpoints reordered for ``-`` chains).
        Returns every distinct interval the chains support -- more than one means
        the source interval is multi-mapping.
        """
        if end <= start:
            raise LiftoverError(f"empty/invalid interval [{start}, {end})")
        last = end - 1
        results: list[tuple[str, int, int, str]] = []
        for chain in self._candidate_chains(contig, start):
            first_lift = chain.lift_point(start)
            last_lift = chain.lift_point(last)
            if first_lift is None or last_lift is None:
                continue  # interval straddles a gap within this chain
            f_chrom, f_pos, strand = first_lift
            l_chrom, l_pos, _ = last_lift
            if f_chrom != l_chrom:
                continue  # split across query contigs -> not a clean lift
            lo, hi = (f_pos, l_pos) if f_pos <= l_pos else (l_pos, f_pos)
            interval = (f_chrom, lo, hi + 1, strand)
            if interval not in results:
                results.append(interval)
        return results


def _parse_chain_stream(handle) -> list[_Chain]:
    chains: list[_Chain] = []
    header: list[str] | None = None
    blocks: list[_Block] = []
    t_pos = q_pos = 0
    for raw in handle:
        line = raw.strip()
        if not line:
            if header is not None:
                chains.append(_finalize_chain(header, blocks))
                header = None
                blocks = []
            continue
        if line.startswith("chain"):
            if header is not None:  # a chain with no trailing blank line
                chains.append(_finalize_chain(header, blocks))
                blocks = []
            header = line.split()
            fields = header
            # chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
            t_pos = int(fields[5])
            q_pos = int(fields[10])
            continue
        if header is None:
            raise LiftoverError(f"chain data line before any chain header: {line!r}")
        parts = line.split()
        size = int(parts[0])
        blocks.append(_Block(t=t_pos, q=q_pos, size=size))
        if len(parts) >= 3:  # more blocks follow
            dt, dq = int(parts[1]), int(parts[2])
            t_pos += size + dt
            q_pos += size + dq
        # a lone `size` line closes the block list of this chain
    if header is not None:
        chains.append(_finalize_chain(header, blocks))
    return chains


def _finalize_chain(header: list[str], blocks: list[_Block]) -> _Chain:
    if len(header) < 12:
        raise LiftoverError(f"malformed chain header: {header!r}")
    q_strand = header[9]
    if q_strand not in {"+", "-"}:
        raise LiftoverError(f"unsupported query strand {q_strand!r}")
    if header[4] != "+":
        raise LiftoverError(f"unsupported target strand {header[4]!r} (must be +)")
    if not blocks:
        raise LiftoverError(f"chain {header!r} has no alignment blocks")
    return _Chain(
        t_name=header[2],
        t_size=int(header[3]),
        t_start=int(header[5]),
        t_end=int(header[6]),
        q_name=header[7],
        q_size=int(header[8]),
        q_strand=q_strand,
        q_start=int(header[10]),
        q_end=int(header[11]),
        blocks=tuple(blocks),
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_chain_file(path: str | Path) -> ChainMap:
    """Parse a UCSC chain file (plain or gzipped) into a :class:`ChainMap`.

    The file's SHA-256 is computed and carried on the map so every lifted frame
    can record exactly which chain produced its coordinates.
    """
    path = Path(path)
    if not path.exists():
        raise LiftoverError(
            f"chain file not found: {path}. Stage the UCSC chain from "
            f"{UCSC_HG19_TO_HG38_CHAIN_URL}"
        )
    sha = _sha256(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        chains = _parse_chain_stream(handle)
    if not chains:
        raise LiftoverError(f"no chains parsed from {path}")
    return ChainMap(chains, path=path, sha256=sha)


@dataclass(frozen=True, slots=True)
class LiftResult:
    """Outcome of lifting an eCLIP peak frame to GRCh38."""

    frame: encode.AssemblyDataFrame
    n_input: int
    n_lifted: int
    n_unmapped: int
    n_multimapped: int

    @property
    def mapped_fraction(self) -> float:
        return self.n_lifted / self.n_input if self.n_input else 0.0


def lift_eclip_peaks(
    frame: pd.DataFrame,
    chain: ChainMap,
    *,
    chrom_col: str = "chrom",
    start_col: str = "start",
    end_col: str = "end",
    strand_col: str = "strand",
) -> LiftResult:
    """Lift an hg19 eCLIP peak frame to GRCh38 and re-stamp its provenance.

    The input must be an hg19 frame as produced by
    :func:`feng.rbp.encode.read_eclip_peaks` (``assembly == 'hg19'`` and
    ``needs_liftover`` True); lifting a non-hg19 or already-lifted frame is
    refused so the fence cannot be bypassed. Unmapped and multi-mapped peaks are
    dropped and counted. The returned frame is an
    :class:`~feng.rbp.encode.AssemblyDataFrame` stamped ``assembly='GRCh38'``
    (column and ``attrs``), ``target_assembly='GRCh38'``, ``needs_liftover``
    False, with the chain path + SHA-256 recorded -- exactly what
    :func:`feng.rbp.encode.assert_same_assembly` needs to accept a join against a
    GRCh38 coordinate frame.
    """
    for col in (chrom_col, start_col, end_col):
        if col not in frame.columns:
            raise LiftoverError(f"input frame missing required column {col!r}")
    _require_unlifted_hg19(frame)

    lifted_rows: list[dict] = []
    n_unmapped = 0
    n_multimapped = 0
    has_strand = strand_col in frame.columns
    for record in frame.to_dict("records"):
        start = record[start_col]
        end = record[end_col]
        if pd.isna(start) or pd.isna(end):
            n_unmapped += 1
            continue
        results = chain.lift_interval(str(record[chrom_col]), int(start), int(end))
        if not results:
            n_unmapped += 1
            continue
        if len(results) > 1:
            n_multimapped += 1
            continue
        new_chrom, new_start, new_end, new_strand = results[0]
        row = dict(record)
        row[chrom_col] = new_chrom
        row[start_col] = new_start
        row[end_col] = new_end
        if has_strand:
            row[strand_col] = _combine_strand(record[strand_col], new_strand)
        lifted_rows.append(row)

    lifted = encode.AssemblyDataFrame(pd.DataFrame(lifted_rows, columns=list(frame.columns)))
    lifted["source_assembly"] = encode.SOURCE_ASSEMBLY
    lifted["assembly"] = encode.TARGET_ASSEMBLY
    lifted["target_assembly"] = encode.TARGET_ASSEMBLY
    lifted["needs_liftover"] = False
    lifted["liftover_chain_path"] = str(chain.path)
    lifted["liftover_chain_sha256"] = chain.sha256
    lifted.attrs["assembly"] = encode.TARGET_ASSEMBLY
    lifted.attrs["needs_liftover"] = False
    lifted.attrs["source_assembly"] = encode.SOURCE_ASSEMBLY
    lifted.attrs["liftover_chain_sha256"] = chain.sha256

    return LiftResult(
        frame=lifted,
        n_input=len(frame),
        n_lifted=len(lifted_rows),
        n_unmapped=n_unmapped,
        n_multimapped=n_multimapped,
    )


def _require_unlifted_hg19(frame: pd.DataFrame) -> None:
    """Refuse to lift anything that is not a fresh hg19 source frame."""
    declared = frame.attrs.get("assembly")
    if "assembly" in frame.columns:
        values = set(frame["assembly"].dropna().unique())
        if values and values != {encode.SOURCE_ASSEMBLY}:
            raise LiftoverError(
                f"refusing to lift a frame that is not uniformly hg19: {values}"
            )
    if declared is not None and declared != encode.SOURCE_ASSEMBLY:
        raise LiftoverError(
            f"refusing to lift a frame declared assembly={declared!r}; "
            f"only {encode.SOURCE_ASSEMBLY!r} sources are liftable"
        )


def _combine_strand(feature_strand: object, chain_strand: str) -> str:
    """Flip the feature strand when the chain maps to the query minus strand."""
    fs = str(feature_strand)
    if chain_strand == "+":
        return fs
    flip = {"+": "-", "-": "+"}
    return flip.get(fs, fs)
