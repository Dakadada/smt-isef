"""Frozen sequence scores for the Arm D locked nusinersen benchmark.

Two scores rank every exon (positives and matched backgrounds alike):

* BASELINE -- ``min 18-nt Hamming distance`` of the exon's transcript-strand
  ``+/- window`` sequence to nusinersen's reverse complement. A pure
  sequence-distance model: the lowest Hamming distance over every 18-nt window
  is the best complementarity-proxy match. Lower = a closer canonical
  nusinersen off-target site.
* TEST -- best (most negative) IntaRNA interaction energy of nusinersen against
  the same ``+/- window`` target, accessibility-aware (IntaRNA's RNAplfold
  accessibility model). More negative = a more favourable predicted duplex.

Both scores are deterministic functions of genomic sequence and the fixed
nusinersen oligo. They are therefore *inputs* that get frozen at preregistration
-- computing them does not read, and cannot reveal, the sealed treatment ΔPSI or
the benchmark outcome. This module never touches a treatment alignment or a
contrast table, and never computes a rank, percentile, ΔR, or permutation
p-value (that is :mod:`feng.bench.permute`, guarded shut until Phase 9).

Sequence extraction reuses the same pysam-preferred / samtools-fallback,
fail-closed contract as :mod:`feng.bench.strata` GC attachment: if any required
base cannot be read, no scores are returned.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable
import argparse
import shutil
import subprocess
import tempfile

import pandas as pd

from feng.bench.strata import GenomeFastaError, _resolve_fasta_contig

# --- nusinersen oligo --------------------------------------------------------

# Nusinersen (Spinraza) is an 18-nt 2'-MOE phosphorothioate ASO targeting ISS-N1
# in SMN2 intron 7. Written 5'->3' as its DNA-equivalent sequence. The reverse
# complement is the sense/target site the oligo is complementary to; the
# baseline score measures how closely an exon window resembles that target site.
NUSINERSEN_5P_TO_3P = "TCACTTTCATAATGCTGG"
OLIGO_LENGTH = len(NUSINERSEN_5P_TO_3P)

_COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

DEFAULT_WINDOW = 100  # +/- flank added to each exon before scoring


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


NUSINERSEN_RC = reverse_complement(NUSINERSEN_5P_TO_3P)


# --- sequence extraction -----------------------------------------------------


@dataclass(frozen=True, slots=True)
class WindowSpec:
    """A window to score: transcript-strand sequence around one exon."""

    exon_key: str
    chrom: str
    start_0base: int  # window start, 0-based half-open (clamped at 0)
    end: int          # window end, 0-based half-open
    strand: str


def exon_window_specs(
    events: pd.DataFrame, window: int = DEFAULT_WINDOW
) -> list[WindowSpec]:
    """Build +/- ``window`` window specs from exon-identity columns."""
    required = {"exon_key", "chrom", "exon_start_0base", "exon_end", "strand"}
    missing = required.difference(events.columns)
    if missing:
        raise ValueError(f"events lack columns for windowing: {sorted(missing)}")
    specs = []
    for row in events.itertuples(index=False):
        start = max(0, int(row.exon_start_0base) - window)
        end = int(row.exon_end) + window
        specs.append(
            WindowSpec(
                exon_key=str(row.exon_key),
                chrom=str(row.chrom),
                start_0base=start,
                end=end,
                strand=str(row.strand),
            )
        )
    return specs


def _fetch_transcript_sequences(
    specs: list[WindowSpec], genome_fasta: Path, *, samtools: Path | None = None
) -> dict[str, str]:
    """Fetch each window on the transcript strand, fail-closed.

    A '-' strand exon's transcript-strand sequence is the reverse complement of
    the genomic window, so the returned sequence is always 5'->3' along the
    transcript -- the strand the ASO binds.
    """
    genome_fasta = Path(genome_fasta)
    if not genome_fasta.is_file() or genome_fasta.stat().st_size == 0:
        raise GenomeFastaError(f"reference FASTA is missing or empty: {genome_fasta}")
    try:
        return _fetch_with_pysam(specs, genome_fasta)
    except GenomeFastaError:
        executable = Path(samtools) if samtools is not None else None
        if executable is None or not executable.is_file():
            raise
        return _fetch_with_samtools(specs, genome_fasta, executable)


def _orient(sequence: str, strand: str) -> str:
    sequence = sequence.upper()
    if not sequence:
        raise GenomeFastaError("reference FASTA returned an empty window sequence")
    return reverse_complement(sequence) if strand == "-" else sequence


def _fetch_with_pysam(specs: list[WindowSpec], genome_fasta: Path) -> dict[str, str]:
    try:
        import pysam
    except ImportError as exc:  # pragma: no cover - depends on build environment
        raise GenomeFastaError("pysam is unavailable") from exc
    out: dict[str, str] = {}
    try:
        with pysam.FastaFile(str(genome_fasta)) as fasta:
            available = set(fasta.references)
            for spec in specs:
                contig = _resolve_fasta_contig(spec.chrom, available)
                seq = fasta.fetch(contig, spec.start_0base, spec.end)
                out[spec.exon_key] = _orient(seq, spec.strand)
    except GenomeFastaError:
        raise
    except Exception as exc:  # pragma: no cover - IO failure path
        raise GenomeFastaError(f"failed to open/read FASTA {genome_fasta}") from exc
    return out


def _fetch_with_samtools(
    specs: list[WindowSpec], genome_fasta: Path, samtools: Path
) -> dict[str, str]:
    with tempfile.TemporaryDirectory(prefix="feng-scores-") as tmp:
        tmpdir = Path(tmp)
        fai = tmpdir / "reference.fai"
        try:
            subprocess.run(
                [str(samtools), "faidx", "--fai-idx", str(fai), str(genome_fasta)],
                check=True, capture_output=True, text=True,
            )
        except (OSError, subprocess.CalledProcessError) as exc:
            stderr = getattr(exc, "stderr", "")
            raise GenomeFastaError(f"failed to index FASTA {genome_fasta}: {stderr}") from exc
        available = {
            line.split("\t", 1)[0]
            for line in fai.read_text(encoding="utf-8").splitlines()
            if line
        }
        regions = tmpdir / "regions.txt"
        regions.write_text(
            "".join(
                f"{_resolve_fasta_contig(s.chrom, available)}:{s.start_0base + 1}-{s.end}\n"
                for s in specs
            ),
            encoding="utf-8",
        )
        try:
            result = subprocess.run(
                [str(samtools), "faidx", "--fai-idx", str(fai), "-r", str(regions), str(genome_fasta)],
                check=True, capture_output=True, text=True,
            )
        except (OSError, subprocess.CalledProcessError) as exc:
            stderr = getattr(exc, "stderr", "")
            raise GenomeFastaError(f"failed to read FASTA {genome_fasta}: {stderr}") from exc
    sequences: list[str] = []
    current: list[str] = []
    for line in result.stdout.splitlines():
        if line.startswith(">"):
            if current:
                sequences.append("".join(current))
                current = []
        else:
            current.append(line.strip())
    if current:
        sequences.append("".join(current))
    if len(sequences) != len(specs):
        raise GenomeFastaError(
            f"FASTA returned {len(sequences)} windows; expected {len(specs)}"
        )
    return {spec.exon_key: _orient(seq, spec.strand) for spec, seq in zip(specs, sequences)}


# --- baseline: min 18-nt Hamming distance ------------------------------------


def min_hamming_to_target(sequence: str, target: str = NUSINERSEN_RC) -> int:
    """Minimum Hamming distance between ``target`` and any equal-length window.

    Only fully spanning windows composed of ACGT bases are compared; a window
    containing an ``N`` (or any non-ACGT base) is skipped so an ambiguous base
    cannot masquerade as a match. Returns ``len(target)`` when no comparable
    window exists (the worst possible distance, i.e. no evidence of a site).
    """
    seq = sequence.upper()
    n = len(target)
    if len(seq) < n:
        return n
    best = n
    valid = set("ACGT")
    for i in range(len(seq) - n + 1):
        window = seq[i : i + n]
        if not set(window) <= valid:
            continue
        dist = sum(1 for a, b in zip(window, target) if a != b)
        if dist < best:
            best = dist
            if best == 0:
                break
    return best


def baseline_hamming_scores(
    events: pd.DataFrame,
    genome_fasta: Path,
    *,
    window: int = DEFAULT_WINDOW,
    samtools: Path | None = None,
    target: str = NUSINERSEN_RC,
) -> pd.DataFrame:
    """Attach ``baseline_min_hamming`` for every exon (fail-closed on sequence)."""
    specs = exon_window_specs(events, window=window)
    seqs = _fetch_transcript_sequences(specs, genome_fasta, samtools=samtools)
    scores = []
    for spec in specs:
        seq = seqs[spec.exon_key]
        scores.append(min_hamming_to_target(seq, target))
    out = events.copy()
    out["baseline_min_hamming"] = scores
    return out


# --- test: IntaRNA interaction energy ----------------------------------------


class IntaRNAError(RuntimeError):
    """Raised when the IntaRNA backend is missing or fails."""


def find_intarna(explicit: Path | None = None) -> Path | None:
    """Locate an IntaRNA executable, preferring an explicit path then PATH."""
    if explicit is not None:
        p = Path(explicit)
        return p if p.is_file() else None
    found = shutil.which("IntaRNA")
    return Path(found) if found else None


def intarna_min_energy(
    target_sequence: str,
    *,
    intarna: Path,
    query: str = NUSINERSEN_5P_TO_3P,
    extra_args: Iterable[str] = (),
) -> float:
    """Return the best (most negative) IntaRNA interaction energy in kcal/mol.

    IntaRNA slides the query over the target internally and reports the minimum
    free energy of the best predicted duplex; accessibility is IntaRNA's default
    RNAplfold-based model. Returns ``0.0`` when IntaRNA reports no interaction
    (its convention for "no favourable duplex"). Raises :class:`IntaRNAError`
    on any backend failure -- energies are never fabricated.
    """
    intarna = Path(intarna)
    if not intarna.is_file():
        raise IntaRNAError(f"IntaRNA executable not found: {intarna}")
    query_rna = query.upper().replace("T", "U")
    target_rna = target_sequence.upper().replace("T", "U")
    cmd = [
        str(intarna),
        "-q", query_rna,
        "-t", target_rna,
        "--outMode", "C",          # CSV
        "--outCsvCols", "E",       # energy only
        "-n", "1",                 # best interaction only
        *extra_args,
    ]
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        stderr = getattr(exc, "stderr", "")
        raise IntaRNAError(f"IntaRNA failed: {stderr.strip() or exc}") from exc
    lines = [ln for ln in result.stdout.splitlines() if ln.strip()]
    # First line is the header ("E"); a missing data line means no interaction.
    if len(lines) < 2:
        return 0.0
    try:
        return float(lines[1].split(";")[0])
    except (ValueError, IndexError) as exc:
        raise IntaRNAError(f"could not parse IntaRNA energy from: {result.stdout!r}") from exc


def intarna_energy_scores(
    events: pd.DataFrame,
    genome_fasta: Path,
    *,
    intarna: Path,
    window: int = DEFAULT_WINDOW,
    samtools: Path | None = None,
    query: str = NUSINERSEN_5P_TO_3P,
) -> pd.DataFrame:
    """Attach ``intarna_min_energy`` for every exon.

    WARNING -- one IntaRNA call per exon over a full universe is a 100-400
    CPU-hour sweep. Run it under ``feng/scripts/bench_intarna_scores.sbatch``,
    not inline. This helper exists so the sbatch and the tiny unit fixture share
    one code path.
    """
    exe = find_intarna(intarna)
    if exe is None:
        raise IntaRNAError(f"IntaRNA executable not available at {intarna}")
    specs = exon_window_specs(events, window=window)
    seqs = _fetch_transcript_sequences(specs, genome_fasta, samtools=samtools)
    energies = [intarna_min_energy(seqs[s.exon_key], intarna=exe, query=query) for s in specs]
    out = events.copy()
    out["intarna_min_energy"] = energies
    return out


# --- CLI (baseline only; IntaRNA sweep goes through the sbatch) ---------------


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Compute frozen Arm D sequence scores.")
    parser.add_argument("--events", type=Path, required=True, help="parquet with exon-identity columns")
    parser.add_argument("--genome-fasta", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--window", type=int, default=DEFAULT_WINDOW)
    parser.add_argument("--samtools", type=Path, default=None)
    parser.add_argument("--score", choices=["baseline", "intarna"], default="baseline")
    parser.add_argument("--intarna", type=Path, default=None)
    parser.add_argument("--n-shards", type=int, default=1, help="split the input into N strided shards")
    parser.add_argument("--shard-index", type=int, default=0, help="0-based shard to process")
    args = parser.parse_args(argv)

    if args.n_shards < 1 or not (0 <= args.shard_index < args.n_shards):
        parser.error("require --n-shards >= 1 and 0 <= --shard-index < --n-shards")

    events = pd.read_parquet(args.events)
    if args.n_shards > 1:
        # Strided slice keeps IntaRNA cost balanced across shards; row order is
        # irrelevant because each exon is scored independently.
        events = events.iloc[args.shard_index :: args.n_shards].copy()
    if args.score == "baseline":
        scored = baseline_hamming_scores(
            events, args.genome_fasta, window=args.window, samtools=args.samtools
        )
    else:
        if args.intarna is None:
            parser.error("--intarna is required for --score intarna")
        scored = intarna_energy_scores(
            events, args.genome_fasta, intarna=args.intarna, window=args.window, samtools=args.samtools
        )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    scored.to_parquet(args.out, index=False)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
