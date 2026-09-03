"""Symmetric nusinersen length-titration ladder builder (Arm E / StericRuler).

For each even length ``L`` in ``[min_len, max_len]`` this builds the ``L``-nt
antisense oligo obtained by trimming (``L < 18``) or extending (``L > 18``) the
true nusinersen 18-mer symmetrically -- one nucleotide per end -- around its
ISS-N1 target window. The 18-mer sits at the centre of the ladder, so the whole
series is a genuine truncation/extension series of the real drug rather than a
shifted reporter-construct offset.

Every sequence is fetched from the reference genome at the confirmed anchor and
derived by reverse complement; no ASO sequence is hardcoded. Fail-closed: the
``L = 18`` rung must reproduce nusinersen exactly (both the ASO and its target),
otherwise :func:`build_ladder` raises -- a wrong frame cannot silently propagate
to the other rungs.

Freeze-safe: this module reads only reference genome sequence. It never touches a
treatment alignment, a contrast table, or the sealed Phase-9 ΔPSI, and computes
no rank, percentile, or benchmark outcome. Every emitted row carries honesty
state ``candidate``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import shutil
import sys

from feng.bench.scores import (
    NUSINERSEN_5P_TO_3P,
    NUSINERSEN_RC,
    OLIGO_LENGTH,
    WindowSpec,
    _fetch_transcript_sequences,
    reverse_complement,
)
from feng.bench.strata import GenomeFastaError
from feng.honesty.states import ResultState

# --- confirmed ISS-N1 nusinersen anchor (GRCh38 v50, + strand) ---------------
# 1-based inclusive; verified by extraction in
# feng/artifacts/phase2_variant_recon/README.md (lines 167-177). This 18-nt
# window is the L = 18 target; its reverse complement is nusinersen itself.
ANCHOR_CHROM = "chr5"
ANCHOR_START_1BASED = 70076584
ANCHOR_END_1BASED = 70076601
ANCHOR_STRAND = "+"

_CANDIDATE = ResultState.CANDIDATE.value

# The anchor window must be exactly the nusinersen oligo length, or the whole
# symmetric-trim scheme is mis-framed.
assert ANCHOR_END_1BASED - ANCHOR_START_1BASED + 1 == OLIGO_LENGTH


@dataclass(frozen=True, slots=True)
class LadderRung:
    """One length rung: coordinates + derived target and ASO sequences."""

    length: int
    chrom: str
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive
    strand: str
    target_5to3: str
    aso_5to3: str
    state: str

    def as_row(self) -> dict:
        return asdict(self)


def _resolve_samtools(samtools: Path | None) -> Path | None:
    """Prefer an explicit path, then ``$PATH``, then the env sibling binary."""
    if samtools is not None:
        return Path(samtools)
    found = shutil.which("samtools")
    if found:
        return Path(found)
    sibling = Path(sys.executable).resolve().parent / "samtools"
    return sibling if sibling.is_file() else None


def _window_for_length(length: int) -> tuple[int, int]:
    """1-based inclusive (start, end) for an even ``length`` rung.

    ``half`` is positive when trimming (``length < 18``) and negative when
    extending (``length > 18``); both bounds move symmetrically.
    """
    half = (OLIGO_LENGTH - length) // 2
    return ANCHOR_START_1BASED + half, ANCHOR_END_1BASED - half


def _ladder_lengths(min_len: int, max_len: int) -> list[int]:
    if min_len < 2 or max_len < min_len:
        raise ValueError(
            f"invalid ladder bounds: min_len={min_len}, max_len={max_len}"
        )
    if min_len % 2 or max_len % 2:
        raise ValueError(
            f"ladder lengths must be even; got min_len={min_len}, max_len={max_len}"
        )
    if not (min_len <= OLIGO_LENGTH <= max_len):
        raise ValueError(
            "ladder range must include L=18 so the nusinersen anchor invariant "
            f"is checkable; got [{min_len}, {max_len}]"
        )
    return list(range(min_len, max_len + 1, 2))


def build_ladder(
    *,
    genome_fasta: str | Path,
    min_len: int = 2,
    max_len: int = 30,
    samtools: Path | None = None,
) -> list[dict]:
    """Build the symmetric length ladder as a list of tabular row dicts.

    One row per even length in ``[min_len, max_len]`` with keys ``length,
    chrom, start, end, strand, target_5to3, aso_5to3, state``. Coordinates are
    1-based inclusive; ``target_5to3`` is the + strand genomic window and
    ``aso_5to3`` its reverse complement. Raises if the ``L = 18`` rung does not
    reproduce nusinersen exactly, or if any window cannot be read.
    """
    lengths = _ladder_lengths(min_len, max_len)
    resolved_samtools = _resolve_samtools(samtools)

    specs: list[WindowSpec] = []
    for length in lengths:
        start_1based, end_1based = _window_for_length(length)
        if end_1based - start_1based + 1 != length:
            raise AssertionError(
                f"window arithmetic error at L={length}: "
                f"[{start_1based}, {end_1based}]"
            )
        specs.append(
            WindowSpec(
                exon_key=f"E_L{length:02d}",
                chrom=ANCHOR_CHROM,
                start_0base=start_1based - 1,  # 0-based half-open for the fetch
                end=end_1based,
                strand=ANCHOR_STRAND,
            )
        )

    sequences = _fetch_transcript_sequences(
        specs, Path(genome_fasta), samtools=resolved_samtools
    )

    rungs: list[LadderRung] = []
    for length, spec in zip(lengths, specs):
        target = sequences[spec.exon_key]
        if len(target) != length:
            raise GenomeFastaError(
                f"L={length} window returned {len(target)} nt; expected {length}"
            )
        rungs.append(
            LadderRung(
                length=length,
                chrom=ANCHOR_CHROM,
                start=spec.start_0base + 1,
                end=spec.end,
                strand=ANCHOR_STRAND,
                target_5to3=target,
                aso_5to3=reverse_complement(target),
                state=_CANDIDATE,
            )
        )

    _assert_nusinersen_anchor(rungs)
    return [rung.as_row() for rung in rungs]


def _assert_nusinersen_anchor(rungs: list[LadderRung]) -> None:
    """Fail-closed: the L=18 rung must be nusinersen and its RC target."""
    anchor = next((r for r in rungs if r.length == OLIGO_LENGTH), None)
    if anchor is None:  # pragma: no cover - guarded by _ladder_lengths
        raise ValueError("ladder is missing the L=18 anchor rung")
    if anchor.aso_5to3 != NUSINERSEN_5P_TO_3P or anchor.target_5to3 != NUSINERSEN_RC:
        raise ValueError(
            "L=18 rung does not reproduce nusinersen -- frame is wrong: "
            f"aso={anchor.aso_5to3!r} target={anchor.target_5to3!r}"
        )
