"""Per-mer CPU scorers for the Arm E length-titration ladder (StericRuler).

This layer takes the rungs from :func:`feng.titration.ladder.build_ladder` and
attaches, for each length, a panel of sequence-only biophysical scores computed
with the same external tools Arm D already trusts (IntaRNA + the ViennaRNA
suite) plus the in-repo seed-and-extend off-target scanner. Every score is
parsed from the real tool's stdout -- nothing is fabricated; a missing binary
yields ``None`` plus an availability flag, never a fake number.

Scores per mer:

* ``intarna_dg`` / ``intarna_no_duplex`` -- on-target IntaRNA interaction energy
  (query = the mer's ASO, target = its target window). IntaRNA reports ``0.0``
  for "no favourable duplex"; short mers hit this, so ``intarna_no_duplex``
  flags it and the ``0.0`` must never be read as a weak-but-real energy.
* ``rnaduplex_dg`` -- ViennaRNA ``RNAduplex`` ASO:target hybridisation ΔG.
* ``rnacofold_dg`` -- ViennaRNA ``RNAcofold`` ASO&target co-fold ΔG.
* ``rnaup_dg`` -- ViennaRNA ``RNAup`` accessibility-corrected binding energy
  (the total in its ``(-x = -y + ...)`` decomposition): a better on-target
  proxy than IntaRNA's default because it charges the cost of opening target
  structure.
* ``rnaplfold_accessibility`` -- mean per-base unpaired probability of the
  target bases under ``RNAplfold``, computed with the target embedded in
  genomic flank so the score reflects the ISS-N1 site's accessibility in
  pre-mRNA context rather than an isolated L-nt fragment folding on itself
  (0..1; higher = more accessible). Requires the genome FASTA; ``None`` when
  it is not supplied.
* ``rnafold_mfe`` -- ViennaRNA ``RNAfold`` MFE of the ASO alone (self-structure;
  a strongly self-paired ASO is a worse binder).
* ``off_target_hits`` -- count of near-complement sites for the ASO in a
  reference, via :func:`feng.scan.seed_and_extend`.

Freeze-safe: reads only reference/ASO sequence and public tool output; touches
no treatment alignment, contrast, rank, or sealed ΔPSI. Enriched rows keep the
``candidate`` honesty state of the ladder rows they extend.
"""

from __future__ import annotations

from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
from typing import Mapping, Sequence

from feng.bench.scores import (
    IntaRNAError,
    WindowSpec,
    _fetch_transcript_sequences,
    find_intarna,
    intarna_min_energy,
    reverse_complement,
)
from feng.bench.strata import GenomeFastaError
from feng.honesty.states import ResultState
from feng.scan.seed_extend import HybridizationHit, seed_and_extend

_CANDIDATE = ResultState.CANDIDATE.value

# A fully bracketed free energy, e.g. a parenthesised negative float at the end
# of a Vienna structure line. Vienna prints it as the last such token.
_PAREN_ENERGY = re.compile(r"\(\s*(-?\d+(?:\.\d+)?)\s*\)")
# RNAup decomposes its energy as "(total = bind + open ...)"; the first number
# (before the "=") is the total accessibility-corrected binding energy.
_RNAUP_TOTAL = re.compile(r"\(\s*(-?\d+(?:\.\d+)?)\s*=")

# The Vienna binaries this module knows how to drive.
VIENNA_TOOLS = ("RNAduplex", "RNAcofold", "RNAup", "RNAplfold", "RNAfold")

# Accessibility defaults (candidate, pre-registration-eligible). The target is
# embedded in ``ACCESSIBILITY_PAD`` nt of genomic flank on each side -- >= the
# averaging window ``W`` -- so its bases see real pre-mRNA context instead of
# folding on themselves. ``W`` (RNAplfold's averaging window) and ``L`` (maximum
# base-pair span, bounded below ``W`` to the range where local folding is
# reliable) follow common RNAplfold accessibility practice.
ACCESSIBILITY_PAD = 70
ACCESSIBILITY_WINDOW = 70
ACCESSIBILITY_SPAN = 40


class ViennaError(RuntimeError):
    """Raised when a ViennaRNA backend is missing or emits unparseable output."""


def _to_rna(seq: str) -> str:
    """Uppercase and map DNA ``T`` to RNA ``U`` for the ViennaRNA/IntaRNA tools."""
    return seq.upper().replace("T", "U")


def _resolve_tool(tool: str, bindir: Path | None = None) -> Path | None:
    """Locate an executable: explicit ``bindir`` first, then ``$PATH``, then the
    sibling of the running interpreter (the feng-env bin, which is not
    necessarily on ``$PATH`` under a bare ``env/bin/python`` invocation)."""
    if bindir is not None:
        candidate = Path(bindir) / tool
        return candidate if candidate.is_file() else None
    found = shutil.which(tool)
    if found:
        return Path(found)
    sibling = Path(sys.executable).resolve().parent / tool
    return sibling if sibling.is_file() else None


def find_vienna(tool: str, bindir: Path | None = None) -> Path | None:
    """Locate a ViennaRNA executable (see :func:`_resolve_tool`)."""
    return _resolve_tool(tool, bindir)


def resolve_intarna(intarna: Path | None = None) -> Path | None:
    """Locate IntaRNA with the interpreter-sibling fallback ``find_intarna``
    lacks, so it resolves in feng-env without an activated ``$PATH``."""
    if intarna is not None:
        exe = find_intarna(intarna)
        if exe is not None:
            return exe
    return _resolve_tool("IntaRNA")


def _resolve_samtools(samtools: Path | None = None) -> Path | None:
    """Prefer an explicit path, then ``$PATH``, then the env sibling binary
    (the fetch fallback used when pysam is absent)."""
    if samtools is not None:
        return Path(samtools)
    found = shutil.which("samtools")
    if found:
        return Path(found)
    sibling = Path(sys.executable).resolve().parent / "samtools"
    return sibling if sibling.is_file() else None


def _run(cmd: Sequence[str], *, stdin: str | None = None, cwd: Path | None = None) -> str:
    """Run a subprocess with the never-fabricate discipline of ``scores.py``."""
    try:
        result = subprocess.run(
            [str(c) for c in cmd],
            input=stdin,
            check=True,
            capture_output=True,
            text=True,
            cwd=str(cwd) if cwd is not None else None,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        stderr = getattr(exc, "stderr", "")
        raise ViennaError(f"{cmd[0]} failed: {stderr.strip() or exc}") from exc
    return result.stdout


def _last_paren_energy(text: str, *, tool: str) -> float:
    matches = _PAREN_ENERGY.findall(text)
    if not matches:
        raise ViennaError(f"could not parse energy from {tool} output: {text!r}")
    return float(matches[-1])


# --- A. on-target IntaRNA ----------------------------------------------------


def intarna_on_target(aso_5to3: str, target_5to3: str, *, intarna: Path) -> tuple[float, bool]:
    """Return ``(intarna_dg, no_duplex)`` for the ASO against its target window.

    ``no_duplex`` is ``True`` exactly when IntaRNA reports ``0.0`` (its "no
    favourable duplex" convention) -- the caller must not read that ``0.0`` as a
    weak real energy.
    """
    dg = intarna_min_energy(target_5to3, intarna=intarna, query=aso_5to3)
    return dg, dg == 0.0


# --- B. duplex ΔG (RNAduplex / RNAcofold) ------------------------------------


def rnaduplex_dg(aso_5to3: str, target_5to3: str, *, rnaduplex: Path) -> float:
    """ViennaRNA ``RNAduplex`` hybridisation ΔG (kcal/mol) for ASO:target."""
    stdin = f"{_to_rna(aso_5to3)}\n{_to_rna(target_5to3)}\n"
    out = _run([rnaduplex], stdin=stdin)
    return _last_paren_energy(out, tool="RNAduplex")


def rnacofold_dg(aso_5to3: str, target_5to3: str, *, rnacofold: Path) -> float:
    """ViennaRNA ``RNAcofold`` co-fold ΔG (kcal/mol) for ASO&target."""
    stdin = f"{_to_rna(aso_5to3)}&{_to_rna(target_5to3)}\n"
    out = _run([rnacofold, "--noPS"], stdin=stdin)
    return _last_paren_energy(out, tool="RNAcofold")


# --- C. accessibility (RNAup / RNAplfold) ------------------------------------


def rnaup_dg(aso_5to3: str, target_5to3: str, *, rnaup: Path) -> float:
    """ViennaRNA ``RNAup`` accessibility-corrected binding energy (kcal/mol).

    Parses the total from RNAup's ``(-total = -bind + open ...)`` line. Runs in
    a temp dir because RNAup writes an ``RNA_w*_u*.out`` side file.
    """
    stdin = f"{_to_rna(aso_5to3)}\n{_to_rna(target_5to3)}\n"
    with tempfile.TemporaryDirectory() as tmp:
        out = _run([rnaup, "-b"], stdin=stdin, cwd=Path(tmp))
    match = _RNAUP_TOTAL.search(out)
    if match is None:
        # RNAup with no favourable interaction can omit the decomposition.
        raise ViennaError(f"could not parse RNAup total energy from: {out!r}")
    return float(match.group(1))


def mean_unpaired_over_window(
    padded_5to3: str,
    target_offset: int,
    target_len: int,
    *,
    rnaplfold: Path,
    window: int = ACCESSIBILITY_WINDOW,
    span: int = ACCESSIBILITY_SPAN,
) -> float:
    """Mean per-base unpaired probability over the target bases of a padded seq.

    ``padded_5to3`` is the target embedded in genomic flank; the target occupies
    the half-open slice ``[target_offset, target_offset + target_len)`` (0-based).
    RNAplfold folds the whole padded sequence -- so target bases can pair into the
    flank -- and we average the single-base (``l=1``) unpaired probabilities over
    only the target positions, never the pad. Fail-closed: a missing ``_lunp`` or
    a target-row count that does not match ``target_len`` raises.
    """
    if target_len <= 0 or target_offset < 0 or target_offset + target_len > len(padded_5to3):
        raise ViennaError("target slice falls outside the padded sequence")
    rna = _to_rna(padded_5to3)
    eff_window = min(window, len(rna))
    eff_span = min(span, eff_window)
    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        _run(
            [rnaplfold, "-W", str(eff_window), "-L", str(eff_span), "-u", "1"],
            stdin=f">win\n{rna}\n",
            cwd=tmp_path,
        )
        lunp = tmp_path / "win_lunp"
        if not lunp.is_file():
            raise ViennaError(f"RNAplfold produced no _lunp file for {rna!r}")
        lo = target_offset + 1  # _lunp positions are 1-based
        hi = target_offset + target_len
        probs: list[float] = []
        for line in lunp.read_text().splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) < 2 or fields[1] == "NA":
                continue
            pos = int(fields[0])
            if lo <= pos <= hi:
                probs.append(float(fields[1]))
    if len(probs) != target_len:
        raise ViennaError(
            f"RNAplfold returned {len(probs)} usable target rows; expected {target_len}"
        )
    return sum(probs) / len(probs)


def rnaplfold_accessibility(
    rung: Mapping,
    *,
    rnaplfold: Path,
    genome_fasta: str | Path,
    samtools: Path | None = None,
    pad: int = ACCESSIBILITY_PAD,
    window: int = ACCESSIBILITY_WINDOW,
    span: int = ACCESSIBILITY_SPAN,
) -> float:
    """Mean per-base unpaired probability of the target IN GENOMIC CONTEXT (0..1).

    The rung's target window is re-fetched with ``pad`` nt of genomic flank on
    each side and folded by RNAplfold, so accessibility reflects whether the
    genomic ISS-N1 site is open in pre-mRNA context rather than an isolated
    fragment trivially folding on itself. Reuses the ladder's proven
    pysam->samtools fetch. Fail-closed: a missing FASTA, a padded window that
    runs off the contig, or a frame/strand mismatch (the rung target not sitting
    at the expected offset) raises -- the caller maps that to ``None``, never a
    fabricated accessibility.
    """
    chrom = str(rung["chrom"])
    start_1based = int(rung["start"])
    end_1based = int(rung["end"])
    strand = str(rung["strand"])
    target = str(rung["target_5to3"]).upper()
    target_len = end_1based - start_1based + 1
    if target_len != len(target):
        raise ViennaError(
            f"rung target length {len(target)} != coordinate span {target_len}"
        )

    fetch_start_0base = (start_1based - 1) - pad
    if fetch_start_0base < 0:
        raise ViennaError("padded window runs off the chromosome start")
    spec = WindowSpec(
        exon_key="acc_window",
        chrom=chrom,
        start_0base=fetch_start_0base,
        end=end_1based + pad,
        strand=strand,
    )
    padded = _fetch_transcript_sequences(
        [spec], Path(genome_fasta), samtools=_resolve_samtools(samtools)
    )["acc_window"].upper()

    # With a symmetric pad the transcript-strand target sits at offset ``pad``
    # regardless of strand: on '-' the flanks reverse-complement-swap but keep
    # equal length. Verify it before trusting the offset -- a fail-closed frame
    # check that also catches any boundary clamping.
    if len(padded) != target_len + 2 * pad:
        raise ViennaError(
            f"padded window is {len(padded)} nt; expected {target_len + 2 * pad}"
        )
    if padded[pad : pad + target_len] != target:
        raise ViennaError(
            "padded window does not contain the rung target at the expected offset"
        )

    return mean_unpaired_over_window(
        padded, pad, target_len, rnaplfold=rnaplfold, window=window, span=span
    )


# --- D. ASO self-structure (RNAfold) -----------------------------------------


def rnafold_mfe(aso_5to3: str, *, rnafold: Path) -> float:
    """ViennaRNA ``RNAfold`` minimum-free-energy of the ASO alone (kcal/mol)."""
    with tempfile.TemporaryDirectory() as tmp:
        out = _run([rnafold, "--noPS"], stdin=f"{_to_rna(aso_5to3)}\n", cwd=Path(tmp))
    return _last_paren_energy(out, tool="RNAfold")


# --- E. off-target specificity (seed-and-extend) -----------------------------


def off_target_hits(
    aso_5to3: str,
    reference: str | Mapping[str, str],
    *,
    max_mismatch: int | None = None,
    max_edit: int | None = None,
) -> tuple[HybridizationHit, ...]:
    """Near-complement hits of the ASO in ``reference``.

    An ASO hybridises to a sense window equal to its reverse complement (i.e. the
    target orientation), so the scan query is ``reverse_complement(aso)`` -- the
    same convention Arm C's exon scan uses. Returns the raw hit tuple; callers
    take ``len(...)`` for a burden count. Exactly one of ``max_mismatch`` /
    ``max_edit`` selects Hamming vs edit mode (see :func:`seed_and_extend`).
    """
    query = reverse_complement(aso_5to3)
    return seed_and_extend(
        query, reference, max_mismatch=max_mismatch, max_edit=max_edit
    )


# --- assembler ---------------------------------------------------------------

_SCORE_KEYS = (
    "intarna_dg",
    "intarna_no_duplex",
    "rnaduplex_dg",
    "rnacofold_dg",
    "rnaup_dg",
    "rnaplfold_accessibility",
    "rnafold_mfe",
    "off_target_hits",
)


def evaluate_ladder(
    rungs: Sequence[Mapping],
    *,
    intarna: Path | None = None,
    vienna_bindir: Path | None = None,
    genome_fasta: str | Path | None = None,
    samtools: Path | None = None,
    off_target_reference: str | Mapping[str, str] | None = None,
    off_target_max_edit: int | None = 5,
    off_target_max_mismatch: int | None = None,
    accessibility_pad: int = ACCESSIBILITY_PAD,
    accessibility_window: int = ACCESSIBILITY_WINDOW,
    accessibility_span: int = ACCESSIBILITY_SPAN,
) -> list[dict]:
    """Attach the per-mer CPU score panel to each ladder rung.

    Every input rung dict (from :func:`build_ladder`) is copied and extended with
    the score columns in ``_SCORE_KEYS``. A backend that is unavailable leaves its
    column ``None`` for every rung rather than aborting the sweep -- degradation
    is silent-but-visible via the ``None`` value, never a fabricated score. When
    ``off_target_reference`` is ``None`` the off-target column stays ``None``
    (the genome/transcriptome scan is an opt-in integration path, not run by
    default). ``rnaplfold_accessibility`` requires ``genome_fasta`` (the target is
    scored in genomic context); without it that column stays ``None`` rather than
    falling back to the misleading isolated-window fold. Enriched rows retain the
    ``candidate`` honesty state.
    """
    intarna_exe = resolve_intarna(intarna)
    vienna = {tool: find_vienna(tool, vienna_bindir) for tool in VIENNA_TOOLS}

    enriched: list[dict] = []
    for rung in rungs:
        aso = rung["aso_5to3"]
        target = rung["target_5to3"]
        row = dict(rung)
        row["state"] = row.get("state", _CANDIDATE)

        # A. IntaRNA on-target
        if intarna_exe is not None:
            try:
                dg, no_duplex = intarna_on_target(aso, target, intarna=intarna_exe)
            except IntaRNAError:
                dg, no_duplex = None, None
            row["intarna_dg"], row["intarna_no_duplex"] = dg, no_duplex
        else:
            row["intarna_dg"], row["intarna_no_duplex"] = None, None

        # B/C/D. ViennaRNA panel
        row["rnaduplex_dg"] = _safe(
            vienna["RNAduplex"], lambda exe: rnaduplex_dg(aso, target, rnaduplex=exe)
        )
        row["rnacofold_dg"] = _safe(
            vienna["RNAcofold"], lambda exe: rnacofold_dg(aso, target, rnacofold=exe)
        )
        row["rnaup_dg"] = _safe(
            vienna["RNAup"], lambda exe: rnaup_dg(aso, target, rnaup=exe)
        )
        if genome_fasta is not None and vienna["RNAplfold"] is not None:
            try:
                row["rnaplfold_accessibility"] = rnaplfold_accessibility(
                    rung,
                    rnaplfold=vienna["RNAplfold"],
                    genome_fasta=genome_fasta,
                    samtools=samtools,
                    pad=accessibility_pad,
                    window=accessibility_window,
                    span=accessibility_span,
                )
            except (ViennaError, GenomeFastaError):
                row["rnaplfold_accessibility"] = None
        else:
            row["rnaplfold_accessibility"] = None
        row["rnafold_mfe"] = _safe(
            vienna["RNAfold"], lambda exe: rnafold_mfe(aso, rnafold=exe)
        )

        # E. off-target burden (opt-in)
        if off_target_reference is not None:
            hits = off_target_hits(
                aso,
                off_target_reference,
                max_mismatch=off_target_max_mismatch,
                max_edit=off_target_max_edit,
            )
            row["off_target_hits"] = len(hits)
        else:
            row["off_target_hits"] = None

        enriched.append(row)
    return enriched


def _safe(exe: Path | None, call):
    """Run ``call(exe)`` if the binary exists; on absence or failure return None."""
    if exe is None:
        return None
    try:
        return call(exe)
    except ViennaError:
        return None
