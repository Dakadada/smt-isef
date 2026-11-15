"""Two raw hybridization-energy models; no wet-lab or derived efficacy claim."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import math
import re
import subprocess
import tempfile
from typing import Mapping

from feng.titration.evaluate import intarna_on_target

from .sequences import to_rna


_DUPLEX_LINE = re.compile(
    r"^\s*(?P<structure>\S+).*\(\s*(?P<energy>-?\d+(?:\.\d+)?)\s*\)\s*$"
)


class DuplexError(RuntimeError):
    """RNAduplex was unavailable, failed, or returned an unparseable result."""


def intarna_many_targets(
    aso_5to3: str,
    targets: Mapping[str, str],
    *,
    intarna: Path,
    threads: int = 1,
) -> dict[str, tuple[float, bool]]:
    """Score isolated targets in one IntaRNA call using wrapper-equivalent settings.

    The ASO remains IntaRNA's query (sequence 2) and every genomic target
    variant is target sequence 1, matching :func:`intarna_on_target`. IntaRNA
    omits target rows having no favourable interaction; those rows retain its
    documented ``0.0`` plus an explicit no-duplex flag.
    """
    executable = Path(intarna)
    if not executable.is_file():
        raise DuplexError(f"IntaRNA executable not found: {executable}")
    if not targets:
        raise ValueError("targets must not be empty")
    if threads < 1 or threads > 64:
        raise ValueError("threads must be in [1, 64]")
    for target_id in targets:
        if not target_id or any(char in target_id for char in ";\t\r\n "):
            raise ValueError(f"invalid target identifier: {target_id!r}")

    with tempfile.TemporaryDirectory(prefix="item3-intarna-") as tmp:
        target_fasta = Path(tmp) / "targets.fa"
        with target_fasta.open("w") as handle:
            for target_id, sequence in targets.items():
                handle.write(f">{target_id}\n{to_rna(sequence)}\n")
        cmd = [
            str(executable),
            "-q",
            to_rna(aso_5to3),
            "-t",
            str(target_fasta),
            "--outMode",
            "C",
            "--outCsvCols",
            "id1,E",
            "-n",
            "1",
            "--threads",
            str(threads),
        ]
        try:
            completed = subprocess.run(cmd, text=True, capture_output=True, check=True)
        except (OSError, subprocess.CalledProcessError) as exc:
            stderr = getattr(exc, "stderr", "")
            raise DuplexError(f"IntaRNA failed: {stderr.strip() or exc}") from exc

    result = {target_id: (0.0, True) for target_id in targets}
    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    if not lines or lines[0].split(";") != ["id1", "E"]:
        raise DuplexError(f"unexpected IntaRNA CSV header: {completed.stdout!r}")
    for line in lines[1:]:
        fields = line.split(";")
        if len(fields) != 2 or fields[0] not in result:
            raise DuplexError(f"unexpected IntaRNA CSV row: {line!r}")
        if result[fields[0]] != (0.0, True):
            raise DuplexError(f"duplicate IntaRNA target row: {fields[0]!r}")
        try:
            energy = float(fields[1])
        except ValueError as exc:
            raise DuplexError(f"invalid IntaRNA energy row: {line!r}") from exc
        result[fields[0]] = (energy, energy == 0.0)
    return result


@dataclass(frozen=True, slots=True)
class DuplexScore:
    structure: str
    energy_kcal_mol: float
    paired_bases: int
    raw_reference_dg_kcal_mol: float
    salt_correction_kcal_mol: float
    temperature_c: float = 37.0
    salt_molar: float | None = None
    reference_salt_molar: float = 1.021
    chemistry: str = "RNA:RNA"


@dataclass(frozen=True, slots=True)
class PairScore:
    intarna_dg_kcal_mol: float
    intarna_no_duplex: bool
    nn_dg_kcal_mol: float
    nn_raw_reference_dg_kcal_mol: float
    nn_salt_correction_kcal_mol: float
    nn_structure: str
    nn_paired_bases: int
    nn_salt_molar: float
    temperature_c: float


def parse_rnaduplex(stdout: str) -> DuplexScore:
    lines = [line for line in stdout.splitlines() if line.strip()]
    if not lines:
        raise DuplexError("could not parse empty RNAduplex output")
    match = _DUPLEX_LINE.match(lines[-1])
    if match is None:
        raise DuplexError(f"could not parse RNAduplex output: {stdout!r}")
    structure = match.group("structure")
    return DuplexScore(
        structure=structure,
        energy_kcal_mol=float(match.group("energy")),
        paired_bases=structure.count("("),
        raw_reference_dg_kcal_mol=float(match.group("energy")),
        salt_correction_kcal_mol=0.0,
    )


def monovalent_salt_correction_kcal_mol(
    paired_bases: int,
    salt_molar: float,
    *,
    temperature_c: float = 37.0,
    reference_salt_molar: float = 1.021,
) -> float:
    """Return the SantaLucia entropy-derived ΔG salt adjustment.

    ViennaRNA 2.7.2 reports its default parameter reference as 1.021 M but its
    RNAduplex MFE executable does not change output when only ``--salt`` is
    varied. We therefore apply ΔS = 0.368*(N-1)*ln([M+]/[M+]ref) cal/K/mol and
    convert -TΔS to kcal/mol. ``N`` is the number of formed base pairs parsed
    from the predicted duplex, which is explicit for partial scramble duplexes.
    """
    if paired_bases < 1:
        raise ValueError("paired_bases must be positive")
    if salt_molar <= 0 or reference_salt_molar <= 0:
        raise ValueError("salt concentrations must be positive")
    temperature_k = temperature_c + 273.15
    delta_s = (
        0.368
        * (paired_bases - 1)
        * math.log(salt_molar / reference_salt_molar)
    )
    return -temperature_k * delta_s / 1000.0


def effective_monovalent_from_mixed_salt(
    monovalent_molar: float,
    divalent_molar: float,
) -> float:
    """Return the von Ahsen (2001) Na+-equivalent for a mixed mono/divalent buffer.

    [Na+]_eq (mM) = [monovalent] (mM) + 120 * sqrt([divalent] (mM)); returned in M.
    This folds a divalent (Mg2+) contribution into the single effective monovalent
    concentration consumed by ``monovalent_salt_correction_kcal_mol``. It is a
    DNA-derived empirical equivalence applied here to an RNA:RNA duplex, so the
    resulting ΔG remains a heuristic salt-adjusted value, not a validated RNA/TBE
    measurement (same caveat class as the monovalent correction itself).
    """
    if monovalent_molar < 0 or divalent_molar < 0:
        raise ValueError("salt concentrations must be non-negative")
    if monovalent_molar <= 0 and divalent_molar <= 0:
        raise ValueError("at least one salt concentration must be positive")
    monovalent_mm = monovalent_molar * 1000.0
    divalent_mm = divalent_molar * 1000.0
    na_eq_mm = monovalent_mm + 120.0 * math.sqrt(divalent_mm)
    return na_eq_mm / 1000.0


def rnaduplex_score(
    aso_5to3: str,
    target_5to3: str,
    *,
    rnaduplex: Path,
    salt_molar: float,
    temperature_c: float = 37.0,
    parameter_file: Path | str | None = None,
    chemistry: str = "RNA:RNA",
) -> DuplexScore:
    """Run ViennaRNA's Turner-NN duplex model at explicit temperature and salt.

    ``parameter_file`` makes chemistry parameters replaceable. The bundled
    default is RNA:RNA; a future DNA:RNA rerun must provide a validated hybrid
    parameter file rather than silently treating Vienna's DNA preset as hybrid.
    """
    if chemistry not in {"RNA:RNA", "DNA:RNA"}:
        raise ValueError(f"unsupported chemistry: {chemistry}")
    if chemistry == "DNA:RNA" and parameter_file is None:
        raise ValueError("DNA:RNA requires an explicit validated hybrid parameter file")
    if not (0 < salt_molar <= 2):
        raise ValueError("salt_molar must be in (0, 2]")
    executable = Path(rnaduplex)
    if not executable.is_file():
        raise DuplexError(f"RNAduplex executable not found: {executable}")
    cmd = [str(executable), "--temp", str(temperature_c)]
    if parameter_file is not None:
        cmd.extend(["--paramFile", str(parameter_file)])
    try:
        completed = subprocess.run(
            cmd,
            input=f"{to_rna(aso_5to3)}\n{to_rna(target_5to3)}\n",
            text=True,
            capture_output=True,
            check=True,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        stderr = getattr(exc, "stderr", "")
        raise DuplexError(f"RNAduplex failed: {stderr.strip() or exc}") from exc
    parsed = parse_rnaduplex(completed.stdout)
    correction = monovalent_salt_correction_kcal_mol(
        parsed.paired_bases, salt_molar, temperature_c=temperature_c
    )
    return DuplexScore(
        structure=parsed.structure,
        energy_kcal_mol=parsed.energy_kcal_mol + correction,
        paired_bases=parsed.paired_bases,
        raw_reference_dg_kcal_mol=parsed.energy_kcal_mol,
        salt_correction_kcal_mol=correction,
        temperature_c=temperature_c,
        salt_molar=salt_molar,
        chemistry=chemistry,
    )


def score_pair(
    aso_5to3: str,
    target_5to3: str,
    *,
    intarna: Path,
    rnaduplex: Path,
    salt_molar: float,
    temperature_c: float = 37.0,
    parameter_file: Path | str | None = None,
    chemistry: str = "RNA:RNA",
) -> PairScore:
    if temperature_c != 37.0:
        raise ValueError(
            "the reused intarna_on_target wrapper is fixed to IntaRNA's 37 C default"
        )
    intarna_dg, no_duplex = intarna_on_target(
        aso_5to3, target_5to3, intarna=Path(intarna)
    )
    nn = rnaduplex_score(
        aso_5to3,
        target_5to3,
        rnaduplex=rnaduplex,
        salt_molar=salt_molar,
        temperature_c=temperature_c,
        parameter_file=parameter_file,
        chemistry=chemistry,
    )
    return PairScore(
        intarna_dg_kcal_mol=intarna_dg,
        intarna_no_duplex=no_duplex,
        nn_dg_kcal_mol=nn.energy_kcal_mol,
        nn_raw_reference_dg_kcal_mol=nn.raw_reference_dg_kcal_mol,
        nn_salt_correction_kcal_mol=nn.salt_correction_kcal_mol,
        nn_structure=nn.structure,
        nn_paired_bases=nn.paired_bases,
        nn_salt_molar=salt_molar,
        temperature_c=temperature_c,
    )
