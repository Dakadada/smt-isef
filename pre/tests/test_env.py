"""A0 environment smoke test.

Verifies every external dependency the thermodynamic core and steric lane rely
on is importable and functional. Run this before anything else; a silent
version mismatch here invalidates every downstream number.
"""

import shutil
import subprocess

import pytest


def test_viennarna_imports_and_folds():
    import RNA

    # Well-characterized hairpin; MFE must be negative and structure balanced.
    seq = "GGGGAAAACCCC"
    structure, mfe = RNA.fold(seq)
    assert len(structure) == len(seq)
    assert structure.count("(") == structure.count(")")
    assert mfe < 0


def test_viennarna_partition_function():
    import RNA

    fc = RNA.fold_compound("GGGGAAAACCCC")
    fc.pf()
    # Ensemble free energy must be at or below MFE.
    _, mfe = RNA.fold("GGGGAAAACCCC")
    _, efe = fc.pf()
    assert efe <= mfe + 1e-6


def test_viennarna_hard_constraints_available():
    """A2 depends on hard constraints to impose an effector footprint."""
    import RNA

    fc = RNA.fold_compound("GGGGAAAACCCC")
    # Force position 1 (1-indexed) to stay unpaired.
    fc.hc_add_up(1, RNA.CONSTRAINT_CONTEXT_ALL_LOOPS)
    structure, _ = fc.mfe()
    assert structure[0] == "."


def test_rnaplfold_binary_present():
    assert shutil.which("RNAplfold"), "RNAplfold needed for A2 accessibility"


def test_rnaup_binary_present():
    assert shutil.which("RNAup"), "RNAup needed for A2 accessibility"


def test_intarna_runs():
    exe = shutil.which("IntaRNA")
    assert exe, "IntaRNA needed for A2 hybridization energies"
    result = subprocess.run(
        [exe, "-t", "GGGGAAAACCCCGGGGAAAACCCC", "-q", "GGGGUUUUCCCC", "--outMode=C"],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stderr


def test_biopython_imports():
    from Bio.Seq import Seq

    assert str(Seq("ATGC").reverse_complement()) == "GCAT"


def test_freesasa_imports():
    """A4 steric shadow uses rolling-probe SASA on cryo-EM coordinates."""
    import freesasa

    assert hasattr(freesasa, "calc")


def test_gemmi_imports():
    """A4 loads Cas13 ternary complex structures."""
    import gemmi

    assert hasattr(gemmi, "read_structure")


def test_scientific_stack():
    import numpy
    import pandas
    import scipy

    assert numpy.__version__
    assert pandas.__version__
    assert scipy.__version__


@pytest.mark.parametrize("mod", ["RNA", "Bio", "freesasa", "gemmi", "numpy", "scipy", "pandas"])
def test_version_recorded(mod, capsys):
    """Print versions so the pytest log doubles as a provenance record."""
    import importlib

    m = importlib.import_module(mod)
    version = getattr(m, "__version__", None) or getattr(m, "__RNAVERSION__", "unknown")
    with capsys.disabled():
        print(f"  {mod}: {version}")
