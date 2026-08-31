"""Absence and failure must never share a return value.

EVERY silent failure this project has catalogued has the same shape: a code
path that returns a PLAUSIBLE VALUE where it should have raised. The value then
flows downstream and gets reported as a result.

  - ISS-N1 vs GC-core footprint mixup      -> a confident inverted sign
  - UNC13A exon table truncated 2 of 44    -> printed its own limitation, did
                                              not enforce it, wrong gene length
  - SMN1/SMN2 aligner tie-break            -> PSI = 1.000 in every sample
  - P hallmarks scanned one strand         -> an inverted arm guaranteed to null
  - prereg panel pseudo-replicated         -> one site measured six times
  - pose grader summed altlocs             -> W53 as 39 atoms vs a 32-atom SMILES
  - ligand identity matched by name        -> every real comparison REFUSED
  - ligand size guard never fired          -> strontium ions scored as ligands
  - IntaRNA absent from PATH               -> EVERY oligo "forms no duplex"

The last one is the sharpest: ``goldilocks_offtarget.feasible`` returns None to
mean "no duplex forms", and a missing binary made it return None for every
candidate. A wet-lab panel built on that would have reported that no ASO works.

These tests pin the distinction where it is load-bearing. They are cheap and
they are meant to grow: when a new silent failure is found, it gets a case here.
"""

from __future__ import annotations

import pytest


def test_feasible_raises_on_missing_intarna_rather_than_returning_none(monkeypatch):
    """A missing binary must not be indistinguishable from an infeasible oligo."""
    import shutil

    from smt.validate import goldilocks_offtarget as g

    real = shutil.which

    def fake(name, *a, **k):
        return None if name == "IntaRNA" else real(name, *a, **k)

    monkeypatch.setattr(g.shutil if hasattr(g, "shutil") else shutil, "which", fake)
    monkeypatch.setattr(shutil, "which", fake)

    w = g.Winner(length=6, start=32059, worst_fom=1.0,
                 target_rna="CUGCCA", aso="UGGCAG", edge_excluded=True)
    with pytest.raises(RuntimeError, match="IntaRNA"):
        g.feasible(w)


def test_feasible_returns_a_real_energy_when_the_tool_is_present():
    """The mirror: the guard must not simply always raise.

    UGGCAG over the GC core binds at -9.09 kcal/mol -- the value recorded in
    goldilocks_offtarget's own docstring. If this drifts, the feasibility
    pre-gate is no longer measuring what it claims.
    """
    import shutil

    if shutil.which("IntaRNA") is None:
        pytest.skip("IntaRNA not on PATH in this environment")

    from smt.validate.goldilocks_offtarget import Winner, feasible

    w = Winner(length=6, start=32059, worst_fom=1.0,
               target_rna="CUGCCA", aso="UGGCAG", edge_excluded=True)
    out = feasible(w)
    assert out is not None
    assert out[0] == pytest.approx(-9.09, abs=0.05)


def test_smn2_panel_is_never_silently_empty():
    """A panel of zero candidates must be impossible to mistake for a result.

    The first run of smn2_panel emitted '0 candidates, 36 rejected' purely
    because IntaRNA was off PATH. Shipping that would have told a wet lab that
    every structure-ranked ASO was infeasible.
    """
    import shutil

    if shutil.which("IntaRNA") is None:
        pytest.skip("IntaRNA not on PATH in this environment")

    from smt.validate.smn2_panel import build

    panel, rejected = build(top=4)
    assert panel, "an empty panel with a working IntaRNA is a bug, not a finding"
    assert all(r["dG_hybrid"] < 0 for r in panel)


def test_a_5ss_refuses_a_footprint_overlapping_its_own_readout():
    """Forcing an occupied site UNPAIRED reads as trivially accessible.

    The framework models INDIRECT structural effects only, so a footprint that
    overlaps its readout must raise rather than return a flattering number.
    """
    from smt.seq.construct import SS5_EXON7
    from smt.seq.reference import Interval
    from smt.thermo.constrained import a_5ss

    with pytest.raises(Exception):
        a_5ss(footprint=Interval(SS5_EXON7.start, SS5_EXON7.end))
