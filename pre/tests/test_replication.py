"""Tests for Arm V: Ottesen 2021 independent replication.

The load-bearing risk in this arm is NOT the statistics -- it is the inverted
sign convention. The published dPSI is control - treated while every other
dPSI in this project is treated - control, so a silent failure here would flip
every conclusion while still producing plausible numbers. That is the same
class as P's single-strand scan and R3's SMN1/SMN2 tie-break, both of which
produced clean, publishable-looking, wrong answers.

Covers:
* the sign convention is RE-DERIVED from the pinned workbook, not restated:
  control - treated must reproduce the published column and treated - control
  must NOT, checked as a contrast so a symmetric bug cannot pass both;
* ``load_events`` returns PROJECT convention -- verified by recomputing
  treated - control from the raw PSI columns for real rows and matching;
* ``verify_sign_convention`` actually RAISES on a mis-signed workbook: a
  grader that cannot fail is not a grader, so the failure path is exercised
  on a synthetic sheet rather than assumed;
* the semantic cross-check, independent of the arithmetic one: in the
  exon-INCLUSION sheet the ASO must on balance INCREASE inclusion under
  project convention, and in the exon-SKIPPING sheet decrease it. If both the
  arithmetic and the semantics agree, the convention is not a coincidence;
* the retention ratio is invariant to the sign convention, which is why V2
  survives even if the flip were wrong;
* the design facts the pre-registration was written against still hold, so
  the prereg cannot silently stop matching the file it describes.
"""

from __future__ import annotations

import pytest

openpyxl = pytest.importorskip("openpyxl")

from smt.sieve.ottesen import (
    DEFAULT_XLSX,
    EVENT_SHEETS,
    SignConventionError,
    best_per_gene,
    load_events,
    verify_sign_convention,
)

pytestmark = pytest.mark.skipif(
    not DEFAULT_XLSX.exists(), reason="Ottesen 2021 workbook not on disk"
)


def _raw_rows():
    """Raw (control, hia_psi, published_dpsi) triples straight from the file."""
    wb = openpyxl.load_workbook(DEFAULT_XLSX, read_only=True)
    out = []
    for sheet in EVENT_SHEETS:
        ws = wb[sheet]
        it = ws.iter_rows(values_only=True)
        idx = {str(h): i for i, h in enumerate(next(it))}
        for row in it:
            try:
                out.append((float(row[idx["10mer PSI"]]),
                            float(row[idx["HiA PSI"]]),
                            float(row[idx["HiA ∆PSI"]])))
            except (TypeError, ValueError):
                continue
    return out


def test_published_dpsi_is_control_minus_treated_and_not_the_reverse():
    """Both directions checked, so a symmetric bug cannot satisfy the test."""
    rows = _raw_rows()
    assert len(rows) > 1000

    ct = sum(1 for c, h, d in rows if abs((c - h) - d) < 1e-3)
    tc = sum(1 for c, h, d in rows if abs((h - c) - d) < 1e-3)

    assert ct == len(rows), "control - treated must reproduce every row"
    assert tc / len(rows) < 0.01, "treated - control must NOT reproduce the column"


def test_load_events_returns_project_convention():
    """dpsi_hia must equal treated - control, recomputed from the PSI columns.

    Compared POSITIONALLY rather than by Ensembl ID: a gene contributes several
    rows (one per alternative exon), so keying on the ID silently collapses
    them. Zipping also asserts that load_events preserves file order.
    """
    events = load_events()

    expected = []
    wb = openpyxl.load_workbook(DEFAULT_XLSX, read_only=True)
    for sheet in EVENT_SHEETS:
        ws = wb[sheet]
        it = ws.iter_rows(values_only=True)
        idx = {str(h): i for i, h in enumerate(next(it))}
        for row in it:
            if row[idx["Gene name"]] is None:
                continue
            try:
                expected.append(float(row[idx["HiA PSI"]]) - float(row[idx["10mer PSI"]]))
            except (TypeError, ValueError):
                expected.append(None)

    assert len(expected) == len(events)

    checked = 0
    for ev, want in zip(events, expected):
        if want is None or ev.dpsi_hia is None:
            continue
        assert ev.dpsi_hia == pytest.approx(want, abs=1e-3)
        checked += 1
    assert checked > 1000


def test_verify_sign_convention_raises_on_a_mis_signed_workbook(tmp_path):
    """The guard must actually fire. A grader that cannot fail is not a grader."""
    header = ["Ensembl ID", "Gene name", "Chromosome", "Strand",
              "Alt exon start", "Alt exon end", "Up exon start", "Up exon end",
              "Down exon start", "Down exon end", "FDR", "10mer PSI",
              "HiA PSI", "LoA PSI", "HiA \u2206PSI", "LoA \u2206PSI"]
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    for sheet in EVENT_SHEETS:
        ws = wb.create_sheet(sheet)
        ws.append(header)
        for i in range(50):
            ctl, hia = 0.2, 0.8
            # treated - control: the WRONG convention for this file.
            ws.append([f"ENSG{i}", f"G{i}", "chr1", "+", 1, 2, 3, 4, 5, 6,
                       0.01, ctl, hia, 0.5, hia - ctl, 0.1])
    path = tmp_path / "bad.xlsx"
    wb.save(path)

    with pytest.raises(SignConventionError):
        verify_sign_convention(path)


def test_verify_sign_convention_accepts_a_correctly_signed_workbook(tmp_path):
    """The mirror of the above, so the guard is not simply always raising."""
    header = ["Ensembl ID", "Gene name", "Chromosome", "Strand",
              "Alt exon start", "Alt exon end", "Up exon start", "Up exon end",
              "Down exon start", "Down exon end", "FDR", "10mer PSI",
              "HiA PSI", "LoA PSI", "HiA \u2206PSI", "LoA \u2206PSI"]
    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    for sheet in EVENT_SHEETS:
        ws = wb.create_sheet(sheet)
        ws.append(header)
        for i in range(50):
            ctl, hia = 0.2, 0.8
            ws.append([f"ENSG{i}", f"G{i}", "chr1", "+", 1, 2, 3, 4, 5, 6,
                       0.01, ctl, hia, 0.5, ctl - hia, 0.1])
    path = tmp_path / "good.xlsx"
    wb.save(path)

    n, frac = verify_sign_convention(path)
    assert n == 200 and frac == 1.0


def test_semantics_agree_with_the_arithmetic():
    """Independent of the arithmetic check: EIN must rise, ESK must fall."""
    events = [e for e in load_events() if e.dpsi_hia is not None]
    ein = [e.dpsi_hia for e in events if e.etype == "EIN"]
    esk = [e.dpsi_hia for e in events if e.etype == "ESK"]
    assert ein and esk

    # Under project convention the exon-INCLUSION sheet must be mostly positive
    # (the ASO increased inclusion) and the skipping sheet mostly negative.
    assert sum(1 for v in ein if v > 0) / len(ein) > 0.9
    assert sum(1 for v in esk if v < 0) / len(esk) > 0.9


def test_retention_ratio_is_sign_convention_invariant():
    """V2's statistic survives even if the flip were wrong -- by construction."""
    events = [e for e in load_events()
              if e.dpsi_hia and e.dpsi_loa is not None][:200]
    assert events
    for ev in events:
        flipped = (-ev.dpsi_loa) / (-ev.dpsi_hia)
        assert ev.retention == pytest.approx(flipped, rel=1e-12)


def test_design_facts_the_preregistration_was_written_against():
    """docs/preregistration.md quotes these. They must not drift silently."""
    events = load_events()
    both = best_per_gene(events)
    assert len(events) == 4300
    assert len(both) == 2644

    n_pairs, frac = verify_sign_convention()
    assert n_pairs == 4300
    assert frac == 1.0
