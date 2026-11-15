"""Synthetic fixtures only: no laboratory workbook is read by tests."""
import importlib.util
from pathlib import Path
import sys

import pytest
from openpyxl import Workbook

SPEC = importlib.util.spec_from_file_location(
    "exon3_parse", Path(__file__).parents[1] / "exon3_parse.py")
mod = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = mod
SPEC.loader.exec_module(mod)


@pytest.mark.parametrize("i,e,expected", [(3, 1, 75), (0, 4, None),
    (4, 0, None), (None, 4, None), (4, None, None), (0, 0, None),
    (-1, 4, None), (float('nan'), 4, None), ('-', 4, None)])
def test_signal_fraction_preserves_undetected_and_invalid(i, e, expected):
    assert mod.signal_fraction(i, e) == expected


def fixture(tmp_path, missing=False, duplicate=False, bad_reference=False):
    w = Workbook()
    w.remove(w.active)
    for name in ['Read me', 'Summary', 'Quantification', 'Raw 7-7', 'Raw 7-9', 'Methods']:
        w.create_sheet(name)
    for name in ['Raw 7-7', 'Raw 7-9']:
        w[name].append(['title'])
        w[name].cell(4, 1, 'Lane #')
    raw = w['Raw 7-7']
    for r, values in enumerate([[5,1,'ASO A (nt 2-19)',1,.4,560,30,75,'-'],
                              [5,2,'ASO A (nt 2-19)',1,.6,290,10,25,'-']], 5):
        if missing and r == 6:
            continue
        for c, v in enumerate(values, 1): raw.cell(r,c,v)
    if duplicate:
        for c in range(1,10): raw.cell(7,c,raw.cell(5,c).value)
    q = w['Quantification']
    vals = ['7-7',5,'ASO A (nt 2-19)',"='Raw 7-7'!G5",
            0 if missing else "='Raw 7-7'!G6", '=D5+E5',
            '=100*E5/F5' if bad_reference else '=100*D5/F5', '=100*E5/F5']
    for c,v in enumerate(vals,1): q.cell(5,c,v)
    s = w['Summary']
    for c,v in enumerate([5,'7-7','ASO A (nt 2-19)',2,19,'=Quantification!G5'],1):
        s.cell(5,c,v)
    w['Methods']['A5']='Endpoint'
    w['Methods']['B5']='Synthetic metadata'
    p=tmp_path/'synthetic.xlsx'
    w.save(p)
    return p


def test_raw_recomputation_and_no_cached_value(tmp_path):
    result=mod.parse_workbook(fixture(tmp_path))
    row=result['positional'][0]
    assert row['signal_fraction_pct']==75
    assert row['nt_start']==2 and row['nt_end']==19
    assert row['biological_replicate_id'] is None
    assert row['technical_replicate_id'] is None
    assert row['expected_included_bp'] is None
    assert row['band_size_difference_proxy']==270
    assert result['metadata']['Methods!B5']=='Synthetic metadata'
    assert any(a['cached_status']=='unavailable' for a in result['psi_audit'])
    assert result['comparability']=='not comparable'


def test_missing_band_not_workbook_zero(tmp_path):
    result=mod.parse_workbook(fixture(tmp_path,missing=True))
    row=result['positional'][0]
    assert row['excluded_volume'] is None
    assert row['signal_fraction_pct'] is None
    assert row['excluded_status']=='absent_from_export'
    assert result['psi_audit'][0]['formula_value']==100
    assert result['psi_audit'][0]['independent_status']=='undefined_missing_or_undetected_band'


def test_duplicate_band_fails_closed(tmp_path):
    with pytest.raises(ValueError, match='duplicate'):
        mod.parse_workbook(fixture(tmp_path,duplicate=True))


def test_wrong_formula_detected_independently(tmp_path):
    result=mod.parse_workbook(fixture(tmp_path,bad_reference=True))
    assert result['positional'][0]['signal_fraction_pct']==75
    assert result['psi_audit'][0]['formula_status']=='discrepancy'


def test_stale_cache_reported(tmp_path):
    from zipfile import ZipFile, ZIP_DEFLATED
    p=fixture(tmp_path)
    with ZipFile(p) as z: parts={n:z.read(n) for n in z.namelist()}
    import re
    parts['xl/worksheets/sheet3.xml'], count = re.subn(
        rb'(<f>100\*D5/F5</f>)<v(?:\s*/>|></v>)',
        rb'\1<v>12</v>', parts['xl/worksheets/sheet3.xml'])
    assert count == 1
    with ZipFile(p,'w',ZIP_DEFLATED) as z:
        for n,v in parts.items(): z.writestr(n,v)
    result=mod.parse_workbook(p)
    assert result['psi_audit'][0]['cached_value']==12
    assert result['psi_audit'][0]['cached_status']=='discrepancy'


def test_output_schema_covers_every_value_and_source_unchanged(tmp_path):
    import hashlib
    p=fixture(tmp_path)
    before=hashlib.sha256(p.read_bytes()).hexdigest()
    result=mod.parse_workbook(p)
    mod.write_outputs(result,tmp_path/'out')
    assert hashlib.sha256(p.read_bytes()).hexdigest()==before
    assert set(result['positional'][0])==set(result['evidence_schema']['positional'])
    assert set(result['psi_audit'][0])==set(result['evidence_schema']['psi_audit'])
    assert (tmp_path/'out'/'positional.png').stat().st_size>0
    assert (tmp_path/'out'/'positional.tsv').exists()


def test_position_and_migration_conflicts_reported(tmp_path):
    from openpyxl import load_workbook
    p=fixture(tmp_path)
    w=load_workbook(p)
    w['Summary']['E5']=20
    w['Raw 7-7']['E6']=.3
    w.save(p)
    result=mod.parse_workbook(p)
    assert any('position_label_discrepancy' in x for x in result['issues'])
    assert any('position_length_discrepancy' in x for x in result['issues'])
    assert any('band_order_discrepancy' in x for x in result['issues'])


@pytest.mark.parametrize('formula', ['=SUM(D5:E5)', '=G5', '=1/0'])
def test_unsupported_cyclic_or_divzero_formula_has_no_score(tmp_path, formula):
    from openpyxl import load_workbook
    p=fixture(tmp_path)
    w=load_workbook(p)
    w['Quantification']['G5']=formula
    w.save(p)
    result=mod.parse_workbook(p)
    assert result['psi_audit'][0]['formula_value'] is None
    assert result['psi_audit'][0]['formula_error']
    assert result['positional'][0]['signal_fraction_pct']==75


def test_exon7_rejected_before_opening():
    with pytest.raises(ValueError, match='excluded'):
        mod.parse_workbook('/does/not/exist/SMN2_Exon7_7_27+28_Complete_Data.xlsx')
