"""Read-only, exon-3-specific workbook audit. No inference of missing bands.

All signal fractions are conditional on the workbook's unverified band roles.
Run with the feng-env interpreter; see README.md for the command.
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import math
import operator
from pathlib import Path
import re

from openpyxl import load_workbook


def finite_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def signal_fraction(included, excluded):
    """Undetected/zero bands have no usable quantitative inclusion estimate."""
    if not all(finite_number(v) and v > 0 for v in (included, excluded)):
        return None
    return 100 * included / (included + excluded)


class FormulaReader:
    """Evaluate only arithmetic and cell references; never execute workbook code."""

    reference = re.compile(r"(?:(?:'([^']+)'|([A-Za-z_][A-Za-z_0-9 ]*))!)?(\$?[A-Z]{1,3}\$?\d+)")

    def __init__(self, workbook):
        self.workbook = workbook

    def value(self, sheet, coordinate, active=None):
        active = set() if active is None else active
        key = (sheet, coordinate)
        if key in active:
            raise ValueError('cyclic formula reference')
        value = self.workbook[sheet][coordinate].value
        if not isinstance(value, str) or not value.startswith('='):
            if not finite_number(value):
                raise ValueError(f'non-numeric cell {sheet}!{coordinate}')
            return value
        active = active | {key}
        refs = {}

        def replace(match):
            name = f'ref_{len(refs)}'
            refs[name] = self.value(match[1] or match[2] or sheet,
                                    match[3].replace('$', ''), active)
            return name

        expression = self.reference.sub(replace, value[1:])
        binary = {ast.Add: operator.add, ast.Sub: operator.sub,
                  ast.Mult: operator.mul, ast.Div: operator.truediv}

        def evaluate(node):
            if isinstance(node, ast.Constant) and finite_number(node.value):
                return node.value
            if isinstance(node, ast.Name) and node.id in refs:
                return refs[node.id]
            if isinstance(node, ast.BinOp) and type(node.op) in binary:
                return binary[type(node.op)](evaluate(node.left), evaluate(node.right))
            if isinstance(node, ast.UnaryOp) and isinstance(node.op, (ast.UAdd, ast.USub)):
                return evaluate(node.operand) * (-1 if isinstance(node.op, ast.USub) else 1)
            raise ValueError('unsupported formula syntax')

        answer = evaluate(ast.parse(expression, mode='eval').body)
        if not finite_number(answer):
            raise ValueError('nonfinite formula result')
        return answer


def status(expected, actual, tolerance=1e-6):
    if expected is None:
        return 'undefined_missing_or_undetected_band'
    if not finite_number(actual):
        return 'unavailable'
    return 'match' if abs(expected - actual) <= tolerance else 'discrepancy'


def parse_workbook(path):
    path = Path(path)
    # Deliberately refuse the prohibited workbook, even on accidental invocation.
    if 'exon7' in path.name.lower().replace('_', ''):
        raise ValueError('exon-7 workbook is excluded')
    original = path.read_bytes()
    workbook = load_workbook(path, data_only=False)
    cached = load_workbook(path, data_only=True)
    required = {'Read me', 'Summary', 'Quantification', 'Raw 7-7', 'Raw 7-9', 'Methods'}
    if set(workbook.sheetnames) != required:
        raise ValueError('unexpected workbook sheets; inspect schema before parsing')
    reader = FormulaReader(workbook)
    cells, metadata = [], {}
    for sheet in workbook:
        for row in sheet:
            for cell in row:
                if cell.value is None:
                    continue
                cells.append({'sheet': sheet.title, 'cell': cell.coordinate,
                              'value': cell.value,
                              'cached_value': cached[sheet.title][cell.coordinate].value,
                              'evidence_class': 'observed_workbook_cell'})
                if sheet.title in ('Read me', 'Methods'):
                    metadata[f'{sheet.title}!{cell.coordinate}'] = cell.value
    bands = {}
    raw_count = {}
    for sheet_name in ('Raw 7-7', 'Raw 7-9'):
        sheet = workbook[sheet_name]
        gel = sheet_name.removeprefix('Raw ')
        raw_count[sheet_name] = 0
        for r in range(5, sheet.max_row + 1):
            lane, band = sheet.cell(r, 1).value, sheet.cell(r, 2).value
            if lane is None:
                continue
            raw_count[sheet_name] += 1
            if lane == 1:  # retain ladder in inventory, exclude from fractions
                continue
            key = (gel, lane)
            entries = bands.setdefault(key, {})
            if band in entries:
                raise ValueError(f'duplicate band {key}, {band}')
            entries[band] = {'volume': sheet.cell(r, 7).value,
                             'mw': sheet.cell(r, 6).value,
                             'rf': sheet.cell(r, 5).value,
                             'ratio': sheet.cell(r, 8).value,
                             'label': sheet.cell(r, 3).value,
                             'source': f'{sheet_name}!A{r}:I{r}'}
            if band not in (1, 2):
                raise ValueError(f'unexpected sample band {key}, {band}')

    summary = {}
    s = workbook['Summary']
    for r in range(5, s.max_row + 1):
        if not isinstance(s.cell(r, 1).value, int):
            continue
        key = (s.cell(r, 2).value, s.cell(r, 1).value)
        if key in summary:
            raise ValueError(f'duplicate summary lane {key}')
        summary[key] = r
    q = workbook['Quantification']
    positional, audit, issues = [], [], []
    seen = set()
    for r in range(5, q.max_row + 1):
        gel, lane, condition = [q.cell(r, c).value for c in (1, 2, 3)]
        if gel is None:
            continue
        key = (gel, lane)
        if key in seen or key not in bands or key not in summary:
            raise ValueError(f'duplicate or unmatched quantification lane {key}')
        seen.add(key)
        entries = bands[key]
        inc, exc = entries.get(1, {}), entries.get(2, {})
        i, e = inc.get('volume'), exc.get('volume')
        independent = signal_fraction(i, e)
        sr = summary[key]
        start, end = s.cell(sr, 4).value, s.cell(sr, 5).value
        position_match = re.search(r'\((?:nt\s*)?(\d+)-(\d+)\)', str(condition))
        if position_match and (start, end) != tuple(map(int, position_match.groups())):
            issues.append(f'position_label_discrepancy: Quantification!C{r} vs Summary!D{sr}:E{sr}')
        if start is not None and (not isinstance(start, int) or not isinstance(end, int)
                                  or end - start + 1 != 18):
            issues.append(f'position_length_discrepancy: Summary!D{sr}:E{sr}')
        for band in entries.values():
            # Raw B label omits "nt"; compare parsed interval rather than string.
            raw_match = re.search(r'\((?:nt\s*)?(\d+)-(\d+)\)', str(band['label']))
            if raw_match and position_match and raw_match.groups() != position_match.groups():
                issues.append(f'raw_position_label_discrepancy: {band["source"]}')
        imw, emw = inc.get('mw'), exc.get('mw')
        row = dict(gel=gel, lane=lane, condition=condition, nt_start=start, nt_end=end,
                   target_length_from_labels=None if start is None else end-start+1,
                   included_volume=i, excluded_volume=e,
                   included_status='absent_from_export' if not inc else ('reported_positive' if finite_number(i) and i > 0 else 'undetected_or_invalid'),
                   excluded_status='absent_from_export' if not exc else ('reported_positive' if finite_number(e) and e > 0 else 'undetected_or_invalid'),
                   signal_fraction_pct=independent,
                   included_size_proxy=imw, excluded_size_proxy=emw,
                   band_size_difference_proxy=imw-emw if all(finite_number(v) for v in (imw, emw)) else None,
                   expected_included_bp=None, expected_excluded_bp=None,
                   biological_replicate_id=None, technical_replicate_id=None,
                   experimental_unit_id=f'{gel}:lane{lane}',
                   band_assignment='supplied_unverified',
                   differing_product_lengths='needs_human_confirmation',
                   on_gel_control_available=gel == '7-7',
                   source_included=inc.get('source'), source_excluded=exc.get('source'),
                   source_position=f'Summary!D{sr}:E{sr}')
        positional.append(row)
        if exc and inc and not (finite_number(inc['rf']) and finite_number(exc['rf'])
                                and inc['rf'] < exc['rf'] and finite_number(imw)
                                and finite_number(emw) and imw > emw):
            issues.append(f'band_order_discrepancy: {gel} lane {lane}')
        for sheet, coordinate in [('Quantification', f'G{r}'), ('Summary', f'F{sr}')]:
            try:
                formula_value = reader.value(sheet, coordinate)
                formula_error = None
            except (ValueError, KeyError, SyntaxError, ZeroDivisionError) as error:
                formula_value, formula_error = None, str(error)
            cache = cached[sheet][coordinate].value
            audit.append(dict(gel=gel, lane=lane, source=f'{sheet}!{coordinate}',
                              formula=workbook[sheet][coordinate].value,
                              formula_value=formula_value, formula_error=formula_error,
                              cached_value=cache, independent_value=independent,
                              independent_status='defined_conditional_signal_fraction' if independent is not None else 'undefined_missing_or_undetected_band',
                              formula_status=status(independent, formula_value),
                              cached_status=status(independent, cache) if independent is not None else ('unavailable' if cache is None else 'workbook_only_missing_band'),
                              cache_vs_formula_status=status(formula_value, cache),
                              export_ratio=inc.get('ratio'),
                              export_status=status(independent, inc.get('ratio'), .000500001),
                              difference_formula_pp=None if independent is None or formula_value is None else formula_value-independent,
                              difference_cached_pp=None if independent is None or not finite_number(cache) else cache-independent))
        # Audit quantification volume references independently of their downstream PSI.
        for col, expected in [('D', i), ('E', e)]:
            try:
                actual = reader.value('Quantification', f'{col}{r}')
            except (ValueError, KeyError, SyntaxError, ZeroDivisionError):
                actual = None
            if expected is None:
                issues.append(f'missing_band_workbook_value: Quantification!{col}{r}={actual}; independent=None')
            elif status(expected, actual) != 'match':
                issues.append(f'volume_reference_discrepancy: Quantification!{col}{r}')
    if seen != set(bands) or seen != set(summary):
        raise ValueError('raw/quantification/summary lane sets differ')

    schemas = {}
    observed = {'gel', 'lane', 'condition', 'nt_start', 'nt_end', 'included_volume',
                'excluded_volume', 'formula', 'cached_value', 'export_ratio'}
    missing = {'expected_included_bp', 'expected_excluded_bp',
               'biological_replicate_id', 'technical_replicate_id'}
    for name, records in [('positional', positional), ('psi_audit', audit)]:
        schemas[name] = {key: ('missing_needs_human_confirmation' if key in missing else
                               'proxy_workbook_MW_not_verified_product_size' if 'proxy' in key else
                               'observed_workbook_unverified' if key in observed else
                               'source_locator' if key.startswith('source') else
                               'derived_conditional_on_supplied_band_roles' if key in ('signal_fraction_pct', 'independent_value') else
                               'derived_audit_or_policy') for key in records[0]} if records else {}
    needs = [
        'ASO A-K sequences, chemistry, 5-to-3 orientation and strand; labels alone do not establish binding',
        'Full SMN2Sup construct sequence/version and exon-3 coordinate origin; verify labeled intervals against sequence',
        'Primer sequences, binding sites and expected inclusion/skipping amplicon sizes; construct specificity unverified',
        'Band identities via junction sequencing or product-size anchors; relative migration order is only internally consistent',
        'Differing product lengths: within-lane size separation and between-gel MW shifts; calibration vs alternate products unresolved',
        'Original gel images and PDF exports; workbook transcription and background/saturation/linear range unverified',
        'Detection limits for absent/weak bands; absent export row cannot establish zero transcript',
        'Biological versus technical experimental units, sample lineage, cell system, timing and replicate IDs; n=1 is workbook-reported only',
        'Matched on-gel controls for gel 7-9; no cross-gel normalized effect estimated',
        'Genome/transcript/construct alignment and assembly before any computational-feature comparison',
    ]
    return dict(source={'path': str(path.resolve()), 'bytes': len(original),
                        'sha256': hashlib.sha256(original).hexdigest(),
                        'sheet_names': workbook.sheetnames,
                        'evidence_class': 'observed_file_and_derived_hash'},
                metadata=metadata, cells=cells, raw_row_counts=raw_count,
                positional=positional, psi_audit=audit, issues=issues,
                needs_human_confirmation=needs, comparability='not comparable',
                evidence_schema=schemas,
                evidence_conventions={
                    'metadata': 'observed_workbook_statement_not_independent_verification',
                    'cells': 'observed_workbook_cells_including_formulas_and_cached_values',
                    'raw_row_counts': 'derived_count_including_ladder',
                    'issues': 'derived_audit_finding',
                    'needs_human_confirmation': 'missing_or_unverified_evidence',
                    'comparability': 'policy_no_verified_system_coordinate_compatibility',
                    'null': 'missing_or_undefined_never_imputed; see status columns',
                    'schema': 'column classes apply to every value in JSON and TSV, including nulls',
                    'figures': 'plots of positional columns; same evidence classes and units'})


def write_outputs(result, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output/'audit.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    for name in ('positional', 'psi_audit'):
        with (output/f'{name}.tsv').open('w', newline='') as handle:
            records = result[name]
            writer = csv.DictWriter(handle, fieldnames=list(result['evidence_schema'][name]), delimiter='\t')
            writer.writeheader()
            writer.writerows(records)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharey=True)
    for ax, gel in zip(axes, ('7-7', '7-9')):
        for row in result['positional']:
            if row['gel'] != gel or row['nt_start'] is None or row['signal_fraction_pct'] is None:
                continue
            y = row['signal_fraction_pct']
            ax.plot([row['nt_start'], row['nt_end']], [y, y], marker='|')
            ax.annotate(row['condition'].split(' (')[0], ((row['nt_start']+row['nt_end'])/2, y),
                        xytext=(0, 5), textcoords='offset points', ha='center', fontsize=8)
        ax.set(xlim=(0, 175), ylim=(0, 105), xlabel='Workbook exon-3 nt labels (unverified)',
               title=f'Gel {gel}'+(' — no on-gel controls' if gel == '7-9' else ''))
    axes[0].set_ylabel('Conditional Band 1 signal fraction (%)')
    fig.suptitle('Exploratory only: supplied band identities; n=1 reported per ASO')
    fig.text(.5, .01, 'No length correction; missing bands undefined; size shifts unresolved; no prospective/model comparison.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .05, 1, .93))
    fig.savefig(output/'positional.png', dpi=160)
    fig.savefig(output/'positional.svg')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('workbook', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = parse_workbook(args.workbook)
    write_outputs(result, args.output)
    print(json.dumps({'lanes': len(result['positional']), 'PSI_cells_audited': len(result['psi_audit']),
                      'issues': result['issues'], 'source': result['source']}, indent=2))


if __name__ == '__main__':
    main()
