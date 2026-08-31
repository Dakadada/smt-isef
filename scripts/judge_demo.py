#!/usr/bin/env python3
"""Build a fail-closed, offline judge artifact from real Stage-A results."""
from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import math
import os
import shutil
import tempfile
from pathlib import Path

BANNER = "PROSPECTIVE PREDICTION — NOT EXPERIMENTALLY VALIDATED"
HEADLINE = ("A compute-native SMN2 pipeline recovered the expected splice direction across five deposited "
            "steric intervention arms—ASO1, F18MOE, F20PMO, F18OMe, and dCas13d D[333,357]—with replicate observations "
            "within each arm; genuine scramble-MOE and "
            "non-targeting-guide null arms stayed near zero; and the sequence-committed SMN2 target yielded "
            "frozen prospective designs for follow-up testing.")
INTERVENTIONS = {"ASO1", "F18MOE", "F20PMO", "F18OMe", "dCas13d_D[333,357]", "scramble_MOE", "non_targeting_guide"}
EXPECTED = {"ASO1": "increase", "F18MOE": "increase", "F20PMO": "increase", "F18OMe": "increase",
            "dCas13d_D[333,357]": "decrease",
            "scramble_MOE": "no_change", "non_targeting_guide": "no_change"}


def sha256(data: bytes) -> str: return hashlib.sha256(data).hexdigest()


def rows_from(data: bytes):
    return list(csv.DictReader(data.decode("utf-8").splitlines(), delimiter="\t"))


def validate_panel(data: bytes):
    panel = json.loads(data)
    if (panel.get("schema_version") != "prospective-design-panel-report/2"
            or panel.get("prediction_status") != "prospective-not-experimentally-validated"
            or not panel.get("banner", "").startswith(BANNER)):
        raise ValueError("panel schema/status/banner is not prospective")
    designs = [x for x in panel.get("designs", []) if x.get("target_id") == "SMN2_exon7_inclusion"]
    if not designs: raise ValueError("panel contains zero real SMN2 designs")
    for design in designs:
        if not design.get("banner", "").startswith(BANNER): raise ValueError("design banner missing")
        if not design.get("candidates"): raise ValueError("SMN2 design contains no candidates")
    backends = []
    for design in designs:
        backend = design.get("provenance", {}).get("run_manifest", {}).get("backend", {})
        qualified_name = backend.get("qualified_name")
        if not qualified_name: raise ValueError("SMN2 design lacks folding-backend identity")
        backends.append(qualified_name)
    return panel, designs, tuple(sorted(set(backends)))


def select_evidence(tidy, effects):
    selected = [row for row in tidy if row.get("event") == "SMN2_exon7_inclusion"
                and row.get("intervention") in INTERVENTIONS and row.get("direction") != "MISSING"]
    if not selected: raise ValueError("no eligible SMN2 evidence rows")
    present = {row["intervention"] for row in selected}
    missing = INTERVENTIONS - present
    if missing: raise ValueError(f"required evidence arms absent: {sorted(missing)}")
    arm_counts = {intervention: 0 for intervention in INTERVENTIONS}
    for row in selected: arm_counts[row["intervention"]] += 1
    under_replicated = sorted(intervention for intervention, count in arm_counts.items() if count < 2)
    if under_replicated:
        raise ValueError(f"arm has fewer than 2 replicate observations: {under_replicated[0]}")
    null_effects = {}
    for row in selected:
        if EXPECTED[row["intervention"]] != "no_change": continue
        try:
            psi = float(row["psi"])
            control_mean_psi = float(row["control_mean_psi"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError("genuine null arm PSI values must be numeric") from exc
        if not (math.isfinite(psi) and math.isfinite(control_mean_psi)):
            raise ValueError("genuine null arm PSI values must be finite")
        key = (row.get("family"), row["intervention"])
        null_effects.setdefault(key, []).append(psi - control_mean_psi)
    for key, replicate_effects in null_effects.items():
        arm_effect = sum(replicate_effects) / len(replicate_effects)
        if abs(arm_effect) > 2:
            raise ValueError(f"genuine null arm exceeds 2 dPSI: {key}")
    effect_index = {}
    for effect in effects:
        key = (effect.get("family"), effect.get("sample"), effect.get("intervention"))
        if key in effect_index: raise ValueError(f"duplicate sample-effect identity: {key}")
        effect_index[key] = effect
    evidence = []
    for row in selected:
        expected = EXPECTED[row["intervention"]]; magnitude = float(row["magnitude"])
        if not math.isfinite(magnitude): raise ValueError("selected evidence magnitude must be finite")
        if expected == "increase" and magnitude <= 0: raise ValueError("positive-direction evidence mismatch")
        if expected == "decrease" and magnitude >= 0: raise ValueError("negative-direction evidence mismatch")
        if row.get("coverage_pass") != "true": raise ValueError("selected evidence lacks coverage")
        effect = effect_index.get((row["family"], row["sample"], row["intervention"]))
        if (effect is None or not effect.get("run_accession")
                or int(effect["informative_reads"]) < 20):
            raise ValueError("selected evidence lacks run identity/informative reads")
        evidence.append({**row, "run_accession": effect["run_accession"],
                         "informative_reads": int(effect["informative_reads"])})
    return evidence


def build(run_dir: Path, panel_path: Path, out: Path, timestamp: str, allow_incomplete: bool):
    complete = (run_dir / "STAGE_A_COMPLETE").exists()
    if not complete and not allow_incomplete:
        raise RuntimeError("Stage A is incomplete; pass --allow-incomplete-snapshot only for a labeled snapshot")
    tidy_path = run_dir / "work/harmonized/tidy_splicing.tsv"
    effects_path = run_dir / "work/psi/sample_effects.tsv"
    tidy_data, effects_data, panel_data = tidy_path.read_bytes(), effects_path.read_bytes(), panel_path.read_bytes()
    panel, designs, backends = validate_panel(panel_data)
    evidence = select_evidence(rows_from(tidy_data), rows_from(effects_data))
    excluded = [row for row in rows_from(tidy_data) if row.get("event") == "SMN2_exon7_inclusion"
                and row.get("direction") == "MISSING"]
    if complete and any(row.get("intervention") in INTERVENTIONS for row in excluded):
        raise ValueError("Stage A completion marker conflicts with MISSING SMN2 evidence rows")
    state = "COMPLETE STAGE-A RUN" if complete else "INCOMPLETE STAGE-A SNAPSHOT"
    rendered_headline = HEADLINE if complete else f"{state}: {HEADLINE} Unresolved arms, including scramble_PMO, are excluded."
    payload = {"schema_version":"judge-demo-evidence/1", "stage_a_complete":complete,
               "snapshot_label":state, "headline":rendered_headline, "observations":evidence,
               "excluded_missing_rows":excluded, "tidy_sha256":sha256(tidy_data),
               "sample_effects_sha256":sha256(effects_data)}
    candidates = [c for d in designs for c in d["candidates"] if c.get("pareto")]
    backend_label = ", ".join(backends)
    is_stub = any("stub" in name.lower() for name in backends)
    scoring_notice = (f"Scoring backend: {backend_label}. "
        + ("Deterministic placeholder only; no real thermodynamic folding or biological validation."
           if is_stub else "See frozen panel provenance for backend identity and validation scope."))
    cards = "".join(f"<article><h3>{html.escape(c['candidate']['candidate_id'])}</h3>"
        f"<strong>{BANNER}</strong><p>Target: {html.escape(c['candidate']['target_sequence'])}; "
        f"effector: {html.escape(c['candidate']['effector_class'])}; rank: {c.get('rank')}</p></article>" for c in candidates)
    points = "".join(f"<tr><td>{html.escape(x['intervention'])}</td><td>{html.escape(x['sample'])}</td>"
        f"<td>{html.escape(x['magnitude'])}</td><td>{x['informative_reads']}</td></tr>" for x in evidence)
    document = f"""<!doctype html><html lang=en><meta charset=utf-8><title>SMN2 judge demo</title>
<style>body{{font:16px system-ui;max-width:1050px;margin:auto;padding:2rem}} .warning{{background:#fff3cd;padding:1rem}} article{{border:1px solid #999;padding:1rem;margin:1rem 0}} table{{border-collapse:collapse}}td,th{{border:1px solid #aaa;padding:.4rem}}</style>
<h1>SMN2 raw-read evidence → prospective design</h1><p class=warning><b>{state}</b></p>
<h2>Headline</h2><p>{html.escape(rendered_headline)}</p><h2>Real observations (dPSI percentage points)</h2>
<table><tr><th>Intervention</th><th>Sample</th><th>dPSI</th><th>Informative reads</th></tr>{points}</table>
<p>Only scramble-MOE and non-targeting-guide rows are interpreted as genuine null-intervention controls. Untreated ASO1 replicate scatter is not presented as specificity evidence. Missing rows, including scramble_PMO, are excluded and retained in evidence.json.</p>
<h2>Frozen prospective SMN2 designs</h2><p class=warning>{html.escape(scoring_notice)}</p>{cards}<h2>Claim boundary</h2><p>Not a powered general benchmark, cross-locus validation, clinical claim, validated molecule, or new mechanism.</p></html>"""
    provenance = {"schema_version":"judge-demo-provenance/1", "build_timestamp":timestamp,
        "panel_sha256":sha256(panel_data), "panel_manifest_sha256":panel["panel_manifest_sha256"],
        "evidence_sha256":sha256(json.dumps(payload,sort_keys=True,separators=(",",":")).encode()),
        "renderer_sha256":sha256(Path(__file__).read_bytes())}
    out.parent.mkdir(parents=True, exist_ok=True)
    temp = Path(tempfile.mkdtemp(prefix=f".{out.name}.", dir=out.parent))
    try:
        (temp/"index.html").write_text(document, encoding="utf-8")
        (temp/"evidence.json").write_text(json.dumps(payload,indent=2,sort_keys=True)+"\n")
        (temp/"provenance.json").write_text(json.dumps(provenance,indent=2,sort_keys=True)+"\n")
        (temp/"prospective-design-panel.json").write_bytes(panel_data)
        (temp/"README.txt").write_text(f"{state}\n{rendered_headline}\n")
        if out.exists(): raise FileExistsError(f"refusing to overwrite {out}")
        os.replace(temp, out)
    except Exception:
        shutil.rmtree(temp, ignore_errors=True); raise
    return out


def main():
    p=argparse.ArgumentParser(); p.add_argument("--run-dir",type=Path,required=True); p.add_argument("--panel",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True); p.add_argument("--build-timestamp",required=True); p.add_argument("--allow-incomplete-snapshot",action="store_true")
    a=p.parse_args(); print(build(a.run_dir,a.panel,a.out,a.build_timestamp,a.allow_incomplete_snapshot))


if __name__ == "__main__": main()
