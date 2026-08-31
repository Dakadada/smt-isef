import json
import tempfile
import unittest
from datetime import datetime, timezone
from pathlib import Path

from design_panel import render_json, run_panel
from rawread_reanalysis.common import write_tsv
from scripts.judge_demo import BANNER, INTERVENTIONS, build


class JudgeDemoTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.root = Path(self.temp.name)
        self.run = self.root/"run"; (self.run/"work/harmonized").mkdir(parents=True); (self.run/"work/psi").mkdir()
        self.panel = self.root/"panel.json"
        self.panel.write_text(render_json(run_panel(clock=lambda: datetime(2026,8,30,tzinfo=timezone.utc))))
        base={"event":"SMN2_exon7_inclusion","coverage_pass":"true","psi":"80","control_mean_psi":"60","pair_set_id":"p","approximation_label":"targeted"}
        specs=[("F1","S1","ASO1","treated","increase","25.4"),("F2","S2","F18MOE","treated","increase","14.3"),
               ("F3","S3","dCas13d_D[333,357]","treated","decrease","-21.8"),("F2","S4","scramble_MOE","control","no_change","0"),
               ("F3","S5","non_targeting_guide","control","no_change","0.58"),("F4","S6","F20PMO","treated","increase","27.87"),
               ("F4","S7","F18OMe","treated","increase","24.08"),("F4","S8","scramble_PMO","control","MISSING","")]
        self.rows=[{**base,"family":f,"sample":s,"intervention":i,"condition":c,"direction":d,"magnitude":m} for f,s,i,c,d,m in specs]
        for row in self.rows:
            if row["direction"] == "no_change": row["psi"] = row["control_mean_psi"]
        qualifying = {row["intervention"]: row for row in self.rows if row["intervention"] in INTERVENTIONS}
        self.rows.extend({**qualifying[intervention], "sample":f"S{n}"}
                         for n, intervention in enumerate(sorted(INTERVENTIONS), start=10))
        write_tsv(self.run/"work/harmonized/tidy_splicing.tsv",self.rows,self.rows[0].keys())
        effects=[{"family":r["family"],"sample":r["sample"],"intervention":r["intervention"],"run_accession":f"R{n}","informative_reads":"30"} for n,r in enumerate(self.rows)]
        write_tsv(self.run/"work/psi/sample_effects.tsv",effects,effects[0].keys())

    def tearDown(self): self.temp.cleanup()

    def test_incomplete_snapshot_is_labeled_and_every_design_is_bannered(self):
        out=build(self.run,self.panel,self.root/"demo","2026-08-30T00:00:00Z",True)
        page=(out/"index.html").read_text(); evidence=json.loads((out/"evidence.json").read_text())
        self.assertIn("INCOMPLETE STAGE-A SNAPSHOT",page)
        self.assertIn("scramble_PMO",json.dumps(evidence["excluded_missing_rows"]))
        self.assertGreaterEqual(page.count(BANNER),1)
        self.assertNotIn("untreated",page)
        self.assertIn("Deterministic placeholder only",page)
        self.assertIn("replicate observations within each arm",page)

    def test_incomplete_run_fails_without_explicit_snapshot_mode(self):
        with self.assertRaisesRegex(RuntimeError,"Stage A is incomplete"):
            build(self.run,self.panel,self.root/"demo","2026-08-30T00:00:00Z",False)

    def test_missing_genuine_null_arm_fails_closed(self):
        rows=[r for r in self.rows if r["intervention"]!="scramble_MOE"]
        write_tsv(self.run/"work/harmonized/tidy_splicing.tsv",rows,rows[0].keys())
        with self.assertRaisesRegex(ValueError,"required evidence arms absent"):
            build(self.run,self.panel,self.root/"demo","2026-08-30T00:00:00Z",True)

    def test_arm_with_fewer_than_two_replicates_fails_closed(self):
        f18moe_rows=[row for row in self.rows if row["intervention"] == "F18MOE"]
        self.rows.remove(f18moe_rows[-1])
        write_tsv(self.run/"work/harmonized/tidy_splicing.tsv",self.rows,self.rows[0].keys())

        with self.assertRaisesRegex(ValueError,"arm has fewer than 2 replicate observations: F18MOE"):
            build(self.run,self.panel,self.root/"demo","2026-08-30T00:00:00Z",True)

    def test_completion_marker_allows_out_of_panel_missing_rows(self):
        (self.run/"STAGE_A_COMPLETE").touch()

        out=build(self.run,self.panel,self.root/"demo","2026-08-30T00:00:00Z",False)

        evidence=json.loads((out/"evidence.json").read_text())
        self.assertTrue(evidence["stage_a_complete"])
        self.assertIn("scramble_PMO",json.dumps(evidence["excluded_missing_rows"]))

    def test_completion_marker_with_required_missing_row_fails_closed(self):
        (self.run/"STAGE_A_COMPLETE").touch()
        aso1=next(row for row in self.rows if row["intervention"] == "ASO1")
        self.rows.append({**aso1, "sample":"S9", "direction":"MISSING", "magnitude":"", "coverage_pass":"false"})
        write_tsv(self.run/"work/harmonized/tidy_splicing.tsv",self.rows,self.rows[0].keys())

        with self.assertRaisesRegex(ValueError,"completion marker conflicts"):
            build(self.run,self.panel,self.root/"demo","2026-08-30T00:00:00Z",False)

    def test_nonfinite_magnitude_fails_closed(self):
        self.rows[0]["magnitude"]="nan"
        write_tsv(self.run/"work/harmonized/tidy_splicing.tsv",self.rows,self.rows[0].keys())
        with self.assertRaisesRegex(ValueError,"magnitude must be finite"):
            build(self.run,self.panel,self.root/"demo","2026-08-30T00:00:00Z",True)

    def test_duplicate_sample_effect_identity_fails_closed(self):
        effects_path=self.run/"work/psi/sample_effects.tsv"
        with effects_path.open() as handle:
            effects=list(__import__("csv").DictReader(handle,delimiter="\t"))
        effects.append(dict(effects[0]))
        write_tsv(effects_path,effects,effects[0].keys())
        with self.assertRaisesRegex(ValueError,"duplicate sample-effect identity"):
            build(self.run,self.panel,self.root/"demo","2026-08-30T00:00:00Z",True)

    def test_replicate_null_noise_uses_arm_level_psi_effect(self):
        null = next(r for r in self.rows if r["intervention"] == "scramble_MOE")
        null.update(psi="61", magnitude="4")
        replicate = {**null, "sample":"S7", "psi":"59", "magnitude":"-4"}
        self.rows.append(replicate)
        write_tsv(self.run/"work/harmonized/tidy_splicing.tsv",self.rows,self.rows[0].keys())
        effects_path=self.run/"work/psi/sample_effects.tsv"
        with effects_path.open() as handle:
            effects=list(__import__("csv").DictReader(handle,delimiter="\t"))
        effects.append({"family":replicate["family"],"sample":replicate["sample"],
                        "intervention":replicate["intervention"],"run_accession":"R7",
                        "informative_reads":"30"})
        write_tsv(effects_path,effects,effects[0].keys())

        out=build(self.run,self.panel,self.root/"demo","2026-08-30T00:00:00Z",True)

        self.assertTrue((out/"evidence.json").exists())

    def test_genuine_arm_level_null_violation_fails_closed(self):
        for row in self.rows:
            if row["intervention"] == "scramble_MOE":
                row.update(psi="63", magnitude="0")
        write_tsv(self.run/"work/harmonized/tidy_splicing.tsv",self.rows,self.rows[0].keys())

        with self.assertRaisesRegex(ValueError,"genuine null arm exceeds 2 dPSI"):
            build(self.run,self.panel,self.root/"demo","2026-08-30T00:00:00Z",True)


if __name__ == "__main__": unittest.main()
