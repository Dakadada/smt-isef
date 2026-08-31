import csv
import json
import subprocess
import sys
import tempfile
import unittest
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

HEADER = ["study_accession", "run_accession", "sample_accession", "secondary_sample_accession",
          "sample_alias", "sample_title", "experiment_alias", "experiment_title", "fastq_ftp",
          "fastq_md5", "fastq_bytes", "library_layout", "library_strategy"]
RUNNABLE = "GSE307431_PRJNA1321660,GSE167762_PRJNA705187,PRJNA1052769_dCas13d"
DCAS_RUNS = [
    ("SRR27213209", "endogenous_RNAseq_D333-357_Rep1"),
    ("SRR27213208", "endogenous_RNAseq_D333-357_Rep2"),
    ("SRR27213210", "endogenous_RNAseq_non-targeting_Rep2"),
    ("SRR27213211", "endogenous_RNAseq_non-targeting_Rep1"),
]


def write_project(metadata_dir, project, rows):
    with (metadata_dir / f"{project}.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=HEADER, delimiter="\t")
        writer.writeheader()
        for run_accession, alias, title, strategy in rows:
            writer.writerow({"study_accession": project, "run_accession": run_accession,
                             "sample_accession": alias, "sample_alias": alias, "sample_title": title,
                             "fastq_ftp": f"ftp/x/{run_accession}.fastq.gz", "fastq_md5": "0" * 32,
                             "fastq_bytes": "100", "library_layout": "SINGLE", "library_strategy": strategy})


def make_metadata(base, ambiguous_geo=False):
    metadata_dir = base / "ena"
    metadata_dir.mkdir()
    gsm1 = [f"GSM{x}" for x in range(9224134, 9224154)]
    if ambiguous_geo:
        gsm1[-1] = "UNMAPPED_SAMPLE"
    write_project(metadata_dir, "PRJNA1321660",
                  [(f"SMN_A_{i}", alias, alias, "RNA-Seq") for i, alias in enumerate(gsm1)])
    aliases = [f"GSM{x}" for x in range(6063682, 6063694)] + [f"AUX{x}" for x in range(16)]
    write_project(metadata_dir, "PRJNA705187",
                  [(f"SMN_B_{i}", alias, alias, "RNA-Seq" if alias.startswith("GSM") else "ChIP-Seq")
                   for i, alias in enumerate(aliases)])
    write_project(metadata_dir, "PRJNA1052769",
                  [(run, title, title, "RNA-Seq") for run, title in DCAS_RUNS] +
                  [(f"SCREEN{i}", f"SCREEN{i}", "pooled screen", "AMPLICON") for i in range(12)])
    write_project(metadata_dir, "PRJEB13980",
                  [(f"DMD_RUN{i}", f"DMD{i}", "U7E53_OFFTARGET_ANALYSIS", "RNA-Seq") for i in range(12)])
    return metadata_dir


def run_builder(base, metadata_dir, *extra):
    out = base / "out"
    proc = subprocess.run(
        [sys.executable, "scripts/real_run/build_manifest.py", "--out", str(out),
         "--metadata-dir", str(metadata_dir), *extra],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return proc, out


class RealRunManifestTest(unittest.TestCase):
    def test_subset_includes_only_requested_families_and_defers_dmd(self):
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            proc, out = run_builder(base, make_metadata(base), "--families", RUNNABLE)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            cfg = json.loads((out / "stageA.targeted.json").read_text())
            self.assertEqual(set(cfg["studies"]), set(RUNNABLE.split(",")))
            self.assertEqual({event["family"] for event in cfg["events"]}, set(cfg["studies"]))
            self.assertFalse((out / "PRJEB13980_U7.frozen.tsv").exists())
            with (out / "all_runs.frozen.tsv").open(newline="") as handle:
                families = {row["family"] for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(families, set(cfg["studies"]))

    def test_smn2_subset_has_hashed_exact_coverage_event_table(self):
        from rawread_reanalysis.common import read_tsv, validate_config

        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            proc, out = run_builder(base, make_metadata(base), "--families", RUNNABLE)
            self.assertEqual(proc.returncode, 0, proc.stderr)
            config_path = out / "stageA.targeted.json"
            cfg = validate_config(config_path)
            event_path = out / "event_junctions.GRCh38.run.tsv"
            self.assertEqual(cfg["alignment"]["tier"], "targeted")
            self.assertEqual({row["event"] for row in read_tsv(event_path)}, {"SMN2_exon7_inclusion"})
            self.assertEqual({event["id"] for event in cfg["events"]}, {"SMN2_exon7_inclusion"})
            self.assertNotIn(b"DMD_exon53_inclusion", event_path.read_bytes())
            self.assertEqual(cfg["asset_sha256"]["event_junctions"], hashlib.sha256(event_path.read_bytes()).hexdigest())
            master_lines = Path("config/event_junctions.GRCh38.tsv").read_bytes().splitlines(keepends=True)
            expected = master_lines[0] + next(line for line in master_lines[1:] if line.startswith(b"SMN2_exon7_inclusion\t"))
            self.assertEqual(event_path.read_bytes(), expected)

    def test_embedded_dcas13d_mapping_has_two_treated_and_two_control(self):
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            proc, out = run_builder(base, make_metadata(base), "--families", "PRJNA1052769_dCas13d")
            self.assertEqual(proc.returncode, 0, proc.stderr)
            cfg = json.loads((out / "stageA.targeted.json").read_text())
            sample_map = cfg["studies"]["PRJNA1052769_dCas13d"]["sample_map"]
            self.assertEqual(sorted(sample_map), sorted(run for run, _ in DCAS_RUNS))
            conditions = [row["condition"] for row in sample_map.values()]
            self.assertEqual(conditions.count("treated"), 2)
            self.assertEqual(conditions.count("control"), 2)
            self.assertTrue(all(row["contrast"] == "targeting_vs_nontargeting" for row in sample_map.values()))
            self.assertTrue(all(row["pair_set_id"] == "targeting_vs_nontargeting" for row in sample_map.values()))
            expected_covariates = {"effector": "dCas13d", "system": "HEK293T"}
            self.assertTrue(all(row["batch_covariates"] == expected_covariates for row in sample_map.values()))

    def test_ambiguous_arm_in_selected_family_still_fails_closed(self):
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            proc, _ = run_builder(base, make_metadata(base, ambiguous_geo=True),
                                  "--families", "GSE307431_PRJNA1321660")
            self.assertEqual(proc.returncode, 2)
            self.assertIn("has no verified arm mapping", proc.stderr)

    def test_default_all_families_behavior_still_requires_dmd_mapping(self):
        with tempfile.TemporaryDirectory() as td:
            base = Path(td)
            proc, _ = run_builder(base, make_metadata(base))
            self.assertEqual(proc.returncode, 2)
            self.assertIn("PRJEB13980_U7", proc.stderr)


if __name__ == "__main__":
    unittest.main()
