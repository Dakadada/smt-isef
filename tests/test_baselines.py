import csv
import hashlib
import tempfile
import unittest
from pathlib import Path

from baselines.backends import FrozenBackendContract, StubAccessibilityBackend, StubCisBackend, VerifiedCisBackend
from baselines.evaluation import BenchmarkHarness, EvaluationConfig, _cluster_uncertainty, load_harmonized
from baselines.features import (FEATURE_SETS, FeatureBuilder, InterventionInput, SpliceSite,
                                hybridization_dg, motif_features, occlusion_importance)
from baselines.models import CapacityMatchedLinear


def label(i, observed=0.0, family=None):
    family = family or f"f{i}"
    return {"family":family,"study":f"study{i}","locus":f"locus{i}","sample":(f"sample{i}",),
            "intervention":f"int{i}","chemistry":f"chem{i}","replicate_clusters":(f"rep{i}",),
            "event":"event","observed":observed,"magnitude_eligible":False}


class BaselineHarnessTest(unittest.TestCase):
    def item(self, i=0):
        seq="ACGTGAAGAAACGTGGAGGATTTGGGACGT"; start,end=5,11
        target=seq[start:end]; oligo=target.translate(str.maketrans("ACGT","TGCA"))[::-1]
        return InterventionInput(f"f{i}",f"study{i}",f"sample{i}",f"int{i}","event",f"locus{i}",seq,start,end,
                                 (SpliceSite(4,"donor",1),SpliceSite(16,"acceptor",-1)),oligo,"moe",dose=i+1)

    def test_strong_features_and_missing_oligo_fail_closed(self):
        builder=FeatureBuilder([StubCisBackend("spliceai"),StubCisBackend("pangolin")],StubAccessibilityBackend(),2)
        features,metadata=builder.build(self.item())
        self.assertTrue(all(any(k==s or k.startswith(s) or k.endswith(s) for k in features for s in selectors) for selectors in FEATURE_SETS.values()))
        self.assertTrue(any(k.startswith("motif_ESE_") for k in features)); self.assertTrue(any("splice_distance_donor" in k for k in features))
        self.assertTrue(all(v["contrast_kind"].endswith("not-native-treatment-output") for v in metadata.values()))
        with self.assertRaisesRegex(ValueError,"required"): hybridization_dg("AAAA","","moe")
        self.assertNotEqual(hybridization_dg("GAAGAA","TTCTTC","moe",.15),hybridization_dg("GAAGAA","TTCTTC","dna",1.0))

    def test_motif_scan_is_overlapping(self):
        values=motif_features("TAGGGTAGGGT")
        self.assertLess(values["motif_ESS_hnRNPA1_density"],0)
        self.assertIn("motif_ESE_SF2ASF_max",values)

    def test_identity_components_block_every_leakage_dimension(self):
        for dimension in ("study","locus","intervention","chemistry","sample"):
            labels={}; features={}
            for i in range(5):
                row=label(i,float(i*5-10)); labels[row["intervention"]]=row; features[row["intervention"]]={"affinity_dg":float(i)}
            keys=list(labels); value=labels[keys[0]][dimension]
            labels[keys[1]][dimension]=value
            result=BenchmarkHarness(EvaluationConfig(max_effective_df=3,ridge_grid=(10,100))).run(labels,features,physics_effective_degrees_of_freedom=2)
            self.assertEqual(result["affinity_dg"]["n_identity_components"],4)
            self.assertTrue(all(f["identity_disjoint"] for f in result["affinity_dg"]["folds"]))

    def test_missing_identity_fails_closed(self):
        rows={f"int{i}":label(i) for i in range(4)}; rows["int0"]["study"]=""
        with self.assertRaisesRegex(ValueError,"identity study"):
            BenchmarkHarness().run(rows,{k:{"affinity_dg":1} for k in rows},physics_effective_degrees_of_freedom=2)

    def test_ridge_edf_matches_analytic_one_predictor_result(self):
        rows=[{"x":-1.0},{"x":1.0}]; model=CapacityMatchedLinear(3,2).fit(rows,[-1,1])
        self.assertAlmostEqual(model.effective_degrees_of_freedom,1+2/(2+2),places=10)

    def test_channel_mapping_prevents_donor_acceptor_cancellation(self):
        class Opposing:
            name="opposing"
            def score_channels(self,seq):
                trigger=1.0 if seq[0]=="C" else 0.0
                return {"donor":[trigger]*len(seq),"acceptor":[-trigger]*len(seq)}
        result=occlusion_importance("AAAA",0,1,(SpliceSite(1,"donor",1),SpliceSite(2,"acceptor",-1)),Opposing())
        self.assertGreater(result["site_occlusion"][0],0); self.assertGreater(result["site_occlusion"][1],0)
        self.assertAlmostEqual(result["occlusion_importance"],sum(result["site_occlusion"]))

    def test_verified_backend_contract_and_golden_placeholder(self):
        digest="a"*64; golden=(("official-case-1",(.1,.2)),)
        contract=FrozenBackendContract("spliceai","1.3.1",digest,10000,"N","0-based-sequence-position","transcript-5to3","GRCh38-GENCODE-v44",("donor","acceptor"),golden)
        scorer=lambda seq,**kwargs:{"donor":[.1]*len(seq),"acceptor":[.2]*len(seq)}
        with self.assertRaisesRegex(ValueError,"golden"): VerifiedCisBackend(contract,scorer,loaded_weight_sha256=digest,reference_convention="GRCh38-GENCODE-v44",golden_observations=())
        backend=VerifiedCisBackend(contract,scorer,loaded_weight_sha256=digest,reference_convention="GRCh38-GENCODE-v44",golden_observations=golden)
        self.assertEqual(set(backend.score_channels("AC")),{"donor","acceptor"})

    def test_non_lowerable_gates(self):
        with self.assertRaisesRegex(ValueError,"cannot be lower"): EvaluationConfig(min_ci_families=7)
        labels={f"int{i}":label(i,10.0) for i in range(7)}; features={k:{"affinity_dg":float(i)} for i,k in enumerate(labels)}
        result=BenchmarkHarness(EvaluationConfig(max_effective_df=3,ridge_grid=(10,100))).run(labels,features,physics_effective_degrees_of_freedom=2)
        self.assertIsNone(result["affinity_dg"]["cluster_uncertainty"]); self.assertIsNone(result["affinity_dg"]["mae"])

    def test_contrast_provenance_and_theta_schema(self):
        fields=["cluster_family","cluster_study","cluster_locus","cluster_biological_sample","cluster_intervention","cluster_chemistry","cluster_replicate","event","primary_event","condition","magnitude","magnitude_eligible","coverage_pass","contrast_type","control_cluster","matched_control_count","sign_convention","units","contrast_schema","theta_evidence_type","theta_input_hash","theta_protocol_hash"]
        row=dict.fromkeys(fields,""); row.update(cluster_family="f1",cluster_study="st1",cluster_locus="l1",cluster_biological_sample="s1",cluster_intervention="i1",cluster_chemistry="moe",cluster_replicate="r1",event="e1",primary_event="true",condition="treated",magnitude="8",magnitude_eligible="false",coverage_pass="true",contrast_type="treated-minus-control",control_cluster="c1",matched_control_count="1",sign_convention="treated-minus-control",units="delta-PSI-points",contrast_schema="harmonized-treated-minus-matched-control/v1",theta_evidence_type="none")
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/"labels.tsv"
            def write():
                with path.open("w",newline="") as h: w=csv.DictWriter(h,fields,delimiter="\t"); w.writeheader(); w.writerow(row)
            write(); self.assertEqual(load_harmonized(path)["i1"]["observed"],8)
            row["matched_control_count"]="2"; write()
            with self.assertRaisesRegex(ValueError,"contrast provenance"): load_harmonized(path)
            row["matched_control_count"]="1"; row["magnitude_eligible"]="true"; write()
            with self.assertRaisesRegex(ValueError,"theta-evidence"): load_harmonized(path)

    def test_bootstrap_reports_real_support_and_small_cluster_sensitivity(self):
        rows=[]
        for family in range(8):
            for intervention in range(3): rows.append({"family":f"f{family}","intervention":f"i{family}-{intervention}","observed":10.0,"predicted":10.0})
        report=_cluster_uncertainty(rows,5,1,draws=200)
        self.assertIn("bootstrap_unique_family_support",report); self.assertIn("bootstrap_effective_clusters",report)
        self.assertIn("leave_one_family_out_sensitivity",report); self.assertLessEqual(report["bootstrap_unique_family_support"]["max"],8)

    def test_signal_recovery_proves_each_baseline_is_not_strawman(self):
        signals={"affinity_dg":"affinity_dg","accessibility":"accessibility","motif_model":"motif_ESE_SF2ASF_mean","splice_site_distance":"splice_distance_donor_include","kmer":"kmer_AAA","cis_occlusion":"spliceai_donor_4_occlusion","nonstructural_stack":"dose"}
        for baseline,key in signals.items():
            labels={}; features={}
            for i in range(12):
                x=(i-5.5)/5.5; row=label(i,20*x); labels[row["intervention"]]=row
                features[row["intervention"]]={key:x}
            result=BenchmarkHarness(EvaluationConfig(max_effective_df=3,ridge_grid=(.1,1,10,100))).run(labels,features,physics_effective_degrees_of_freedom=2)
            self.assertGreaterEqual(result[baseline]["directional_concordance"],.75,baseline)
            self.assertTrue(all(f["baseline_effective_degrees_of_freedom"]<=3 for f in result[baseline]["folds"]))


if __name__=="__main__": unittest.main()
