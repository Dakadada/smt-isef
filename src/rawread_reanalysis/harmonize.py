from pathlib import Path

from .common import Checkpoints, read_tsv, validate_config, write_tsv

FIELDS = ["family", "sample", "intervention", "event", "condition", "direction", "magnitude", "therapeutic_polarity_multiplier", "therapeutic_effect", "pair_set_id", "magnitude_eligible",
          "psi", "control_mean_psi", "coverage_pass", "batch_covariates", "cluster_family", "cluster_biological_sample",
          "cluster_intervention", "primary_event", "quantification_tier", "approximation_label", "confound_flags"]


def run(config_path, effects_path=None):
    cfg = validate_config(config_path); work = Path(cfg["paths"]["work_dir"])
    effects_path = Path(effects_path or work / "psi/sample_effects.tsv"); output = work / "harmonized/tidy_splicing.tsv"
    cp = Checkpoints(cfg["paths"]["checkpoints"], "harmonize", config_path,cfg)
    if cp.valid("all", [output],[effects_path]): return output
    events = {e["id"]: e for e in cfg["events"]}; rows = []
    for row in read_tsv(effects_path):
        # Magnitude is observed and retained, but eligibility explicitly prevents conditional magnitude evaluation without theta.
        rows.append({"family": row["family"], "sample": row["sample"], "intervention": row["intervention"], "event": row["event"],
            "condition": row["condition"], "direction": row["direction"], "magnitude": row["dpsi"], "therapeutic_polarity_multiplier":row["therapeutic_polarity_multiplier"], "therapeutic_effect":row["therapeutic_effect"], "pair_set_id":row["pair_set_id"],
            "magnitude_eligible": row["magnitude_eligible"], "psi": row["psi"], "control_mean_psi": row["control_mean_psi"],
            "coverage_pass": row["coverage_pass"], "batch_covariates": row["batch_covariates"], "cluster_family": row["family"],
            "cluster_biological_sample": f"{row['family']}::{row['sample']}", "cluster_intervention": f"{row['family']}::{row['intervention']}",
            "primary_event": str(events[row["event"]].get("primary", True)).lower(), "quantification_tier": row["quantification_tier"],
            "approximation_label": row["approximation_label"], "confound_flags": row["confound_flags"]})
    write_tsv(output, rows, FIELDS); cp.complete("all", [output],inputs=[effects_path], hierarchy="family>biological_sample>intervention>event", rescaled=False)
    return output
