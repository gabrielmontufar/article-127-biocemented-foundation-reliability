from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


def read_csv(name: str) -> pd.DataFrame:
    path = DATA / name
    if not path.exists():
        raise FileNotFoundError(f"Missing required data file: {path}")
    return pd.read_csv(path)


def build_parameter_evidence_matrix() -> pd.DataFrame:
    rows = [
        {
            "ae_issue": "cohesion_and_friction_degradation_parameters",
            "model_item": "cb0_prime, Delta_tan_phi, eta(t)",
            "published_data_anchor": "Kulkarni et al. 2021 plate-load/UCS ranges; de Rezende et al. 2022 triaxial bounds; Gao et al. 2019 biocemented sand mechanics",
            "calibration_or_constraint": "Bound initial strength increment with PLT/UCS/triaxial component evidence; do not fit cohesion and friction from one composite gain simultaneously.",
            "manuscript_claim_allowed": "component-constrained pre-design parameterization",
            "manuscript_claim_blocked": "universal calibrated footing-scale cohesion/friction degradation law",
            "required_text_boundary": "State that final design requires site-specific triaxial, direct-shear, plate-load, or monitoring update.",
        },
        {
            "ae_issue": "degradation_rate_and_residual_cementation",
            "model_item": "lambda_e, eta_r",
            "published_data_anchor": "Ahenkorah et al. 2023 durability tests; Sharma et al. 2021 retention data already digitized in supplement",
            "calibration_or_constraint": "Use exposure-clock degradation priors and hold-out retention checks; convert to service time only through an explicit environmental-history assumption.",
            "manuscript_claim_allowed": "data-updatable degradation clock and sensitivity envelope",
            "manuscript_claim_blocked": "universal annual degradation rate",
            "required_text_boundary": "Separate exposure-cycle calibration from calendar-year prediction.",
        },
        {
            "ae_issue": "treated_zone_contribution_function",
            "model_item": "chi_H, kH, w(z), kz",
            "published_data_anchor": "Kulkarni et al. 2021 PLT treated/untreated response; independent limit-analysis/OpenSeesPy/CalculiX virtual footing audits; cemented-layer D/B analogs only as bounding support",
            "calibration_or_constraint": "Use a capped monotonic H/B screening rule and report it as virtual/mechanistic support, not physical MICP depth calibration.",
            "manuscript_claim_allowed": "bounded treated-depth screening function",
            "manuscript_claim_blocked": "experimentally calibrated treated-depth contribution function",
            "required_text_boundary": "Require project PLT/FEM calibration for thin treated crusts and final design.",
        },
        {
            "ae_issue": "localization_fragility_factor",
            "model_item": "chi_s, r_s",
            "published_data_anchor": "No direct post-peak/localization data in current package",
            "calibration_or_constraint": "Use chi_s only as sensitivity or conservative cap; base-case decisions must be reported with chi_s=1 or a declared sensitivity range.",
            "manuscript_claim_allowed": "fragility sensitivity and falsification flag",
            "manuscript_claim_blocked": "identified localization-fragility parameter",
            "required_text_boundary": "Do not present chi_s as calibrated unless post-peak data or elastoplastic FEM back-analysis is added.",
        },
        {
            "ae_issue": "direct_experimental_validation",
            "model_item": "V4/V5 validation level",
            "published_data_anchor": "Published component data plus virtual degraded-footing comparisons",
            "calibration_or_constraint": "Declare V4 component-constrained indirect validation; reserve V5 for direct degraded MICP footing/plate tests.",
            "manuscript_claim_allowed": "falsifiable reliability pre-design and updating workflow",
            "manuscript_claim_blocked": "direct physical degraded-footing validation",
            "required_text_boundary": "All figures that forecast 50-year behavior must be labelled as calibrated-component scenario predictions.",
        },
    ]
    return pd.DataFrame(rows)


def build_calibrated_degradation_priors() -> pd.DataFrame:
    metrics = read_csv("dataset_level_degradation_metrics.csv")
    params = read_csv("dataset_level_degradation_fitted_parameters.csv")
    param_map = {row["parameter"]: row["value"] for _, row in params.iterrows()}
    val = metrics[metrics["Role"].str.lower() == "validation"].copy()
    return pd.DataFrame(
        [
            {
                "source": "Ahenkorah et al. 2023 / Sharma et al. 2021 durability evidence envelope",
                "calibrated_parameter": "eta_r",
                "value": param_map.get("eta_r"),
                "units": "-",
                "use_in_manuscript": "residual cementation prior/sensitivity, not universal residual",
                "holdout_metric": f"validation MAPE={float(val['MAPE'].iloc[0]):.2f}%"
                if not val.empty
                else "not available",
            },
            {
                "source": "Ahenkorah et al. 2023 / Sharma et al. 2021 durability evidence envelope",
                "calibrated_parameter": "lambda_E",
                "value": param_map.get("lambda_E"),
                "units": "per exposure cycle",
                "use_in_manuscript": "exposure-clock degradation prior; convert to years only with environmental history",
                "holdout_metric": f"validation final-retention error={float(val['Final-retention error (%)'].iloc[0]):.2f}%"
                if not val.empty
                else "not available",
            },
        ]
    )


def build_component_calibration_metrics() -> pd.DataFrame:
    load = read_csv("dataset_level_load_settlement_metrics.csv")
    degr = read_csv("dataset_level_degradation_metrics.csv")
    external = read_csv("external_quantitative_validation_round3.csv")
    rows = []
    for _, row in load.iterrows():
        rows.append(
            {
                "component": "load_settlement",
                "source": row["Source"],
                "role": row["Role"],
                "metric_1": f"NRMSE={float(row['NRMSE']):.3f}",
                "metric_2": f"MAPE={float(row['MAPE']):.2f}%",
                "interpretation": row["Interpretation"],
            }
        )
    for _, row in degr.iterrows():
        rows.append(
            {
                "component": "degradation_retention",
                "source": row["Source"],
                "role": row["Role"],
                "metric_1": f"RMSE={float(row['RMSE']):.3f}",
                "metric_2": f"MAPE={float(row['MAPE']):.2f}%",
                "interpretation": row["Interpretation"],
            }
        )
    for _, row in external.iterrows():
        rows.append(
            {
                "component": "external_range_check",
                "source": row["source"],
                "role": "range_check",
                "metric_1": row["observed_range"],
                "metric_2": row["model_range"],
                "interpretation": row["comparison"],
            }
        )
    return pd.DataFrame(rows)


def build_model_hierarchy_ablation() -> pd.DataFrame:
    model = read_csv("v4_model_m0_m8_comparison.csv")
    validation = read_csv("v4_validation_scale_v0_v6.csv")
    v4_status = validation.loc[validation["level"] == "V4", "manuscript_status"].iloc[0]
    v5_status = validation.loc[validation["level"] == "V5", "manuscript_status"].iloc[0]
    model["resubmission_interpretation"] = model["model"].map(
        {
            "M0": "untreated baseline comparator",
            "M1": "permanent-credit model used to show over-credit risk",
            "M3/M4": "degradation-aware closed-form comparator",
            "M7": "capped treated-depth comparator",
            "M8": f"preferred V4 indirect model; {v4_status}; V5 boundary: {v5_status}",
        }
    )
    return model


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    outputs = {
        "ae_parameter_evidence_matrix.csv": build_parameter_evidence_matrix(),
        "calibrated_degradation_priors.csv": build_calibrated_degradation_priors(),
        "component_calibration_metrics.csv": build_component_calibration_metrics(),
        "model_hierarchy_ablation.csv": build_model_hierarchy_ablation(),
    }
    for filename, frame in outputs.items():
        frame.to_csv(DATA / filename, index=False)

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(ROOT),
        "outputs": sorted(outputs),
        "validation_boundary": "V4 component-constrained indirect validation; V5 direct degraded-footing validation is not claimed.",
    }
    (DATA / "ae_resubmission_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
