from __future__ import annotations

import csv
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"


AE_FINDINGS = [
    {
        "ae_issue": "cohesion_and_friction_degradation_parameters",
        "editorial_risk": "The proposed c and phi degradation parameters can be read as unsupported calibration.",
        "current_local_evidence": "benchmark_parameter_audit_strengthened.csv flags cb0 as bounded and Delta tan phi as not directly isolated; v4_component_dataset_inventory.csv lists component data only.",
        "support_level": "PARTIAL_COMPONENT_SUPPORT_NOT_DIRECT_CALIBRATION",
        "claim_allowed": "bounded pre-design parameterization informed by component evidence",
        "claim_blocked": "calibrated footing-scale cohesion/friction degradation law",
        "minimum_fix": "Use a conservative envelope, report Delta tan phi = 0 sensitivity, and require project-specific triaxial/direct-shear or plate-load updating before final design.",
    },
    {
        "ae_issue": "treated_zone_contribution_function",
        "editorial_risk": "The H/B, kH, kz and treated-zone weighting can look speculative without a treated-depth dataset.",
        "current_local_evidence": "v4_closed_form_capped_fem_comparison.csv and virtual footing audits support a cap/sensitivity rule, not experimental H/B calibration.",
        "support_level": "VIRTUAL_MECHANISTIC_SUPPORT_NOT_EXPERIMENTAL_CALIBRATION",
        "claim_allowed": "screening cap and sensitivity rule for pre-design",
        "claim_blocked": "experimentally calibrated treated-zone contribution function",
        "minimum_fix": "Keep H/B caps mandatory and require FEM or plate-load calibration for thin treated crusts and final design.",
    },
    {
        "ae_issue": "localization_fragility_factor",
        "editorial_risk": "The chi_s/localization factor may mask brittle post-peak behavior without direct evidence.",
        "current_local_evidence": "benchmark_parameter_audit_strengthened.csv classifies chi_s, rs as sensitivity-envelope parameters with no calibration.",
        "support_level": "SENSITIVITY_ONLY_BLOCKED_FOR_CALIBRATION",
        "claim_allowed": "conservative fragility sensitivity and falsification flag",
        "claim_blocked": "identified localization-fragility parameter",
        "minimum_fix": "State chi_s is not identified; use it only for sensitivity unless post-peak/localization data or calibrated elastoplastic FEM are added.",
    },
    {
        "ae_issue": "direct_experimental_validation",
        "editorial_risk": "The paper may still lack direct experimental data showing degraded MICP shallow-footing performance over time.",
        "current_local_evidence": "v4_validation_scale_v0_v6.csv states V4 coupled indirect validation, not V5 direct degraded-footing validation.",
        "support_level": "V4_ONLY_V5_BLOCKED",
        "claim_allowed": "component-constrained pre-design screening workflow",
        "claim_blocked": "direct degraded-footing validation and universal design model",
        "minimum_fix": "Resubmit only as falsifiable pre-design method or add direct degraded MICP plate/footing tests.",
    },
]


NOVELTY_ROWS = [
    {
        "prior_family": "MICP biogeochemical/multiphysics FEM design models",
        "representative_sources": "Fauriel and Laloui 2012; Wang and Nackenhorst 2022; Bosch et al. 2024",
        "what_prior_already_covers": "MICP process, treatment evolution, cementation state and/or shallow-foundation strengthening mechanics.",
        "article_127_differentiator": "service-life reliability pre-design with degradation discount, beta(t), pf(t), inspection or retreatment triggers.",
        "evidence_in_package": "time-dependent reliability scripts, Monte Carlo/FOSM checks, V4 validation architecture",
        "novelty_verdict": "PARTIAL_PASS_IF_FRAMED_AS_RELIABILITY_PREDESIGN",
        "claim_adjustment": "Do not claim a new MICP constitutive law; claim a reliability pre-design decision layer around existing MICP evidence.",
    },
    {
        "prior_family": "MICP plate-load and footing-scale improvement tests",
        "representative_sources": "Kulkarni et al. 2021 and related footing-scale MICP tests",
        "what_prior_already_covers": "Initial MICP bearing-capacity and settlement improvement under plate/footing loading.",
        "article_127_differentiator": "discounting initial improvement over service life and identifying when permanent-credit designs become unsafe.",
        "evidence_in_package": "Kulkarni-derived trend anchors and service/load-settlement calibration files",
        "novelty_verdict": "PARTIAL_PASS_ONLY_WITH_CLEAR_INITIAL_VS_DEGRADED_DECISION_GAIN",
        "claim_adjustment": "Do not claim new footing test evidence; claim decision impact of not treating initial MICP credit as permanent.",
    },
    {
        "prior_family": "MICP durability and cyclic/environmental degradation studies",
        "representative_sources": "Sharma et al. 2021; Ahenkorah et al. 2023; durability/aging literature",
        "what_prior_already_covers": "Strength, stiffness, CaCO3 retention and mass-loss degradation under wetting-drying, freezing-thawing, aging or temperature.",
        "article_127_differentiator": "translating durability retention uncertainty into foundation reliability and pre-design decisions.",
        "evidence_in_package": "degradation calibration files plus Ahenkorah 2023 open-access 40-UCS-test source as external support",
        "novelty_verdict": "PASS_AS_TRANSLATION_TO_FOUNDATION_RELIABILITY_NOT_AS_DURABILITY_LAW",
        "claim_adjustment": "Do not claim universal annual degradation law; claim a data-updatable retention-to-reliability workflow.",
    },
    {
        "prior_family": "Reliability-based shallow-foundation and deterioration models",
        "representative_sources": "general geotechnical RBD and deterioration/service-life reliability literature",
        "what_prior_already_covers": "Uncertain loads and soil parameters, reliability indices, service-life and deterioration concepts.",
        "article_127_differentiator": "MICP-specific internal cementation state and measurable-update route for biocemented shallow foundations.",
        "evidence_in_package": "Bayesian updating and parameter-inference mapping supplements",
        "novelty_verdict": "PARTIAL_PASS_IF_MICP_MEASURABLE_UPDATE_ROUTE_IS_PRIMARY",
        "claim_adjustment": "Position as MICP-specific pre-design/update workflow, not as general reliability theory.",
    },
    {
        "prior_family": "Direct degraded-footing validation studies",
        "representative_sources": "No local accessible direct degraded MICP footing validation located",
        "what_prior_already_covers": "Potentially decisive V5 evidence if located or produced.",
        "article_127_differentiator": "none until direct degraded-footing data are obtained",
        "evidence_in_package": "explicit V4/V5 boundary",
        "novelty_verdict": "FAIL_FOR_DIRECT_VALIDATED_MODEL_CLAIM",
        "claim_adjustment": "Block any claim of direct degraded-footing validation.",
    },
]


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    ae_path = DATA / "ae_rejection_evidence_gate.csv"
    novelty_path = DATA / "q1_novelty_positioning_gate.csv"
    write_csv(ae_path, AE_FINDINGS)
    write_csv(novelty_path, NOVELTY_ROWS)

    summary = {
        "status": "CONDITIONAL_RESUBMISSION_ONLY_NOT_Q1_HIGH_CLOSED",
        "novelty_status": "NOVELTY_DEFENDABLE_ONLY_AS_MICP_RELIABILITY_PREDESIGN_AND_UPDATE_WORKFLOW",
        "ae_evidence_status": "AE_CORE_DATA_GAP_PARTIALLY_REDUCED_BUT_NOT_CLOSED",
        "q1_high_status": "NOT_Q1_HIGH_CLOSED",
        "strongest_novelty_claim": (
            "a falsifiable reliability pre-design and updating workflow that translates MICP component durability "
            "and initial footing evidence into beta(t), pf(t), service-pressure limits, inspection and retreatment decisions"
        ),
        "claims_to_remove_or_downgrade": [
            "universal calibrated degradation law for MICP-treated shallow foundations",
            "direct experimental validation of degraded MICP footing performance",
            "identified localization fragility factor",
            "experimentally calibrated treated-zone contribution function",
            "fully supported cohesion/friction degradation split",
        ],
        "minimum_editorial_reframe": (
            "resubmit as a V4 component-constrained pre-design method, not as a V5 experimentally validated design model"
        ),
        "external_support_added_to_gate": {
            "source": "Ahenkorah et al. 2023, Journal of Rock Mechanics and Geotechnical Engineering",
            "doi": "10.1016/j.jrmge.2022.08.007",
            "basis": "open-access 40-test UCS durability study for MICP/EICP treated sands; supports degradation concern but not footing-scale calibration",
        },
        "ae_issue_count": len(AE_FINDINGS),
        "ae_issues_closed_for_q1_high": 0,
        "ae_issues_reduced": 4,
        "novelty_rows": len(NOVELTY_ROWS),
        "novelty_pass_or_partial_rows": sum("PASS" in row["novelty_verdict"] for row in NOVELTY_ROWS),
        "required_before_q1_high": [
            "direct degraded MICP footing or plate-load data, or equivalent project-specific calibration",
            "triaxial/direct-shear data to separate cb0 from Delta tan phi or a conservative Delta tan phi = 0 design path",
            "treated-depth data or calibrated FEM/limit-analysis protocol for kH/kz",
            "post-peak/localization evidence or explicit chi_s sensitivity-only status",
        ],
    }
    (DATA / "ae_rejection_novelty_evidence_gate_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
