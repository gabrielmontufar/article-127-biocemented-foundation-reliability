from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
SUPP = ROOT / "04 Supplemental data and code" / "Supplementary files"
DATA = SUPP / "data"
CODE = SUPP / "code"


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    CODE.mkdir(parents=True, exist_ok=True)

    identifiability = pd.DataFrame(
        [
            {
                "Parameter": "eta(t,z)",
                "Physical role": "Cementation state and its temporal/depth variation",
                "Recommended data": "CaCO3, Vs, UPV, resistivity, image-based cement maps, treatment-delivery records, durability-retention data",
                "Identifiability level": "High if measured independently",
                "Confounded with": "cb0', Delta tan phi, lambda_e",
                "Recommended treatment": "Constrain before fitting footing response",
            },
            {
                "Parameter": "cb0'",
                "Physical role": "Initial cohesion-like cementation increment",
                "Recommended data": "UCS, triaxial compression, direct shear, treated and untreated element tests",
                "Identifiability level": "Moderate to high with element tests",
                "Confounded with": "eta, chiH, chis",
                "Recommended treatment": "Estimate from element tests or use informative priors",
            },
            {
                "Parameter": "Delta tan phi",
                "Physical role": "Friction-angle enhancement expressed in tangent space",
                "Recommended data": "Drained triaxial or direct shear tests comparing treated and untreated sand",
                "Identifiability level": "Moderate with stress-path data",
                "Confounded with": "cb0', bearing-factor model error",
                "Recommended treatment": "Estimate from treated-minus-untreated friction response",
            },
            {
                "Parameter": "kH",
                "Physical role": "Depth scale controlling mobilized treated-zone contribution to ultimate capacity",
                "Recommended data": "Plate-load tests, centrifuge tests, limit analysis, or FEM with varied H/B",
                "Identifiability level": "Low unless multiple treated depths are available",
                "Confounded with": "chis, cb0', eta0",
                "Recommended treatment": "Calibrate only from depth series; otherwise report sensitivity envelope",
            },
            {
                "Parameter": "kE",
                "Physical role": "Depth scale controlling treated-zone contribution to serviceability stiffness",
                "Recommended data": "Load-settlement curves, secant stiffness, initial stiffness, stiffness profiles with varied H/B",
                "Identifiability level": "Low unless settlement data span multiple H/B values",
                "Confounded with": "aE, stiffness calibration, eta",
                "Recommended treatment": "Fit from settlement/stiffness data; otherwise fix and bracket",
            },
            {
                "Parameter": "chimin, rs",
                "Physical role": "Localization or brittle-response cap on cemented strength contribution",
                "Recommended data": "Post-peak triaxial/direct shear, localization measurements, centrifuge tests, plate-load failure modes, calibrated elastoplastic FEM",
                "Identifiability level": "Low without post-peak or localization-sensitive data",
                "Confounded with": "kH, cb0', eta0",
                "Recommended treatment": "Use conservative envelope unless post-peak data exist",
            },
        ]
    )

    workflow = pd.DataFrame(
        [
            {
                "Stage": "1",
                "Data source": "CaCO3, Vs, UPV, resistivity, delivery records, durability-retention tests",
                "Parameters estimated": "eta(t,z), eta0, eta_r, lambda_e",
                "Parameters fixed": "cb0', Delta tan phi, kH, kE, chimin, rs",
                "Output checked": "Cementation retention and spatial/depth trends",
                "Risk if skipped": "Footing data may absorb cementation uncertainty into strength or modifier terms",
            },
            {
                "Stage": "2",
                "Data source": "UCS, drained triaxial, direct shear",
                "Parameters estimated": "cb0', Delta tan phi",
                "Parameters fixed": "kH, kE, chimin, rs",
                "Output checked": "Treated-minus-untreated strength increments",
                "Risk if skipped": "Cohesion and friction increments become non-unique in footing calibration",
            },
            {
                "Stage": "3",
                "Data source": "Plate-load, centrifuge, FEM, or limit analysis with varied H/B",
                "Parameters estimated": "kH or chiH(H/B)",
                "Parameters fixed": "eta(t,z), cb0', Delta tan phi, kE, chimin, rs",
                "Output checked": "Ultimate capacity trend with treated depth",
                "Risk if skipped": "Depth benefit may be overassigned to strength or localization terms",
            },
            {
                "Stage": "4",
                "Data source": "Load-settlement curves or stiffness profiles with varied H/B",
                "Parameters estimated": "kE or chiE(H/B)",
                "Parameters fixed": "eta(t,z), cb0', Delta tan phi, kH, chimin, rs",
                "Output checked": "Initial/secant stiffness and settlement index",
                "Risk if skipped": "Serviceability response may be fitted by ultimate-capacity parameters",
            },
            {
                "Stage": "5",
                "Data source": "Post-peak, softening, localization, centrifuge, or calibrated elastoplastic FEM",
                "Parameters estimated": "chimin, rs",
                "Parameters fixed": "eta(t,z), cb0', Delta tan phi, kH, kE",
                "Output checked": "Peak/post-peak ratio and localized failure mode",
                "Risk if skipped": "Brittleness cap remains non-unique and should not be fitted",
            },
            {
                "Stage": "6",
                "Data source": "Independent footing geometry, treatment depth, or exposure history",
                "Parameters estimated": "None, validation only",
                "Parameters fixed": "All calibrated parameters",
                "Output checked": "qu, settlement, beta50, pf50",
                "Risk if skipped": "Calibration may reproduce one curve but fail outside the fitted case",
            },
        ]
    )

    identifiability.to_csv(DATA / "modifier_calibration_identifiability_table.csv", index=False)
    workflow.to_csv(DATA / "modifier_calibration_workflow_table.csv", index=False)

    readme = """# Modifier calibration and identifiability supplement

This supplement documents the calibration role of chiH, chiE and chis in the
screening model. The modifiers are not universal material constants. They are
project-calibrated or sensitivity-controlled factors used to separate element
scale cementation measurements from footing-scale mobilization, serviceability
stiffness, and localization or brittle response.

The recommended use is staged calibration:

1. Constrain eta(t,z), eta0, eta_r and lambda_e from independent cementation
   and durability observations.
2. Estimate cb0' and Delta tan phi from element tests.
3. Estimate kH only from bearing-capacity data with more than one treated depth.
4. Estimate kE only from settlement or stiffness data, preferably across more
   than one treated depth.
5. Estimate chimin and rs only if post-peak or localization-sensitive data are
   available; otherwise use conservative sensitivity envelopes.
6. Validate on data not used during calibration.

The companion script modifier_identifiability_analysis.py generates local
sensitivity, correlation, and profile-objective diagnostics that show why a
single load-settlement curve cannot identify all modifier and strength
parameters simultaneously.
"""
    (SUPP / "README_modifier_calibration_identifiability.md").write_text(readme, encoding="utf-8")
    print("Wrote modifier calibration guidance tables and README.")


if __name__ == "__main__":
    main()


