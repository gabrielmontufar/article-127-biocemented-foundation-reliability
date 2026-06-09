from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
SUPP = ROOT / "04 Supplemental data and code" / "Supplementary files"
DATA = SUPP / "data"
FIGS = SUPP / "figures"


def treatment_schedule_to_caco3(ninj, curea, cca, qinj, tcure, temp_c, activity):
    ceq = np.minimum(curea, cca) / 0.5
    ftemp = np.exp(-((temp_c - 30.0) / 18.0) ** 2)
    fflow = np.exp(-0.35 * np.maximum(qinj / 1.0 - 1.0, 0.0))
    feff = np.clip(activity * ftemp * fflow * (1.0 - np.exp(-tcure / 24.0)), 0.0, 1.0)
    ccmax = 12.0
    return ccmax * (1.0 - np.exp(-0.19 * ninj * ceq * feff))


def eta_from_caco3(cc, cv=0.25, a_eta=0.42, alpha_u=0.85):
    eta_raw = 1.0 - np.exp(-a_eta * cc)
    uniformity = np.exp(-alpha_u * cv)
    return np.clip(eta_raw * uniformity, 0.0, 1.0)


def cb0_from_caco3(cc, cbmax=85.0, bc=0.22):
    return cbmax * (1.0 - np.exp(-bc * cc))


def dtan_from_caco3(cc, aphi=0.22, bphi=0.16, fdensity=1.0, fdistribution=1.0):
    return aphi * (1.0 - np.exp(-bphi * cc)) * fdensity * fdistribution


def lambda_environment(temp_c, salinity, ph, si_calcite, hydraulic_i, wd_freq, kdamage):
    lambda_ref = 0.018
    ftemp = 1.08 ** ((temp_c - 25.0) / 10.0)
    fchem = 1.0 + 0.90 * np.maximum(0.0, -si_calcite) + 0.28 * np.maximum(0.0, 7.5 - ph)
    fsal = 1.0 + 0.18 * (salinity / 35.0) ** 0.8
    fflow = 1.0 + 0.22 * (hydraulic_i / 0.5) ** 0.7
    fwd = 1.0 + 0.06 * wd_freq
    ffrac = 1.0 + 0.35 * kdamage
    return lambda_ref * ftemp * fchem * fsal * fflow * fwd * ffrac


def write_tables() -> None:
    mapping = pd.DataFrame(
        [
            ["eta0", "CaCO3 content, spatial CaCO3 distribution, Vs/UPV, CPT, injection records", "eta0(Cc) times uniformity factor Uc", "Carbonate profiles, wave velocity, CPT/CPTu", "Carbonate morphology and spatial variability", "Confounded with cb0' if only capacity is used"],
            ["cb0'", "UCS, drained triaxial, direct shear, CaCO3", "xi_c UCS increment or saturating cb0'(Cc)", "Element tests on treated and untreated sand", "UCS-to-c' conversion and fabric", "Confounded with eta0 and chi_s"],
            ["Delta tan phi", "Drained triaxial/direct shear, dilatancy, density, CaCO3", "tan(phi'_treated)-tan(phi'_untreated)", "Treated vs untreated shear tests", "Density, fabric, dilation and confining stress", "Confounded with cb0' in footing tests"],
            ["lambda_e", "Wetting-drying, salinity, pH, calcite saturation, flushing, temperature, durability retention", "lambda_ref Fchem Fflow FWD Ftemp Fsalinity", "Retention curves for UCS, UPV, Vs, mass, stiffness", "Separating chemical, hydraulic and cyclic mechanisms", "Confounded with eta_r without long-duration data"],
        ],
        columns=["Parameter", "Measurable inputs", "Suggested mapping", "Calibration data", "Main uncertainty", "Identifiability warning"],
    )
    mapping.to_csv(DATA / "micp_measurable_inputs_mapping.csv", index=False)

    literature = pd.DataFrame(
        [
            ["Al Qabany and Soga (2013)", "chemical treatment concentration", "Treatment chemistry changes carbonate precipitation pattern and engineering properties", "eta0, cb0', Delta tan phi", "Supports using reagent concentration and CaCO3 distribution as priors"],
            ["Harkes et al. (2010)", "bacterial activity fixation and distribution", "Spatially variable activity controls carbonate precipitation for ground reinforcement", "eta0 and Uc", "Motivates uniformity factor and injection-record constraints"],
            ["Gomez et al. (2015); van Paassen et al. (2010)", "field/large-scale bio-cementation response", "Field delivery and monitoring affect treatment uniformity and scale-up", "eta0, lambda_e, kH, kz", "Links measurable treatment delivery to reliability inputs"],
            ["Montoya and DeJong (2015)", "stress-strain response of cemented sands", "Cementation changes strength and stiffness response", "cb0', Delta tan phi", "Supports element-test separation of cohesion-like and frictional increments"],
            ["Choi et al. (2020)", "compiled CaCO3, UCS, Mohr-Coulomb parameters, permeability", "Engineering properties are correlated with carbonate content but remain material dependent", "eta0, cb0', Delta tan phi", "Supports calibrable, non-universal maps"],
            ["Terzis and Laloui (2018)", "3-D micro-architecture and mechanical response", "Calcite bond architecture controls mechanical response", "Uc, chi_s, fracture/damage modifiers", "Supports morphology/uniformity factors"],
            ["Sharma et al. (2021)", "UCS, STS, UPV, shear wave velocity, durability", "Strength, stiffness and durability metrics respond to treatment cycles and exposure", "cb0', lambda_e", "Supports durability-retention calibration"],
            ["Konstantinou et al. (2023)", "fracture toughness and strength of bio-cemented sand", "Fracture properties vary with degree of cementation", "chi_s and lambda_e modifiers", "Supports fracture/bond-damage path for fragility and degradation"],
            ["Lin et al. (2023)", "cementation stress characteristic curve", "Cementation state can be related to stress-level dependent response", "eta0, cb0'", "Supports using cementation-state priors rather than raw CaCO3 alone"],
        ],
        columns=["Study", "Measured input", "Reported MICP effect", "Parameter informed", "How used here"],
    )
    literature.to_csv(DATA / "micp_literature_parameter_links.csv", index=False)


def make_figures() -> None:
    cc = np.linspace(0, 12, 200)
    fig, ax1 = plt.subplots(figsize=(7.2, 4.5), dpi=180)
    for cv, ls in [(0.10, "-"), (0.25, "--"), (0.45, ":")]:
        ax1.plot(cc, eta_from_caco3(cc, cv=cv), ls, lw=2.0, label=f"eta0, CV_CaCO3={cv:.2f}")
    ax1.set_xlabel("Calcium carbonate content, Cc (%)")
    ax1.set_ylabel("Initial effective cementation state, eta0")
    ax2 = ax1.twinx()
    ax2.plot(cc, cb0_from_caco3(cc), color="#b23a48", lw=2.0, label="cb0' saturating map")
    ax2.fill_between(cc, 0.75 * cb0_from_caco3(cc), 1.25 * cb0_from_caco3(cc), color="#b23a48", alpha=0.14, label="illustrative uncertainty band")
    ax2.set_ylabel("Initial cohesion-like increment cb0' (kPa)")
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, frameon=False, fontsize=7, loc="lower right")
    ax1.grid(True, alpha=0.25)
    ax1.set_title("Example mapping from CaCO3 to eta0 and cb0'")
    fig.tight_layout()
    fig.savefig(FIGS / "figure_caco3_to_eta_cb0.png")
    plt.close(fig)

    x = np.linspace(0, 1, 80)
    cases = pd.DataFrame(
        {
            "driver": ["pH/SI chemistry", "salinity", "hydraulic flushing", "wetting-drying frequency", "temperature", "fracture damage"],
            "lambda_e": [
                lambda_environment(25, 35, 6.7, -0.7, 0.5, 0, 0),
                lambda_environment(25, 70, 7.5, 0.0, 0.5, 0, 0),
                lambda_environment(25, 35, 7.5, 0.0, 1.2, 0, 0),
                lambda_environment(25, 35, 7.5, 0.0, 0.5, 12, 0),
                lambda_environment(40, 35, 7.5, 0.0, 0.5, 0, 0),
                lambda_environment(25, 35, 7.5, 0.0, 0.5, 0, 0.8),
            ],
        }
    )
    cases.to_csv(DATA / "micp_environment_to_lambda_cases.csv", index=False)
    fig, ax = plt.subplots(figsize=(7.2, 4.3), dpi=180)
    ax.barh(cases["driver"], cases["lambda_e"], color="#4c78a8")
    ax.axvline(0.018, color="#333333", lw=1.0, ls="--", label="lambda_ref")
    ax.set_xlabel("Effective degradation rate lambda_e (1/yr)")
    ax.set_title("Environmental acceleration factors for lambda_e")
    ax.grid(True, axis="x", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "figure_environment_to_lambda.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 4.4), dpi=180)
    ax.axis("off")
    boxes = [
        (0.03, 0.62, "Treatment schedule\nNinj, urea, Ca, flow,\ncuring, temperature,\nactivity"),
        (0.27, 0.62, "Measured state\nCaCO3(z), uniformity,\nVs/UPV, CPT"),
        (0.51, 0.62, "Element response\nUCS, triaxial,\ndirect shear"),
        (0.75, 0.62, "Durability\npH, salinity,\nflushing, WD,\ntemperature"),
        (0.18, 0.20, "eta0, cb0',\nDelta tan phi"),
        (0.50, 0.20, "lambda_e,\neta_r, chi_s"),
        (0.74, 0.20, "Eq. (8)\nreliability"),
    ]
    for x0, y0, text in boxes:
        rect = plt.Rectangle((x0, y0), 0.19, 0.20, ec="#333333", fc="#edf4fb", lw=1.0)
        ax.add_patch(rect)
        ax.text(x0 + 0.095, y0 + 0.10, text, ha="center", va="center", fontsize=7.3)
    arrows = [((0.22, 0.72), (0.27, 0.72)), ((0.46, 0.72), (0.51, 0.72)), ((0.70, 0.72), (0.75, 0.72)), ((0.36, 0.62), (0.28, 0.40)), ((0.60, 0.62), (0.28, 0.40)), ((0.84, 0.62), (0.59, 0.40)), ((0.37, 0.30), (0.50, 0.30)), ((0.69, 0.30), (0.74, 0.30))]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", color="#333333", lw=1.0))
    ax.set_title("Inference workflow from measurable MICP inputs to reliability parameters", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "figure_parameter_inference_workflow.png")
    plt.close(fig)


def bayesian_demo() -> None:
    rng = np.random.default_rng(127)
    ninj = np.array([2, 4, 6, 8, 10, 12])
    cc_true = treatment_schedule_to_caco3(ninj, 0.5, 0.5, 0.8, 48, 30, 0.90)
    cc_obs = cc_true + rng.normal(0.0, 0.35, len(ninj))
    eta_obs = eta_from_caco3(cc_obs, cv=0.22) + rng.normal(0.0, 0.025, len(ninj))
    ucs_obs = 110.0 + 5.8 * cb0_from_caco3(cc_obs) + rng.normal(0.0, 35.0, len(ninj))
    demo = pd.DataFrame({"Ninj": ninj, "CaCO3_percent_synthetic": cc_obs, "eta0_indicator_synthetic": eta_obs, "UCS_kPa_synthetic": ucs_obs})
    demo.to_csv(DATA / "micp_bayesian_parameter_update_demo.csv", index=False)
    summary = pd.DataFrame(
        [
            ["eta0", float(np.mean(eta_obs[-2:])), "updated from synthetic CaCO3 and eta indicator", "demonstration only"],
            ["cb0'", float(np.mean(cb0_from_caco3(cc_obs[-2:]))), "updated from synthetic CaCO3-strength map", "demonstration only"],
            ["Delta tan phi", float(np.mean(dtan_from_caco3(cc_obs[-2:]))), "updated from synthetic tangent-space map", "demonstration only"],
            ["lambda_e", float(lambda_environment(30, 35, 7.2, -0.2, 0.6, 8, 0.2)), "updated from synthetic environmental descriptors", "demonstration only"],
        ],
        columns=["Parameter", "Posterior_demo_mean", "Update basis", "Caveat"],
    )
    summary.to_csv(DATA / "micp_parameter_inference_summary.csv", index=False)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    write_tables()
    make_figures()
    bayesian_demo()
    readme = """# MICP parameter inference supplement

This supplement provides calibrable maps from measurable MICP treatment,
strength, stiffness, and durability indicators to eta0, cb0', Delta tan phi,
and lambda_e. The numerical demonstration is synthetic and is not presented as
a universal calibration.
"""
    (SUPP / "README_micp_parameter_inference.md").write_text(readme, encoding="utf-8")
    print(DATA / "micp_parameter_inference_summary.csv")


if __name__ == "__main__":
    main()


