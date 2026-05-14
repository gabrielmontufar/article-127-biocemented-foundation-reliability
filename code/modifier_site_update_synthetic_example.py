from __future__ import annotations

from pathlib import Path
from math import erf

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT
DATA = SUPP / "data"
FIGS = SUPP / "figures"


def chi_h(h_over_b: np.ndarray | float, kh: float) -> np.ndarray | float:
    return 1.0 - np.exp(-np.asarray(h_over_b) / kh)


def weighting(z_over_b: np.ndarray, kz: float) -> np.ndarray:
    raw = np.exp(-z_over_b / kz)
    area = np.trapezoid(raw, z_over_b)
    return raw / area


def normal_pf(beta: float) -> float:
    return 0.5 * (1.0 + erf(-beta / np.sqrt(2.0)))


def calibrate_kh_from_plate_load() -> dict[str, float]:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    h_over_b = np.array([0.50, 0.75, 1.00, 1.50, 2.00])
    kh_true = 0.72
    bcr_max_gain = 0.54
    deterministic_error = np.array([-0.015, 0.010, -0.004, 0.008, -0.006])
    bcr_obs = 1.0 + bcr_max_gain * chi_h(h_over_b, kh_true) + deterministic_error

    grid = np.linspace(0.25, 1.80, 400)
    rows = []
    for kh in grid:
        pred = 1.0 + bcr_max_gain * chi_h(h_over_b, kh)
        residual = pred - bcr_obs
        rows.append(
            {
                "kH": kh,
                "objective": float(np.mean(residual**2)),
                "rmse_bcr": float(np.sqrt(np.mean(residual**2))),
                "mape_percent": float(np.mean(np.abs(residual / bcr_obs)) * 100.0),
                "max_abs_error_bcr": float(np.max(np.abs(residual))),
                "bias_bcr": float(np.mean(residual)),
            }
        )
    objective = pd.DataFrame(rows)
    best = objective.loc[objective["objective"].idxmin()].to_dict()
    best_kh = float(best["kH"])

    fit = pd.DataFrame(
        {
            "H_over_B": h_over_b,
            "BCR_observed_synthetic": bcr_obs,
            "BCR_fitted": 1.0 + bcr_max_gain * chi_h(h_over_b, best_kh),
            "BCR_default_kH_0p70": 1.0 + bcr_max_gain * chi_h(h_over_b, 0.70),
        }
    )
    objective.to_csv(DATA / "kH_profile_objective.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.0, 4.4), dpi=180)
    x = np.linspace(0.0, 2.2, 200)
    ax.plot(x, 1.0 + bcr_max_gain * chi_h(x, best_kh), color="#1f77b4", lw=2.2, label=f"fitted kH = {best_kh:.2f}")
    ax.plot(x, 1.0 + bcr_max_gain * chi_h(x, 0.40), "--", color="#777777", lw=1.4, label="low kH = 0.40")
    ax.plot(x, 1.0 + bcr_max_gain * chi_h(x, 1.20), ":", color="#444444", lw=1.8, label="high kH = 1.20")
    ax.scatter(h_over_b, bcr_obs, s=42, color="#b23a48", zorder=5, label="synthetic PLT depth series")
    ax.set_xlabel("Treated depth, H/B")
    ax.set_ylabel("Bearing-capacity ratio, BCR")
    ax.set_title("Depth-benefit calibration for kH")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "figure_kH_depth_benefit_fit.png")
    plt.close(fig)

    return {
        "kH_fit": best_kh,
        "kH_rmse_bcr": float(best["rmse_bcr"]),
        "kH_mape_percent": float(best["mape_percent"]),
        "kH_max_abs_error_bcr": float(best["max_abs_error_bcr"]),
        "kH_bias_bcr": float(best["bias_bcr"]),
    }


def calibrate_kz_from_cpt_profile() -> dict[str, float]:
    z_over_b = np.linspace(0.0, 3.0, 91)
    ir_cpt = 1.12 * np.exp(-z_over_b / 0.78) + 0.14 * np.exp(-((z_over_b - 1.35) / 0.36) ** 2)
    ir_cpt += 0.025 * np.sin(4.0 * z_over_b)
    target_kz = 0.88
    observed_effective_ir = float(np.trapezoid(weighting(z_over_b, target_kz) * ir_cpt, z_over_b) - 0.010)

    grid = np.linspace(0.25, 1.90, 400)
    rows = []
    for kz in grid:
        w = weighting(z_over_b, kz)
        effective = float(np.trapezoid(w * ir_cpt, z_over_b))
        rows.append(
            {
                "kz": kz,
                "effective_IR_CPT": effective,
                "target_effective_IR_CPT": observed_effective_ir,
                "objective": (effective - observed_effective_ir) ** 2,
            }
        )
    objective = pd.DataFrame(rows)
    best = objective.loc[objective["objective"].idxmin()].to_dict()
    best_kz = float(best["kz"])
    objective.to_csv(DATA / "kz_cpt_profile_objective.csv", index=False)

    fig, ax1 = plt.subplots(figsize=(7.0, 4.6), dpi=180)
    ax1.plot(ir_cpt, z_over_b, color="#b23a48", lw=2.2, label="synthetic CPT improvement ratio")
    ax1.invert_yaxis()
    ax1.set_xlabel("CPT improvement ratio, IR_CPT")
    ax1.set_ylabel("Depth, z/B")
    ax1.grid(True, alpha=0.22)
    ax2 = ax1.twiny()
    for kz, style, label in [(0.40, "--", "w(z), kz = 0.40"), (best_kz, "-", f"w(z), fitted kz = {best_kz:.2f}"), (1.50, ":", "w(z), kz = 1.50")]:
        ax2.plot(weighting(z_over_b, kz), z_over_b, style, lw=1.8, label=label)
    ax2.set_xlabel("Normalized influence weight, w(z)")
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax2.legend(lines + lines2, labels + labels2, frameon=False, fontsize=8, loc="lower right")
    ax1.set_title("CPT-profile updating for kz")
    fig.tight_layout()
    fig.savefig(FIGS / "figure_kz_cpt_weighting_fit.png")
    plt.close(fig)

    return {
        "kz_fit": best_kz,
        "kz_target_effective_ir": observed_effective_ir,
        "kz_fitted_effective_ir": float(best["effective_IR_CPT"]),
        "kz_objective": float(best["objective"]),
    }


def reliability_from_modifiers(kh: float, kz: float, chi_min: float, rs: float) -> dict[str, float]:
    eta50 = 0.25 + (0.92 - 0.25) * np.exp(-0.028 * 50.0)
    h_over_b = 1.0
    profile_factor = 0.92 + 0.13 * (0.85 / kz) ** 0.22
    eta_eff = eta50 * profile_factor
    chis = max(chi_min, 1.0 - rs * eta_eff**2)
    ch = float(chi_h(h_over_b, kh))
    qu = 760.0 + (2.10 * 105.0 * eta_eff**1.15 + 880.0 * 0.12 * eta_eff**0.85) * ch * chis
    beta = (qu - 950.0) / (0.13 * qu)
    return {
        "eta_eff_50yr": eta_eff,
        "chiH": ch,
        "chis_50yr": chis,
        "qu50_kPa": qu,
        "beta50": beta,
        "pf50": normal_pf(beta),
    }


def modifier_update_and_sensitivity(kh_fit: float, kz_fit: float) -> dict[str, float]:
    rows = []
    for kh in [0.40, kh_fit, 1.20]:
        for kz in [0.40, kz_fit, 1.50]:
            for chi_min in [0.65, 0.75, 0.85]:
                for rs in [0.10, 0.18, 0.25]:
                    out = reliability_from_modifiers(kh, kz, chi_min, rs)
                    rows.append({"kH": kh, "kz": kz, "chi_min": chi_min, "rs": rs, **out})
    sensitivity = pd.DataFrame(rows)
    sensitivity.to_csv(DATA / "chis_sensitivity_summary.csv", index=False)

    base = reliability_from_modifiers(kh_fit, kz_fit, 0.75, 0.18)
    summary = pd.DataFrame(
        [
            {
                "item": "synthetic_PLT_depth_series",
                "calibrated_parameter": "kH",
                "calibrated_value": kh_fit,
                "calibration_role": "Controls saturation of treated-depth benefit chi_H",
                "data_used": "Synthetic BCR values for H/B = 0.5, 0.75, 1.0, 1.5, and 2.0",
            },
            {
                "item": "synthetic_CPT_profile",
                "calibrated_parameter": "kz",
                "calibrated_value": kz_fit,
                "calibration_role": "Controls vertical weighting of CPT-derived improvement profile",
                "data_used": "Synthetic before-after CPT improvement profile and effective settlement improvement target",
            },
            {
                "item": "localization_sensitivity",
                "calibrated_parameter": "chi_min, rs",
                "calibrated_value": "not uniquely calibrated in this synthetic pre-peak example",
                "calibration_role": "Bounds localized or brittle mobilization loss through chi_s",
                "data_used": "Sensitivity envelope because post-peak or localization observations are absent",
            },
            {
                "item": "baseline_updated_prediction",
                "calibrated_parameter": "beta50",
                "calibrated_value": base["beta50"],
                "calibration_role": "Reliability outcome after staged modifier update",
                "data_used": "kH and kz fits with chi_min = 0.75 and rs = 0.18",
            },
        ]
    )
    summary.to_csv(DATA / "kH_kz_chis_calibration_summary.csv", index=False)

    eta = np.linspace(0.20, 1.00, 160)
    fig, ax = plt.subplots(figsize=(7.0, 4.4), dpi=180)
    for chi_min in [0.65, 0.75, 0.85]:
        for rs, style in [(0.10, "--"), (0.18, "-"), (0.25, ":")]:
            chis = np.maximum(chi_min, 1.0 - rs * eta**2)
            ax.plot(eta, chis, style, lw=1.4, label=f"chi_min={chi_min:.2f}, rs={rs:.2f}")
    ax.set_xlabel("Effective cementation state, eta_eff")
    ax.set_ylabel("Localization modifier, chi_s")
    ax.set_title("Fragility/localization sensitivity of chi_s")
    ax.set_ylim(0.60, 1.02)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=6.8, ncol=3)
    fig.tight_layout()
    fig.savefig(FIGS / "figure_chis_localization_sensitivity.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.2, 4.4), dpi=180)
    ax.axis("off")
    boxes = [
        (0.05, 0.64, "Untreated PLT\nfit phi_m, qu,m"),
        (0.28, 0.64, "Element tests\ncb0, Delta tan phi"),
        (0.52, 0.64, "CPT/Vs/UPV profile\nupdate kz"),
        (0.75, 0.64, "PLT depth series\nupdate kH"),
        (0.28, 0.24, "Post-peak/localization\nconstrain chi_s"),
        (0.55, 0.24, "Joint audit\nJ, RMSE, bias"),
        (0.78, 0.24, "Reliability screen\nbeta(t), pf(t)"),
    ]
    for x, y, text in boxes:
        rect = plt.Rectangle((x, y), 0.18, 0.18, ec="#333333", fc="#edf4fb", lw=1.1)
        ax.add_patch(rect)
        ax.text(x + 0.09, y + 0.09, text, ha="center", va="center", fontsize=8)
    arrows = [
        ((0.23, 0.73), (0.28, 0.73)),
        ((0.46, 0.73), (0.52, 0.73)),
        ((0.70, 0.73), (0.75, 0.73)),
        ((0.84, 0.64), (0.64, 0.42)),
        ((0.37, 0.64), (0.37, 0.42)),
        ((0.46, 0.33), (0.55, 0.33)),
        ((0.73, 0.33), (0.78, 0.33)),
    ]
    for start, end in arrows:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", lw=1.1, color="#333333"))
    ax.set_title("Staged site-updating workflow for kH, kz, and chi_s", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "figure_modifier_update_workflow.png")
    plt.close(fig)

    return base


def write_readme() -> None:
    text = """# Practical updating of kH, kz, and chi_s

This supplement provides a synthetic, reproducible site-updating example for
the footing-scale modifiers kH, kz, and chi_s. The data are not field
measurements; they are deliberately labeled synthetic and are used only to
demonstrate the calibration workflow.

The staged workflow is:

1. Fit untreated plate-load response before any treated-zone modifier is tuned.
2. Bound cb0' and Delta tan phi from element tests or conservative priors.
3. Update kz from depth-resolved CPT/CPTu, Vs, UPV, or carbonate profiles.
4. Update kH from a plate-load or numerical depth-benefit series.
5. Keep chi_s as a sensitivity parameter unless post-peak, localization, or
   calibrated elastoplastic numerical evidence is available.

Generated outputs include kH_profile_objective.csv,
kz_cpt_profile_objective.csv, chis_sensitivity_summary.csv,
kH_kz_chis_calibration_summary.csv, and four PNG figures documenting the
depth-benefit fit, CPT weighting fit, localization sensitivity, and staged
updating workflow.
"""
    (SUPP / "README_kH_kz_chis_site_updating.md").write_text(text, encoding="utf-8")


def main() -> None:
    kh = calibrate_kh_from_plate_load()
    kz = calibrate_kz_from_cpt_profile()
    base = modifier_update_and_sensitivity(kh["kH_fit"], kz["kz_fit"])
    write_readme()
    print("kH fit:", kh)
    print("kz fit:", kz)
    print("baseline reliability:", base)


if __name__ == "__main__":
    main()
