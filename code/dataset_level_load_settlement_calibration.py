from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT / "04 Supplemental data and code" / "Supplementary files"
DATA = SUPP / "data"
FIGS = SUPP / "figures"
CODE = SUPP / "code"


SOURCE = "Kulkarni et al. (2021), Figure 6 and Table 6; DOI:10.5614/j.eng.technol.sci.2021.53.6.2"


def pixel_to_data(x_px: float, y_px: float) -> tuple[float, float]:
    # Calibration from Figure 6 crop:
    # x=131.8 px -> 0 kPa, x=926.0 px -> 270 kPa;
    # y=166.6 px -> 0 mm, y=650.6 px -> 40 mm settlement.
    q = (x_px - 131.8) * 270.0 / (926.0 - 131.8)
    s = (y_px - 166.6) * 40.0 / (650.6 - 166.6)
    return max(q, 0.0), max(s, 0.0)


def digitized_points() -> pd.DataFrame:
    untreated_px = [
        (131.8, 166.6),
        (140.5, 193.4),
        (149.8, 208.9),
        (168.5, 223.8),
        (185.5, 259.5),
        (200.4, 291.3),
        (217.0, 331.9),
        (231.0, 387.5),
        (246.1, 436.1),
        (267.2, 488.1),
    ]
    treated_px = [
        (131.8, 166.6),
        (171.8, 178.1),
        (209.9, 183.6),
        (287.3, 192.6),
        (365.2, 205.1),
        (443.3, 209.9),
        (521.1, 231.7),
        (597.8, 275.8),
        (677.5, 308.1),
        (753.9, 369.1),
        (830.3, 438.3),
        (910.2, 511.5),
    ]
    rows = []
    for treatment, pts in [("untreated", untreated_px), ("treated", treated_px)]:
        for x, y in pts:
            q, s = pixel_to_data(x, y)
            rows.append(
                {
                    "settlement_mm": s,
                    "pressure_kpa": q,
                    "geometry": "square_120mm_x_120mm",
                    "treatment": treatment,
                    "source": SOURCE,
                    "digitization_note": "Manual pixel digitization from the open PDF crop of Figure 6.",
                }
            )
    return pd.DataFrame(rows)


def hyperbolic_q(s: np.ndarray, k0: float, qu: float) -> np.ndarray:
    return k0 * s / (1.0 + k0 * s / qu)


def fit_k0(df: pd.DataFrame, qu_fixed: float, role: str | None = None) -> float:
    fit_df = df[df["pressure_kpa"] <= qu_fixed + 1e-9].copy()
    if role is not None:
        fit_df = fit_df[fit_df["role"] == role].copy()
    s = fit_df["settlement_mm"].to_numpy()
    q = fit_df["pressure_kpa"].to_numpy()

    def residual(p):
        return (hyperbolic_q(s, np.exp(p[0]), qu_fixed) - q) / qu_fixed

    res = least_squares(residual, x0=np.log([12.0]), bounds=(np.log([0.1]), np.log([1000.0])))
    return float(np.exp(res.x[0]))


def metrics(obs: pd.DataFrame, k0: float, qu: float, role: str, treatment: str, params: str) -> dict[str, float | str]:
    eval_df = obs[(obs["pressure_kpa"] <= qu + 1e-9) & (obs["role"] == role)].copy()
    s = eval_df["settlement_mm"].to_numpy()
    q_obs = eval_df["pressure_kpa"].to_numpy()
    q_pred = hyperbolic_q(s, k0, qu)
    err = q_pred - q_obs
    nonzero = q_obs >= 5.0
    denom = np.maximum(np.abs(q_obs[nonzero]), 1.0)
    rmse = float(np.sqrt(np.mean(err**2)))
    return {
        "Source": "Kulkarni et al. 2021",
        "Geometry": "square 120 mm x 120 mm",
        "Role": role,
        "Treatment": treatment,
        "Calibrated parameters": params,
        "NRMSE": rmse / qu,
        "MAPE": float(np.mean(np.abs(err[nonzero]) / denom) * 100.0),
        "qu error (%)": 0.0,
        "SRF error (%)": np.nan,
        "Bias (kPa)": float(np.mean(err)),
        "Interpretation": "Pre-ultimate branch defined by Table 6 qu.",
    }


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    CODE.mkdir(parents=True, exist_ok=True)

    df = digitized_points()
    for treatment in ["untreated", "treated"]:
        mask = df["treatment"] == treatment
        idx = df[mask].sort_values("settlement_mm").index.to_list()
        # Use the early branch for calibration and reserve later pre-ultimate points
        # for validation. Points beyond Table 6 qu are retained in the CSV but
        # not used in the hyperbolic fit because qu is a tangent-defined limit.
        n_cal = 5 if treatment == "untreated" else 6
        for n, row_idx in enumerate(idx):
            df.loc[row_idx, "role"] = "calibration" if n < n_cal else "validation"
    df.to_csv(DATA / "kulkarni2021_digitized_plate_load_data.csv", index=False)
    for treatment in ["untreated", "treated"]:
        df[df["treatment"] == treatment].to_csv(
            DATA / f"kulkarni2021_square_120mm_{treatment}.csv", index=False
        )

    table6 = pd.DataFrame(
        [
            ["50 diameter", "circular", 19.0, 13.5, 93.0, 7.1, 4.89, 0.526, 53.76],
            ["100 diameter", "circular", 42.0, 9.2, 133.0, 5.4, 3.16, 0.587, 53.76],
            ["120 diameter", "circular", 29.0, 7.2, 170.0, 3.6, 5.86, 0.500, 53.76],
            ["50 x 50", "square", 25.0, 17.4, 111.5, 7.3, 4.46, 0.420, 34.90],
            ["100 x 100", "square", 49.0, 16.2, 166.0, 5.3, 3.38, 0.327, 34.90],
            ["120 x 120", "square", 41.0, 14.6, 177.0, 4.4, 4.31, 0.300, 34.90],
        ],
        columns=[
            "plate_size",
            "geometry",
            "qu_untreated_kpa",
            "settlement_untreated_mm",
            "qu_treated_kpa",
            "settlement_treated_mm",
            "BCR",
            "SRF",
            "settlement_reduction_percent",
        ],
    )
    table6["source"] = SOURCE
    table6.to_csv(DATA / "kulkarni2021_table6_bcr_srf_summary.csv", index=False)

    unt = df[df["treatment"] == "untreated"]
    tr = df[df["treatment"] == "treated"]
    qu_u = 41.0
    qu_t = 177.0
    km = fit_k0(unt, qu_u, "calibration")
    k0_t = fit_k0(tr, qu_t, "calibration")
    aE_equiv = k0_t / km - 1.0
    srf_pred = qu_t / k0_t / (qu_u / km)
    # Table 6 defines SRF as delta'_u/delta_u; for 120 mm square it is 0.30.
    srf_obs = 0.30

    met = [
        metrics(unt, km, qu_u, "calibration", "untreated", f"km={km:.3f} kPa/mm; qu fixed from Table 6"),
        metrics(unt, km, qu_u, "validation", "untreated", f"km={km:.3f} kPa/mm; qu fixed from Table 6"),
        metrics(
            tr,
            k0_t,
            qu_t,
            "calibration",
            "treated",
            f"k0={k0_t:.3f} kPa/mm; aE_equiv={aE_equiv:.3f}; qu fixed from Table 6",
        ),
        metrics(
            tr,
            k0_t,
            qu_t,
            "validation",
            "treated",
            f"k0={k0_t:.3f} kPa/mm; aE_equiv={aE_equiv:.3f}; qu fixed from Table 6",
        ),
    ]
    for row in met:
        if row["Treatment"] == "treated":
            row["SRF error (%)"] = float((srf_pred - srf_obs) / srf_obs * 100.0)
    pd.DataFrame(met).to_csv(DATA / "dataset_level_load_settlement_metrics.csv", index=False)

    param = pd.DataFrame(
        [
            {"parameter": "km", "value": km, "units": "kPa/mm", "meaning": "Untreated matrix initial stiffness"},
            {"parameter": "k0_treated", "value": k0_t, "units": "kPa/mm", "meaning": "Treated initial stiffness"},
            {"parameter": "aE_equiv", "value": aE_equiv, "units": "-", "meaning": "Equivalent stiffness multiplier for eta_eff*chi_E=1"},
            {"parameter": "qu_untreated", "value": qu_u, "units": "kPa", "meaning": "Table 6 ultimate pressure for square 120 mm untreated"},
            {"parameter": "qu_treated", "value": qu_t, "units": "kPa", "meaning": "Table 6 ultimate pressure for square 120 mm treated"},
        ]
    )
    param.to_csv(DATA / "dataset_level_load_settlement_fitted_parameters.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.0, 4.8), dpi=180)
    for treatment, color, marker, qu, k0 in [
        ("untreated", "#1f5a7a", "D", qu_u, km),
        ("treated", "#b43b3b", "s", qu_t, k0_t),
    ]:
        obs = df[df["treatment"] == treatment]
        for role, face in [("calibration", color), ("validation", "none")]:
            part = obs[obs["role"] == role]
            ax.scatter(
                part["settlement_mm"],
                part["pressure_kpa"],
                facecolors=face,
                edgecolors=color,
                marker=marker,
                label=f"{treatment} {role}",
            )
        smax = obs.loc[obs["pressure_kpa"] <= qu + 1e-9, "settlement_mm"].max()
        ss = np.linspace(0, smax, 150)
        ax.plot(ss, hyperbolic_q(ss, k0, qu), color=color, lw=2, label=f"{treatment} fitted")
        ax.axhline(qu, color=color, lw=0.8, ls=":", alpha=0.6)
    ax.set_xlabel("Settlement, s (mm)")
    ax.set_ylabel("Pressure, q (kPa)")
    ax.set_title("Dataset-level load-settlement calibration")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "figure_dataset_load_settlement_validation.png", bbox_inches="tight")
    plt.close(fig)

    (CODE / "dataset_level_load_settlement_calibration.py").write_text(Path(__file__).read_text(encoding="utf-8"), encoding="utf-8")
    print(DATA / "kulkarni2021_digitized_plate_load_data.csv")
    print(DATA / "dataset_level_load_settlement_metrics.csv")
    print(FIGS / "figure_dataset_load_settlement_validation.png")


if __name__ == "__main__":
    main()
