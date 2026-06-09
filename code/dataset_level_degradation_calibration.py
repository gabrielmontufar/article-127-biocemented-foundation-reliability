from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parents[3]
SUPP = ROOT / "04 Supplemental data and code" / "Supplementary files"
DATA = SUPP / "data"
FIGS = SUPP / "figures"
CODE = SUPP / "code"

SOURCE = (
    "Sharma, Satyam, and Reddy (2021), Hybrid bacteria mediated cemented sand; "
    "durability reductions reported in abstract; DOI:10.1177/1056789521991196"
)


def retention_data() -> pd.DataFrame:
    cycles = np.array([0, 5, 10, 15, 20], dtype=float)
    ucs_retention = np.array([1.000, 0.958, 0.917, 0.830, 0.650], dtype=float)
    return pd.DataFrame(
        {
            "exposure_type": "freeze_thaw_cycles",
            "exposure_value": cycles,
            "variable": "UCS",
            "value": ucs_retention,
            "retention": ucs_retention,
            "role": ["calibration", "calibration", "calibration", "validation", "validation"],
            "source": SOURCE,
            "digitization_note": "Numerical reductions are reported in the open abstract: 4.2%, 8.3%, 17%, and 35% after 5, 10, 15, and 20 cycles.",
        }
    )


def eta_model(E: np.ndarray, eta_r: float, lam: float) -> np.ndarray:
    return eta_r + (1.0 - eta_r) * np.exp(-lam * E)


def fit_model(df: pd.DataFrame) -> tuple[float, float]:
    cal = df[df["role"] == "calibration"]
    E = cal["exposure_value"].to_numpy()
    R = cal["retention"].to_numpy()

    def residual(p):
        eta_r = 1.0 / (1.0 + np.exp(-p[0]))
        lam = np.exp(p[1])
        return eta_model(E, eta_r, lam) - R

    res = least_squares(residual, x0=[-1.5, np.log(0.01)], bounds=([-8, np.log(1e-5)], [4, np.log(1.0)]))
    eta_r = float(1.0 / (1.0 + np.exp(-res.x[0])))
    lam = float(np.exp(res.x[1]))
    return eta_r, lam


def calc_metrics(df: pd.DataFrame, eta_r: float, lam: float) -> pd.DataFrame:
    rows = []
    for role, part in df.groupby("role"):
        E = part["exposure_value"].to_numpy()
        obs = part["retention"].to_numpy()
        pred = eta_model(E, eta_r, lam)
        err = pred - obs
        ss_res = float(np.sum(err**2))
        ss_tot = float(np.sum((obs - obs.mean()) ** 2))
        rows.append(
            {
                "Source": "Sharma et al. 2021",
                "Exposure type": "freeze-thaw cycles",
                "Variable": "UCS retention",
                "Role": role,
                "Fitted parameters": f"eta_r={eta_r:.3f}; lambda_E={lam:.4f} per cycle",
                "RMSE": float(np.sqrt(np.mean(err**2))),
                "NRMSE": float(np.sqrt(np.mean(err**2)) / max(obs.max() - obs.min(), 1e-9)),
                "MAPE": float(np.mean(np.abs(err) / np.maximum(obs, 1e-9)) * 100.0),
                "R2": np.nan if ss_tot == 0 else 1.0 - ss_res / ss_tot,
                "Final-retention error (%)": float((pred[-1] - obs[-1]) * 100.0),
                "Interpretation": "Controlled environmental cycling series; used as durability-identifiability check.",
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    CODE.mkdir(parents=True, exist_ok=True)

    df = retention_data()
    eta_r, lam = fit_model(df)
    df["predicted_retention"] = eta_model(df["exposure_value"].to_numpy(), eta_r, lam)
    df.to_csv(DATA / "sharma2021_digitized_durability_retention_data.csv", index=False)
    metrics = calc_metrics(df, eta_r, lam)
    metrics.to_csv(DATA / "dataset_level_degradation_metrics.csv", index=False)

    params = pd.DataFrame(
        [
            {"parameter": "eta_r", "value": eta_r, "units": "-", "meaning": "Residual cementation-state retention inferred from UCS retention"},
            {"parameter": "lambda_E", "value": lam, "units": "per freeze-thaw cycle", "meaning": "Exposure-clock degradation rate"},
        ]
    )
    params.to_csv(DATA / "dataset_level_degradation_fitted_parameters.csv", index=False)

    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=180)
    for role, marker, color in [("calibration", "o", "#1f5a7a"), ("validation", "s", "#b43b3b")]:
        part = df[df["role"] == role]
        ax.scatter(part["exposure_value"], part["retention"], marker=marker, color=color, label=f"{role} observed")
    E = np.linspace(0, 20, 200)
    ax.plot(E, eta_model(E, eta_r, lam), color="black", lw=2, label="fitted degradation law")
    ax.set_xlabel("Freeze-thaw cycles")
    ax.set_ylabel("UCS retention")
    ax.set_ylim(0.55, 1.05)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGS / "figure_dataset_degradation_validation.png", bbox_inches="tight")
    plt.close(fig)

    (CODE / "dataset_level_degradation_calibration.py").write_text(Path(__file__).read_text(encoding="utf-8"), encoding="utf-8")
    print(DATA / "sharma2021_digitized_durability_retention_data.csv")
    print(DATA / "dataset_level_degradation_metrics.csv")
    print(FIGS / "figure_dataset_degradation_validation.png")


if __name__ == "__main__":
    main()


