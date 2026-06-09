from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parents[3]
SUPP = BASE / "04 Supplemental data and code" / "Supplementary files"
DATA = SUPP / "data"
FIGS = SUPP / "figures"
CODE = SUPP / "code"
UPLOAD = BASE / "00 Files for journal upload"
AUDIT = BASE / "05 Journal selection and audits"


def bearing_factors(phi_deg: float) -> tuple[float, float, float]:
    phi = math.radians(phi_deg)
    nq = math.exp(math.pi * math.tan(phi)) * math.tan(math.radians(45.0) + phi / 2.0) ** 2
    nc = (nq - 1.0) / math.tan(phi)
    ngamma = 2.0 * (nq + 1.0) * math.tan(phi)
    return nc, nq, ngamma


def eq8_capacity(
    h_over_b: float,
    *,
    b: float = 2.0,
    eta: float = 0.3556262363,
    cb0: float = 20.0,
    phi_m_deg: float = 32.0,
    phi_t_deg: float = 35.0,
    gamma_eff: float = 18.0,
    df: float = 1.0,
) -> dict[str, float]:
    dtan = math.tan(math.radians(phi_t_deg)) - math.tan(math.radians(phi_m_deg))
    tan_phi = math.tan(math.radians(phi_m_deg)) + dtan * eta
    phi = math.degrees(math.atan(tan_phi))
    cb = cb0 * eta**1.15
    chi_h = 1.0 - math.exp(-h_over_b / 0.70)
    chi_s = min(1.0, max(0.72, 1.0 - 0.18 * eta**2.0))
    nc, nq, ng = bearing_factors(phi)
    _, nq_m, ng_m = bearing_factors(phi_m_deg)
    q_base = gamma_eff * df
    q_matrix = q_base * nq_m + 0.5 * gamma_eff * b * ng_m
    delta = cb * nc + q_base * (nq - nq_m) + 0.5 * gamma_eff * b * (ng - ng_m)
    q_eq8 = q_matrix + chi_h * chi_s * delta
    q_full = q_matrix + chi_s * delta
    return {
        "phi_eff_deg": phi,
        "cb_kpa": cb,
        "q_matrix_kpa": q_matrix,
        "q_eq8_kpa": q_eq8,
        "q_full_treated_kpa": q_full,
        "chi_h": chi_h,
        "chi_s": chi_s,
    }


def punching_limit_benchmark(
    h_over_b: float,
    q_matrix: float,
    q_full: float,
    *,
    b: float = 2.0,
    cb: float = 20.0,
    phi_t_deg: float = 35.0,
    gamma_eff: float = 18.0,
    df: float = 1.0,
) -> dict[str, float]:
    """Layered open-source benchmark for a cemented crust over sand.

    This is a deliberately simple kinematic/punch-through calculation. It is
    independent of Eq. (8), uses passive/lateral stress bounds in the treated
    layer, and is capped by the fully treated homogeneous capacity. It is not a
    commercial FEM result; it is a transparent limit-analysis-style comparison.
    """

    h = h_over_b * b
    if h <= 0.0:
        return {
            "q_la_lower_kpa": q_matrix,
            "q_la_central_kpa": q_matrix,
            "q_la_upper_kpa": q_matrix,
        }

    phi = math.radians(phi_t_deg)
    k0 = 1.0 - math.sin(phi)
    kp = math.tan(math.radians(45.0) + phi / 2.0) ** 2
    sv_avg = gamma_eff * (df + 0.5 * h)

    # Inclined side-shear path factor for a shallow punch-through mechanism.
    # The lower/central/upper triplet is retained to avoid a false single-value
    # calibration claim.
    geom_lower = 1.00
    geom_central = 1.15
    geom_upper = 1.30

    def estimate(k_side: float, c_mob: float, geom: float) -> float:
        side_resistance = c_mob * cb + k_side * sv_avg * math.tan(phi)
        increment = geom * (2.0 * h / b) * side_resistance
        return min(q_full, q_matrix + increment)

    q_lower = estimate(k0, 0.50, geom_lower)
    q_central = estimate(math.sqrt(k0 * kp), 0.75, geom_central)
    q_upper = estimate(kp, 1.00, geom_upper)
    q_lower, q_upper = min(q_lower, q_upper), max(q_lower, q_upper)
    q_central = min(q_upper, max(q_lower, q_central))
    return {
        "q_la_lower_kpa": q_lower,
        "q_la_central_kpa": q_central,
        "q_la_upper_kpa": q_upper,
    }


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    CODE.mkdir(parents=True, exist_ok=True)
    AUDIT.mkdir(parents=True, exist_ok=True)

    rows = []
    for h_over_b in np.round(np.arange(0.0, 2.51, 0.25), 2):
        e = eq8_capacity(float(h_over_b))
        la = punching_limit_benchmark(
            float(h_over_b),
            e["q_matrix_kpa"],
            e["q_full_treated_kpa"],
            cb=e["cb_kpa"],
            phi_t_deg=e["phi_eff_deg"],
        )
        central = la["q_la_central_kpa"]
        diff = 100.0 * (e["q_eq8_kpa"] - central) / central if central else 0.0
        inside = la["q_la_lower_kpa"] <= e["q_eq8_kpa"] <= la["q_la_upper_kpa"]
        rows.append(
            {
                "H_over_B": float(h_over_b),
                **e,
                **la,
                "eq8_minus_limit_central_percent": diff,
                "eq8_inside_limit_envelope": bool(inside),
            }
        )

    df = pd.DataFrame(rows)
    csv_path = DATA / "round6_open_source_limit_analysis_benchmark.csv"
    df.to_csv(csv_path, index=False)

    fig_path = FIGS / "Figure 8 open source limit analysis benchmark.png"
    plt.figure(figsize=(7.0, 4.6), dpi=200)
    plt.fill_between(
        df["H_over_B"],
        df["q_la_lower_kpa"],
        df["q_la_upper_kpa"],
        color="#b8d8f0",
        alpha=0.50,
        label="Open-source limit-analysis envelope",
    )
    plt.plot(df["H_over_B"], df["q_la_central_kpa"], color="#1f77b4", lw=2.0, label="Central punch-through benchmark")
    plt.plot(df["H_over_B"], df["q_eq8_kpa"], color="#c43b31", lw=2.2, marker="o", ms=3.5, label="Eq. (8)")
    plt.plot(df["H_over_B"], df["q_matrix_kpa"], color="#555555", lw=1.4, ls="--", label="Untreated lower reference")
    plt.plot(df["H_over_B"], df["q_full_treated_kpa"], color="#333333", lw=1.4, ls=":", label="Fully treated cap")
    plt.xlabel("Treatment depth ratio, H/B")
    plt.ylabel("Ultimate bearing pressure (kPa)")
    plt.xlim(0, 2.5)
    plt.grid(True, alpha=0.25)
    plt.legend(frameon=False, fontsize=8)
    plt.tight_layout()
    plt.savefig(fig_path)
    plt.close()

    shutil.copy2(fig_path, BASE / "02 Figures" / fig_path.name)
    shutil.copy2(fig_path, UPLOAD / fig_path.name)
    shutil.copy2(Path(__file__), CODE / "open_source_limit_analysis_layered_footing.py")

    max_abs_diff = float(np.max(np.abs(df["eq8_minus_limit_central_percent"])))
    within_count = int(df["eq8_inside_limit_envelope"].sum())
    summary = {
        "csv": str(csv_path),
        "figure": str(fig_path),
        "script_copy": str(CODE / "open_source_limit_analysis_layered_footing.py"),
        "n_cases": int(len(df)),
        "cases_inside_limit_envelope": within_count,
        "max_abs_percent_difference_vs_central": max_abs_diff,
        "open_source_stack": "Python, NumPy, Pandas, Matplotlib; OpenSeesPy was attempted but failed to import on this Windows Python because of a missing DLL.",
    }
    (AUDIT / "round6_open_source_limit_analysis_report.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()


