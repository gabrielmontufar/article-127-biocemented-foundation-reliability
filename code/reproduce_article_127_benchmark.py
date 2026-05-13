from __future__ import annotations

import json
import math
from pathlib import Path
from statistics import NormalDist

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data"
OUT_DIR = ROOT / "reproduced_outputs"
NORM = NormalDist()


def bearing_factors(phi_deg):
    phi = np.radians(phi_deg)
    n_q = np.exp(np.pi * np.tan(phi)) * np.tan(np.radians(45.0) + phi / 2.0) ** 2
    n_c = (n_q - 1.0) / np.tan(phi)
    n_gamma = 2.0 * (n_q + 1.0) * np.tan(phi)
    return n_c, n_q, n_gamma


def capacity(b, h_over_b, eta0, eta_r, lam, cb0, phi_m_deg, delta_tan_phi, gamma_eff, d_f, t_year):
    eta = eta_r + (eta0 - eta_r) * np.exp(-lam * t_year)
    tan_phi = np.tan(np.radians(phi_m_deg)) + delta_tan_phi * eta
    phi_deg = np.degrees(np.arctan(tan_phi))
    cb = cb0 * eta**1.15
    chi_h = 1.0 - np.exp(-h_over_b / 0.70)
    chi_s = np.clip(1.0 - 0.18 * eta**2.0, 0.72, 1.0)
    n_c, n_q, n_g = bearing_factors(phi_deg)
    _, n_q_m, n_g_m = bearing_factors(phi_m_deg)
    q_base = gamma_eff * d_f
    q_matrix = q_base * n_q_m + 0.5 * gamma_eff * b * n_g_m
    delta_cement = cb * n_c + q_base * (n_q - n_q_m) + 0.5 * gamma_eff * b * (n_g - n_g_m)
    q_ult = q_matrix + chi_h * chi_s * delta_cement
    return eta, phi_deg, cb, q_ult


def settlement_index(q_service, q_ult, eta, h_over_b):
    return 0.18 * q_service / (q_ult * (1.0 + 4.2 * eta * (1.0 - np.exp(-h_over_b / 0.55))))


def main():
    OUT_DIR.mkdir(exist_ok=True)
    years = np.linspace(0, 50, 101)
    dtan = math.tan(math.radians(35.0)) - math.tan(math.radians(32.0))
    rows = []
    for t in years:
        eta, phi, cb, qu = capacity(2.0, 1.0, 1.0, 0.22, 0.035, 20.0, 32.0, dtan, 18.0, 1.0, t)
        _, _, _, qu_matrix = capacity(2.0, 0.0, 1.0, 0.22, 0.035, 20.0, 32.0, dtan, 18.0, 1.0, t)
        _, _, _, qu_nd = capacity(2.0, 1.0, 1.0, 1.0, 0.0, 20.0, 32.0, dtan, 18.0, 1.0, t)
        rows.append({
            "time_year": t,
            "eta": float(eta),
            "phi_deg": float(phi),
            "c_b_kpa": float(cb),
            "q_u_bdc_kpa": float(qu),
            "q_u_matrix_kpa": float(qu_matrix),
            "q_u_non_degrading_kpa": float(qu_nd),
            "factor_of_safety": float(qu / 900.0),
            "settlement_index": float(settlement_index(900.0, qu, eta, 1.0)),
        })
    det = pd.DataFrame(rows)
    det.to_csv(OUT_DIR / "deterministic_time_history_reproduced.csv", index=False)

    if (DATA_DIR / "deterministic_time_history.csv").exists():
        ref = pd.read_csv(DATA_DIR / "deterministic_time_history.csv")
        merged = det.merge(ref, on="time_year", suffixes=("_new", "_reference"))
        max_abs = {}
        for col in ["q_u_bdc_kpa", "factor_of_safety", "settlement_index"]:
            max_abs[col] = float(np.max(np.abs(merged[f"{col}_new"] - merged[f"{col}_reference"])))
    else:
        max_abs = {"warning": "reference CSV not found"}

    summary = {"status": "reproduced", "max_absolute_difference_against_reference": max_abs}
    (OUT_DIR / "reproduction_check.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
