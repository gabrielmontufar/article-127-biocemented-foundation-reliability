from __future__ import annotations

from pathlib import Path
from statistics import NormalDist

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT
DATA = SUPP / "data"
FIGS = SUPP / "figures"
CODE = SUPP / "code"
NORM = NormalDist()
Q_COV = 0.08
SI_THRESHOLD = 0.045


def lognormal_params(mean: np.ndarray | float, cov: np.ndarray | float):
    cov = np.asarray(cov)
    sigma = np.sqrt(np.log1p(cov**2))
    mu = np.log(mean) - 0.5 * sigma**2
    return mu, sigma


def lognormal_beta(mean_r: float, cov_r: float, mean_q: float, cov_q: float) -> tuple[float, float]:
    mu_r, sig_r = lognormal_params(mean_r, cov_r)
    mu_q, sig_q = lognormal_params(mean_q, cov_q)
    beta = float((mu_r - mu_q) / np.sqrt(sig_r**2 + sig_q**2))
    return beta, float(NORM.cdf(-beta))


def service_exceedance(mean_r: float, cov_r: float, mean_q: float, cov_q: float, stiffness_den: float) -> float:
    mu_r, sig_r = lognormal_params(mean_r, cov_r)
    mu_q, sig_q = lognormal_params(mean_q, cov_q)
    threshold = np.log(SI_THRESHOLD * stiffness_den / 0.18)
    z = (mu_q - mu_r - threshold) / np.sqrt(sig_r**2 + sig_q**2)
    return float(NORM.cdf(z))


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    det = pd.read_csv(DATA / "deterministic_time_history.csv")
    mc = pd.read_csv(DATA / "monte_carlo_reliability_time_history.csv")
    df = det.merge(mc, on="time_year", how="inner")
    pressures = [250.0, 500.0, 650.0, 900.0]
    years_table = [0.0, 10.0, 30.0, 50.0]

    rows = []
    for q in pressures:
        pf_life = None
        beta_life = None
        for _, r in df.iterrows():
            eta = float(r["eta"])
            stiffness_den = 1.0 + 4.2 * eta * (1.0 - np.exp(-1.0 / 0.55))
            beta, pf = lognormal_beta(float(r["mean_q_u_kpa"]), float(r["cov_q_u"]), q, Q_COV)
            p_service = service_exceedance(float(r["mean_q_u_kpa"]), float(r["cov_q_u"]), q, Q_COV, stiffness_den)
            if float(r["time_year"]) == 50.0:
                pf_life = pf
                beta_life = beta

            rows.append(
                {
                    "qservice_mean_kpa": q,
                    "time_year": float(r["time_year"]),
                    "qu_deterministic_kpa": float(r["q_u_bdc_kpa"]),
                    "qu_mc_mean_kpa": float(r["mean_q_u_kpa"]),
                    "cov_qu": float(r["cov_q_u"]),
                    "FS": float(r["q_u_bdc_kpa"] / q),
                    "sI": float(0.18 * q / (float(r["q_u_bdc_kpa"]) * stiffness_den)),
                    "pf_U": pf,
                    "beta_U": beta,
                    "p_service_index_gt_0_045": p_service,
                }
            )
        for row in rows:
            if row["qservice_mean_kpa"] == q:
                row["pf_life"] = pf_life
                row["beta_life"] = beta_life

    results = pd.DataFrame(rows)
    results.to_csv(DATA / "service_pressure_comparison_results.csv", index=False)
    lower = results[(results["qservice_mean_kpa"] == 500.0) & (results["time_year"].isin(years_table))].copy()
    lower.to_csv(DATA / "lower_demand_service_results.csv", index=False)

    final = results[results["time_year"] == 50.0].copy()
    final["Case"] = final["qservice_mean_kpa"].map(
        {
            250.0: "Low service-demand illustration",
            500.0: "Lower-demand service case",
            650.0: "Moderate-to-high service case",
            900.0: "High-demand stress test",
        }
    )
    final["Interpretation"] = final["qservice_mean_kpa"].map(
        {
            250.0: "Ultimate reliability remains high; serviceability proxy is low.",
            500.0: "Degradation reduces reserve but does not control ultimate reliability at 50 years.",
            650.0: "Reliability margin remains positive but serviceability becomes more relevant.",
            900.0: "Stress test exposes sensitivity to permanent initial cementation credit.",
        }
    )
    final[
        [
            "Case",
            "qservice_mean_kpa",
            "qu_deterministic_kpa",
            "FS",
            "pf_U",
            "beta_U",
            "p_service_index_gt_0_045",
            "Interpretation",
        ]
    ].to_csv(DATA / "lower_demand_maintext_table.csv", index=False)

    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=180)
    for q, style in [(250.0, ":"), (500.0, "-"), (650.0, "-."), (900.0, "--")]:
        sub = results[results["qservice_mean_kpa"] == q]
        ax.plot(sub["time_year"], sub["beta_U"], style, lw=2.0, label=f"qservice = {int(q)} kPa")
    ax.axhline(0.0, color="#333333", lw=0.9)
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Ultimate reliability index, beta_U(t)")
    ax.set_title("Reliability trajectories across service-pressure levels")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "figure_service_pressure_beta_comparison.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=180)
    for q, style in [(250.0, ":"), (500.0, "-"), (650.0, "-."), (900.0, "--")]:
        sub = results[results["qservice_mean_kpa"] == q]
        ax.semilogy(sub["time_year"], sub["pf_U"], style, lw=2.0, label=f"qservice = {int(q)} kPa")
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Ultimate failure probability, pf_U(t)")
    ax.set_title("Failure-probability trajectories across service-pressure levels")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "figure_service_pressure_pf_comparison.png")
    plt.close(fig)

    readme = """# Lower-demand service benchmark

This script keeps the original high-demand stress test and adds a lower-demand
service comparison for the main text. The same deterministic capacity history,
Monte Carlo mean capacity, capacity COV, resistance uncertainty, load COV, and
serviceability proxy are used. Only the service-pressure mean is changed.

The 500 kPa case is a practitioner-oriented lower-demand illustration, not a
universal allowable pressure for shallow foundations on sand.
"""
    (SUPP / "README_lower_demand_service_benchmark.md").write_text(readme, encoding="utf-8")
    print(DATA / "lower_demand_service_results.csv")
    print(DATA / "service_pressure_comparison_results.csv")
    print(FIGS / "figure_service_pressure_beta_comparison.png")
    print(FIGS / "figure_service_pressure_pf_comparison.png")


if __name__ == "__main__":
    main()
