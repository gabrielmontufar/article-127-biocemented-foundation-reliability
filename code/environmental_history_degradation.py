from __future__ import annotations

from math import erf, log, sqrt
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from environmental_scenarios_generator import build_histories, scenario_table


ROOT = Path(__file__).resolve().parents[3]
SUPP = ROOT / "04 Supplemental data and code" / "Supplementary files"
DATA = SUPP / "data"
FIGS = SUPP / "figures"
CODE = SUPP / "code"


def norm_cdf(x: np.ndarray | float) -> np.ndarray | float:
    return 0.5 * (1.0 + np.vectorize(erf)(np.asarray(x) / sqrt(2.0)))


def reliability_from_eta(eta: np.ndarray, rng_seed: int = 127) -> pd.DataFrame:
    q_matrix = 760.0
    cb0 = 105.0
    dtanphi = 0.12
    h_over_b = 1.0
    kH = 0.70
    chimin = 0.55
    rs = 0.22
    q_service = 950.0
    cov_qu = 0.13
    chi_h = 1.0 - np.exp(-h_over_b / kH)
    chi_s = np.maximum(chimin, 1.0 - rs * eta**2)
    qu = q_matrix + 2.10 * cb0 * eta**1.15 * chi_h * chi_s + 880.0 * dtanphi * eta**0.85 * chi_h * chi_s
    sigma_r = cov_qu
    sigma_q = 0.08
    mu_r = -0.5 * sigma_r**2
    mu_q = -0.5 * sigma_q**2
    zfail = (np.log(q_service / qu) + mu_q - mu_r) / np.sqrt(sigma_r**2 + sigma_q**2)
    pf = norm_cdf(zfail)
    beta = -zfail
    fs = qu / q_service

    rng = np.random.default_rng(rng_seed)
    n = 120_000
    theta_r = rng.lognormal(mean=mu_r, sigma=sigma_r, size=n)
    theta_q = rng.lognormal(mean=mu_q, sigma=sigma_q, size=n)
    fail = (qu[:, None] * theta_r[None, :]) <= (q_service * theta_q[None, :])
    pf_life = float(np.any(fail, axis=0).mean())
    return pd.DataFrame({"qu": qu, "FS": fs, "beta": beta, "pf": pf, "pf_life": pf_life})


def run_environmental_model(hist: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    eta0 = 0.92
    etar = 0.25
    lambda_base = 0.028
    kappa_eta = lambda_base
    all_rows = []
    summary_rows = []
    for sid, grp in hist.groupby("scenario"):
        grp = grp.sort_values("time_year").reset_index(drop=True)
        t = grp["time_year"].to_numpy()
        senv = grp["S_env"].to_numpy()
        eta = np.empty_like(t)
        eenv = np.zeros_like(t)
        eta[0] = eta0
        for i in range(1, len(t)):
            dt = t[i] - t[i - 1]
            dE = 0.5 * (senv[i] + senv[i - 1]) * dt
            eenv[i] = eenv[i - 1] + dE
            eta[i] = etar + (eta[i - 1] - etar) * np.exp(-kappa_eta * dE)
        rel = reliability_from_eta(eta, rng_seed=127 + ord(sid))
        out = pd.concat([grp, pd.DataFrame({"E_env": eenv, "eta_eff": eta}), rel], axis=1)
        all_rows.append(out)
        eta_t = float(eta[-1])
        lambda_eff = -log((eta_t - etar) / (eta0 - etar)) / t[-1]
        summary_rows.append(
            {
                "Scenario": sid,
                "eta_eff_50": eta_t,
                "lambda_eff": lambda_eff,
                "qu50": float(rel["qu"].iloc[-1]),
                "FS50": float(rel["FS"].iloc[-1]),
                "pf50": float(rel["pf"].iloc[-1]),
                "beta50": float(rel["beta"].iloc[-1]),
                "pf_life": float(rel["pf_life"].iloc[-1]),
                "Interpretation": interpretation(sid),
            }
        )
    return pd.concat(all_rows, ignore_index=True), pd.DataFrame(summary_rows)


def interpretation(sid: str) -> str:
    return {
        "A": "Mild stationary exposure close to reduced base case.",
        "B": "Aggressive chemistry accelerates cementation loss.",
        "C": "Flushing pulses increase cumulative exposure during high-gradient periods.",
        "D": "Seasonal wetting-drying creates gradual cyclic exposure accumulation.",
        "E": "Combined exposure gives the largest degradation and reliability loss.",
        "F": "Mitigation reduces equivalent degradation rate and reliability loss.",
    }[sid]


def driver_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Driver": "Groundwater chemistry",
                "Physical meaning": "Aggressiveness of pore water toward carbonate bonds",
                "Possible measured inputs": "pH, calcite saturation index, alkalinity, ionic composition, salinity",
                "Screening severity function": "S_chem=max(0,-SI_calcite) or bounded pH proxy",
                "Calibration data": "Durability tests under controlled chemistry, UCS/UPV/Vs retention",
                "Limitation": "Does not replace PHREEQC or reactive-transport modeling",
            },
            {
                "Driver": "Hydraulic flushing",
                "Physical meaning": "Transport and renewal of aggressive pore water through treated sand",
                "Possible measured inputs": "Hydraulic gradient, seepage velocity, Darcy flux, pore volumes flushed",
                "Screening severity function": "S_flow=min(Smax,(i/i_ref)^a_i) or pore-volume exposure",
                "Calibration data": "Permeability/flushing tests and retention after known pore volumes",
                "Limitation": "Needs hydraulic boundary conditions and spatial flow paths",
            },
            {
                "Driver": "Wetting-drying",
                "Physical meaning": "Damage from saturation/suction cycling and chemical renewal",
                "Possible measured inputs": "Cycle count, Delta Sr, suction amplitude, exposure duration",
                "Screening severity function": "S_WD=f_WD(Delta Sr/Delta Sr_ref)^a_WD or event sum",
                "Calibration data": "Wetting-drying series with UCS, UPV, Vs, stiffness or mass retention",
                "Limitation": "Cycle amplitude and boundary drainage must be controlled",
            },
            {
                "Driver": "Combined exposure",
                "Physical meaning": "Aggregate environmental severity used by the reliability screen",
                "Possible measured inputs": "Chemistry, flow and wetting-drying records over time",
                "Screening severity function": "S_env=w_chem S_chem+w_flow S_flow+w_WD S_WD",
                "Calibration data": "Mechanism-isolated tests followed by independent combined validation",
                "Limitation": "Weights and kappa_eta are not separately identifiable from one endpoint",
            },
        ]
    )


def make_figures(hist: pd.DataFrame, time_series: pd.DataFrame, summary: pd.DataFrame) -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    labels = dict(zip(scenario_table()["Scenario"], scenario_table()["Name"]))

    fig, axes = plt.subplots(4, 1, figsize=(8.0, 8.8), sharex=True, constrained_layout=True)
    for sid, grp in hist.groupby("scenario"):
        for ax, col in zip(axes, ["S_chem", "S_flow", "S_WD", "S_env"]):
            ax.plot(grp["time_year"], grp[col], label=f"{sid}: {labels[sid]}", lw=1.6)
            ax.set_ylabel(col)
            ax.grid(alpha=0.25)
    axes[0].legend(ncol=2, fontsize=7)
    axes[-1].set_xlabel("Time (yr)")
    fig.savefig(FIGS / "figure_environmental_histories.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 4.8), constrained_layout=True)
    for sid, grp in time_series.groupby("scenario"):
        ax.plot(grp["time_year"], grp["eta_eff"], label=f"{sid}: {labels[sid]}", lw=2)
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("eta_eff(t)")
    ax.set_title("Cementation degradation under environmental histories")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7)
    fig.savefig(FIGS / "figure_eta_environmental_scenarios.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 4.8), constrained_layout=True)
    for sid, grp in time_series.groupby("scenario"):
        ax.plot(grp["time_year"], grp["beta"], label=f"{sid}: {labels[sid]}", lw=2)
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("beta_U(t)")
    ax.set_title("Reliability response to environmental histories")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7)
    fig.savefig(FIGS / "figure_beta_environmental_scenarios.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 4.5), constrained_layout=True)
    ax.bar(summary["Scenario"], summary["lambda_eff"], color="#607d8b")
    ax.axhline(0.028, color="crimson", ls="--", lw=1.5, label="Base lambda_e")
    ax.set_xlabel("Scenario")
    ax.set_ylabel("Equivalent lambda_eff (1/yr)")
    ax.set_title("Equivalent degradation-rate comparison")
    ax.legend()
    fig.savefig(FIGS / "figure_lambda_eff_comparison.png", dpi=300)
    plt.close(fig)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    CODE.mkdir(parents=True, exist_ok=True)
    t = np.linspace(0, 50, 201)
    scen = scenario_table()
    hist = build_histories(t)
    time_series, summary = run_environmental_model(hist)
    drivers = driver_table()
    scen.to_csv(DATA / "environmental_scenarios_table.csv", index=False)
    hist.to_csv(DATA / "environmental_histories.csv", index=False)
    time_series.to_csv(DATA / "environmental_scenario_reliability.csv", index=False)
    summary.to_csv(DATA / "environmental_degradation_summary.csv", index=False)
    summary[["Scenario", "lambda_eff"]].to_csv(DATA / "environmental_lambda_equivalent.csv", index=False)
    drivers.to_csv(DATA / "environmental_driver_screening_proxies.csv", index=False)
    make_figures(hist, time_series, summary)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()


