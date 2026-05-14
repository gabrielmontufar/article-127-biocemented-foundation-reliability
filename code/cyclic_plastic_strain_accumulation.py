from __future__ import annotations

from math import erf, sqrt
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cyclic_loading_scenarios import build_loading_blocks, scenario_metadata


ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT
DATA = SUPP / "data"
FIGS = SUPP / "figures"
CODE = SUPP / "code"


PARAMS = {
    "Acyc": 0.0060,
    "CDR0": 0.08,
    "mcyc": 2.15,
    "Nref": 10_000.0,
    "rcyc": 0.42,
    "eps_p_max": 0.090,
    "ap": 8.0,
    "bp": 5.5,
    "eta0": 0.92,
    "etar": 0.25,
    "lambda_e": 0.028,
    "cb0": 105.0,
    "dtanphi": 0.12,
    "kH": 0.70,
    "chimin": 0.55,
    "rs": 0.22,
    "q_service": 950.0,
}


def norm_cdf(x):
    return 0.5 * (1.0 + np.vectorize(erf)(np.asarray(x) / sqrt(2.0)))


def parameter_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            ["Acyc", "Plastic-strain accumulation coefficient", PARAMS["Acyc"], "0.002-0.015", "Cyclic triaxial, cyclic direct shear, repeated plate load", "Confounded with ap if only final capacity loss is observed"],
            ["CDR0", "Cyclic demand threshold", PARAMS["CDR0"], "0.03-0.15", "Low-amplitude cyclic tests", "Hard to identify without below-threshold tests"],
            ["mcyc", "Amplitude sensitivity exponent", PARAMS["mcyc"], "1.0-4.0", "Multiple cyclic amplitudes", "Confounded with CDR0 if amplitude range is narrow"],
            ["rcyc", "Cycle-count accumulation exponent", PARAMS["rcyc"], "0.20-0.70", "Multiple cycle counts at fixed amplitude", "Confounded with Acyc if only one cycle count is used"],
            ["ap", "Cohesion damage exponent", PARAMS["ap"], "3-14", "Post-cyclic cohesion/UCS or triaxial strength loss", "Confounded with Acyc and eta without independent strain data"],
            ["bp", "Friction-enhancement damage exponent", PARAMS["bp"], "2-10", "Post-cyclic drained shear response", "Difficult to separate from ap using capacity alone"],
            ["eps_p_max", "Bounded accumulated plastic-strain cap", PARAMS["eps_p_max"], "0.03-0.15", "Large-strain cyclic or post-peak tests", "Acts as a numerical/physical limiter, not a universal constant"],
        ],
        columns=["Parameter", "Meaning", "Baseline value", "Sensitivity range", "Calibration data", "Identifiability issue"],
    )


def eta_at(t):
    return PARAMS["etar"] + (PARAMS["eta0"] - PARAMS["etar"]) * np.exp(-PARAMS["lambda_e"] * t)


def softened_eps(raw):
    emax = PARAMS["eps_p_max"]
    return emax * (1.0 - np.exp(-raw / emax))


def capacity_from_state(t, eps_p):
    eta = eta_at(t)
    chi_h = 1.0 - np.exp(-1.0 / PARAMS["kH"])
    chi_s = np.maximum(PARAMS["chimin"], 1.0 - PARAMS["rs"] * eta**2)
    dcp = np.exp(-PARAMS["ap"] * eps_p)
    dphi = np.exp(-PARAMS["bp"] * eps_p)
    cb_eff = PARAMS["cb0"] * eta**1.15 * dcp
    dtan_eff = PARAMS["dtanphi"] * eta**0.85 * dphi
    qu = 760.0 + 2.10 * cb_eff * chi_h * chi_s + 880.0 * dtan_eff * chi_h * chi_s
    return eta, cb_eff, dtan_eff, qu


def reliability(qu, rng_seed=127):
    q = PARAMS["q_service"]
    sigma_r = 0.13
    sigma_q = 0.08
    mu_r = -0.5 * sigma_r**2
    mu_q = -0.5 * sigma_q**2
    zfail = (np.log(q / qu) + mu_q - mu_r) / np.sqrt(sigma_r**2 + sigma_q**2)
    pf = norm_cdf(zfail)
    beta = -zfail
    fs = qu / q
    rng = np.random.default_rng(rng_seed)
    n = 120_000
    theta_r = rng.lognormal(mean=mu_r, sigma=sigma_r, size=n)
    theta_q = rng.lognormal(mean=mu_q, sigma=sigma_q, size=n)
    fail = (qu[:, None] * theta_r[None, :]) <= (q * theta_q[None, :])
    pf_life = float(np.any(fail, axis=0).mean())
    return fs, pf, beta, pf_life


def run_cyclic_model(hist: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    summary = []
    for sid, grp in hist.groupby("scenario"):
        grp = grp.sort_values("time_year").reset_index(drop=True)
        t = grp["time_year"].to_numpy()
        eps_raw = np.zeros_like(t)
        eps_eff = np.zeros_like(t)
        for i in range(1, len(t)):
            dt = max(t[i] - t[i - 1], 1e-9)
            cdr = grp["CDR"].iloc[i]
            cycles = grp["cycles"].iloc[i] * dt
            renv = grp["R_env"].iloc[i]
            dep = PARAMS["Acyc"] * max(0.0, cdr - PARAMS["CDR0"]) ** PARAMS["mcyc"]
            dep *= max(cycles / PARAMS["Nref"], 0.0) ** PARAMS["rcyc"] * renv
            eps_raw[i] = eps_raw[i - 1] + dep
            eps_eff[i] = softened_eps(eps_raw[i])
        eta, cb_eff, dtan_eff, qu = capacity_from_state(t, eps_eff)
        fs, pf, beta, pf_life = reliability(qu, rng_seed=127 + sum(ord(c) for c in sid))
        no_eps_qu = capacity_from_state(t, np.zeros_like(t))[3]
        out = grp.copy()
        out["epsilon_p_raw"] = eps_raw
        out["epsilon_p_eff"] = eps_eff
        out["eta_eff"] = eta
        out["cb_eff"] = cb_eff
        out["Delta_tan_phi_eff"] = dtan_eff
        out["Dcp"] = np.exp(-PARAMS["ap"] * eps_eff)
        out["Dphip"] = np.exp(-PARAMS["bp"] * eps_eff)
        out["qu"] = qu
        out["qu_no_cyclic_damage"] = no_eps_qu
        out["FS"] = fs
        out["pf"] = pf
        out["beta"] = beta
        out["pf_life"] = pf_life
        rows.append(out)
        final = out.iloc[-1]
        base_final = capacity_from_state(np.array([t[-1]]), np.array([0.0]))
        cb_base = base_final[1][0]
        dtan_base = base_final[2][0]
        summary.append(
            {
                "Scenario": sid,
                "epsilon_p_eff_50": final["epsilon_p_eff"],
                "cb_reduction_percent": 100.0 * (1.0 - final["cb_eff"] / cb_base),
                "Delta_tan_phi_reduction_percent": 100.0 * (1.0 - final["Delta_tan_phi_eff"] / dtan_base),
                "qu50": final["qu"],
                "FS50": final["FS"],
                "pf50": final["pf"],
                "beta50": final["beta"],
                "pf_life": final["pf_life"],
                "Interpretation": interpretation(sid),
            }
        )
    return pd.concat(rows, ignore_index=True), pd.DataFrame(summary)


def interpretation(sid: str) -> str:
    return {
        "S0": "Static case; epsilon_p remains zero and the base benchmark is recovered.",
        "C1": "Low-amplitude cycling remains near the threshold and causes minor damage.",
        "C2": "Moderate repeated loading accumulates gradual mechanical degradation.",
        "C3": "High-amplitude intermittent blocks create event-driven damage.",
        "C4": "Variable blocks demonstrate path-dependent accumulation.",
        "C5": "Environmental amplification increases cyclic plastic-strain damage.",
    }[sid]


def sensitivity_table(hist: pd.DataFrame) -> pd.DataFrame:
    rows = []
    base = dict(PARAMS)
    for pname, factors in {"Acyc": [0.5, 1.0, 2.0], "CDR0": [0.75, 1.0, 1.25], "mcyc": [0.75, 1.0, 1.25], "rcyc": [0.75, 1.0, 1.25], "ap": [0.5, 1.0, 1.5], "bp": [0.5, 1.0, 1.5]}.items():
        original = PARAMS[pname]
        for fac in factors:
            PARAMS[pname] = original * fac
            _, s = run_cyclic_model(hist)
            c5 = s[s["Scenario"] == "C5"].iloc[0]
            rows.append({"parameter": pname, "factor": fac, "C5_epsilon_p_eff_50": c5["epsilon_p_eff_50"], "C5_qu50": c5["qu50"], "C5_beta50": c5["beta50"]})
        PARAMS[pname] = base[pname]
    return pd.DataFrame(rows)


def make_figures(hist: pd.DataFrame, ts: pd.DataFrame) -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(3, 1, figsize=(8.0, 7.2), sharex=True, constrained_layout=True)
    for sid, grp in hist.groupby("scenario"):
        axes[0].plot(grp["time_year"], grp["CDR"], label=sid)
        axes[1].plot(grp["time_year"], grp["cycles"])
        axes[2].plot(grp["time_year"], grp["qmax_over_qu"])
    axes[0].set_ylabel("CDR")
    axes[1].set_ylabel("cycles/block")
    axes[2].set_ylabel("qmax/qu")
    axes[2].set_xlabel("Time (yr)")
    axes[0].legend(ncol=6, fontsize=7)
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.savefig(FIGS / "figure_cyclic_loading_histories.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 4.8), constrained_layout=True)
    for sid, grp in ts.groupby("scenario"):
        ax.plot(grp["time_year"], grp["epsilon_p_eff"], label=sid, lw=2)
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("epsilon_p_eff(t)")
    ax.set_title("Accumulated plastic-strain damage")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.savefig(FIGS / "figure_plastic_strain_accumulation.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 4.8), constrained_layout=True)
    for sid, grp in ts.groupby("scenario"):
        ax.plot(grp["time_year"], grp["qu"], label=sid, lw=2)
    ax.plot(ts[ts["scenario"] == "S0"]["time_year"], ts[ts["scenario"] == "S0"]["qu_no_cyclic_damage"], "k--", lw=1.5, label="no cyclic damage")
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("qu(t) (kPa)")
    ax.set_title("Effect of cyclic damage on bearing capacity")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    fig.savefig(FIGS / "figure_cyclic_capacity_degradation.png", dpi=300)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.6, 4.8), constrained_layout=True)
    for sid, grp in ts.groupby("scenario"):
        ax.plot(grp["time_year"], grp["beta"], label=sid, lw=2)
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("beta_U(t)")
    ax.set_title("Reliability under repeated loading")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.savefig(FIGS / "figure_cyclic_beta_degradation.png", dpi=300)
    plt.close(fig)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    CODE.mkdir(parents=True, exist_ok=True)
    t = np.linspace(0, 50, 101)
    meta = scenario_metadata()
    hist = build_loading_blocks(t)
    ts, summary = run_cyclic_model(hist)
    sens = sensitivity_table(hist)
    meta.to_csv(DATA / "cyclic_loading_scenarios_table.csv", index=False)
    parameter_table().to_csv(DATA / "cyclic_accumulation_parameters.csv", index=False)
    hist.to_csv(DATA / "cyclic_loading_histories.csv", index=False)
    ts.to_csv(DATA / "cyclic_scenario_reliability.csv", index=False)
    summary.to_csv(DATA / "cyclic_damage_summary.csv", index=False)
    sens.to_csv(DATA / "cyclic_parameter_sensitivity.csv", index=False)
    make_figures(hist, ts)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
