from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm


ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT / "04 Supplemental data and code" / "Supplementary files"
DATA = SUPP / "data"
FIGS = SUPP / "figures"
CODE = SUPP / "code"


SEED = 1272030
N = 120000
YEARS = np.arange(0.0, 51.0, 1.0)
ETA_R = 0.22
LAMBDA = 0.035
Q_M = 961.0
Q_GAIN = 889.0
Q_SERVICE = 900.0
SI_THRESHOLD = 0.045


def lognormal_params(mean: float, cov: float) -> tuple[float, float]:
    sig = np.sqrt(np.log(1 + cov**2))
    mu = np.log(mean) - 0.5 * sig**2
    return mu, sig


def sample_base(rng: np.random.Generator):
    mu_r, sig_r = lognormal_params(1.0, 0.10)
    mu_q, sig_q = lognormal_params(1.0, 0.08)
    mu_lam, sig_lam = lognormal_params(LAMBDA, 0.35)
    theta_r = np.exp(mu_r + sig_r * rng.standard_normal(N))
    theta_q = np.exp(mu_q + sig_q * rng.standard_normal(N))
    lam = np.exp(mu_lam + sig_lam * rng.standard_normal(N))
    eta0 = np.clip(rng.normal(1.0, 0.08, N), 0.65, 1.15)
    return theta_r, theta_q, lam, eta0


def scenario_catalog() -> pd.DataFrame:
    rows = [
        ["M", "monotonic exponential clock", "stationary q_service", "none", True, "terminal equivalence check"],
        ["E1", "episodic wetting-drying exposure increments", "stationary q_service", "none", True, "clustered exposure increases outcrossing rate"],
        ["E2", "short chemical pulse exposure", "stationary q_service", "none", True, "localized high-degradation pulse"],
        ["E3", "monotonic degradation", "seasonal and temporary load excursions", "none", False, "load variability can cause temporary first crossings"],
        ["E4", "monotonic degradation with partial recovery", "stationary q_service", "years 20 and 35", False, "inspection and retreatment may restore margin after earlier crossings"],
        ["E5", "combined pulses and partial recovery", "seasonal and event load excursions", "year 30", False, "combined non-monotonic exposure and loading"],
    ]
    return pd.DataFrame(rows, columns=["Scenario", "Degradation history", "Loading history", "Retreatment", "Monotonic?", "Purpose"])


def exposure_increment(scenario: str, years: np.ndarray) -> np.ndarray:
    dE = np.ones((len(years) - 1,))
    mids = 0.5 * (years[:-1] + years[1:])
    if scenario == "E1":
        dE += np.where(((mids >= 8) & (mids <= 12)) | ((mids >= 24) & (mids <= 28)) | ((mids >= 40) & (mids <= 44)), 2.5, 0.0)
    elif scenario == "E2":
        dE += np.where(((mids >= 18) & (mids <= 20)) | ((mids >= 37) & (mids <= 38)), 5.0, 0.0)
    elif scenario == "E5":
        dE += np.where(((mids >= 10) & (mids <= 14)) | ((mids >= 32) & (mids <= 34)), 3.0, 0.0)
    return dE


def q_service_history(scenario: str, years: np.ndarray) -> np.ndarray:
    q = np.full_like(years, Q_SERVICE, dtype=float)
    if scenario in {"E3", "E5"}:
        q += 70.0 * np.sin(2 * np.pi * years / 8.0)
        q += np.where(((years >= 15) & (years <= 18)) | ((years >= 34) & (years <= 36)), 160.0, 0.0)
    return q


def evolve_eta(scenario: str, years: np.ndarray, lam: np.ndarray, eta0: np.ndarray) -> np.ndarray:
    eta = np.zeros((len(eta0), len(years)))
    eta[:, 0] = eta0
    dE = exposure_increment(scenario, years)
    for j in range(1, len(years)):
        eta[:, j] = ETA_R + (eta[:, j - 1] - ETA_R) * np.exp(-lam * dE[j - 1])
        if scenario == "E4" and years[j] in [20, 35]:
            eta[:, j] = np.minimum(eta0, eta[:, j] + 0.45 * (eta0 - eta[:, j]))
        if scenario == "E5" and years[j] == 30:
            eta[:, j] = np.minimum(eta0, eta[:, j] + 0.35 * (eta0 - eta[:, j]))
    return eta


def compute_paths(scenario: str):
    rng = np.random.default_rng(SEED + sum(ord(c) for c in scenario))
    theta_r, theta_q, lam, eta0 = sample_base(rng)
    eta = evolve_eta(scenario, YEARS, lam, eta0)
    q_u = Q_M + Q_GAIN * np.clip(eta, 0.0, 1.2) ** 1.15
    qserv = q_service_history(scenario, YEARS)
    g = theta_r[:, None] * q_u - theta_q[:, None] * qserv[None, :]
    si = 0.18 * qserv[None, :] / (q_u * (1 + 4.2 * eta))
    return g, eta, qserv, si


def analyze_scenario(scenario: str):
    g, eta, qserv, si = compute_paths(scenario)
    fail_grid = g <= 0.0
    terminal = fail_grid[:, -1]
    life = fail_grid.any(axis=1)
    first_idx = np.where(life, fail_grid.argmax(axis=1), -1)
    first_time = np.where(life, YEARS[first_idx], np.nan)
    beta_t = -norm.ppf(np.clip(terminal.mean(), 0.5 / N, 1 - 0.5 / N))
    beta_l = -norm.ppf(np.clip(life.mean(), 0.5 / N, 1 - 0.5 / N))

    hazards = []
    survival_mask = np.ones(N, dtype=bool)
    for j in range(1, len(YEARS)):
        down = survival_mask & (g[:, j - 1] > 0.0) & (g[:, j] <= 0.0)
        n_j = int(survival_mask.sum())
        d_j = int(down.sum())
        h = d_j / n_j if n_j else 0.0
        hazards.append(
            {
                "Scenario": scenario,
                "year": YEARS[j],
                "n_survivors": n_j,
                "d_new_crossings": d_j,
                "hazard": h,
                "nu_minus_per_year": h / (YEARS[j] - YEARS[j - 1]),
            }
        )
        survival_mask &= ~down
    hazard_df = pd.DataFrame(hazards)
    s_disc = float(np.prod(1 - hazard_df["hazard"].to_numpy()))
    pf_hazard = 1.0 - s_disc
    beta_h = -norm.ppf(np.clip(pf_hazard, 0.5 / N, 1 - 0.5 / N))
    return {
        "Scenario": scenario,
        "pf_terminal": float(terminal.mean()),
        "beta_terminal": float(beta_t),
        "pf_life_first_passage": float(life.mean()),
        "beta_life": float(beta_l),
        "pf_life_hazard": float(pf_hazard),
        "beta_hazard": float(beta_h),
        "max_hazard_year": float(hazard_df.loc[hazard_df["hazard"].idxmax(), "year"]),
        "max_hazard": float(hazard_df["hazard"].max()),
        "difference_vs_terminal": float(life.mean() - terminal.mean()),
        "p_serviceability_index_gt_0p045": float((si > SI_THRESHOLD).any(axis=1).mean()),
    }, hazard_df, first_time, g


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    CODE.mkdir(parents=True, exist_ok=True)

    scenario_catalog().to_csv(DATA / "outcrossing_environment_scenarios.csv", index=False)
    rows, hazards, first_rows, sample_paths = [], [], [], {}
    for scenario in ["M", "E1", "E2", "E3", "E4", "E5"]:
        row, hdf, first_time, g = analyze_scenario(scenario)
        rows.append(row)
        hazards.append(hdf)
        idx = np.where(~np.isnan(first_time))[0][:200]
        first_rows.extend({"Scenario": scenario, "first_passage_time_year": float(t)} for t in first_time[idx])
        sample_paths[scenario] = g[:30, :]
    summary = pd.DataFrame(rows)
    summary.to_csv(DATA / "outcrossing_hazard_summary.csv", index=False)
    hazard_df = pd.concat(hazards, ignore_index=True)
    hazard_df.to_csv(DATA / "outcrossing_hazard_time_series.csv", index=False)
    pd.DataFrame(first_rows).to_csv(DATA / "first_passage_times.csv", index=False)

    m = summary[summary["Scenario"] == "M"].iloc[0]
    fig, ax = plt.subplots(figsize=(6.2, 4.0), dpi=180)
    ax.bar(["terminal pf(T)", "life first passage", "hazard survival"], [m.pf_terminal, m.pf_life_first_passage, m.pf_life_hazard], color=["#1f5a7a", "#4d7f36", "#b43b3b"])
    ax.set_ylabel("Failure probability")
    ax.set_title("Monotonic terminal equivalence check")
    fig.tight_layout()
    fig.savefig(FIGS / "figure_monotonic_terminal_equivalence.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=180)
    for k in range(10):
        ax.plot(YEARS, sample_paths["E5"][k], lw=0.9, alpha=0.8)
    ax.axhline(0, color="black", ls="--", lw=1)
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Limit-state margin g_U(t)")
    ax.set_title("Non-monotonic sample reliability paths")
    fig.tight_layout()
    fig.savefig(FIGS / "figure_nonmonotonic_sample_paths.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=180)
    for scenario in ["M", "E1", "E2", "E3", "E4", "E5"]:
        h = hazard_df[hazard_df["Scenario"] == scenario]
        ax.plot(h["year"], h["hazard"], label=scenario)
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Discrete hazard h_j")
    ax.set_title("Hazard / outcrossing rate over service life")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, ncol=3)
    fig.tight_layout()
    fig.savefig(FIGS / "figure_hazard_rate_over_time.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=180)
    x = np.arange(len(summary))
    w = 0.25
    ax.bar(x - w, summary["pf_terminal"], width=w, label="terminal pf(T)")
    ax.bar(x, summary["pf_life_first_passage"], width=w, label="first passage")
    ax.bar(x + w, summary["pf_life_hazard"], width=w, label="hazard survival")
    ax.set_xticks(x, summary["Scenario"])
    ax.set_ylabel("Failure probability")
    ax.set_title("Terminal versus lifetime probability")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGS / "figure_lifetime_pf_terminal_vs_outcrossing.png", bbox_inches="tight")
    plt.close(fig)

    (CODE / "lifetime_outcrossing_hazard.py").write_text(Path(__file__).read_text(encoding="utf-8"), encoding="utf-8")
    (CODE / "nonmonotonic_environment_scenarios.py").write_text(Path(__file__).read_text(encoding="utf-8"), encoding="utf-8")
    print(DATA / "outcrossing_hazard_summary.csv")
    print(FIGS / "figure_lifetime_pf_terminal_vs_outcrossing.png")


if __name__ == "__main__":
    main()
