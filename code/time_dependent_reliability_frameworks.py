from __future__ import annotations

from pathlib import Path
from statistics import NormalDist

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT
DATA = SUPP / "data"
FIGS = SUPP / "figures"
NORM = NormalDist()
N = 50_000
YEARS = np.arange(0.0, 51.0, 1.0)
SEED = 127


def lognormal_params(mean: float, cov: float) -> tuple[float, float]:
    sigma = np.sqrt(np.log1p(cov**2))
    mu = np.log(mean) - 0.5 * sigma**2
    return float(mu), float(sigma)


def beta_from_pf(pf: float) -> float:
    pf = min(max(float(pf), 0.5 / N), 1.0 - 0.5 / N)
    return -NORM.inv_cdf(pf)


def qu_from_eta(eta: np.ndarray) -> np.ndarray:
    q_matrix = 961.0457
    gain0 = 888.6384
    return q_matrix + gain0 * (0.62 * eta**1.15 + 0.38 * eta**0.85) * (1.0 - 0.18 * eta**2)


def sample_common(rng: np.random.Generator):
    mu_r, sig_r = lognormal_params(1.0, 0.10)
    mu_q, sig_q = lognormal_params(1.0, 0.08)
    mu_lam, sig_lam = lognormal_params(0.035, 0.35)
    theta_r = np.exp(mu_r + sig_r * rng.standard_normal(N))
    theta_q = np.exp(mu_q + sig_q * rng.standard_normal(N))
    lam = np.exp(mu_lam + sig_lam * rng.standard_normal(N))
    eta0 = np.clip(rng.normal(1.0, 0.08, N), 0.65, 1.15)
    return theta_r, theta_q, lam, eta0


def scenario_paths(scenario: str, rng: np.random.Generator):
    theta_r, theta_q, lam, eta0 = sample_common(rng)
    eta_r = 0.22
    eta = np.zeros((N, len(YEARS)))
    q = np.zeros((N, len(YEARS)))
    states = np.full((min(N, 300), len(YEARS)), "mild", dtype=object)

    if scenario == "T":
        for j, t in enumerate(YEARS):
            eta[:, j] = eta_r + (eta0 - eta_r) * np.exp(-lam * t)
            q[:, j] = 900.0
    elif scenario == "O":
        for j, t in enumerate(YEARS):
            pulse = 0.055 * np.exp(-0.5 * ((t - 18.0) / 3.0) ** 2) + 0.030 * np.exp(-0.5 * ((t - 34.0) / 4.0) ** 2)
            recovery = 0.22 if 28.0 <= t <= 36.0 else 0.0
            eta[:, j] = eta_r + (eta0 - eta_r) * np.exp(-(lam + pulse) * t)
            if recovery:
                eta[:, j] = np.minimum(eta0, eta[:, j] + recovery * (eta0 - eta[:, j]))
            q[:, j] = 830.0 + 170.0 * np.exp(-0.5 * ((t - 18.0) / 2.0) ** 2) + 110.0 * np.exp(-0.5 * ((t - 41.0) / 2.5) ** 2)
    elif scenario == "H":
        for j, t in enumerate(YEARS):
            event = 0.065 if 12 <= t <= 16 or 31 <= t <= 33 else 0.0
            eta[:, j] = eta_r + (eta0 - eta_r) * np.exp(-(lam + event) * t)
            q[:, j] = 850.0 + 180.0 * (12 <= t <= 16) + 120.0 * (31 <= t <= 33)
    elif scenario == "M":
        names = np.array(["mild", "wetting_drying", "chemical_pulse", "flushing", "mitigated"], dtype=object)
        pmat = np.array(
            [
                [0.72, 0.14, 0.03, 0.06, 0.05],
                [0.34, 0.42, 0.08, 0.10, 0.06],
                [0.46, 0.10, 0.22, 0.12, 0.10],
                [0.36, 0.14, 0.08, 0.34, 0.08],
                [0.54, 0.08, 0.02, 0.04, 0.32],
            ]
        )
        severity = np.array([0.000, 0.020, 0.065, 0.035, -0.010])
        state = np.zeros(N, dtype=int)
        eta[:, 0] = eta0
        q[:, 0] = 860.0
        states[:, 0] = names[state[: len(states)]]
        for j in range(1, len(YEARS)):
            u = rng.random(N)
            cs = np.cumsum(pmat[state], axis=1)
            state = (u[:, None] > cs).sum(axis=1)
            dt_lam = np.maximum(lam + severity[state], -0.015)
            eta[:, j] = eta_r + (eta[:, j - 1] - eta_r) * np.exp(-dt_lam)
            eta[:, j] = np.minimum(eta0, np.maximum(eta_r, eta[:, j]))
            q[:, j] = 850.0 + 140.0 * (state == 2) + 90.0 * (state == 3)
            states[:, j] = names[state[: len(states)]]
    elif scenario == "S":
        y = np.zeros(N)
        eta[:, 0] = eta0
        q[:, 0] = 860.0
        for j in range(1, len(YEARS)):
            y = 0.82 * y + 0.32 * rng.standard_normal(N)
            lam_t = lam * np.exp(0.45 * y)
            seasonal = 0.07 * np.sin(2.0 * np.pi * YEARS[j] / 8.0)
            events = (rng.random(N) < 0.035) * rng.lognormal(mean=np.log(0.18), sigma=0.35, size=N)
            q[:, j] = 850.0 * (1.0 + seasonal + events)
            eta[:, j] = eta_r + (eta[:, j - 1] - eta_r) * np.exp(-lam_t)
            eta[:, j] = np.minimum(eta0, np.maximum(eta_r, eta[:, j]))
    else:
        raise ValueError(scenario)

    qu = qu_from_eta(eta)
    g = theta_r[:, None] * qu - theta_q[:, None] * q
    return eta, q, g, states


def summarize_scenario(scenario: str, eta: np.ndarray, q: np.ndarray, g: np.ndarray, states: np.ndarray):
    fail = g <= 0.0
    terminal = fail[:, -1]
    ever = fail.any(axis=1)
    first_idx = np.argmax(fail, axis=1)
    first_time = np.where(ever, YEARS[first_idx], np.nan)
    pf_t = float(np.mean(terminal))
    pf_life = float(np.mean(ever))

    survived = np.ones(N, dtype=bool)
    hazard_rows = []
    for j, t in enumerate(YEARS):
        at_risk = int(survived.sum())
        new_fail = fail[:, j] & survived
        d_j = int(new_fail.sum())
        h_j = d_j / at_risk if at_risk else 0.0
        survived[new_fail] = False
        hazard_rows.append(
            {
                "scenario": scenario,
                "time_year": t,
                "new_first_crossings": d_j,
                "survivors_at_start": at_risk,
                "hazard_hj": h_j,
                "survival_S": float(survived.mean()),
                "pf_life": float(1.0 - survived.mean()),
                "pf_terminal_at_t": float(np.mean(fail[:, j])),
            }
        )
    haz = pd.DataFrame(hazard_rows)
    max_row = haz.loc[haz["hazard_hj"].idxmax()]
    difference = pf_life - pf_t
    markov_states = "mild; wetting_drying; chemical_pulse; flushing; mitigated" if scenario == "M" else "not used"
    return {
        "Scenario": scenario,
        "pf_T": pf_t,
        "beta_T": beta_from_pf(pf_t),
        "pf_life_first_passage": pf_life,
        "beta_life": beta_from_pf(pf_life),
        "max_hazard_year": float(max_row["time_year"]),
        "Markov states used": markov_states,
        "Difference vs terminal": difference,
        "Interpretation": {
            "T": "Terminal and first-passage estimates coincide for monotonic degradation and stationary loading.",
            "O": "Temporary exposure/load pulses create lifetime crossings not visible from the final year alone.",
            "H": "Discrete interval hazards identify the years where new first crossings concentrate.",
            "M": "State transitions translate environmental histories into probabilistic reliability paths.",
            "S": "Autocorrelated degradation and load fluctuations require path-wise lifetime reliability.",
        }[scenario],
    }, haz, first_time


def make_figures(all_results: pd.DataFrame, hazard: pd.DataFrame, path_data: dict):
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=180)
    for scenario, label, style in [("T", "terminal monotonic", "-"), ("O", "non-monotonic outcrossing", "--"), ("M", "Markov environment", "-."), ("S", "stochastic process", ":")]:
        sub = hazard[hazard["scenario"] == scenario]
        ax.plot(sub["time_year"], sub["pf_terminal_at_t"], style, lw=2.0, label=f"{label}: pf(t)")
        if scenario in ["T", "O"]:
            ax.plot(sub["time_year"], sub["pf_life"], style, lw=1.2, alpha=0.65, label=f"{label}: pf_life(t)")
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Failure probability")
    ax.set_title("Terminal versus first-passage reliability")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(FIGS / "figure_terminal_vs_first_passage_reliability.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=180)
    for scenario, color in [("O", "#1f77b4"), ("M", "#2ca02c"), ("S", "#d62728")]:
        g = path_data[scenario]["g"]
        for idx in np.linspace(0, g.shape[0] - 1, 12, dtype=int):
            ax.plot(YEARS, g[idx] / 1000.0, lw=0.8, alpha=0.35, color=color)
    ax.axhline(0, color="#333333", lw=1.0)
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Limit-state margin g_U(t) (MPa-equivalent kPa/1000)")
    ax.set_title("Sample stochastic reliability paths")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(FIGS / "figure_sample_stochastic_reliability_paths.png")
    plt.close(fig)

    sub = hazard[hazard["scenario"] == "H"]
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 6.4), dpi=180, sharex=True)
    axes[0].plot(sub["time_year"], sub["hazard_hj"], color="#b23a48", lw=2.0)
    axes[0].set_ylabel("h_j")
    axes[1].plot(sub["time_year"], sub["survival_S"], color="#1f77b4", lw=2.0)
    axes[1].set_ylabel("S(t)")
    axes[2].plot(sub["time_year"], sub["pf_life"], color="#2ca02c", lw=2.0)
    axes[2].set_ylabel("pf_life(t)")
    axes[2].set_xlabel("Time (yr)")
    for ax in axes:
        ax.grid(True, alpha=0.25)
    axes[0].set_title("Discrete hazard and survival curves")
    fig.tight_layout()
    fig.savefig(FIGS / "figure_discrete_hazard_survival_curves.png")
    plt.close(fig)

    states = path_data["M"]["states"][:40]
    state_map = {"mild": 0, "wetting_drying": 1, "chemical_pulse": 2, "flushing": 3, "mitigated": 4}
    z = np.vectorize(state_map.get)(states)
    sub = hazard[hazard["scenario"] == "M"]
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.2), dpi=180, sharex=True)
    im = axes[0].imshow(z, aspect="auto", interpolation="nearest", cmap="viridis", extent=[YEARS[0], YEARS[-1], z.shape[0], 0])
    axes[0].set_ylabel("sample path")
    axes[0].set_title("Markov environmental state paths and reliability response")
    cbar = fig.colorbar(im, ax=axes[0], ticks=list(state_map.values()))
    cbar.ax.set_yticklabels(list(state_map.keys()), fontsize=7)
    axes[1].plot(sub["time_year"], sub["pf_life"], lw=2.0, color="#b23a48", label="pf_life")
    axes[1].plot(sub["time_year"], sub["pf_terminal_at_t"], lw=1.6, color="#1f77b4", label="pf(t)")
    axes[1].set_xlabel("Time (yr)")
    axes[1].set_ylabel("Probability")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGS / "figure_markov_environment_reliability_response.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=180)
    ax.axis("off")
    boxes = [
        (0.04, 0.64, "Monotonic degradation\nstationary loading", "Terminal\npf(T)"),
        (0.30, 0.64, "Non-monotonic loading\nor exposure", "First passage\nmin g(t)"),
        (0.56, 0.64, "Simulated paths\nor inspection intervals", "Discrete hazard\nh_j, S(t)"),
        (0.17, 0.20, "Discrete exposure states", "Markov\nstate transitions"),
        (0.47, 0.20, "Autocorrelated loads\nor degradation", "Stochastic process\npaths of g(t)"),
    ]
    for x, y, trigger, method in boxes:
        rect = plt.Rectangle((x, y), 0.22, 0.20, ec="#333333", fc="#edf4fb", lw=1.1)
        ax.add_patch(rect)
        ax.text(x + 0.11, y + 0.13, trigger, ha="center", va="center", fontsize=7.5)
        ax.text(x + 0.11, y + 0.045, method, ha="center", va="center", fontsize=8.0, weight="bold")
    for start, end in [((0.26, 0.74), (0.30, 0.74)), ((0.52, 0.74), (0.56, 0.74)), ((0.41, 0.64), (0.27, 0.40)), ((0.67, 0.64), (0.58, 0.40))]:
        ax.annotate("", xy=end, xytext=start, arrowprops=dict(arrowstyle="->", color="#333333", lw=1.0))
    ax.set_title("Selecting a time-dependent reliability framework", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIGS / "figure_tdr_framework_decision_diagram.png")
    plt.close(fig)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    scenario_defs = pd.DataFrame(
        [
            ["T", "monotonic exponential degradation", "stationary service demand", "terminal and first-passage equivalence", "synthetic methodological benchmark", "verify reduced monotonic benchmark"],
            ["O", "exposure pulses with partial recovery", "temporary service-load pulses", "first-passage/outcrossing", "synthetic methodological benchmark", "show pf_life can exceed pf(T)"],
            ["H", "event blocks", "event blocks", "discrete hazard", "synthetic methodological benchmark", "estimate h_j, S(t), and cumulative failure"],
            ["M", "Markov states: mild, wetting_drying, chemical_pulse, flushing, mitigated", "state-conditioned demand", "Markov environmental reliability", "synthetic methodological benchmark", "demonstrate state-transition exposure modeling"],
            ["S", "autocorrelated degradation multiplier", "seasonal and random event demand", "stochastic-process reliability", "synthetic methodological benchmark", "simulate random time-correlated margin paths"],
        ],
        columns=["Scenario", "Environmental history", "Loading history", "Reliability framework", "Synthetic or data-based", "Purpose"],
    )
    scenario_defs.to_csv(DATA / "tdr_scenario_definitions.csv", index=False)

    summaries = []
    hazard_tables = []
    first_times = []
    path_data = {}
    for scenario in ["T", "O", "H", "M", "S"]:
        rng = np.random.default_rng(SEED + ord(scenario))
        eta, q, g, states = scenario_paths(scenario, rng)
        summary, haz, first_time = summarize_scenario(scenario, eta, q, g, states)
        summaries.append(summary)
        hazard_tables.append(haz)
        path_data[scenario] = {"eta": eta[:300], "q": q[:300], "g": g[:300], "states": states}
        first_times.append(pd.DataFrame({"scenario": scenario, "first_passage_time_year": first_time}))
        if scenario == "M":
            pd.DataFrame(states).to_csv(DATA / "markov_environment_paths.csv", index=False)
        if scenario == "S":
            pd.DataFrame(g[:300], columns=[f"year_{int(y)}" for y in YEARS]).to_csv(DATA / "stochastic_process_paths.csv", index=False)

    summary_df = pd.DataFrame(summaries)
    hazard_df = pd.concat(hazard_tables, ignore_index=True)
    first_df = pd.concat(first_times, ignore_index=True)
    summary_df.to_csv(DATA / "tdr_framework_comparison_summary.csv", index=False)
    hazard_df.to_csv(DATA / "hazard_survival_table.csv", index=False)
    first_df.to_csv(DATA / "first_passage_times.csv", index=False)
    make_figures(summary_df, hazard_df, path_data)

    readme = """# Time-dependent reliability framework supplement

This supplement places the screening model within terminal, first-passage,
outcrossing, discrete-hazard, Markov environmental, and stochastic-process
reliability formulations. The scenarios are synthetic methodological examples,
not site-calibrated environmental histories.
"""
    (SUPP / "README_time_dependent_reliability_frameworks.md").write_text(readme, encoding="utf-8")
    print(DATA / "tdr_framework_comparison_summary.csv")
    print(DATA / "hazard_survival_table.csv")


if __name__ == "__main__":
    main()
