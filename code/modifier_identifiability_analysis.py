from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from math import erf


ROOT = Path(__file__).resolve().parents[3]
SUPP = ROOT / "04 Supplemental data and code" / "Supplementary files"
DATA = SUPP / "data"
FIGS = SUPP / "figures"
CODE = SUPP / "code"

PARAMS = ["eta0", "etar", "lambda_e", "cb0", "dtanphi", "kH", "kE", "chimin", "rs"]
BOUNDS = {
    "eta0": (0.50, 1.00),
    "etar": (0.05, 0.70),
    "lambda_e": (0.002, 0.080),
    "cb0": (20.0, 220.0),
    "dtanphi": (0.01, 0.35),
    "kH": (0.20, 1.80),
    "kE": (0.20, 1.80),
    "chimin": (0.35, 0.95),
    "rs": (0.02, 0.70),
}


def nominal() -> dict[str, float]:
    return {
        "eta0": 0.92,
        "etar": 0.25,
        "lambda_e": 0.028,
        "cb0": 105.0,
        "dtanphi": 0.12,
        "kH": 0.70,
        "kE": 0.55,
        "chimin": 0.55,
        "rs": 0.22,
    }


def vec_to_par(x: np.ndarray, names: list[str]) -> dict[str, float]:
    p = nominal()
    p.update({n: float(v) for n, v in zip(names, x)})
    return p


def eta(p: dict[str, float], t: float) -> float:
    return p["etar"] + (p["eta0"] - p["etar"]) * np.exp(-p["lambda_e"] * t)


def chi_h(p: dict[str, float], h_over_b: float) -> float:
    return 1.0 - np.exp(-h_over_b / p["kH"])


def chi_e(p: dict[str, float], h_over_b: float) -> float:
    return 1.0 - np.exp(-h_over_b / p["kE"])


def chi_s(p: dict[str, float], et: float) -> float:
    return max(p["chimin"], 1.0 - p["rs"] * et**2)


def model_outputs(p: dict[str, float], h_over_b: float = 1.0, t: float = 50.0) -> dict[str, float]:
    et = eta(p, t)
    ch = chi_h(p, h_over_b)
    ce = chi_e(p, h_over_b)
    cs = chi_s(p, et)
    q_matrix = 760.0
    cohesion_gain = 2.10 * p["cb0"] * et**1.15 * ch * cs
    friction_gain = 880.0 * p["dtanphi"] * et**0.85 * ch * cs
    qu = q_matrix + cohesion_gain + friction_gain
    stiffness_ratio = 1.0 + 0.72 * ce * (p["cb0"] / 105.0) ** 0.35 * et**0.70
    service_index = 0.060 / stiffness_ratio
    demand = 950.0
    cov_qu = 0.13
    beta = (qu - demand) / (cov_qu * qu)
    pf = 0.5 * (1.0 + erf((-beta) / np.sqrt(2.0)))
    return {"qu50": qu, "sI50": service_index, "beta50": beta, "pf50": pf, "stiffness50": stiffness_ratio}


def rich_observations(p: dict[str, float]) -> tuple[list[dict], np.ndarray, np.ndarray]:
    obs = []
    for t in [0, 10, 25, 50]:
        obs.append({"type": "eta", "t": t, "h": 1.0, "sigma": 0.04})
    for t in [0, 50]:
        obs.append({"type": "element_strength", "t": t, "h": 1.0, "sigma": 12.0})
    for h in [0.5, 1.0, 1.5, 2.0]:
        obs.append({"type": "qu", "t": 50, "h": h, "sigma": 25.0})
    for h in [0.5, 1.0, 1.5, 2.0]:
        obs.append({"type": "stiffness", "t": 50, "h": h, "sigma": 0.08})
    for t in [0, 50]:
        obs.append({"type": "post_peak_ratio", "t": t, "h": 1.0, "sigma": 0.05})
    y = np.array([predict_observation(p, o) for o in obs])
    sigma = np.array([o["sigma"] for o in obs])
    return obs, y, sigma


def sparse_observations(p: dict[str, float]) -> tuple[list[dict], np.ndarray, np.ndarray]:
    obs = [
        {"type": "qu", "t": 50, "h": 1.0, "sigma": 25.0},
        {"type": "stiffness", "t": 50, "h": 1.0, "sigma": 0.08},
    ]
    y = np.array([predict_observation(p, o) for o in obs])
    sigma = np.array([o["sigma"] for o in obs])
    return obs, y, sigma


def predict_observation(p: dict[str, float], obs: dict) -> float:
    et = eta(p, obs["t"])
    if obs["type"] == "eta":
        return et
    if obs["type"] == "element_strength":
        return 70.0 + 1.85 * p["cb0"] * et**1.15 + 680.0 * p["dtanphi"] * et**0.85
    if obs["type"] == "qu":
        return model_outputs(p, obs["h"], obs["t"])["qu50"]
    if obs["type"] == "stiffness":
        return model_outputs(p, obs["h"], obs["t"])["stiffness50"]
    if obs["type"] == "post_peak_ratio":
        return chi_s(p, et)
    raise ValueError(obs["type"])


def residuals_for(names: list[str], x: np.ndarray, obs: list[dict], y: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    p = vec_to_par(x, names)
    pred = np.array([predict_observation(p, o) for o in obs])
    return (pred - y) / sigma


def finite_jacobian(names: list[str], x: np.ndarray, obs: list[dict], y: np.ndarray, sigma: np.ndarray) -> np.ndarray:
    base = residuals_for(names, x, obs, y, sigma)
    jac = np.zeros((len(base), len(names)))
    for j, name in enumerate(names):
        step = 1e-5 * max(abs(x[j]), 1.0)
        xp = x.copy()
        xm = x.copy()
        xp[j] = min(BOUNDS[name][1], x[j] + step)
        xm[j] = max(BOUNDS[name][0], x[j] - step)
        if xp[j] == xm[j]:
            continue
        jac[:, j] = (residuals_for(names, xp, obs, y, sigma) - residuals_for(names, xm, obs, y, sigma)) / (xp[j] - xm[j])
    return jac


def fit_dataset(obs: list[dict], y: np.ndarray, sigma: np.ndarray, names: list[str]) -> tuple[np.ndarray, np.ndarray]:
    x0 = np.array([nominal()[n] for n in names])
    # The diagnostic observations are generated from the nominal benchmark.
    # The fit therefore has a known optimum at x0; a finite-difference Jacobian
    # at that point is enough for local covariance and correlation diagnostics.
    return x0, finite_jacobian(names, x0, obs, y, sigma)


def covariance_from_jac(jac: np.ndarray) -> np.ndarray:
    jtj = jac.T @ jac
    return np.linalg.pinv(jtj)


def corr_from_cov(cov: np.ndarray) -> np.ndarray:
    d = np.sqrt(np.maximum(np.diag(cov), 1e-30))
    return cov / np.outer(d, d)


def local_sensitivity(p: dict[str, float]) -> pd.DataFrame:
    base = model_outputs(p)
    rows = []
    for name in PARAMS:
        step = 1e-4 * max(abs(p[name]), 1.0)
        p_hi = dict(p)
        p_lo = dict(p)
        p_hi[name] = min(BOUNDS[name][1], p[name] + step)
        p_lo[name] = max(BOUNDS[name][0], p[name] - step)
        hi = model_outputs(p_hi)
        lo = model_outputs(p_lo)
        for out in ["qu50", "sI50", "beta50", "pf50"]:
            if hi[out] > 0 and lo[out] > 0 and p_hi[name] > 0 and p_lo[name] > 0:
                sens = (np.log(hi[out]) - np.log(lo[out])) / (np.log(p_hi[name]) - np.log(p_lo[name]))
            else:
                sens = np.nan
            rows.append({"parameter": name, "output": out, "normalized_sensitivity": sens, "baseline_output": base[out]})
    return pd.DataFrame(rows)


def profile_objectives(obs: list[dict], y: np.ndarray, sigma: np.ndarray) -> pd.DataFrame:
    rows = []
    rng = np.random.default_rng(127)
    for fixed in ["kH", "kE", "rs", "chimin"]:
        lo, hi = BOUNDS[fixed]
        grid = np.linspace(lo, hi, 31)
        free = [n for n in PARAMS if n != fixed]
        free_lo = np.array([BOUNDS[n][0] for n in free])
        free_hi = np.array([BOUNDS[n][1] for n in free])
        nominal_free = np.array([nominal()[n] for n in free])
        unit_samples = rng.random((1800, len(free)))
        samples = free_lo + unit_samples * (free_hi - free_lo)
        samples = np.vstack([nominal_free, samples])
        for val in grid:
            best = np.inf
            # Random bounded reoptimization is sufficient here because this is
            # a screening identifiability diagnostic, not a project calibration.
            for x in samples:
                p = vec_to_par(x, free)
                p[fixed] = float(val)
                pred = np.array([predict_observation(p, o) for o in obs])
                obj = float(np.sum(((pred - y) / sigma) ** 2))
                if obj < best:
                    best = obj
            rows.append({"profile_parameter": fixed, "fixed_value": val, "objective": best})
    df = pd.DataFrame(rows)
    df["delta_objective"] = df.groupby("profile_parameter")["objective"].transform(lambda s: s - s.min())
    return df


def identifiability_summary(rich_corr: np.ndarray, sparse_corr: np.ndarray) -> pd.DataFrame:
    warnings = [
        ("eta0", "cb0", "Capacity-only data can trade initial cementation against cohesion increment."),
        ("etar", "lambda_e", "A single degradation endpoint cannot separate residual state from degradation rate."),
        ("cb0", "dtanphi", "Footing response alone may not separate cohesion-like and frictional increments."),
        ("kH", "rs", "A single treated depth can trade mobilized depth against localization cap."),
        ("kE", "cb0", "One settlement curve can trade stiffness-depth scale against strength/stiffness amplitude."),
    ]
    rows = []
    idx = {n: i for i, n in enumerate(PARAMS)}
    for a, b, msg in warnings:
        rows.append(
            {
                "parameter_pair": f"{a} - {b}",
                "rich_data_abs_corr": abs(rich_corr[idx[a], idx[b]]),
                "sparse_data_abs_corr": abs(sparse_corr[idx[a], idx[b]]),
                "warning": msg,
                "recommended_response": "Use staged calibration or report sensitivity envelope if independent data are unavailable.",
            }
        )
    return pd.DataFrame(rows)


def plot_heatmap(corr: np.ndarray, sens: pd.DataFrame) -> None:
    pivot = sens.pivot(index="parameter", columns="output", values="normalized_sensitivity").loc[PARAMS]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), constrained_layout=True)
    im0 = axes[0].imshow(corr, vmin=-1, vmax=1, cmap="coolwarm")
    axes[0].set_xticks(range(len(PARAMS)), PARAMS, rotation=45, ha="right")
    axes[0].set_yticks(range(len(PARAMS)), PARAMS)
    axes[0].set_title("Approximate parameter correlation")
    fig.colorbar(im0, ax=axes[0], fraction=0.046)
    im1 = axes[1].imshow(pivot.values, cmap="viridis", aspect="auto")
    axes[1].set_xticks(range(len(pivot.columns)), pivot.columns)
    axes[1].set_yticks(range(len(PARAMS)), PARAMS)
    axes[1].set_title("Normalized local sensitivity")
    fig.colorbar(im1, ax=axes[1], fraction=0.046)
    fig.savefig(FIGS / "figure_modifier_identifiability_heatmap.png", dpi=300)
    plt.close(fig)


def plot_profiles(df: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 5.0), constrained_layout=True)
    for name, grp in df.groupby("profile_parameter"):
        x = (grp["fixed_value"] - grp["fixed_value"].min()) / (grp["fixed_value"].max() - grp["fixed_value"].min())
        ax.plot(x, grp["delta_objective"], marker="o", ms=3, label=name)
    ax.axhline(3.84, color="0.4", ls="--", lw=1, label="Delta objective = 3.84")
    ax.set_xlabel("Normalized fixed-parameter grid")
    ax.set_ylabel("Profile objective increase")
    ax.set_title("Profile-objective diagnostics for modifier parameters")
    ax.legend()
    fig.savefig(FIGS / "figure_modifier_profile_objectives.png", dpi=300)
    plt.close(fig)


def plot_output_sensitivity(p: dict[str, float]) -> pd.DataFrame:
    rows = []
    grids = {
        "kH": np.linspace(0.25, 1.60, 60),
        "kE": np.linspace(0.25, 1.60, 60),
        "rs": np.linspace(0.02, 0.60, 60),
    }
    for name, grid in grids.items():
        for val in grid:
            pp = dict(p)
            pp[name] = float(val)
            out = model_outputs(pp)
            rows.append({"varied_parameter": name, "value": val, **out})
    df = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.2), constrained_layout=True)
    outputs = [("qu50", "qu50 (kPa)"), ("sI50", "service index"), ("beta50", "beta50")]
    for ax, (out, label) in zip(axes, outputs):
        for name, grp in df.groupby("varied_parameter"):
            ax.plot(grp["value"], grp[out], label=name)
        ax.set_xlabel("Parameter value")
        ax.set_ylabel(label)
        ax.grid(alpha=0.25)
    axes[0].legend()
    fig.suptitle("Effect of chiH, chiE and chis controls on benchmark outputs")
    fig.savefig(FIGS / "figure_modifier_output_sensitivity.png", dpi=300)
    plt.close(fig)
    return df


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    CODE.mkdir(parents=True, exist_ok=True)
    p0 = nominal()
    rich_obs, rich_y, rich_sigma = rich_observations(p0)
    sparse_obs, sparse_y, sparse_sigma = sparse_observations(p0)

    rich_x, rich_jac = fit_dataset(rich_obs, rich_y, rich_sigma, PARAMS)
    sparse_x, sparse_jac = fit_dataset(sparse_obs, sparse_y, sparse_sigma, PARAMS)
    rich_corr = corr_from_cov(covariance_from_jac(rich_jac))
    sparse_corr = corr_from_cov(covariance_from_jac(sparse_jac))

    sens = local_sensitivity(p0)
    profiles = profile_objectives(rich_obs, rich_y, rich_sigma)
    out_sens = plot_output_sensitivity(p0)
    summary = identifiability_summary(rich_corr, sparse_corr)

    sens.to_csv(DATA / "modifier_sensitivity_matrix.csv", index=False)
    profiles.to_csv(DATA / "modifier_profile_likelihood.csv", index=False)
    summary.to_csv(DATA / "modifier_identifiability_summary.csv", index=False)
    out_sens.to_csv(DATA / "modifier_output_sensitivity_curves.csv", index=False)
    pd.DataFrame(rich_corr, index=PARAMS, columns=PARAMS).to_csv(DATA / "modifier_parameter_correlation_rich_data.csv")
    pd.DataFrame(sparse_corr, index=PARAMS, columns=PARAMS).to_csv(DATA / "modifier_parameter_correlation_sparse_data.csv")
    pd.DataFrame({"parameter": PARAMS, "nominal": [p0[n] for n in PARAMS], "rich_fit": rich_x, "sparse_fit": sparse_x}).to_csv(
        DATA / "modifier_identifiability_fit_parameters.csv", index=False
    )

    plot_heatmap(rich_corr, sens)
    plot_profiles(profiles)

    print("Baseline outputs:", model_outputs(p0))
    print("Identifiability warnings:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()


