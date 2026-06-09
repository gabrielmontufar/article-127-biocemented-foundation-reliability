from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.linalg import cholesky
from scipy.stats import norm


ROOT = Path(__file__).resolve().parents[3]
SUPP = ROOT / "04 Supplemental data and code" / "Supplementary files"
DATA = SUPP / "data"
FIGS = SUPP / "figures"
CODE = SUPP / "code"


@dataclass(frozen=True)
class RFParams:
    seed: int = 1272026
    n_realizations: int = 5000
    B: float = 2.0
    Df: float = 1.0
    gamma: float = 18.0
    phi_m_deg: float = 32.0
    cb0_mean: float = 20.0
    delta_tan_phi_mean: float = 0.0753
    eta0_mean: float = 0.95
    eta0_cov: float = 0.22
    eta_r: float = 0.22
    lambda_mean: float = 0.035
    lambda_cov: float = 0.35
    t_year: float = 50.0
    q_service: float = 900.0
    theta_R_mean: float = 1.0
    theta_R_cov: float = 0.10
    theta_Q_mean: float = 1.0
    theta_Q_cov: float = 0.08
    chi_min: float = 0.70
    rs: float = 0.18
    nx: int = 11
    ny: int = 11
    nz: int = 9
    zp_over_B: float = 3.0
    kz: float = 1.15


def bearing_factors(phi_deg: np.ndarray | float):
    phi = np.deg2rad(phi_deg)
    nq = np.exp(np.pi * np.tan(phi)) * np.tan(np.pi / 4.0 + phi / 2.0) ** 2
    nc = (nq - 1.0) / np.tan(phi)
    ng = 2.0 * (nq + 1.0) * np.tan(phi)
    return nc, nq, ng


def lognormal_params(mean: float, cov: float) -> tuple[float, float]:
    sigma = np.sqrt(np.log(1.0 + cov**2))
    mu = np.log(mean) - 0.5 * sigma**2
    return mu, sigma


def logitnormal_params(mean: float, cov: float) -> tuple[float, float]:
    # Moment-matched by deterministic approximation; adequate for bounded
    # sensitivity benchmark and reproducible with fixed seed.
    m = np.clip(mean, 1e-4, 1 - 1e-4)
    sd = min(cov * mean, 0.30)
    sigma = sd / (m * (1 - m))
    mu = np.log(m / (1 - m))
    return mu, sigma


def exp_corr_1d(coords: np.ndarray, theta: float) -> np.ndarray:
    if theta <= 1e-12:
        return np.eye(len(coords))
    d = np.abs(coords[:, None] - coords[None, :])
    c = np.exp(-d / theta)
    c += np.eye(len(coords)) * 1e-8
    return c


def correlated_fields(rng: np.random.Generator, p: RFParams, theta_h_over_B: float, theta_v_over_B: float):
    x = np.linspace(-p.B / 2.0, p.B / 2.0, p.nx)
    y = np.linspace(-p.B / 2.0, p.B / 2.0, p.ny)
    z = np.linspace(0.0, p.zp_over_B * p.B, p.nz)
    lx = cholesky(exp_corr_1d(x, theta_h_over_B * p.B), lower=True)
    ly = cholesky(exp_corr_1d(y, theta_h_over_B * p.B), lower=True)
    lz = cholesky(exp_corr_1d(z, theta_v_over_B * p.B), lower=True)

    raw = rng.standard_normal((p.n_realizations, p.nz, p.ny, p.nx))
    f = np.einsum("ab,nbcd->nacd", lz, raw)
    f = np.einsum("ab,ncbd->ncad", ly, f)
    f = np.einsum("ab,ncdb->ncda", lx, f)
    f = (f - f.mean(axis=(1, 2, 3), keepdims=True)) / f.std(axis=(1, 2, 3), keepdims=True)

    raw2 = rng.standard_normal((p.n_realizations, p.nz, p.ny, p.nx))
    g = np.einsum("ab,nbcd->nacd", lz, raw2)
    g = np.einsum("ab,ncbd->ncad", ly, g)
    g = np.einsum("ab,ncdb->ncda", lx, g)
    g = (g - g.mean(axis=(1, 2, 3), keepdims=True)) / g.std(axis=(1, 2, 3), keepdims=True)

    mu_eta, sig_eta = logitnormal_params(p.eta0_mean, p.eta0_cov)
    eta0 = 1.0 / (1.0 + np.exp(-(mu_eta + sig_eta * f)))
    mu_lam, sig_lam = lognormal_params(p.lambda_mean, p.lambda_cov)
    lam = np.exp(mu_lam + sig_lam * g)
    return x, y, z, eta0, lam


def weights(p: RFParams):
    z = np.linspace(0.0, p.zp_over_B * p.B, p.nz)
    wz = np.exp(-z / (p.kz * p.B))
    w = np.ones((p.nz, p.ny, p.nx))
    w *= wz[:, None, None]
    w /= w.sum()
    return w


def qu_matrix(p: RFParams) -> float:
    nc_m, nq_m, ng_m = bearing_factors(p.phi_m_deg)
    return float(p.gamma * p.Df * nq_m + 0.5 * p.gamma * p.B * ng_m)


def qu_random_field(p: RFParams, eta0: np.ndarray, lam: np.ndarray) -> np.ndarray:
    eta = p.eta_r + (eta0 - p.eta_r) * np.exp(-lam * p.t_year)
    eta = np.clip(eta, 0.0, 1.2)
    chi_s = np.maximum(p.chi_min, 1.0 - p.rs * eta**2)
    cb = p.cb0_mean * eta**1.15
    tan_phi = np.tan(np.deg2rad(p.phi_m_deg)) + p.delta_tan_phi_mean * eta
    phi = np.rad2deg(np.arctan(tan_phi))
    nc, nq, ng = bearing_factors(phi)
    nc_m, nq_m, ng_m = bearing_factors(p.phi_m_deg)
    local_delta = chi_s * (
        cb * nc
        + p.gamma * p.Df * (nq - nq_m)
        + 0.5 * p.gamma * p.B * (ng - ng_m)
    )
    inc = np.sum(local_delta * weights(p)[None, :, :, :], axis=(1, 2, 3))
    return qu_matrix(p) + inc


def reliability_metrics(q: np.ndarray, p: RFParams, rng: np.random.Generator, averaged_beta: float) -> dict[str, float]:
    mu_r, sig_r = lognormal_params(p.theta_R_mean, p.theta_R_cov)
    mu_q, sig_q = lognormal_params(p.theta_Q_mean, p.theta_Q_cov)
    theta_r = np.exp(mu_r + sig_r * rng.standard_normal(len(q)))
    theta_q = np.exp(mu_q + sig_q * rng.standard_normal(len(q)))
    g = theta_r * q - theta_q * p.q_service
    pf = float(np.mean(g <= 0.0))
    beta = float(-norm.ppf(max(pf, 0.5 / len(q)))) if pf > 0 else float(-norm.ppf(0.5 / len(q)))
    si = 0.18 * p.q_service / (q * (1.0 + 4.2 * np.maximum((q - qu_matrix(p)) / 400.0, 0.0)))
    return {
        "mean_qu50": float(np.mean(q)),
        "cov_qu50": float(np.std(q, ddof=1) / np.mean(q)),
        "q5_qu50": float(np.quantile(q, 0.05)),
        "q50_qu50": float(np.quantile(q, 0.50)),
        "q95_qu50": float(np.quantile(q, 0.95)),
        "pf50": pf,
        "beta50": beta,
        "p_serviceability_index_gt_0p045": float(np.mean(si > 0.045)),
        "delta_beta_vs_averaged": beta - averaged_beta,
    }


def run_case(case: str, theta_h: float, theta_v: float, p: RFParams, averaged_beta: float):
    case_offset = sum((i + 1) * ord(ch) for i, ch in enumerate(case)) % 100000
    rng = np.random.default_rng(p.seed + case_offset)
    _, _, _, eta0, lam = correlated_fields(rng, p, theta_h, theta_v)
    q = qu_random_field(p, eta0, lam)
    m = reliability_metrics(q, p, rng, averaged_beta)
    m.update({"case": case, "theta_h_over_B": theta_h, "theta_v_over_B": theta_v})
    return m, q, eta0[0], eta0[-1]


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    CODE.mkdir(parents=True, exist_ok=True)
    p = RFParams()

    rng = np.random.default_rng(p.seed)
    eta_mean = p.eta_r + (p.eta0_mean - p.eta_r) * np.exp(-p.lambda_mean * p.t_year)
    eta0_avg = np.full((p.n_realizations, p.nz, p.ny, p.nx), p.eta0_mean)
    lam_avg = np.full_like(eta0_avg, p.lambda_mean)
    q_avg = qu_random_field(p, eta0_avg, lam_avg)
    avg_metrics = reliability_metrics(q_avg, p, rng, 0.0)
    averaged_beta = avg_metrics["beta50"]
    avg_metrics["delta_beta_vs_averaged"] = 0.0
    avg_metrics.update({"case": "A mean-field averaged model", "theta_h_over_B": 0.0, "theta_v_over_B": 0.0})

    cases = [
        ("B vertical variability only", 1e-9, 0.25),
        ("C short horizontal correlation", 0.25, 0.25),
        ("D moderate horizontal correlation", 1.00, 0.25),
        ("E long horizontal correlation", 2.00, 0.25),
        ("F strong anisotropy horizontal", 2.00, 0.10),
        ("G strong anisotropy vertical", 0.25, 1.00),
    ]
    rows = [avg_metrics]
    q_store = {"A averaged": q_avg}
    samples = {}
    for case, th, tv in cases:
        row, q, eta_low, eta_high = run_case(case, th, tv, p, averaged_beta)
        rows.append(row)
        q_store[case] = q
        if case in ["C short horizontal correlation", "E long horizontal correlation"]:
            samples[case] = (eta_low, eta_high)

    summary = pd.DataFrame(rows)
    summary = summary[
        [
            "case",
            "theta_h_over_B",
            "theta_v_over_B",
            "mean_qu50",
            "cov_qu50",
            "q5_qu50",
            "q50_qu50",
            "q95_qu50",
            "pf50",
            "beta50",
            "delta_beta_vs_averaged",
            "p_serviceability_index_gt_0p045",
        ]
    ]
    summary.to_csv(DATA / "spatial_random_field_summary_metrics.csv", index=False)
    summary.to_csv(DATA / "spatial_random_field_sensitivity.csv", index=False)

    params = pd.DataFrame(
        [
            ["eta0_mean", "Mean initial cementation state", p.eta0_mean, "fixed", "bounded logit-normal field"],
            ["eta0_cov", "COV of initial cementation state", p.eta0_cov, "fixed", "spatial variability amplitude"],
            ["lambda_e_mean", "Mean degradation rate", p.lambda_mean, "fixed", "lognormal positive field"],
            ["lambda_e_cov", "COV of degradation rate", p.lambda_cov, "fixed", "spatial variability amplitude"],
            ["theta_h/B", "Horizontal correlation length", "case-specific", "0.25, 1.0, 2.0", "plan correlation"],
            ["theta_v/B", "Vertical correlation length", "case-specific", "0.10, 0.25, 1.0", "depth correlation"],
            ["grid", "Random-field grid", f"{p.nx} x {p.ny} x {p.nz}", "fixed", "treated volume below footing"],
            ["n_realizations", "Monte Carlo random fields", p.n_realizations, "fixed", "reliability simulation"],
            ["seed", "Reproducibility seed", p.seed, "fixed", "random-number control"],
            ["transformation", "Field transform", "logit-normal eta0; lognormal lambda_e", "fixed", "physical bounds"],
        ],
        columns=["Parameter", "Meaning", "Baseline value", "Sensitivity range", "Role"],
    )
    params.to_csv(DATA / "spatial_random_field_parameters.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(8, 3.6), dpi=180)
    for ax, (case, (eta_low, eta_high)) in zip(axes.ravel(), samples.items()):
        im = ax.imshow(eta_low[0], vmin=0, vmax=1, cmap="viridis", origin="lower")
        ax.set_title(case.replace(" correlation", " corr."))
        ax.set_xlabel("x grid")
        ax.set_ylabel("y grid")
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.86, label="eta0")
    fig.savefig(FIGS / "figure_random_field_realizations.png", bbox_inches="tight")
    plt.close(fig)

    rf = summary[
        summary["case"].isin(
            [
                "C short horizontal correlation",
                "D moderate horizontal correlation",
                "E long horizontal correlation",
            ]
        )
    ].sort_values("theta_h_over_B")
    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=180)
    ax.plot(rf["theta_h_over_B"], rf["beta50"], "o-", color="#1f5a7a")
    ax.axhline(averaged_beta, color="black", ls="--", lw=1, label="averaged model")
    ax.set_xlabel("Horizontal correlation length, theta_h/B")
    ax.set_ylabel("50-year reliability index beta")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False)
    fig.savefig(FIGS / "figure_beta_vs_correlation_length.png", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.4, 4.2), dpi=180)
    for label in ["A averaged", "C short horizontal correlation", "E long horizontal correlation"]:
        q = np.sort(q_store[label])
        cdf = np.linspace(1 / len(q), 1.0, len(q))
        ax.plot(q, cdf, label=label)
    ax.set_xlabel("50-year bearing capacity q_u (kPa)")
    ax.set_ylabel("Empirical CDF")
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIGS / "figure_qu_cdf_random_field_vs_average.png", bbox_inches="tight")
    plt.close(fig)

    ths = [0.25, 0.50, 1.00, 2.00]
    tvs = [0.10, 0.25, 0.50, 1.00]
    heat_rows = []
    heat = np.zeros((len(tvs), len(ths)))
    for i, tv in enumerate(tvs):
        for j, th in enumerate(ths):
            row, _, _, _ = run_case(f"grid_{th}_{tv}", th, tv, p, averaged_beta)
            heat[i, j] = row["delta_beta_vs_averaged"]
            heat_rows.append(row)
    pd.DataFrame(heat_rows).to_csv(DATA / "spatial_random_field_correlation_grid.csv", index=False)
    fig, ax = plt.subplots(figsize=(6.0, 4.4), dpi=180)
    im = ax.imshow(heat, cmap="coolwarm", origin="lower", aspect="auto")
    ax.set_xticks(range(len(ths)), [str(v) for v in ths])
    ax.set_yticks(range(len(tvs)), [str(v) for v in tvs])
    ax.set_xlabel("theta_h/B")
    ax.set_ylabel("theta_v/B")
    ax.set_title("Reliability loss: beta_RF - beta_avg")
    fig.colorbar(im, ax=ax, label="Delta beta")
    fig.savefig(FIGS / "figure_spatial_correlation_heatmap.png", bbox_inches="tight")
    plt.close(fig)

    (CODE / "spatial_random_field_reliability.py").write_text(Path(__file__).read_text(encoding="utf-8"), encoding="utf-8")
    (CODE / "spatial_random_field_cementation.py").write_text(Path(__file__).read_text(encoding="utf-8"), encoding="utf-8")
    print(DATA / "spatial_random_field_summary_metrics.csv")
    print(FIGS / "figure_random_field_realizations.png")


if __name__ == "__main__":
    main()


