from __future__ import annotations

from pathlib import Path
from math import erf, lgamma, sqrt

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
SUPP = ROOT / "04 Supplemental data and code" / "Supplementary files"
DATA = SUPP / "data"
FIGS = SUPP / "figures"

RNG = np.random.default_rng(127)
ETA0 = 0.95
LAMBDA_RANGE = (0.001, 0.10)
ETA_R_RANGE = (0.03, 0.65)
Q_UNTREATED = 961.0
Q_SERVICE = 900.0
BETA_TARGET = 3.0
BASELINE_Q_U_50_KPA = 1277.0568589673255
BASELINE_QSERVICE_MAX_BETA3_KPA = 556.2594690942205
SIGMA_LN_R = np.log(BASELINE_Q_U_50_KPA / BASELINE_QSERVICE_MAX_BETA3_KPA) / BETA_TARGET


def retention_model(t: np.ndarray | float, lambda_e: np.ndarray, eta_r: np.ndarray) -> np.ndarray:
    return eta_r + (ETA0 - eta_r) * np.exp(-lambda_e * np.asarray(t))


def normal_pdf(x: np.ndarray | float, loc: np.ndarray | float = 0.0, scale: float = 1.0) -> np.ndarray:
    z = (np.asarray(x) - loc) / scale
    return np.exp(-0.5 * z**2) / (scale * np.sqrt(2.0 * np.pi))


def normal_logpdf(x: np.ndarray | float, loc: np.ndarray | float = 0.0, scale: float = 1.0) -> np.ndarray:
    z = (np.asarray(x) - loc) / scale
    return -0.5 * z**2 - np.log(scale) - 0.5 * np.log(2.0 * np.pi)


def normal_cdf(x: np.ndarray | float) -> np.ndarray:
    return 0.5 * (1.0 + np.vectorize(erf)(np.asarray(x) / sqrt(2.0)))


def lognormal_pdf_raw(x: np.ndarray, mean: float, sigma_g: float) -> np.ndarray:
    sigma = np.log(sigma_g)
    mu = np.log(mean) - 0.5 * sigma**2
    x = np.asarray(x)
    out = np.zeros_like(x, dtype=float)
    mask = x > 0.0
    out[mask] = np.exp(-0.5 * ((np.log(x[mask]) - mu) / sigma) ** 2) / (x[mask] * sigma * np.sqrt(2.0 * np.pi))
    return out


def truncated_lognormal_pdf(x: np.ndarray, mean: float = 0.035, sigma_g: float = 1.75) -> np.ndarray:
    lo, hi = LAMBDA_RANGE
    raw = lognormal_pdf_raw(x, mean=mean, sigma_g=sigma_g)
    dense = np.linspace(lo, hi, 5000)
    z = np.trapezoid(lognormal_pdf_raw(dense, mean=mean, sigma_g=sigma_g), dense)
    return raw / z


def scaled_beta_pdf(x: np.ndarray, a: float = 2.3, b: float = 3.6) -> np.ndarray:
    lo, hi = ETA_R_RANGE
    z = (x - lo) / (hi - lo)
    out = np.zeros_like(z, dtype=float)
    mask = (z > 0.0) & (z < 1.0)
    log_beta_ab = lgamma(a) + lgamma(b) - lgamma(a + b)
    out[mask] = np.exp((a - 1.0) * np.log(z[mask]) + (b - 1.0) * np.log1p(-z[mask]) - log_beta_ab) / (hi - lo)
    return out


def summarize_weighted(values: np.ndarray, weights: np.ndarray) -> tuple[float, float, float, float]:
    order = np.argsort(values)
    xs = values[order]
    ws = weights[order] / np.sum(weights)
    cdf = np.cumsum(ws)
    mean = float(np.sum(values * weights) / np.sum(weights))
    q025 = float(np.interp(0.025, cdf, xs))
    q500 = float(np.interp(0.500, cdf, xs))
    q975 = float(np.interp(0.975, cdf, xs))
    return mean, q025, q500, q975


def reliability_from_eta(eta50: np.ndarray, delta_q_max: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mean_qu = Q_UNTREATED + delta_q_max * eta50
    beta50 = np.log(mean_qu / Q_SERVICE) / SIGMA_LN_R
    pf50 = normal_cdf(-beta50)
    qmax = mean_qu / np.exp(BETA_TARGET * SIGMA_LN_R)
    return beta50, pf50, qmax


def weighted_resample(lam_grid: np.ndarray, eta_grid: np.ndarray, posterior: np.ndarray, n: int = 8000) -> pd.DataFrame:
    flat = posterior.ravel()
    idx = RNG.choice(flat.size, size=n, replace=True, p=flat / flat.sum())
    ii, jj = np.unravel_index(idx, posterior.shape)
    return pd.DataFrame({"lambda_e": lam_grid[ii], "eta_r": eta_grid[jj]})


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    obs = pd.DataFrame(
        [
            ["y15", "Sharma et al. (2021) digitized durability-retention hold-out", 15.0, 0.83, 0.08, "lambda_e and eta_r"],
            ["y20", "Sharma et al. (2021) digitized durability-retention hold-out", 20.0, 0.65, 0.08, "lambda_e and eta_r"],
        ],
        columns=["observation", "source", "equivalent_exposure_age", "retention", "sigma_y", "parameter_informed"],
    )
    obs.to_csv(DATA / "bayesian_update_observations.csv", index=False)

    lam = np.linspace(LAMBDA_RANGE[0], LAMBDA_RANGE[1], 260)
    eta_r = np.linspace(ETA_R_RANGE[0], ETA_R_RANGE[1], 260)
    L, E = np.meshgrid(lam, eta_r, indexing="ij")

    prior = truncated_lognormal_pdf(L) * scaled_beta_pdf(E)
    prior = prior / prior.sum()

    log_like = np.zeros_like(prior)
    for row in obs.itertuples(index=False):
        y_hat = retention_model(row.equivalent_exposure_age, L, E)
        log_like += normal_logpdf(row.retention, loc=y_hat, scale=row.sigma_y)
    like = np.exp(log_like - np.max(log_like))
    posterior = prior * like
    posterior = posterior / posterior.sum()

    prior_lam = prior.sum(axis=1)
    prior_eta = prior.sum(axis=0)
    post_lam = posterior.sum(axis=1)
    post_eta = posterior.sum(axis=0)

    eta50_prior_grid = retention_model(50.0, L, E)
    eta50_post_grid = eta50_prior_grid
    eta50_prior_mean = float(np.sum(eta50_prior_grid * prior) / np.sum(prior))
    delta_q_max = (BASELINE_Q_U_50_KPA - Q_UNTREATED) / eta50_prior_mean
    beta_prior_grid, pf_prior_grid, qmax_prior_grid = reliability_from_eta(eta50_prior_grid, delta_q_max)
    beta_post_grid, pf_post_grid, qmax_post_grid = reliability_from_eta(eta50_post_grid, delta_q_max)

    summary_rows = []
    for name, values, weights_prior, weights_post, unit, consequence in [
        ("lambda_e", L.ravel(), prior.ravel(), posterior.ravel(), "1/yr", "degradation rate"),
        ("eta_r", E.ravel(), prior.ravel(), posterior.ravel(), "-", "residual cementation"),
        ("eta50", eta50_prior_grid.ravel(), prior.ravel(), posterior.ravel(), "-", "50-year retained cementation credit"),
        ("beta50", beta_prior_grid.ravel(), prior.ravel(), posterior.ravel(), "-", "50-year reliability index at qservice=900 kPa"),
        ("pf50", pf_prior_grid.ravel(), prior.ravel(), posterior.ravel(), "-", "50-year failure probability at qservice=900 kPa"),
        ("qservice_max_beta3", qmax_prior_grid.ravel(), prior.ravel(), posterior.ravel(), "kPa", "allowable service pressure for beta>=3"),
    ]:
        pmean, plo, pmed, phi = summarize_weighted(values, weights_prior)
        qmean, qlo, qmed, qhi = summarize_weighted(values, weights_post)
        summary_rows.append(
            {
                "quantity": name,
                "unit": unit,
                "prior_mean": pmean,
                "prior_median": pmed,
                "prior_ci95_low": plo,
                "prior_ci95_high": phi,
                "posterior_mean": qmean,
                "posterior_median": qmed,
                "posterior_ci95_low": qlo,
                "posterior_ci95_high": qhi,
                "change_mean": qmean - pmean,
                "design_consequence": consequence,
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(DATA / "bayesian_update_posterior_summary.csv", index=False)

    prior_samples = weighted_resample(lam, eta_r, prior, n=8000)
    prior_samples["sample_type"] = "prior"
    posterior_samples = weighted_resample(lam, eta_r, posterior, n=8000)
    posterior_samples["sample_type"] = "posterior"
    samples = pd.concat([prior_samples, posterior_samples], ignore_index=True)
    samples["eta50"] = retention_model(50.0, samples["lambda_e"].to_numpy(), samples["eta_r"].to_numpy())
    samples["beta50"], samples["pf50"], samples["qservice_max_beta3"] = reliability_from_eta(samples["eta50"].to_numpy(), delta_q_max)
    samples.to_csv(DATA / "bayesian_update_prior_posterior_samples.csv", index=False)

    times = np.arange(0.0, 50.1, 1.0)
    traj_rows = []
    for sample_type, weights in [("prior", prior), ("posterior", posterior)]:
        for t in times:
            eta_t = retention_model(t, L, E)
            beta_t, pf_t, qmax_t = reliability_from_eta(eta_t, delta_q_max)
            eta_mean, eta_lo, eta_med, eta_hi = summarize_weighted(eta_t.ravel(), weights.ravel())
            beta_mean, beta_lo, beta_med, beta_hi = summarize_weighted(beta_t.ravel(), weights.ravel())
            pf_mean, pf_lo, pf_med, pf_hi = summarize_weighted(pf_t.ravel(), weights.ravel())
            q_mean, q_lo, q_med, q_hi = summarize_weighted(qmax_t.ravel(), weights.ravel())
            traj_rows.append(
                {
                    "sample_type": sample_type,
                    "time_year": t,
                    "eta_mean": eta_mean,
                    "eta_ci95_low": eta_lo,
                    "eta_median": eta_med,
                    "eta_ci95_high": eta_hi,
                    "beta_mean": beta_mean,
                    "beta_ci95_low": beta_lo,
                    "beta_median": beta_med,
                    "beta_ci95_high": beta_hi,
                    "pf_mean": pf_mean,
                    "pf_ci95_low": pf_lo,
                    "pf_median": pf_med,
                    "pf_ci95_high": pf_hi,
                    "qservice_max_beta3_mean": q_mean,
                    "qservice_max_beta3_ci95_low": q_lo,
                    "qservice_max_beta3_median": q_med,
                    "qservice_max_beta3_ci95_high": q_hi,
                }
            )
    traj = pd.DataFrame(traj_rows)
    traj.to_csv(DATA / "bayesian_update_reliability_trajectories.csv", index=False)

    central = summary.set_index("quantity")
    decision = pd.DataFrame(
        [
            [
                "published-retention update",
                f"lambda_e {central.loc['lambda_e','prior_mean']:.3f}->{central.loc['lambda_e','posterior_mean']:.3f}; eta_r {central.loc['eta_r','prior_mean']:.2f}->{central.loc['eta_r','posterior_mean']:.2f}",
                f"{central.loc['beta50','prior_mean']:.2f}->{central.loc['beta50','posterior_mean']:.2f}",
                "qservice=900 kPa fails beta>=3; updated qservice,max is reported for screening",
            ],
            [
                "posterior central design",
                f"eta50 {central.loc['eta50','prior_mean']:.2f}->{central.loc['eta50','posterior_mean']:.2f}",
                f"pf50 {central.loc['pf50','prior_mean']:.3f}->{central.loc['pf50','posterior_mean']:.3f}",
                "reduce service pressure, require project-specific plate/load-settlement evidence, or redesign",
            ],
            [
                "posterior unfavorable tail",
                f"qservice,max 2.5%={central.loc['qservice_max_beta3','posterior_ci95_low']:.0f} kPa",
                f"beta50 2.5%={central.loc['beta50','posterior_ci95_low']:.2f}",
                "inspection/re-treatment trigger if monitoring follows the lower credible envelope",
            ],
        ],
        columns=["new_evidence", "posterior_change", "effect_on_beta50_or_pf50", "decision"],
    )
    decision.to_csv(DATA / "bayesian_update_design_decisions.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(8.2, 6.2), dpi=180)
    axes[0, 0].plot(lam, prior_lam / np.trapezoid(prior_lam, lam), label="prior", lw=1.8)
    axes[0, 0].plot(lam, post_lam / np.trapezoid(post_lam, lam), label="posterior", lw=1.8)
    axes[0, 0].set_xlabel("lambda_e (1/yr)")
    axes[0, 0].set_ylabel("density")
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[0, 0].grid(alpha=0.25)

    axes[0, 1].plot(eta_r, prior_eta / np.trapezoid(prior_eta, eta_r), label="prior", lw=1.8)
    axes[0, 1].plot(eta_r, post_eta / np.trapezoid(post_eta, eta_r), label="posterior", lw=1.8)
    axes[0, 1].set_xlabel("eta_r")
    axes[0, 1].set_ylabel("density")
    axes[0, 1].grid(alpha=0.25)

    for sample_type, color in [("prior", "#4c78a8"), ("posterior", "#b23a48")]:
        subset = traj[traj["sample_type"] == sample_type]
        axes[1, 0].plot(subset["time_year"], subset["eta_mean"], color=color, label=sample_type, lw=1.9)
        axes[1, 0].fill_between(subset["time_year"], subset["eta_ci95_low"], subset["eta_ci95_high"], color=color, alpha=0.12)
        axes[1, 1].plot(subset["time_year"], subset["beta_mean"], color=color, label=sample_type, lw=1.9)
        axes[1, 1].fill_between(subset["time_year"], subset["beta_ci95_low"], subset["beta_ci95_high"], color=color, alpha=0.12)
    axes[1, 0].scatter(obs["equivalent_exposure_age"], obs["retention"], c="#222222", s=24, zorder=5, label="update data")
    axes[1, 0].set_xlabel("equivalent exposure age")
    axes[1, 0].set_ylabel("retention eta(t)")
    axes[1, 0].grid(alpha=0.25)
    axes[1, 0].legend(frameon=False, fontsize=8)
    axes[1, 1].axhline(BETA_TARGET, color="#333333", ls="--", lw=1.0, label="beta target")
    axes[1, 1].set_xlabel("service life (yr)")
    axes[1, 1].set_ylabel("beta(t), qservice=900 kPa")
    axes[1, 1].grid(alpha=0.25)
    axes[1, 1].legend(frameon=False, fontsize=8)
    fig.suptitle("Bayesian updating of MICP degradation parameters and posterior reliability", fontsize=10.5)
    fig.tight_layout()
    fig.savefig(FIGS / "figure_bayesian_updating_reliability.png")
    fig.savefig(FIGS / "figure_bayesian_updating_reliability.svg")
    plt.close(fig)

    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()


