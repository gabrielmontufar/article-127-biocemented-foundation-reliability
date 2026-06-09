from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
SUPP_FIGS = ROOT / "figures"


def external_figures_dir() -> Path | None:
    """Return the outer manuscript figure folder only when the original
    upload-package structure is present.

    The supplementary archive may be extracted as a standalone ZIP. In that
    case the script should not create or reuse a stray ancestor-level
    `02 Figures` directory outside the package; it should save only to the
    supplement's local figures folder.
    """
    if ROOT.parent.name == "04 Supplemental data and code":
        candidate = ROOT.parent.parent / "02 Figures"
        if candidate.is_dir():
            return candidate
    return None


def save_figure(fig: plt.Figure, filename: str) -> None:
    folders = [SUPP_FIGS]
    main_figs = external_figures_dir()
    if main_figs is not None:
        folders.append(main_figs)
    for folder in folders:
        folder.mkdir(parents=True, exist_ok=True)
        fig.savefig(folder / filename, dpi=300, bbox_inches="tight")


def figure_degradation_band() -> None:
    history = pd.read_csv(DATA / "deterministic_time_history.csv")
    priors = pd.read_csv(DATA / "calibrated_degradation_priors.csv")
    t = history["time_year"].to_numpy()
    eta_r = 0.20
    lambdas = [0.012, 0.035, 0.060]
    eta = {lam: eta_r + (1.0 - eta_r) * np.exp(-lam * t) for lam in lambdas}

    fig, ax = plt.subplots(figsize=(6.6, 4.2))
    ax.fill_between(t, eta[0.060], eta[0.012], color="#8fb6c9", alpha=0.35, label="published-component envelope")
    ax.plot(t, eta[0.035], color="#1f4d5f", lw=2.2, label="benchmark scenario")
    ax.scatter([0, 10, 30, 50], np.interp([0, 10, 30, 50], t, history["eta"]), color="#b23a2e", zorder=3, label="reported benchmark years")
    subtitle = "V4 envelope: exposure-clock priors must be updated before final design"
    ax.set_title("AE Fig. 2. Degradation state with component-supported envelope")
    ax.text(0.01, 0.02, subtitle, transform=ax.transAxes, fontsize=8.5, color="#333333")
    ax.set_xlabel("Time (yr)")
    ax.set_ylabel("Active cementation state eta")
    ax.set_ylim(0.0, 1.05)
    ax.grid(True, alpha=0.25)
    ax.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    save_figure(fig, "Figure 2 AE degradation component envelope.png")
    plt.close(fig)


def figure_model_hierarchy() -> None:
    model = pd.read_csv(DATA / "model_hierarchy_ablation.csv")
    order = ["M0", "M1", "M3/M4", "M7", "M8"]
    model = model.set_index("model").loc[order].reset_index()
    colors = ["#8a8f93", "#c47f34", "#4f7d96", "#8268a8", "#2f6f4e"]

    fig, ax = plt.subplots(figsize=(6.8, 4.3))
    ax.bar(model["model"], model["beta50"], color=colors)
    ax.axhline(1.5, color="#8b1e1e", lw=1.5, ls="--", label="illustrative beta target")
    for idx, row in model.iterrows():
        ax.text(idx, row["beta50"] + 0.05, f"{row['beta50']:.2f}", ha="center", va="bottom", fontsize=9)
    ax.set_title("AE Fig. 4. Model hierarchy and permanent-credit overstatement")
    ax.set_ylabel("50-year reliability index beta")
    ax.set_xlabel("Model form")
    ax.set_ylim(min(-1.0, model["beta50"].min() - 0.25), max(3.0, model["beta50"].max() + 0.35))
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    save_figure(fig, "Figure 4 AE model hierarchy ablation.png")
    plt.close(fig)


def figure_design_boundary() -> None:
    design = pd.read_csv(DATA / "parametric_design_matrix.csv")
    subset = design[design["q_service_kpa"].astype(float) == 950.0].copy()
    subset["passes_fs15"] = subset["FS_50yr"].astype(float) >= 1.5

    fig, ax = plt.subplots(figsize=(6.6, 4.4))
    sc = ax.scatter(
        subset["H_over_B"].astype(float),
        subset["lambda_per_year"].astype(float),
        c=subset["FS_50yr"].astype(float),
        s=95,
        cmap="viridis",
        edgecolor=np.where(subset["passes_fs15"], "#111111", "#b23a2e"),
        linewidth=0.9,
    )
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("FS at 50 yr")
    ax.set_title("AE Fig. 5. Design map with V4 screening boundary")
    ax.set_xlabel("Treatment depth H/B")
    ax.set_ylabel("Degradation rate lambda_e (1/yr)")
    ax.grid(True, alpha=0.25)
    ax.text(0.02, 0.03, "Red outline: FS<1.5 at qservice=950 kPa", transform=ax.transAxes, fontsize=8.5)
    fig.tight_layout()
    save_figure(fig, "Figure 5 AE V4 design boundary.png")
    plt.close(fig)


def main() -> None:
    figure_degradation_band()
    figure_model_hierarchy()
    figure_design_boundary()
    print("Generated AE resubmission figures in:")
    print(SUPP_FIGS)
    main_figs = external_figures_dir()
    if main_figs is not None:
        print(main_figs)


if __name__ == "__main__":
    main()
