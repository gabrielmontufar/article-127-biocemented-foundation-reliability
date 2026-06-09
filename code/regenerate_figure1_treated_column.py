from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[3]
FIG = (
    ROOT
    / "04 Supplemental data and code"
    / "Supplementary files"
    / "figures"
    / "Figure 1 footing treated zone and weighting function.png"
)


def main() -> None:
    B = 2.0
    H = 2.4
    ground_y = 0.0
    footing_h = 0.46
    x0 = -B / 2.0
    x1 = B / 2.0

    fig, (ax, axw) = plt.subplots(
        1,
        2,
        figsize=(8.2, 4.0),
        dpi=220,
        gridspec_kw={"width_ratios": [2.25, 0.85], "wspace": 0.28},
    )

    ax.axhline(ground_y, color="black", lw=1.6)
    ax.add_patch(
        Rectangle((x0, ground_y + 0.10), B, footing_h, facecolor="#d9d9d9", edgecolor="black", lw=1.2)
    )
    ax.add_patch(
        Rectangle(
            (x0, ground_y - H),
            B,
            H,
            facecolor="#dbeaf8",
            edgecolor="#4a86b8",
            lw=1.2,
        )
    )

    ax.annotate(
        "",
        xy=(x0, ground_y + footing_h + 0.25),
        xytext=(x1, ground_y + footing_h + 0.25),
        arrowprops=dict(arrowstyle="<->", lw=1.0, color="black"),
    )
    ax.text(0.0, ground_y + footing_h + 0.34, "footing width B", ha="center", va="bottom", fontsize=10)

    ax.annotate(
        "",
        xy=(x0 - 0.28, ground_y),
        xytext=(x0 - 0.28, ground_y - H),
        arrowprops=dict(arrowstyle="<->", lw=1.0, color="black"),
    )
    ax.text(x0 - 0.43, ground_y - H / 2, "treated depth H", rotation=90, ha="center", va="center", fontsize=10)

    ax.text(0.0, ground_y - H - 0.18, "treated column under footing", ha="center", va="top", fontsize=9)
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-H - 0.55, 0.95)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")

    z = np.linspace(0.0, H, 200)
    kz = 1.15
    w = np.exp(-kz * z / B)
    w = w / np.trapezoid(w, z)
    w = w / w.max()

    axw.plot(w, z, color="#b73745", lw=2.2)
    axw.invert_yaxis()
    axw.set_xlabel("normalized\ninfluence w(z)", fontsize=9)
    axw.set_ylabel("depth z", fontsize=9)
    axw.set_xlim(0.0, 1.08)
    axw.set_ylim(H, 0.0)
    axw.grid(True, alpha=0.25)
    axw.tick_params(labelsize=8)
    axw.spines["top"].set_visible(False)
    axw.spines["right"].set_visible(False)

    fig.savefig(FIG, bbox_inches="tight")
    plt.close(fig)
    print(FIG)


if __name__ == "__main__":
    main()


