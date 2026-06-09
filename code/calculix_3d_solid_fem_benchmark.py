from __future__ import annotations

import math
import os
import re
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ARTICLE_ROOT = Path(__file__).resolve().parents[3]
SUPP_ROOT = ARTICLE_ROOT / "04 Supplemental data and code" / "Supplementary files"
DATA_DIR = SUPP_ROOT / "data"
FIG_DIR = SUPP_ROOT / "figures"
CODE_DIR = SUPP_ROOT / "code"
WORK_DIR = Path(__file__).resolve().parent / "calculix_3d_runs"


@dataclass(frozen=True)
class ModelParams:
    footing_width: float = 2.0
    domain_width: float = 8.0
    domain_depth: float = 6.0
    nx: int = 16
    ny: int = 16
    nz: int = 12
    e_untreated: float = 30.0e6
    nu_untreated: float = 0.30
    e_gain_factor: float = 3.0
    q_service: float = 100.0e3
    lambda_e: float = 0.035
    eta_residual: float = 0.22
    service_year: float = 50.0

    @property
    def eta_50(self) -> float:
        return self.eta_residual + (1.0 - self.eta_residual) * math.exp(
            -self.lambda_e * self.service_year
        )

    @property
    def e_treated(self) -> float:
        return self.e_untreated * (1.0 + self.e_gain_factor * self.eta_50)


def find_ccx() -> str:
    ccx = shutil.which("ccx")
    if ccx:
        return ccx

    env_root = Path(os.environ.get("MAMBA_ROOT_PREFIX", "")) / "envs" / "solid3d"
    candidates = [
        env_root / "Library" / "bin" / "ccx.exe",
        Path.home() / "AppData" / "Roaming" / "mamba" / "envs" / "solid3d" / "Library" / "bin" / "ccx.exe",
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate)
    raise FileNotFoundError("ccx.exe was not found in PATH or in the solid3d environment.")


def node_id(i: int, j: int, k: int, nx: int, ny: int) -> int:
    return 1 + i + (nx + 1) * (j + (ny + 1) * k)


def chunks(values: list[int], size: int = 12) -> list[list[int]]:
    return [values[i : i + size] for i in range(0, len(values), size)]


def write_inp(case_dir: Path, h_over_b: float, params: ModelParams) -> tuple[Path, list[int]]:
    case_dir.mkdir(parents=True, exist_ok=True)
    job = case_dir / f"HB_{str(f'{h_over_b:.1f}').replace('.', 'p')}"
    inp = job.with_suffix(".inp")

    xs = np.linspace(-params.domain_width / 2, params.domain_width / 2, params.nx + 1)
    ys = np.linspace(-params.domain_width / 2, params.domain_width / 2, params.ny + 1)
    zs = np.linspace(0.0, params.domain_depth, params.nz + 1)
    dx = xs[1] - xs[0]
    dy = ys[1] - ys[0]

    half_b = params.footing_width / 2.0
    treatment_depth = h_over_b * params.footing_width
    treated_elems: list[int] = []
    untreated_elems: list[int] = []
    bottom_nodes: list[int] = []
    x_min_nodes: list[int] = []
    x_max_nodes: list[int] = []
    y_min_nodes: list[int] = []
    y_max_nodes: list[int] = []
    footing_nodes: list[int] = []

    for k, z in enumerate(zs):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                nid = node_id(i, j, k, params.nx, params.ny)
                if k == 0:
                    bottom_nodes.append(nid)
                if i == 0:
                    x_min_nodes.append(nid)
                if i == params.nx:
                    x_max_nodes.append(nid)
                if j == 0:
                    y_min_nodes.append(nid)
                if j == params.ny:
                    y_max_nodes.append(nid)
                if k == params.nz and abs(x) <= half_b + 1e-9 and abs(y) <= half_b + 1e-9:
                    footing_nodes.append(nid)

    total_force = params.q_service * params.footing_width**2
    nodal_force = -total_force / len(footing_nodes)

    lines: list[str] = []
    lines.append("*HEADING")
    lines.append(
        f"Article 127 CalculiX 3D solid elastic benchmark, H/B={h_over_b:.1f}"
    )
    lines.append("*NODE")
    for k, z in enumerate(zs):
        for j, y in enumerate(ys):
            for i, x in enumerate(xs):
                lines.append(f"{node_id(i, j, k, params.nx, params.ny)}, {x:.6f}, {y:.6f}, {z:.6f}")

    lines.append("*ELEMENT, TYPE=C3D8, ELSET=ALL_ELEMENTS")
    eid = 1
    for k in range(params.nz):
        zc = 0.5 * (zs[k] + zs[k + 1])
        depth_below_surface = params.domain_depth - zc
        for j in range(params.ny):
            yc = 0.5 * (ys[j] + ys[j + 1])
            for i in range(params.nx):
                xc = 0.5 * (xs[i] + xs[i + 1])
                n1 = node_id(i, j, k, params.nx, params.ny)
                n2 = node_id(i + 1, j, k, params.nx, params.ny)
                n3 = node_id(i + 1, j + 1, k, params.nx, params.ny)
                n4 = node_id(i, j + 1, k, params.nx, params.ny)
                n5 = node_id(i, j, k + 1, params.nx, params.ny)
                n6 = node_id(i + 1, j, k + 1, params.nx, params.ny)
                n7 = node_id(i + 1, j + 1, k + 1, params.nx, params.ny)
                n8 = node_id(i, j + 1, k + 1, params.nx, params.ny)
                lines.append(f"{eid}, {n1}, {n2}, {n3}, {n4}, {n5}, {n6}, {n7}, {n8}")

                in_treated_column = abs(xc) <= half_b and abs(yc) <= half_b
                if treatment_depth > 0 and in_treated_column and depth_below_surface <= treatment_depth:
                    treated_elems.append(eid)
                else:
                    untreated_elems.append(eid)
                eid += 1

    def add_elset(name: str, values: list[int]) -> None:
        lines.append(f"*ELSET, ELSET={name}")
        for part in chunks(values, 16):
            lines.append(", ".join(map(str, part)))

    def add_nset(name: str, values: list[int]) -> None:
        lines.append(f"*NSET, NSET={name}")
        for part in chunks(values, 16):
            lines.append(", ".join(map(str, part)))

    add_elset("UNTREATED_ELEMENTS", untreated_elems)
    if treated_elems:
        add_elset("TREATED_ELEMENTS", treated_elems)

    add_nset("BOTTOM", bottom_nodes)
    add_nset("XMIN", x_min_nodes)
    add_nset("XMAX", x_max_nodes)
    add_nset("YMIN", y_min_nodes)
    add_nset("YMAX", y_max_nodes)
    add_nset("FOOTING_NODES", footing_nodes)

    lines.extend(
        [
            "*MATERIAL, NAME=UNTREATED",
            "*ELASTIC",
            f"{params.e_untreated:.6e}, {params.nu_untreated:.6f}",
            "*MATERIAL, NAME=TREATED",
            "*ELASTIC",
            f"{params.e_treated:.6e}, {params.nu_untreated:.6f}",
            "*SOLID SECTION, ELSET=UNTREATED_ELEMENTS, MATERIAL=UNTREATED",
            "",
        ]
    )
    if treated_elems:
        lines.extend(["*SOLID SECTION, ELSET=TREATED_ELEMENTS, MATERIAL=TREATED", ""])

    lines.extend(
        [
            "*BOUNDARY",
            "BOTTOM, 1, 3, 0.0",
            "XMIN, 1, 1, 0.0",
            "XMAX, 1, 1, 0.0",
            "YMIN, 2, 2, 0.0",
            "YMAX, 2, 2, 0.0",
            "*STEP",
            "*STATIC",
        ]
    )
    for nid in footing_nodes:
        lines.append(f"*CLOAD\n{nid}, 3, {nodal_force:.9e}")
    lines.extend(
        [
            "*NODE PRINT, NSET=FOOTING_NODES",
            "U",
            "*END STEP",
        ]
    )

    inp.write_text("\n".join(lines) + "\n", encoding="ascii")
    return job, footing_nodes


def parse_average_uz(dat_file: Path, footing_nodes: list[int]) -> float:
    text = dat_file.read_text(errors="ignore")
    node_values: dict[int, float] = {}
    node_set = set(footing_nodes)
    pattern = re.compile(
        r"^\s*(\d+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s*$"
    )
    for line in text.splitlines():
        match = pattern.match(line)
        if not match:
            continue
        nid = int(match.group(1))
        if nid in node_set:
            node_values[nid] = float(match.group(4))

    if not node_values:
        raise RuntimeError(f"No FOOTING_NODES displacement rows were parsed from {dat_file}.")
    return float(np.mean(list(node_values.values())))


def run_case(h_over_b: float, params: ModelParams, ccx: str) -> dict[str, float]:
    case_dir = WORK_DIR / f"H_over_B_{h_over_b:.1f}"
    job, footing_nodes = write_inp(case_dir, h_over_b, params)
    result = subprocess.run(
        [ccx, job.name],
        cwd=case_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=180,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"CalculiX failed for H/B={h_over_b:.1f} with code {result.returncode}.\n"
            f"{result.stdout}"
        )
    uz = parse_average_uz(job.with_suffix(".dat"), footing_nodes)
    return {
        "H_over_B": h_over_b,
        "treated_depth_m": h_over_b * params.footing_width,
        "treated_E_MPa": params.e_treated / 1.0e6 if h_over_b > 0 else params.e_untreated / 1.0e6,
        "untreated_E_MPa": params.e_untreated / 1.0e6,
        "eta_50": params.eta_50,
        "mean_footing_Uz_m": uz,
        "mean_settlement_mm": -1000.0 * uz,
    }


def make_outputs(rows: list[dict[str, float]], params: ModelParams) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    CODE_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(rows).sort_values("H_over_B")
    s0 = float(df.loc[df["H_over_B"] == 0.0, "mean_settlement_mm"].iloc[0])
    df["fem_stiffness_ratio"] = s0 / df["mean_settlement_mm"]
    df["settlement_reduction_percent"] = 100.0 * (1.0 - df["mean_settlement_mm"] / s0)
    df["screening_proxy_stiffness_ratio"] = 1.0 + params.e_gain_factor * params.eta_50 * (
        1.0 - np.exp(-df["H_over_B"] / 0.55)
    )
    df["proxy_minus_fem_ratio"] = (
        df["screening_proxy_stiffness_ratio"] - df["fem_stiffness_ratio"]
    )

    csv_path = DATA_DIR / "round9_calculix_3d_solid_fem.csv"
    df.to_csv(csv_path, index=False)

    fig_path = FIG_DIR / "Figure 10 CalculiX 3D solid FEM settlement benchmark.png"
    fig, ax1 = plt.subplots(figsize=(7.2, 4.6), dpi=180)
    ax1.plot(
        df["H_over_B"],
        df["fem_stiffness_ratio"],
        "o-",
        color="#1f5a7a",
        label="CalculiX 3D solid FEM",
    )
    ax1.plot(
        df["H_over_B"],
        df["screening_proxy_stiffness_ratio"],
        "s--",
        color="#a33f2f",
        label="Screening stiffness proxy",
    )
    ax1.set_xlabel("Treatment depth ratio, H/B")
    ax1.set_ylabel("Normalized stiffness ratio")
    ax1.grid(True, alpha=0.28)

    ax2 = ax1.twinx()
    line3 = ax2.plot(
        df["H_over_B"],
        df["settlement_reduction_percent"],
        "^-",
        color="#4d7f36",
        label="Settlement reduction",
    )[0]
    ax2.set_ylabel("Settlement reduction (%)")
    ax2.tick_params(axis="y", colors="#4d7f36")
    handles, labels = ax1.get_legend_handles_labels()
    ax1.legend(handles + [line3], labels + [line3.get_label()], frameon=False, loc="upper left")
    fig.tight_layout()
    fig.savefig(fig_path, bbox_inches="tight")
    plt.close(fig)

    shutil.copy2(Path(__file__).resolve(), CODE_DIR / "calculix_3d_solid_fem_benchmark.py")

    report = Path(__file__).resolve().parent / "calculix_3d_solid_fem_report.txt"
    report.write_text(
        "\n".join(
            [
                "CalculiX 3D solid FEM benchmark for Article 127",
                f"CalculiX input/run directory: {WORK_DIR}",
                f"CSV: {csv_path}",
                f"Figure: {fig_path}",
                f"Script copy: {CODE_DIR / 'calculix_3d_solid_fem_benchmark.py'}",
                "",
                "Model scope:",
                "- 3D linear-elastic solid block with C3D8 brick elements.",
                "- Square footing nodal load distributed over the top footing nodes.",
                "- Treated material assigned to a square column beneath the footing for each H/B.",
                "- Output is a serviceability/stiffness benchmark, not a plastic-collapse capacity theorem.",
                "",
                df.to_string(index=False),
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def main() -> None:
    params = ModelParams()
    ccx = find_ccx()
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    for h_over_b in [0.0, 0.5, 1.0, 1.5, 2.0]:
        rows.append(run_case(h_over_b, params, ccx))
    make_outputs(rows, params)
    print("CalculiX 3D solid FEM benchmark completed.")
    print(f"Run directory: {WORK_DIR}")
    print(f"CSV: {DATA_DIR / 'round9_calculix_3d_solid_fem.csv'}")
    print(f"Figure: {FIG_DIR / 'Figure 10 CalculiX 3D solid FEM settlement benchmark.png'}")


if __name__ == "__main__":
    main()


