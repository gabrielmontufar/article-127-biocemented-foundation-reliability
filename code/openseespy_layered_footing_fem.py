from __future__ import annotations

import json
import math
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import openseespy.opensees as ops


BASE = Path(__file__).resolve().parents[3]
SUPP = BASE / "04 Supplemental data and code" / "Supplementary files"
DATA = SUPP / "data"
FIGS = SUPP / "figures"
CODE = SUPP / "code"
AUDIT = BASE / "05 Journal selection and audits"
UPLOAD = BASE / "00 Files for journal upload"


B = 2.0
HALF_WIDTH = B / 2.0
DOMAIN_X = 8.0
DOMAIN_Y = 6.0
NX = 32
NY = 24
THICKNESS = 1.0


def bearing_factors(phi_deg: float) -> tuple[float, float, float]:
    phi = math.radians(phi_deg)
    nq = math.exp(math.pi * math.tan(phi)) * math.tan(math.radians(45.0) + phi / 2.0) ** 2
    nc = (nq - 1.0) / math.tan(phi)
    ngamma = 2.0 * (nq + 1.0) * math.tan(phi)
    return nc, nq, ngamma


def eq8_capacity_50yr(h_over_b: float) -> float:
    eta = 0.22 + (1.0 - 0.22) * math.exp(-0.035 * 50.0)
    cb0 = 20.0
    phi_m = 32.0
    dtan = math.tan(math.radians(35.0)) - math.tan(math.radians(phi_m))
    phi = math.degrees(math.atan(math.tan(math.radians(phi_m)) + dtan * eta))
    cb = cb0 * eta**1.15
    chi_h = 1.0 - math.exp(-h_over_b / 0.70)
    chi_s = min(1.0, max(0.72, 1.0 - 0.18 * eta**2.0))
    nc, nq, ng = bearing_factors(phi)
    _, nq_m, ng_m = bearing_factors(phi_m)
    gamma = 18.0
    df = 1.0
    q0 = gamma * df
    q_matrix = q0 * nq_m + 0.5 * gamma * B * ng_m
    delta = cb * nc + q0 * (nq - nq_m) + 0.5 * gamma * B * (ng - ng_m)
    return q_matrix + chi_h * chi_s * delta


def dp_material(tag: int, *, strength_scale: float, stiffness_scale: float) -> None:
    # Drucker-Prager parameters are used here as a reproducible open-source
    # numerical benchmark. They are not a calibrated constitutive law.
    e_mod = 30_000.0 * stiffness_scale
    nu = 0.30
    bulk = e_mod / (3.0 * (1.0 - 2.0 * nu))
    shear = e_mod / (2.0 * (1.0 + nu))
    sigma_y = 18.0 * strength_scale
    rho = 0.32
    rho_bar = 0.20
    k_inf = bulk
    k0 = bulk
    ops.nDMaterial("DruckerPrager", tag, bulk, shear, sigma_y, rho, rho_bar, k_inf, k0, 0.0, 0.0, 0.0, 0.0, 1.8)


def build_model(h_over_b: float) -> list[int]:
    ops.wipe()
    ops.model("basic", "-ndm", 2, "-ndf", 2)
    dx = DOMAIN_X / NX
    dy = DOMAIN_Y / NY
    node_id = {}
    tag = 1
    for j in range(NY + 1):
        y = -j * dy
        for i in range(NX + 1):
            x = i * dx
            node_id[(i, j)] = tag
            ops.node(tag, x, y)
            tag += 1

    dp_material(1, strength_scale=1.0, stiffness_scale=1.0)
    dp_material(2, strength_scale=1.0 + 0.85 * 0.3556262363, stiffness_scale=1.0 + 1.25 * 0.3556262363)

    etag = 1
    treated_depth = h_over_b * B
    for j in range(NY):
        y_centroid = -(j + 0.5) * dy
        depth = -y_centroid
        mat = 2 if depth <= treated_depth and h_over_b > 0 else 1
        for i in range(NX):
            n1 = node_id[(i, j)]
            n2 = node_id[(i + 1, j)]
            n3 = node_id[(i + 1, j + 1)]
            n4 = node_id[(i, j + 1)]
            ops.element("quad", etag, n1, n2, n3, n4, THICKNESS, "PlaneStrain", mat, 0.0, 1.8, 0.0, 0.0)
            etag += 1

    top_footing = []
    control = node_id[(0, 0)]
    fix_x = set()
    fix_y = set()
    for i in range(NX + 1):
        n = node_id[(i, 0)]
        x = i * dx
        if abs(x) < 1e-9:
            fix_x.add(n)
        if x <= HALF_WIDTH + 1e-9:
            top_footing.append(n)
            fix_x.add(n)
            if n != control:
                ops.equalDOF(control, n, 2)

    for i in range(NX + 1):
        fix_x.add(node_id[(i, NY)])
        fix_y.add(node_id[(i, NY)])
    for j in range(NY + 1):
        fix_x.add(node_id[(NX, j)])
        fix_x.add(node_id[(0, j)])

    for n in set(fix_x) | set(fix_y):
        ops.fix(n, 1 if n in fix_x else 0, 1 if n in fix_y else 0)

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    ops.load(control, 0.0, -1.0)
    ops.constraints("Transformation")
    ops.numberer("RCM")
    ops.system("BandGeneral")
    ops.test("NormDispIncr", 1.0e-6, 40)
    ops.algorithm("Newton")
    ops.integrator("DisplacementControl", control, 2, -0.001)
    ops.analysis("Static")
    return top_footing


def run_case(h_over_b: float) -> dict[str, float]:
    footing_nodes = build_model(h_over_b)
    records = []
    peak = 0.0
    peak_disp = 0.0
    for step in range(180):
        ok = ops.analyze(1)
        if ok != 0:
            # Retry with a smaller step and modified Newton before accepting nonconvergence.
            ops.integrator("DisplacementControl", footing_nodes[0], 2, -0.0004)
            ops.algorithm("ModifiedNewton")
            ok = ops.analyze(1)
        disp = abs(ops.nodeDisp(footing_nodes[0], 2))
        q = abs(ops.getLoadFactor(1)) / (HALF_WIDTH * THICKNESS)
        records.append((disp, q, ok))
        if q > peak:
            peak = q
            peak_disp = disp
        if ok != 0:
            break
        if step > 40 and q < 0.85 * peak:
            break
    ops.wipe()
    return {
        "H_over_B": h_over_b,
        "fem_peak_pressure_model_units": peak,
        "fem_peak_displacement_m": peak_disp,
        "steps_completed": len(records),
        "last_analyze_code": records[-1][2] if records else 999,
    }


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    CODE.mkdir(parents=True, exist_ok=True)
    AUDIT.mkdir(parents=True, exist_ok=True)

    h_values = [0.0, 0.5, 1.0, 1.5, 2.0]
    rows = [run_case(h) for h in h_values]
    df = pd.DataFrame(rows)
    untreated = float(df.loc[df["H_over_B"] == 0.0, "fem_peak_pressure_model_units"].iloc[0])
    df["fem_normalized_to_untreated"] = df["fem_peak_pressure_model_units"] / untreated
    eq8_vals = np.array([eq8_capacity_50yr(h) for h in h_values])
    df["eq8_kpa_50yr"] = eq8_vals
    df["eq8_normalized_to_untreated"] = eq8_vals / eq8_vals[0]
    df["normalized_difference_percent"] = 100.0 * (
        df["eq8_normalized_to_untreated"] - df["fem_normalized_to_untreated"]
    ) / df["fem_normalized_to_untreated"]

    csv = DATA / "round7_openseespy_layered_footing_fem.csv"
    df.to_csv(csv, index=False)

    fig = FIGS / "Figure 9 OpenSeesPy FEM layered footing benchmark.png"
    plt.figure(figsize=(6.8, 4.4), dpi=200)
    plt.plot(df["H_over_B"], df["fem_normalized_to_untreated"], "o-", lw=2.2, label="OpenSeesPy DP FEM")
    plt.plot(df["H_over_B"], df["eq8_normalized_to_untreated"], "s--", lw=2.0, label="Eq. (8), normalized")
    plt.xlabel("Treatment depth ratio, H/B")
    plt.ylabel("Capacity ratio normalized to untreated case")
    plt.grid(True, alpha=0.25)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(fig)
    plt.close()

    shutil.copy2(Path(__file__), CODE / "openseespy_layered_footing_fem.py")
    shutil.copy2(fig, BASE / "02 Figures" / fig.name)
    shutil.copy2(fig, UPLOAD / fig.name)

    summary = {
        "python": "3.12",
        "openseespy_version": ops.version(),
        "csv": str(csv),
        "figure": str(fig),
        "script": str(CODE / "openseespy_layered_footing_fem.py"),
        "cases": rows,
        "max_abs_normalized_difference_percent": float(np.max(np.abs(df["normalized_difference_percent"]))),
    }
    report = AUDIT / "round7_openseespy_fem_report.json"
    report.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()


