from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
SUPP = ROOT / "04 Supplemental data and code" / "Supplementary files"
DATA = SUPP / "data"


def scenario_table() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Scenario": "A",
                "Name": "Stationary mild groundwater",
                "Chemistry": "Near-neutral pH and non-aggressive calcite saturation proxy",
                "Flow regime": "Low hydraulic flushing",
                "Wetting-drying frequency": "No major cycles",
                "Mitigation": "None required in benchmark",
                "Purpose": "Stationary-exposure reference close to base model",
            },
            {
                "Scenario": "B",
                "Name": "Undersaturated groundwater chemistry",
                "Chemistry": "Negative calcite-saturation proxy or low-pH interval",
                "Flow regime": "Low hydraulic flushing",
                "Wetting-drying frequency": "No major cycles",
                "Mitigation": "None",
                "Purpose": "Isolate chemistry-driven acceleration",
            },
            {
                "Scenario": "C",
                "Name": "Hydraulic flushing pulses",
                "Chemistry": "Mild chemistry",
                "Flow regime": "High-gradient flushing pulses during selected intervals",
                "Wetting-drying frequency": "No major cycles",
                "Mitigation": "None",
                "Purpose": "Isolate hydraulic flushing effects",
            },
            {
                "Scenario": "D",
                "Name": "Seasonal wetting-drying",
                "Chemistry": "Mild chemistry",
                "Flow regime": "Moderate flow",
                "Wetting-drying frequency": "Semiannual cycles with moderate saturation amplitude",
                "Mitigation": "None",
                "Purpose": "Represent cyclic environmental exposure",
            },
            {
                "Scenario": "E",
                "Name": "Combined aggressive environment",
                "Chemistry": "Undersaturated or low-pH intervals",
                "Flow regime": "Repeated flushing pulses",
                "Wetting-drying frequency": "Frequent high-amplitude cycles",
                "Mitigation": "None",
                "Purpose": "Stress-test combined exposure severity",
            },
            {
                "Scenario": "F",
                "Name": "Mitigated environment",
                "Chemistry": "Buffered chemistry",
                "Flow regime": "Controlled drainage and reduced flushing",
                "Wetting-drying frequency": "Reduced cycle amplitude",
                "Mitigation": "Drainage and chemistry control",
                "Purpose": "Show effect of environmental mitigation",
            },
        ]
    )


def build_histories(t: np.ndarray) -> pd.DataFrame:
    rows = []
    for sid in ["A", "B", "C", "D", "E", "F"]:
        for ti in t:
            if sid == "A":
                schem = 1.00
                sflow = 1.00
                swd = 1.00
            elif sid == "B":
                schem = 1.45 + (0.55 if 8 <= ti <= 30 else 0.00)
                sflow = 0.75
                swd = 0.65
            elif sid == "C":
                schem = 0.75
                pulse = (10 <= ti <= 16) or (28 <= ti <= 34)
                sflow = 0.75 + (1.25 if pulse else 0.00)
                swd = 0.60
            elif sid == "D":
                schem = 0.75
                sflow = 0.75
                swd = 0.85 + 0.65 * (0.5 + 0.5 * np.sin(2 * np.pi * ti))
            elif sid == "E":
                schem = 1.35 + (0.65 if 5 <= ti <= 38 else 0.00)
                pulse = (7 <= ti <= 14) or (21 <= ti <= 27) or (37 <= ti <= 44)
                sflow = 1.10 + (1.10 if pulse else 0.00)
                swd = 1.15 + 0.75 * (0.5 + 0.5 * np.sin(4 * np.pi * ti))
            else:
                schem = 0.38
                sflow = 0.42
                swd = 0.25 + 0.15 * (0.5 + 0.5 * np.sin(2 * np.pi * ti))

            # Transparent weighted severity for screening. It is not reactive
            # geochemistry; it maps exposure descriptors to a degradation clock.
            senv = 0.45 * schem + 0.35 * min(sflow, 2.0) + 0.20 * swd
            rows.append(
                {
                    "time_year": ti,
                    "scenario": sid,
                    "S_chem": schem,
                    "S_flow": min(sflow, 2.0),
                    "S_WD": swd,
                    "S_env": senv,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    t = np.linspace(0, 50, 201)
    scenario_table().to_csv(DATA / "environmental_scenarios_table.csv", index=False)
    build_histories(t).to_csv(DATA / "environmental_histories.csv", index=False)
    print(DATA / "environmental_scenarios_table.csv")
    print(DATA / "environmental_histories.csv")


if __name__ == "__main__":
    main()


