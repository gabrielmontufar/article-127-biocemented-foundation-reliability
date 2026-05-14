from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT
DATA = SUPP / "data"


def scenario_metadata() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "Scenario": "S0",
                "CDR range": "0.00",
                "Number of cycles": "0",
                "Environmental amplifier": "1.0",
                "Purpose": "Static monotonic benchmark",
                "Expected damage level": "None",
            },
            {
                "Scenario": "C1",
                "CDR range": "0.07-0.10",
                "Number of cycles": "low annual service cycles",
                "Environmental amplifier": "1.0",
                "Purpose": "Low-amplitude service cycling near or below threshold",
                "Expected damage level": "Small",
            },
            {
                "Scenario": "C2",
                "CDR range": "0.12-0.18",
                "Number of cycles": "high cycle count",
                "Environmental amplifier": "1.0",
                "Purpose": "Moderate repeated loading",
                "Expected damage level": "Moderate",
            },
            {
                "Scenario": "C3",
                "CDR range": "0.28-0.38",
                "Number of cycles": "few intermittent event blocks",
                "Environmental amplifier": "1.0",
                "Purpose": "High-amplitude intermittent loading",
                "Expected damage level": "Severe event-driven",
            },
            {
                "Scenario": "C4",
                "CDR range": "0.05-0.34",
                "Number of cycles": "variable block history",
                "Environmental amplifier": "1.0",
                "Purpose": "History-dependent accumulation",
                "Expected damage level": "Variable",
            },
            {
                "Scenario": "C5",
                "CDR range": "0.12-0.28",
                "Number of cycles": "moderate to high cycle count",
                "Environmental amplifier": "1.0-1.8",
                "Purpose": "Combined environment and cyclic loading",
                "Expected damage level": "Amplified",
            },
        ]
    )


def build_loading_blocks(t: np.ndarray) -> pd.DataFrame:
    rows = []
    for sid in ["S0", "C1", "C2", "C3", "C4", "C5"]:
        for ti in t:
            if sid == "S0":
                cdr, cycles, renv = 0.0, 0.0, 1.0
            elif sid == "C1":
                cdr = 0.07 + 0.03 * (0.5 + 0.5 * np.sin(2 * np.pi * ti / 5.0))
                cycles, renv = 2.0e3, 1.0
            elif sid == "C2":
                cdr = 0.14 + 0.04 * (0.5 + 0.5 * np.sin(2 * np.pi * ti / 10.0))
                cycles, renv = 2.0e4, 1.0
            elif sid == "C3":
                event = (9 <= ti <= 11) or (24 <= ti <= 26) or (39 <= ti <= 41)
                cdr = 0.36 if event else 0.04
                cycles = 600.0 if event else 0.0
                renv = 1.0
            elif sid == "C4":
                if ti < 12:
                    cdr, cycles = 0.06, 3.0e3
                elif ti < 24:
                    cdr, cycles = 0.16, 1.2e4
                elif ti < 33:
                    cdr, cycles = 0.30, 2.0e3
                elif ti < 42:
                    cdr, cycles = 0.10, 4.0e3
                else:
                    cdr, cycles = 0.22, 9.0e3
                renv = 1.0
            else:
                cdr = 0.15 + 0.06 * (0.5 + 0.5 * np.sin(2 * np.pi * ti / 7.0))
                wet_period = (8 <= ti <= 18) or (30 <= ti <= 42)
                cdr += 0.07 if wet_period else 0.0
                cycles = 1.5e4 if not wet_period else 2.5e4
                renv = 1.8 if wet_period else 1.15

            rows.append(
                {
                    "time_year": ti,
                    "scenario": sid,
                    "CDR": cdr,
                    "cycles": cycles,
                    "R_env": renv,
                    "qmax_over_qu": 0.62 + 0.5 * cdr,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    t = np.linspace(0, 50, 101)
    scenario_metadata().to_csv(DATA / "cyclic_loading_scenarios_table.csv", index=False)
    build_loading_blocks(t).to_csv(DATA / "cyclic_loading_histories.csv", index=False)
    print(DATA / "cyclic_loading_scenarios_table.csv")
    print(DATA / "cyclic_loading_histories.csv")


if __name__ == "__main__":
    main()
