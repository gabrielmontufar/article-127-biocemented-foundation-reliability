"""Generate external MICP experimental-data inventory for the GTENG-15832 AE response.
This script stores source metadata and summarized, non-proprietary descriptions only.
It does not redistribute raw data from third-party papers or databases.
"""
from pathlib import Path
import csv

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)

rows = [
    ("DS-PRJ-2453", "Nafisi and Montoya 2019; Nafisi et al. 2021", "curated stress-path/shear database", "c-prime/phi-prime component bounds", "No footing/degradation sequence"),
    ("UCS-402-Talamkhani-2023", "Talamkhani 2023", "402-test UCS literature database", "UCS/CaCO3 treatment-factor prior", "No footing capacity"),
    ("UCS-443-Ahmad-2026", "Ahmad et al. 2026", "443-test UCS literature database", "independent UCS driver check", "No direct degradation"),
    ("PLT-Kulkarni-2021", "Kulkarni et al. 2021", "direct MICP PLT/UCS/permeability", "initial bearing and settlement trend", "No long-term degradation"),
    ("BEARING-DURABILITY-Tao-2025a", "Tao et al. 2025a", "bearing-capacity acid/freeze-thaw retention", "eta/lambda degradation bounds", "surface-crust specimen"),
    ("FIELD-MICP-Tao-2025b", "Tao et al. 2025b", "field MICP surface bearing/durability", "field treated-crust contribution", "not footing PLT"),
    ("FIELD-MICP-Meng-2021", "Meng et al. 2021", "field MICP bearing/wind erosion", "field surface strength persistence", "not footing PLT"),
    ("FIELD-DURABILITY-Zhang-2024", "Zhang et al. 2024", "field rainfall durability", "environmental falsification", "slope erosion context"),
    ("FIELD-WEATHERING-Ji-2024", "Ji et al. 2024", "16-month weathering", "residual cementation loss", "slope/drought context"),
    ("FIELD-EICP-Martin-2024", "Martin et al. 2024", "field EICP plate-load analogue", "carbonate-biocementation field bearing/stiffness", "EICP not MICP"),
    ("FOOTING-SLOPE-Pusadkar-2017", "Pusadkar et al. 2017", "strip footing on MICP-treated slope", "footing-scale trend", "slope geometry"),
    ("DISSOLUTION-Ribeiro-Gomez-2023", "Ribeiro and Gomez 2023", "dissolution columns", "eta/lambda degradation prior", "no mechanical footing capacity"),
]

out = DATA / "external_public_data_inventory_minimal_rebuild.csv"
with out.open("w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["source_id", "citation_key", "evidence_type", "supports", "limitation"])
    w.writerows(rows)
print(f"Wrote {out}")
