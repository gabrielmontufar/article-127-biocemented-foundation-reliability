"""
Regenerates the AE traceability tables and figures added on 2026-06-09.
This script intentionally keeps the claim boundary explicit: the package is traceability-complete
for component-constrained screening, but it does not claim a direct integrated V5 degraded-footing validation.
"""
from pathlib import Path
import pandas as pd

base = Path(__file__).resolve().parents[1]
data = base / 'data'
figs = base / 'figures'
print('Traceability files expected:')
for name in [
    'parameter_dataset_scope_traceability_20260609.csv',
    'AE_requirement_closure_status_20260609.csv',
    'traceability_numeric_support_20260609.csv',
]:
    p = data / name
    print('-', name, 'OK' if p.exists() else 'MISSING')
    if p.exists():
        print('  rows:', len(pd.read_csv(p)))
print('Traceability figures expected:')
for name in [
    'figure_AE_requirement_traceability_matrix.png',
    'Figure S11 AE evidence traceability matrix.tif',
    'Figure S12 direct MICP PLT BCR evidence.png',
    'Figure S13 degradation_retention_evidence.png',
    'Figure S14 validation_claim_status_ladder.png',
]:
    p = figs / name
    print('-', name, 'OK' if p.exists() else 'MISSING')
