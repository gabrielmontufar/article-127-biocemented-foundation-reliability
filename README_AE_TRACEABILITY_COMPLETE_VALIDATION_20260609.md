# AE traceability-complete validation update for GTENG-15832

This file documents the final revision made to respond to the AE request for data supporting the degradation parameters, treated-zone contribution function, localization fragility factor, and scenario figures.

## Validation claim now used

The revised package does **not** claim a direct universal V5 degraded-footing validation. It claims a traceability-complete, component-constrained validation architecture for pre-design and local updating. Each contested quantity is assigned one of four statuses:

1. **calibrated/direct**: direct experimental data exist for the relevant response, e.g., MICP plate-load/load-settlement BCR from Kulkarni et al. (2021);
2. **bounded**: experimental evidence supports a range or trend, but not a unique universal constant, e.g., stress-path data for strength increments and durability data for `eta(t)`;
3. **sensitivity/falsification**: evidence supports capping or rejection of a model assumption but not calibration of a parameter, e.g., `chi_s` localization fragility;
4. **remaining gap**: public evidence does not support the requested direct integrated validation, e.g., degraded MICP footing capacity with paired c-prime/phi-prime/stiffness/localization measurements.

## Files added

- `data/parameter_dataset_scope_traceability_20260609.csv`
- `data/AE_requirement_closure_status_20260609.csv`
- `data/traceability_numeric_support_20260609.csv`
- `figures/figure_AE_requirement_traceability_matrix.png`
- `figures/Figure S11 AE evidence traceability matrix.tif`
- `figures/Figure S12 direct MICP PLT BCR evidence.png`
- `figures/Figure S13 degradation_retention_evidence.png`
- `figures/Figure S14 validation_claim_status_ladder.png`
- `code/traceability_figures_and_tables.py`

## Remaining direct V5 gap

The remaining gap is not hidden: no public dataset was found that repeats shallow-foundation or plate-load capacity tests after controlled MICP degradation and also reports paired `c'`, `phi'`, stiffness, CaCO3 redistribution, and post-peak localization. The revised manuscript therefore converts the AE response into a data-traceability argument rather than a universal-law assertion.
