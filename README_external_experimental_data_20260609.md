# External experimental data search added for GTENG-15832 resubmission

Date added: 2026-06-09.

This note records the public-data/evidence search added after the AE critique and expanded after a second search pass for closer bearing-capacity and degradation evidence. No public source was identified that directly reports long-term bearing-capacity degradation of shallow foundations on MICP-treated sand. Therefore, the manuscript remains bounded as V4 component-constrained indirect validation and does not claim V5 direct degraded-footing validation.

The second search pass added closer evidence streams:

- Direct laboratory MICP plate-load/load-settlement evidence from Kulkarni et al. (2021).
- MICP bearing-capacity loss/retention under acid and freeze-thaw exposure from Tao et al. (2025a) and Tao et al. (2025c).
- Field MICP bearing-capacity and crust-performance evidence from Meng et al. (2021) and Tao et al. (2025b).
- Field exposure and durability-limit evidence from Zhang et al. (2024) and Ji et al. (2024).
- Field EICP plate-load evidence from Martin et al. (2024), used only as an analogue for carbonate biocementation, not as MICP calibration.

Added/updated files:

- `data/public_external_experimental_data_inventory_20260609.csv`
- `data/external_dataset_to_AE_requirement_map_20260609.csv`
- `data/direct_bearing_degradation_evidence_20260609.csv`
- `data/v4_component_dataset_inventory.csv`
- `code/external_public_data_inventory.py`

The inventory distinguishes strength databases, PLT/footing tests, durability/dissolution tests, bearing-retention tests, field-scale context sources, and analogue field EICP evidence. The AE map states which model component each source can support and which claim remains unsupported.


## 2026-06-09 traceability-complete revision

This revision converts the external-evidence discussion into a parameter--dataset--claim-boundary traceability package. The purpose is to answer the AE request without overstating the validation level.

New files added:

- `data/parameter_dataset_traceability_matrix_20260609.csv`: maps each AE concern to model component, supporting datasets, evidence scale, allowed claim, disallowed claim, and supplemental files.
- `data/reviewer_request_completion_matrix_20260609.csv`: summarizes which requested items are now supported, bounded, or left outside the claim.
- `data/validation_gap_closure_matrix_20260609.csv`: explicitly states the remaining direct-V5 evidence gap and how the manuscript handles it.
- `public_literature_search_gap_report_20260609.md`: search report documenting the evidence mosaic and the absence of a single public degraded-MICP-footing dataset.
- `figures/figure_traceability_evidence_map.png`: visual map from AE concern to evidence class.
- `figures/figure_reviewer_completion_status.png`: component-by-component support status.
- `figures/figure_kulkarni_bcr_traceability.png`: direct short-term PLT/BCR evidence used for the treated-zone contribution.
- `figures/figure_degradation_retention_traceability.png`: degradation-retention calibration/validation evidence.
- `figures/figure_validation_gap_closure_map.png`: what is supported, bounded, or still a gap.

Editorial boundary: the evidence supports a component-calibrated and uncertainty-bounded screening model. It does not support a direct universal V5 claim for long-term degraded MICP footing performance.

## Traceability-complete update added after final evidence review

A final traceability layer was added to align the revision with the AE request without over-claiming direct V5 validation. The added layer does not state that a universal constant set has been validated for every soil, climate, and geometry. Instead, it documents which model quantities are calibrated, bounded, sensitized, or still non-identified.

New files added in this update:

- `data/parameter_dataset_scope_traceability_20260609.csv`: parameter-dataset-scope matrix linking each AE concern to evidence, numerical anchor, model operation, allowed claim, and disallowed claim.
- `data/AE_requirement_closure_status_20260609.csv`: closure-status table stating what has been substantiated and what remains outside the public evidence envelope.
- `data/traceability_numeric_support_20260609.csv`: numeric support table with PLT BCR values, settlement reduction, degradation-retention values, and field bearing-capacity anchors.
- `figures/figure_AE_requirement_traceability_matrix.png` and `figures/Figure S11 AE evidence traceability matrix.tif`: evidence stream to model-component matrix.
- `figures/Figure S12 direct MICP PLT BCR evidence.png`: direct MICP PLT bearing-capacity ratio evidence.
- `figures/Figure S13 degradation_retention_evidence.png`: normalized durability and bearing-retention evidence used to bound `eta(t)`.
- `figures/Figure S14 validation_claim_status_ladder.png`: explicit classification of calibrated, bounded, sensitivity-only, and remaining-gap components.
- `code/traceability_figures_and_tables.py`: compact verification script for the traceability files.

The remaining direct V5 gap is intentionally preserved: no public integrated dataset was identified that combines degraded MICP footing/plate capacity, independent c-prime/phi-prime evolution, stiffness, and post-peak localization in the same experimental campaign. The manuscript therefore frames the claim as traceability-complete component-constrained validation for pre-design and local updating, not as a direct universal degraded-footing law.


## 2026-06-09 traceability update

The package now includes `parameter_dataset_scope_traceability_20260609.csv`, `traceability_validation_status_20260609.csv`, `ae_requirement_dataset_support_matrix_20260609.csv`, `degradation_retention_anchor_points_20260609.csv`, `remaining_V5_gap_register_20260609.csv`, and Figs. S11-S15. These files implement a calibrated/bounded/sensitized/not-identified evidence classification for the AE-requested components.
