# AE traceability evidence files added on 2026-06-09

This folder adds the data-to-parameter traceability layer requested after the AE concern. The purpose is not to claim direct universal V5 validation. The purpose is to show exactly which component is calibrated, bounded, sensitized, or not identified.

New data files:

- `parameter_dataset_scope_traceability_20260609.csv`: model component -> dataset -> support status -> allowed claim -> remaining gap.
- `traceability_validation_status_20260609.csv`: compact support score and remaining gap by component.
- `ae_requirement_dataset_support_matrix_20260609.csv`: numeric matrix behind Fig. S11 and manuscript Fig. 8.
- `remaining_V5_gap_register_20260609.csv`: direct V5 gaps that remain outside the public-data envelope.
- `degradation_retention_anchor_points_20260609.csv`: extracted or reported retention anchor points used to bound degradation scenarios.

New figures:

- `Figure S11 AE-requested traceability map.png`
- `Figure S12 component validation status ladder.png`
- `Figure S13 normalized PLT load-settlement anchor.png`
- `Figure S14 degradation retention anchors.png`
- `Figure S15 data-to-parameter traceability architecture.png`

Interpretation:

The revised package uses total traceability rather than an over-claim of universal direct validation. It supports a component-calibrated, uncertainty-bounded pre-design workflow. Direct V5 validation of degraded MICP footing performance remains a specific experimental gap because no public dataset was found that combines footing/PLT response, controlled degradation, c-prime/phi-prime measurement, and localization/post-peak observation in a single campaign.
