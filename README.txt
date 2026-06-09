Time-Dependent Reliability of Shallow Foundations on Degrading Biocemented Sand

This package contains the reproducible benchmark data, scripts, figures, and
tables used for the ASCE JGGE resubmission. Scripts resolve paths from their own
location so the package can be copied and rerun outside the original archive.

AE resubmission validation
--------------------------
Run from the supplementary package root:

python code/run_ae_resubmission_validation.py

This command regenerates:

- data/ae_parameter_evidence_matrix.csv
- data/calibrated_degradation_priors.csv
- data/component_calibration_metrics.csv
- data/model_hierarchy_ablation.csv
- data/ae_resubmission_validation_manifest.json

The AE-specific figure files can be regenerated with:

python code/generate_ae_resubmission_figures.py

This writes the AE degradation envelope, model-hierarchy ablation, and V4 design
boundary figures to both the supplement figures folder and the manuscript figure
folder.

The validation boundary is V4 component-constrained indirect validation. Direct
physical degraded-footing validation is not claimed.

Main benchmark reproduction
---------------------------
Run from the supplementary package root:

python code/reproduce_article_127_benchmark.py

The benchmark script writes outputs to reproduced_outputs and does not delete or
overwrite the source data.

Public repository
-----------------
https://github.com/gabrielmontufar/article-127-biocemented-foundation-reliability

Traceability-first validation update (2026-06-09)
-------------------------------------------------
This archive now includes traceability_total_evidence_matrix_20260609.csv, reviewer_requirement_completion_audit_20260609.csv, validation_residual_gap_register_20260609.csv, and three traceability figures. These files map each AE request to evidence, allowed claims, blocked claims, and residual gaps.
