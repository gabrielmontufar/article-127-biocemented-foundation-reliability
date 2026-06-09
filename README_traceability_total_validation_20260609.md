# Traceability-first validation package for AE response (GTENG-15832)

This supplement implements the revised editorial strategy: every reviewer-requested model component is mapped to a dataset class, a model use, an allowed claim, a blocked claim, and a remaining gap.

## What this package proves

The package shows that the resubmitted manuscript is no longer an unsupported speculative model. It provides traceability from experimental or numerical evidence to the degradation parameters, treated-zone contribution, scenario figures, and decision metrics.

## What this package does not claim

The package does not claim a single universal set of constants for every soil, climate, treatment protocol, or foundation geometry. It also does not claim direct V5 validation of long-term degraded MICP-treated shallow footings, because no public integrated dataset was found that combines footing/plate-load response, controlled degradation, c-prime/phi-prime measurement, and localization observation in one experiment.

## New files added in this revision

- `data/traceability_total_evidence_matrix_20260609.csv`: maps AE concerns to model components, evidence sources, allowed claims, blocked claims, and remaining gaps.
- `data/reviewer_requirement_completion_audit_20260609.csv`: gives the practical completion status of each AE request.
- `data/validation_residual_gap_register_20260609.csv`: lists what remains incomplete if the AE requires direct integrated V5 validation.
- `figures/figure_traceability_reviewer_requirement_map.png`: graphical evidence coverage map.
- `figures/figure_claim_status_by_component.png`: claim-status ladder showing calibrated, bounded, scenario-bounded, sensitivity-only, and not-identified components.
- `figures/figure_dataset_to_claim_traceability_chain.png`: dataset-to-claim traceability workflow.
- `code/generate_traceability_validation_figures.py`: script to regenerate these files.

## Editorial position

The strongest defendable claim is: component-calibrated and uncertainty-bounded validation for pre-design, with total traceability from evidence to claim boundary. The only major remaining gap is a direct integrated experiment of degraded MICP-treated shallow footing performance.
