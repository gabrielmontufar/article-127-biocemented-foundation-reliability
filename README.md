# Time-Dependent Reliability of Shallow Foundations on Degrading Biocemented Sand

This repository contains the reproducible benchmark data, scripts, figures, and
tables for the article 127 resubmission package.

## GTENG-15832 traceability update

The repository was updated from the corrected package:

`GTENG-15832_traceability_evidence_upload_package_20260609/GTENG-15832_corrected_package_20260609`

The 2026-06-09 update adds traceability-first validation material requested in
the editorial review cycle, including:

- AE requirement closure and traceability matrices;
- dataset-level degradation and load-settlement calibration evidence;
- external quantitative validation and residual-gap registers;
- traceability, validation-status, and component-architecture figures;
- reproducible scripts for the AE resubmission validation and traceability
  figures.

The validation boundary remains evidence-bounded: the package supports
component-constrained indirect validation and does not claim direct degraded
MICP footing validation.

## Reproduce

Run from the repository root:

```bash
python code/reproduce_article_127_benchmark.py
python code/run_ae_resubmission_validation.py
python code/generate_ae_resubmission_figures.py
python code/regenerate_traceability_evidence_figures.py
```

The scripts write outputs to `reproduced_outputs` and do not delete source data.

## Public Repository

https://github.com/gabrielmontufar/article-127-biocemented-foundation-reliability
