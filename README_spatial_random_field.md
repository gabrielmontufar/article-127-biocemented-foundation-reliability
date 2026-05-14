# Spatial random-field cementation benchmark

This supplement adds a synthetic spatial random-field sensitivity benchmark for Article 127. It is not field validation because project-scale maps of MICP cementation were not available.

## Purpose

The benchmark evaluates whether the mean depth-weighted cementation model can hide spatially coherent weak zones beneath a footing. The random-field model treats:

- `eta_0(x,y,z)` as a bounded logit-normal field.
- `lambda_e(x,y,z)` as a positive lognormal field.
- `eta(t,x,y,z)=eta_r+[eta_0(x,y,z)-eta_r] exp[-lambda_e(x,y,z)t]`.

Local cementation increments are integrated with a normalized 3D influence function over the treated volume below the footing.

## Correlation model

The implemented generator uses separable anisotropic exponential covariance in x, y, and z, which is computationally efficient and reproducible for the 11 x 11 x 9 grid. Horizontal and vertical correlation lengths are reported as `theta_h/B` and `theta_v/B`.

## Main output

The mean-field model gives `beta50 = 3.35`. Spatial random-field cases reduce beta by about `0.20` to `0.63`, with the largest loss in the long-horizontal-correlation case. This supports the manuscript statement that spatial averaging is acceptable only as a screening mean-field approximation; final design should use spatial characterization or conservative lower-tail cementation parameters when coherent under-treated zones are plausible.

## Files

- `code/spatial_random_field_reliability.py`
- `code/spatial_random_field_cementation.py`
- `data/spatial_random_field_parameters.csv`
- `data/spatial_random_field_summary_metrics.csv`
- `data/spatial_random_field_sensitivity.csv`
- `data/spatial_random_field_correlation_grid.csv`
- `figures/figure_random_field_realizations.png`
- `figures/figure_beta_vs_correlation_length.png`
- `figures/figure_qu_cdf_random_field_vs_average.png`
- `figures/figure_spatial_correlation_heatmap.png`
