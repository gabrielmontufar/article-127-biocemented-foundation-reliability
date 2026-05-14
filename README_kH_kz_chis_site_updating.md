# Practical updating of kH, kz, and chi_s

This supplement provides a synthetic, reproducible site-updating example for
the footing-scale modifiers kH, kz, and chi_s. The data are not field
measurements; they are deliberately labeled synthetic and are used only to
demonstrate the calibration workflow.

The staged workflow is:

1. Fit untreated plate-load response before any treated-zone modifier is tuned.
2. Bound cb0' and Delta tan phi from element tests or conservative priors.
3. Update kz from depth-resolved CPT/CPTu, Vs, UPV, or carbonate profiles.
4. Update kH from a plate-load or numerical depth-benefit series.
5. Keep chi_s as a sensitivity parameter unless post-peak, localization, or
   calibrated elastoplastic numerical evidence is available.

Generated outputs include kH_profile_objective.csv,
kz_cpt_profile_objective.csv, chis_sensitivity_summary.csv,
kH_kz_chis_calibration_summary.csv, and four PNG figures documenting the
depth-benefit fit, CPT weighting fit, localization sensitivity, and staged
updating workflow.
