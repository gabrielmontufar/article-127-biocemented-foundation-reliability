# Modifier calibration and identifiability supplement

This supplement documents the calibration role of chiH, chiE and chis in the
screening model. The modifiers are not universal material constants. They are
project-calibrated or sensitivity-controlled factors used to separate element
scale cementation measurements from footing-scale mobilization, serviceability
stiffness, and localization or brittle response.

The recommended use is staged calibration:

1. Constrain eta(t,z), eta0, eta_r and lambda_e from independent cementation
   and durability observations.
2. Estimate cb0' and Delta tan phi from element tests.
3. Estimate kH only from bearing-capacity data with more than one treated depth.
4. Estimate kE only from settlement or stiffness data, preferably across more
   than one treated depth.
5. Estimate chimin and rs only if post-peak or localization-sensitive data are
   available; otherwise use conservative sensitivity envelopes.
6. Validate on data not used during calibration.

The companion script modifier_identifiability_analysis.py generates local
sensitivity, correlation, and profile-objective diagnostics that show why a
single load-settlement curve cannot identify all modifier and strength
parameters simultaneously.
