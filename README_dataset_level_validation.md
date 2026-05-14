# Dataset-level calibration and validation notes

This folder adds limited dataset-level checks for Article 127.

## Load-settlement data

Source: Kulkarni et al. (2021), Journal of Engineering and Technological Sciences, DOI `10.5614/j.eng.technol.sci.2021.53.6.2`.

The open PDF contains Figure 6, a load-settlement curve for the 120 mm x 120 mm square footing in untreated and MICP-treated sand, and Table 6 with bearing-capacity ratio and settlement-reduction summary for circular and square plates. The individual circular plate load-settlement curves are not visible in the open PDF. Therefore:

- `kulkarni2021_digitized_plate_load_data.csv` contains the digitized square-footing curves from Figure 6.
- `kulkarni2021_square_120mm_untreated.csv` and `kulkarni2021_square_120mm_treated.csv` split those digitized curves by treatment.
- `kulkarni2021_table6_bcr_srf_summary.csv` preserves the Table 6 circular/square summary values.
- The hyperbolic fit uses the pre-ultimate branch up to the Table 6 tangent-defined `qu`; post-`qu` plotted points are retained in the raw CSV but not used to calibrate the asymptotic branch.

This is not a universal calibration. It is a transparent dataset-level demonstration that the matrix stiffness and treated stiffness multiplier can be identified from an experimental curve.

## Durability retention data

Source: Sharma, Satyam, and Reddy (2021), DOI `10.1177/1056789521991196`.

The open abstract reports UCS strength reductions of 4.2%, 8.3%, 17%, and 35% after 5, 10, 15, and 20 freeze-thaw cycles for the 1PV-12TC treated samples. These values were converted to UCS retention ratios and fitted with the same residual exponential structure used for the cementation state.

The Geoderma Sharma and Satyam (2021) wetting-drying/ageing article confirms the exposure grid and variables, but no open full-text figure data were available from the institutional record. Therefore this supplement uses the accessible freeze-thaw environmental-cycling series as a durability identifiability check and explicitly avoids presenting it as wetting-drying calibration.

## Generated files

- `dataset_level_load_settlement_calibration.py`
- `dataset_level_degradation_calibration.py`
- `kulkarni2021_digitized_plate_load_data.csv`
- `kulkarni2021_square_120mm_untreated.csv`
- `kulkarni2021_square_120mm_treated.csv`
- `kulkarni2021_table6_bcr_srf_summary.csv`
- `dataset_level_load_settlement_fitted_parameters.csv`
- `dataset_level_load_settlement_metrics.csv`
- `sharma2021_digitized_durability_retention_data.csv`
- `dataset_level_degradation_fitted_parameters.csv`
- `dataset_level_degradation_metrics.csv`
- `figure_dataset_load_settlement_validation.png`
- `figure_dataset_degradation_validation.png`
