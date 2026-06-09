# MICP parameter inference supplement

This supplement provides calibrable maps from measurable MICP treatment,
strength, stiffness, and durability indicators to eta0, cb0', Delta tan phi,
and lambda_e. The numerical demonstration is synthetic and is not presented as
a universal calibration.

AE rejection and Q1 novelty gate:

- `code/ae_rejection_novelty_evidence_gate.py`
- `data/ae_rejection_evidence_gate.csv`
- `data/q1_novelty_positioning_gate.csv`
- `data/ae_rejection_novelty_evidence_gate_summary.json`

This gate keeps the claim bounded to a V4 component-constrained pre-design and
updating workflow. It explicitly blocks V5 direct degraded-footing validation,
universal annual degradation laws, experimentally calibrated treated-depth
functions, and identified localization-fragility claims until direct data are
added.
