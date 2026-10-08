# Pulmonary oxygen reduced-order model

Steady reduced-order resistance network used to separate distal arterial resistance, perfusion heterogeneity, capillary transit time, alveolar ventilation and right-ventricular oxygen supply.

This repository accompanies the English manuscript *Distal arterial resistance, perfusion heterogeneity and blood oxygenation: a steady reduced-order model of the pulmonary circulation* (Hong and Yang, Sichuan University).

## Layout

- `code/build_lung_model.py` builds the geometry and the baseline solver.
- `code/phase2_scans.py` runs the mechanism scans.
- `code/phase3_scans.py` runs the comparisons reported in the manuscript.
- `data/` holds the saved result tables cited in the text.

The generation-12 root count is calibrated so that healthy pulmonary vascular resistance is 1.06 Wood units. That value is a calibration target, not an independent validation.

## Reproduce

The scan scripts import `build_lung_model.py` from the same `code/` directory. Run them with Python 3, NumPy, pandas and Matplotlib.

Result tables already saved:

- `data/goal4_ladder.csv`
- `data/goal9_transit.csv`
- `data/goal10_vq.csv`
- `data/goal11_coronary.csv`
- `data/goal12_exchange.csv`
- `data/goal12_ordering.csv`

No patient data are included.
