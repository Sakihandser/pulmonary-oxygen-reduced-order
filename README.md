# Pulmonary oxygen reduced-order model

Steady reduced-order resistance network used to separate distal arterial resistance, perfusion heterogeneity, capillary transit time, alveolar ventilation and right-ventricular oxygen supply.

This repository accompanies the English manuscript *Perfusion redistribution and oxygen exchange at matched pulmonary vascular resistance in a reduced-order model* (Hong and Yang, Sichuan University). Mingze Hong is the corresponding author.

## Layout

- `code/build_lung_model.py` builds the geometry and the baseline solver.
- `code/phase2_scans.py` runs the mechanism scans.
- `code/phase3_scans.py` runs the comparisons reported in the manuscript.
- `code/revision_scans.py` runs the oxygen-conservation, regional-sensitivity,
  exchange-calibration, convergence, and repeated-sampling checks added after
  review of manuscript v04.
- `tests/test_model.py` contains 11 automated verification tests.
- `data/` holds the saved result tables cited in the text.

The order-12 root count is obtained from an approximate large-vessel plus tree
target.  Solving the complete five-lobe network gives a healthy pulmonary
vascular resistance of 1.05693 Wood units.  This is a calibration result, not
independent validation.  The cardiac-output parameter of 5 L/min is a
reference scale, not a hard upper bound; the solved healthy value is 5.06470
L/min under the prescribed pressure-flow relation.

The parameter historically named `shunt` is now defined as an effective
post-exchange venous-admixture fraction.  It does not create a separately
resolved anatomical bypass.  Pulmonary effective oxygen uptake and systemic
Fick uptake are both multiplied by the non-admixed fraction, so their mass
balance closes for admixture fractions 0, 0.03, 0.10, and 0.20, including a
low-ventilation test.

## Reproduce

Install the dependencies, then run the tests from the repository root:

```text
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
```

The scan scripts import `build_lung_model.py` from the same `code/` directory. Run them with Python 3:

```text
python code/phase2_scans.py --no-plots
python code/phase3_scans.py
python code/revision_scans.py
```

`phase2_scans.py` accepts `--no-plots`. `phase3_scans.py` writes result tables
by default and accepts `--plots` for its legacy figures.  The review-revision
script writes CSV tables and a PNG sensitivity figure to `data/`.

Result tables already saved:

- `data/goal3_pao2.csv`
- `data/goal3_shunt.csv`
- `data/goal4_ladder.csv`
- `data/goal6_rv_power.csv`
- `data/goal9_transit.csv`
- `data/goal10_vq.csv`
- `data/goal11_coronary.csv`
- `data/goal12_exchange.csv`
- `data/goal12_ordering.csv`
- `data/revision_oxygen_conservation.csv`
- `data/revision_contrast_definition.csv`
- `data/revision_core_sensitivity.csv`
- `data/revision_exchange_calibration.csv`
- `data/revision_exchange_model_sensitivity.csv`
- `data/revision_numerical_convergence.csv`
- `data/revision_sampling_stability.csv`
- `data/revision_model_parameters.csv`
- `data/revision_organ_parameters.csv`
- `data/revision_calibration_diagnostics.csv`
- `data/revision_validation_benchmarks.csv`

No patient data are included. The authors declare no competing interests.

## License

Code, tests and the saved result tables are released under the MIT License. See `LICENSE`.

