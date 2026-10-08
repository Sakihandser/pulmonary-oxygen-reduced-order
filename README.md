# Pulmonary oxygen reduced-order model

Steady reduced-order resistance network used to separate distal arterial resistance, perfusion heterogeneity, capillary transit time, alveolar ventilation and right-ventricular oxygen supply.

This repository accompanies the English manuscript *Distal arterial resistance, perfusion heterogeneity and blood oxygenation: a steady reduced-order model of the pulmonary circulation* (Hong and Yang, Sichuan University). Mingze Hong is the corresponding author.

## Layout

- `code/build_lung_model.py` builds the geometry and the baseline solver.
- `code/phase2_scans.py` runs the mechanism scans.
- `code/phase3_scans.py` runs the comparisons reported in the manuscript.
- `tests/test_model.py` contains the seven automated verification tests cited in the manuscript.
- `data/` holds the saved result tables cited in the text.

The generation-12 root count is calibrated so that healthy pulmonary vascular resistance is 1.06 Wood units. That value is a calibration target, not an independent validation.

## Reproduce

Install the dependencies, then run the tests from the repository root:

```text
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
```

The scan scripts import `build_lung_model.py` from the same `code/` directory. Run them with Python 3:

```text
python code/phase2_scans.py --no-plots
python code/phase3_scans.py --no-plots
```

Result tables already saved:

- `data/goal4_ladder.csv`
- `data/goal9_transit.csv`
- `data/goal10_vq.csv`
- `data/goal11_coronary.csv`
- `data/goal12_exchange.csv`
- `data/goal12_ordering.csv`

No patient data are included. The authors declare no competing interests.

## License

Code, tests and the saved result tables are released under the MIT License. See `LICENSE`.
