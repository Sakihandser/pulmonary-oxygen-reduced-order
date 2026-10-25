# Pulmonary oxygen exchange benchmark

Code and numerical data accompanying **Calibration-Dependent Transferability of Pulmonary Oxygen Exchange Models in a Controlled Numerical Benchmark**, by Mingze Hong and Bojun Yang.

The study compares saturation relaxation, content relaxation and constant conductance after common endpoint or inlet-slope calibration. Fitted parameters are frozen for synthetic inlet, flow and alveolar tests. Additional calculations explain the dissolved-oxygen limit and regional allocation sensitivity.

## Reproduce the results

```text
cd numerical
python -m pip install -r requirements-lock.txt
python reproduce.py
```

Use a fresh working copy because the runner writes outputs. After dependency installation it runs offline, checks file integrity, regenerates 53 result CSVs, verifies numerical main Tables 4–8 and exports seven main and two supplementary figures. See [the result-to-code index](numerical/RESULT_REPRODUCTION_INDEX.md) for exact scripts, data files and commands. The locked dependencies describe the tested software environment.

For figure-only export, from `numerical/` run `python sensitivity/reporting/reproduce_figures.py`. To run the tests:

```text
python -m unittest discover -s tests
python -m unittest discover -s aggregate_analysis/tests
python -m unittest discover -s observation_analysis/tests
```

## Interpretation

These are controlled numerical comparisons. Synthetic states and grid counts under assumed bounded errors are not patient observations or clinical discrimination probabilities. CT inspiratory air-volume shares are regional allocation proxies. Gas-derived capacity uses unpaired aggregate observations with unresolved normalization. The external Alvin pressure submodule supplies numerical qualification, not validation of the full physiological model. Source cells, assumptions and provenance are retained with the code.

## License and citation

Original project software is MIT-licensed. The external pressure-kernel comparison retains GPL-3.0 terms in [its license](numerical/sensitivity/LICENSE), with attribution and upstream commit in [its provenance](numerical/sensitivity/external_alvin/PROVENANCE.md). The MIT license does not relicense third-party material. See [data provenance](numerical/LICENSE_DATA.md) and [citation metadata](CITATION.cff). Cite the fixed Git commit used so the calculations can be traced to the same files.
