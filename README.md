# Pulmonary oxygen reduced-order model

Code and saved numerical results for **Calibration-Dependent Transferability of Pulmonary Oxygen Exchange Models in a Controlled Numerical Benchmark**, by Mingze Hong and Bojun Yang. The current study snapshot is **Modelling v20**.

## Current study and reproduction

Start in [studies/modelling-v20](studies/modelling-v20/README.md). Its `numerical/` and `reporting/` directories are byte-preserving extractions of the v20 S2 and S3 submission supplements. The paper should cite this directory at a fixed Git commit, rather than relying on the moving main branch.

```text
cd studies/modelling-v20/numerical
python -m pip install -r requirements-lock.txt
python reproduce.py
```

Run in a fresh working copy because reproduction writes outputs. After dependency installation the calculations run offline. The entry point verifies supplied SHA256 files, regenerates 53 numerical result CSVs, checks numerical main Tables 4–8, and exports seven main figures and two supplementary figures. See [RESULT_REPRODUCTION_INDEX.md](studies/modelling-v20/numerical/RESULT_REPRODUCTION_INDEX.md) for individual scripts, CSVs and commands. The requirements lock records the tested package versions; it does not guarantee availability on every platform.

For figures only:

```text
cd studies/modelling-v20/reporting
python reproduce_figures.py
```

For tests, from `studies/modelling-v20/numerical`:

```text
python -m unittest discover -s tests
python -m unittest discover -s aggregate_analysis/tests
python -m unittest discover -s observation_analysis/tests
```

## What v20 adds

Three exchange laws are matched to a common endpoint or a common inlet content slope; fitted coefficients are frozen before changing states. The results include the 441-state inlet–flow grid for each matching rule, signed state maps with three illustrative output-error bounds, and a dissolved-oxygen ablation explaining the affine limit of saturation/content relaxation. Regional conductance and allocation are studied in 1,925 deterministic states; the historical 128-condition comparison remains available separately.

A clean v20 reproduction recorded 53 CSVs and 76,776 numeric cells with maximum difference zero, exported all nine figures, and passed 33 unit/regression tests. The controlled matrix also reports 17 mathematical checks. These are implementation and reproducibility checks, not independent physiological validation. Run records and exact versions are retained in the numerical supplement.

## Scope and evidence

The states are synthetic. Grid counts under an assumed bounded error are not clinical success probabilities. The CT inspiratory air-volume shares are allocation proxies, not measured capillary blood volumes. Gas-derived capacity uses unpaired aggregate observations with unresolved normalization. The external Alvin pressure submodule supplies limited numerical qualification, not a validation of the full published model or of patient physiology. No identifiable participant records, restricted source data, third-party article full text, internal reviews or manuscript-authoring scripts are included.

## Earlier snapshots

[Modelling v16](studies/modelling-v16/README.md) and the earlier root `code/`, `data/` and `tests/` are retained for traceability. They are historical entry points and do not reproduce the current v20 results. The addition of v20 does not overwrite the v16 study.

## License and citation

Original project software uses the root MIT license. The adapted external pressure-kernel comparison retains GPL-3.0 terms under [numerical/sensitivity/LICENSE](studies/modelling-v20/numerical/sensitivity/LICENSE), with upstream source, commit and attribution in its provenance record. The root MIT license does not relicense that material or third-party data. See [LICENSE_DATA.md](studies/modelling-v20/numerical/LICENSE_DATA.md) and [CITATION.cff](CITATION.cff). Cite the exact Git commit used. No repository DOI is claimed.
