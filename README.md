# Decision Consistency and Mechanistic Ambiguity in Pulmonary Oxygen Models

Research software, numerical data and figure sources for the study by Mingze Hong and Bojun Yang. The current snapshot accompanies the v31 manuscript. It examines when pulmonary oxygen models that fit the same observations agree on a specified decision and when different mechanisms remain compatible with the data.

## Reproduce the numerical results

From a fresh checkout:

```sh
cd numerical
python -m pip install -r requirements-lock.txt
python reproduce.py
```

The runner verifies the release SHA-256 manifest, compares 142 result CSVs, regenerates 25 figures and runs 131 tests across 11 suites. It runs offline after dependency installation and writes generated outputs to the working copy. The [result-to-code index](numerical/RESULT_REPRODUCTION_INDEX.md) links results to their inputs and scripts. The [multi-state analysis protocol](numerical/certainty_extension/PROTOCOL.md) describes the common-model constraints, finite counterexamples and capacity bounds added in this snapshot.

## Regenerate the figures

From a fresh checkout:

```sh
cd figure_sources
python -m pip install -r requirements-lock.txt
python reproduce_final_figures.py
```

This regenerates 25 submission TIFFs in `submission_figures/` from the frozen plotting data. The directory also contains editable SVGs and CSVs for all 12 main tables. See the [figure and table source guide](figure_sources/README_S3.md).

## Data, licensing and citation

The observed-state audit uses attributed published data. Its decision conclusions are conditional on the specified observations, model class and error bounds. See [data provenance](numerical/LICENSE_DATA.md) for the sources and redistribution scope.

Original project software is MIT-licensed. The external pressure-kernel comparison retains its [GPL-3.0 license](numerical/sensitivity/LICENSE) and [upstream attribution](numerical/sensitivity/external_alvin/PROVENANCE.md). When citing the software, include the specific Git commit and use the [citation metadata](CITATION.cff).
