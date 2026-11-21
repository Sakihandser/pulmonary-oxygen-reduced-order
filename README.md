# Arterial Oxygenation in Pulmonary Hypertension: Decision Margins and Mechanistic Ambiguity in Reduced Models

Numerical analysis, attributed study summaries and figure sources for the manuscript by Mingze Hong and Bojun Yang. This release adds external exercise and drug challenges to tests of arterial oxygenation, metabolic demand and perfusion redistribution. The new studies, conditional predictions, comparison rules and failed comparisons are documented in [External validation](numerical/oxygen_decisions/EXTERNAL_VALIDATION.md).

## Reproduce the numerical results

From a fresh checkout:

```sh
cd numerical
python -m pip install -r requirements-lock.txt
python reproduce.py
```

The runner verifies the release SHA-256 manifest, compares 146 prior result CSVs plus nine oxygen-analysis outputs, regenerates the figures and runs 137 unit tests plus additional regression and allocation-bound checks. It runs offline after dependency installation and writes generated outputs to the working copy. The [result-to-code index](numerical/RESULT_REPRODUCTION_INDEX.md) links results to their inputs and scripts. The [multi-state analysis protocol](numerical/certainty_extension/PROTOCOL.md) describes the common-model constraints, finite counterexamples and capacity bounds underlying the synthetic examples.

## Regenerate the figures

From a fresh checkout:

```sh
cd figure_sources
python -m pip install -r requirements-lock.txt
python reproduce_final_figures.py
```

This regenerates 26 submission TIFFs at 600 dpi in `submission_figures/` from the frozen plotting data. The directory also contains editable SVGs and CSVs for all 12 main tables. See the [figure and table source guide](figure_sources/README_S3.md).

## Data, licensing and citation

The observed-state audit uses attributed published data. Its decision conclusions are conditional on the specified observations, model class and error bounds. See [data provenance](numerical/LICENSE_DATA.md) for the sources and redistribution scope.

Original project software is MIT-licensed. The external pressure-kernel comparison retains its [GPL-3.0 license](numerical/sensitivity/LICENSE) and [upstream attribution](numerical/sensitivity/external_alvin/PROVENANCE.md). When citing the software, include the specific Git commit and use the [citation metadata](CITATION.cff).
