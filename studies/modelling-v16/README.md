# Modelling v16 code and saved results

This English snapshot supports *Oxygen exchange model divergence under matched pulmonary vascular resistance and implications for model development*. It contains the article code, scripts, aggregate inputs, 128-condition matrix, previously recorded mathematical checks and v16 reporting summaries.

## Numerical distribution

`numerical/` contains Supplementary Material S2 in English. Numerical results derive from the v14/v15 study; category labels, directories and documentation were translated for v16. Related imports and filters were updated consistently. Numeric CSV cells were compared as exact strings and code arithmetic was compared statically after translating string constants. No model, numerical verification or research test was run during this revision.

Run `python reproduce.py` from `numerical/` after installing its requirements-lock.txt. The entry point invokes `11_diagnostic_reproduction/code/reproduce_diagnostic.py`. Production comparisons reside in controlled_exchange_frozen.py; independently formulated endpoint functions reside in independent_endpoint_functions.py. Earlier phase directories retain required inputs and dependencies under English names.

release_sha256.json lists the current distribution files. translation_provenance.json records source hashes and the language-only transformation. Historical run manifests retain original run hashes and an explanatory note; they do not describe a new run of the English distribution. Keep relative paths intact and reproduce in a disposable copy because outputs can change.

## Reporting directory

`reporting/` contains allocation contrasts, actual-use/unmet-demand tables, capacity shares, derivation provenance and saved CSV sources for Figures 1-4. Its relative-path figure export script writes 600 dpi PNG/TIFF and SVG; it reads saved results and does not run the model. Submission TIFF figures are supplied separately.

## Provenance and scope

PUBLICATION_MANIFEST.json records current file hashes and the source submission archives. Independent physiological oxygenation and PAH validation remain absent. No repository DOI, restricted data, article full text or unrelated subsequent-study participant records are included. Source cells retain study attribution; the MIT software license does not sublicense third-party articles or datasets.
