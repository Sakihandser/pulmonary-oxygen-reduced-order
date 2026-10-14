# Pulmonary oxygen reduced-order model

Code and saved results accompanying the manuscript **Oxygen exchange model divergence under matched pulmonary vascular resistance and implications for model development**, by Mingze Hong and Bojun Yang. The current manuscript is an original methodological research Article prepared for **Modelling**, version 16; it has not been formally submitted or accepted.

The study compares three exchange laws and four constructions under matched total pulmonary vascular resistance. The diagnostic matrix contains 128 conditions crossing construction, regional capacity allocation, hemoglobin, inlet feedback, synthetic post-exchange admixture and perfusion redistribution. It distinguishes changes in exchange variables, endpoint matching and gas-based capacity calibration. This is not a complete exchange-law by calibration factorial design, and the 128 conditions are not patient observations.

## Current manuscript materials

Start in [`studies/modelling-v16/`](studies/modelling-v16/README.md).

- [`numerical/`](studies/modelling-v16/numerical/): exact extracted contents of the frozen numerical supplement, including production code, independent endpoint functions, archived mathematical checks, source-cell provenance, inputs, outputs, and the offline runner.
- [`reporting/`](studies/modelling-v16/reporting/): v16 allocation and systemic-use tables, retained figure-source CSV files, provenance and the portable 600 dpi figure export script.
- [`PUBLICATION_MANIFEST.json`](studies/modelling-v16/PUBLICATION_MANIFEST.json): byte hashes and archive origins for the published snapshot.
- [`CITATION.cff`](CITATION.cff): citation metadata for this code and data snapshot. No manuscript DOI or journal publication is claimed.

## Reproduction

Use a separate working copy because the numerical runner writes outputs. From the repository root:

```text
cd studies/modelling-v16/numerical
python -m pip install -r requirements-lock.txt
python reproduce.py
```

The runner checks the frozen archive hashes, recomputes the original production matrix and independently checks endpoints. The exact locked dependencies preserve the previously recorded environment; this publication update did not rerun the models, tests or numerical verification, and did not establish a fresh-environment reproduction.

For figures only, from the repository root:

```text
cd studies/modelling-v16/reporting
python tools/export_figures_v16.py
```

This reads saved CSV results and writes `Figures_v16/` without recomputing models. Arial is used; font substitution may alter rendered pixels.

## Interpretation and limitations

Matched flow follows a prescribed common pressure-flow closure. Total capillary volume 86 mL and relaxation time 0.30 s are diagnostic assumptions, not independent patient estimates. Oxygen content is mixed by flow; synthetic admixture acts after exchange and does not resolve a separate anatomical bypass. Regional allocation effects depend on capacity reserve. The systemic feedback reports actual use subject to assigned extraction ceilings, so unmet demand is a model diagnostic rather than measured organ hypoxia.

Independent physiological oxygenation and PAH validation remain absent. The future research agenda is not a set of completed outcomes. The repository contains no identifiable participant records, restricted imaging, third-party article full text or PVDOMICS materials. Separate later public-data and coupling studies are not included in this manuscript snapshot or claimed as external validation.

## Earlier studies

The original `code/`, `data/`, `tests/` and root `requirements.txt` remain for historical continuity. They belong to earlier manuscript development and are not the current v16 reproduction entry point. Current numerical files retain their original Chinese phase-directory names and internal v14/v15 labels to preserve relative imports, bytes and provenance. Historical descriptions of archives as unpublished record their preparation state; the current public snapshot is the directory documented above.

## License and citation

Project software and original technical documentation retain the MIT license in [`LICENSE`](LICENSE). Published aggregate facts retain study and source-cell attribution; the software license does not sublicense third-party articles or datasets. See the numerical supplement's data provenance notes. Cite the repository together with the fixed commit used; GitHub publication does not itself provide a repository DOI. See CITATION.cff for authors and software version.
