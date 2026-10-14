# Modelling v16 code and saved results

This snapshot supports *Oxygen exchange model divergence under matched pulmonary vascular resistance and implications for model development*. It publishes the paper's code, scripts, aggregate inputs, 128-condition matrix, previously recorded mathematical checks and v16 reporting summaries. It does not publish the authors' internal reviews, Word-building scripts or unrelated subsequent studies.

## Numerical directory

`numerical/` is the byte-preserving extraction of Supplementary Material S2, whose internal version is v15 and numerical source is v14. Its original `release_sha256.json`, environment lock, provenance notes and archived validation-status fields are unchanged historical records. All original files, including historical figure exports, are retained so `reproduce.py` can verify the original manifest before recomputation. Run from this directory, not from the repository root. Do not replace the original relative phase paths.

The entry point is `reproduce.py`; it invokes `11_第八轮_公开数据与守恒耦合/code/reproduce_diagnostic.py`. Production comparisons reside in `controlled_exchange_frozen.py`, and independently formulated endpoint functions reside in `independent_endpoint_functions.py`. Their dependencies and hydraulic, gas-mapping and source-cell inputs are preserved in the preceding phase directories. See root README for commands. These commands have not been rerun during the v16 preparation/publication update.

## Reporting directory

`reporting/` contains allocation contrasts, actual-use/unmet-demand tables, assigned capacity shares, derivation provenance and CSV sources for Figures 1–4. `tools/export_figures_v16.py` uses paths relative to its extracted reporting directory and exports 600 dpi figures. Its only publication change is portable file location plus a corrected resolution comment; exchange algorithms and archived numerical results were not changed. PNG/TIFF outputs are generated on demand and are not duplicated here.

## Traceability and status

`PUBLICATION_MANIFEST.json` maps every file to its source archive and SHA256. It also records hashes of the original S2 and revised S3 submission ZIPs. Snapshot publication and file-integrity checks do not constitute new numerical or physiological validation. No repository DOI has been assigned. Source cells are attributed aggregate facts; see LICENSE and the numerical data provenance notes for license boundaries.
