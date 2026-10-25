# Modelling v20 code and numerical results

This snapshot supports **Calibration-Dependent Transferability of Pulmonary Oxygen Exchange Models in a Controlled Numerical Benchmark**. It includes the code, aggregate inputs, numerical outputs and figure sources used for the v20 paper.

`numerical/` exactly reproduces the distributed contents of S2; `reporting/` exactly reproduces S3. Their original `release_sha256.json` files are preserved and verified. `PUBLICATION_MANIFEST.json` records the SHA256 and Git blob hashes of all published study files and the two source archive hashes. Manuscripts, TIFF/PNG exports, internal reviews, working files and unpublished subsequent studies are not included; SVG sources and exporters regenerate the images.

## Commands

From `numerical/`, install `requirements-lock.txt` and run `python reproduce.py`. This entry point verifies file hashes before regenerating numerical results and figures. Use a separate working copy. From `reporting/`, `python reproduce_figures.py` regenerates nine figures from saved CSV inputs. See `numerical/RESULT_REPRODUCTION_INDEX.md` for the result-to-code map and the root README for all three test-discovery commands.

The current main result is common calibration followed by frozen prediction: Table 4 and Figure 2. `common_calibration/code/run.py` generates shared-target fits, the inlet–flow comparison and regional allocation scans. `common_calibration/code/mechanism.py` generates the dissolved-oxygen controls and signed state maps. Earlier convention-based comparisons remain distinct and should not be pooled with the matched fits.

## Verification and limits

A clean v20 run regenerated 53 numerical result CSVs and compared 76,776 numeric cells with maximum difference zero. It exported seven main and two supplementary figures. The 33 unit/regression tests and the controlled matrix's 17 numerical checks passed. These records do not establish patient accuracy, a clinical discrimination probability, or cross-platform reproducibility beyond the recorded environment. No clinical observations or additional external physiological validation were added.

Original project software is MIT-licensed; the external pressure-kernel comparison and source retain GPL-3.0 terms in `numerical/sensitivity/LICENSE`. The source author, fixed upstream commit and limits are recorded in `numerical/sensitivity/external_alvin/PROVENANCE.md`. Aggregate source cells retain attribution under `numerical/LICENSE_DATA.md`.
