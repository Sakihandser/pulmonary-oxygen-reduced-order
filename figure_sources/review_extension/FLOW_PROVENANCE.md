# Flow measurement and independence audit

Checked against full source methods on 2026-11-14. Corrects the overbroad sentence that all modern flow values were direct Fick. The distinction does not change any transcribed number.

| Source | Method | Location | Interpretation |
| --- | --- | --- | --- |
| Carta 2022 | Thermodilution | Assessments, right-heart catheterization; https://www.frontiersin.org/journals/medicine/articles/10.3389/fmed.2022.791423/full | Q is not algebraically derived from the arterial oxygen endpoint; published condition means still lack joint individual data. |
| Lichtblau 2021 | Triplicate cold-saline thermodilution; CI=CO/BSA | Right heart catheterization; https://www.frontiersin.org/journals/medicine/articles/10.3389/fmed.2021.681473/full | Same independence distinction. Grouped with Carta conservatively. |
| Haefliger 2025 | Direct Fick | Section 2.4.4, Eq. 3; https://www.mdpi.com/2077-0383/14/23/8407 | Q=VO2/(CaO2-CvO2), with appropriate unit conversion. Current arterial oxygen is an ancestor of Q. No independently measured Q is available in the extracted comparison. |
| Lichtblau 2025 | Direct Fick | https://pmc.ncbi.nlm.nih.gov/articles/PMC12005039/ | Secondary derived-Pv exercise set, excluded from this extension. |
| Melot 1983 | Independence not established from available method description | https://doi.org/10.1378/chest.83.2.203 | Do not label the historical states as independent-flow validation. |

The beta=0 controls cancel current flow algebraically and avoid the Fick-Q input channel. Baseline oxygen is intentionally used for calibration, and current Pv and PaCO2 remain measured boundary conditions. This does not make the analysis prospective or establish causal independence of all physiological inputs. At beta=1, and for an adaptive beta, computational exclusion of the current Pa field is insufficient to establish observational independence if Q is direct Fick.

The early thermodilution family's four nonhyperoxic changes have Q/Q0 between 0.9667 and 1.0000, providing little flow leverage. The modern direct-Fick family spans 1.0000–1.0851. The four historical changes span 1.4000–1.7368. A pooled scaling contrast is consequently not a strong modern independent-flow validation.
