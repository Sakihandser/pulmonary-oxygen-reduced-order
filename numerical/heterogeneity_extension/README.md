# Baseline constrained heterogeneity and recovery bounds

This extension tests the interpretation of the single-exponential flow criterion. It introduces no patient data. Known flow-diffusion and concavity results are used to bound compensation over unresolved exchange-path distributions.

Run `python heterogeneity_extension/code/reproduce_heterogeneity.py` from the full numerical archive. This regenerates five result CSVs and Figure 7, compares 5262 finite numerical cells, and runs ten tests. `PROTOCOL.md` records choices made before calculation.

`matched_models.csv` specifies every path of two finite models with the same baseline effective fraction and the same total path capacity. All path capacities and flow weights remain fixed as total flow increases. They nevertheless have opposite initial arterial-content responses. `bounded_recovery.csv` searches for first upward recovery only in Q/Q0 in [1,2]; NaN denotes no crossing in that interval, not global impossibility. The ideal limit has infinite capacity and is used only as a sharp limiting bound.

`sharp_bounds.csv` gives attainable or limiting bounds conditional on the baseline mean, a common beta, and the stated support of path residuals. These bounds do not additionally constrain total capacity. The support restrictions are illustrative, not estimated physiology. Classification refers to the specified flow, not every lower flow. All bounds are intersected with nonnegative venous content. They are not confidence intervals.

`local_slope_bounds.csv` and `heterogeneous_response.csv` distinguish the response of effective aggregate capacity from the fixed sum of path capacities. Numerical checks use independent venous-balance roots, finite differences, bound-attaining distributions, finite limiting sequences and deterministic random mixtures. They verify mathematics and implementation, not physiological validity.
