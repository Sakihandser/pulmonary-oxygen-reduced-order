# Compatibility of pulmonary oxygen model assumptions

Exploratory protocol written after strict review of the existing model claims and an initial nominal diagnostic. This is not preregistered and is not an independent validation dataset.

Use all existing primary paired PH records: eleven nonhyperoxic changes and one separate pure-oxygen change. Keep the derived-venous exercise series secondary. Do not add participants, infer individual covariance, or treat CTEPH as PAH.

1. Reconstruct the effective content residual s=(Ce-Ca)/(Ce-Cv), using the inherited common Hill relation and gas convention. This is an outcome-based diagnostic, not a predictor.
2. For any positive flow ratio r, fixed weights, common end-capillary content and path exponents beta_i in [0,1] imply s1 in [min(s0,s0^(1/r)), max(s0,s0^(1/r))]. Endpoints allow limiting zero or infinite path capacities. This is a necessary two-state condition, not a full trajectory or mechanism test.
3. Compute distances to compatibility: minimum common absolute residual halfwidth, and minimum common arterial-pressure halfwidth with all other inputs fixed. These are deterministic sensitivity margins, not confidence intervals or clinical measurement tolerances.
4. Repeat the existing 108 shared-input scenarios (P50, Hill exponent, R, missing Hb), changed-state P50 shifts, reported-saturation conversion where present, and an explicitly illustrative +/-10% flow-ratio interval. Do not call a count of scenarios a probability.
5. Quantify the conditional Pa separation between beta=0 and beta=1. Half the separation is only the equal-error budget for disjoint predictions with all other inputs exact. Actual error includes shared gas, baseline and structural uncertainty.
6. Document which required structural quantities are actually observed. A pure PAH multi-state independent-flow cohort with full gas inputs is not supplied by the current collection.

Use the compatibility audit to evaluate whether v29 certificates can be applied. A failed nominal audit can reflect model failure, gas-conversion assumptions, rounding or aggregation. It does not locate the physiological cause. A successful audit does not identify beta or certify the entire structure.
