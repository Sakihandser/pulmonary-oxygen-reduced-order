# Pulmonary hypertension and arterial oxygenation
Run `python oxygen_decisions/verify.py` from a fresh extraction after installing requirements-lock.txt. The root reproduction command regenerates the upstream inputs first. Source observations and conversions are in applicability_extension. Shared model bounds, gas boxes and finite witnesses come from certainty_extension and measurement_extension.

The observational calculation crosses residual halfwidths 0, 0.005, 0.01 and 0.02 with relative flow-ratio halfwidths 0, 0.05 and 0.10. Repeated visits and cohort means retain their source labels. observed.csv reports the residual discrepancy for all combinations and both available gas conversions. The 0.01/0.10 setting is an illustrative sensitivity scenario.

margin_grid.csv reports conservative lower margins under independently bounded observation, demand and target changes at exact query flow. witnesses.csv reports feasible query predictions. The construction and original node weights in measurement_extension allow direct reconstruction of their contents. redistribution.csv holds the two-path allocation example. allocation_bound.csv reports the monotonicity and total-variation upper bound under the additional joint query-input box in Supplement V. A nonpositive exclusion margin labels loss of this bound. The nominal B structure witnesses nonemptiness.

Figure 12 is regenerated in SVG and 600 dpi TIFF. The finite random allocation checks supplement the analytic derivative proof. Outputs characterize conditional model properties and exploratory published PH aggregates.


pah_exercise.csv transcribes published stage medians from Tobita et al., Table 2 (https://doi.org/10.1161/JAHA.122.026890). The three rows are summaries of exercise in PAH. pah_exercise_budget.csv divides the baseline arteriovenous oxygen-content difference by each stage value. This gives the maximum ratio of residual odds consistent with maintaining baseline arterial content under a shared endpoint. Group medians cannot identify individual responses.

Three external result files report model checks in mmHg, percentage points, and L/min. See EXTERNAL_VALIDATION.md for all source tables, model comparators, numerical conventions and independence limitations.
