# Output-specific model selection

Run `python adequacy_extension/code/reproduce_adequacy.py` from the numerical archive root. Seven result CSVs and twelve tests accompany Figure 8. No additional patient data are introduced.

Given a sharp feasible content interval [L,U], an arbitrary point predictor p has exact worst absolute error max(|p-L|,|p-U|). The midpoint is minimax with error (U-L)/2. An endpoint is adequate at tolerance delta if U-L <= delta; the midpoint is adequate if U-L <= 2 delta. If U-L > 2 delta, no point prediction based on that information can certify delta over the entire class. This is an output- and loss-specific statement, not identification of the true physiology.

The baseline-only class assumes fixed path weights, common current Ce and Cv, and content relaxation. Current Ce and consumption can change and must be supplied. The common capacity exponent may be known or lie in a declared interval. Support restrictions must be justified independently. Physiological feasibility is imposed at the queried state; the interval does not constrain an entire unobserved trajectory.

The sharp two-state result assumes a KNOWN common beta>0 and a second state at r=2**(1/beta). With y=sqrt(z), exact s0=E[y^2], s_second=E[y] determine two extremal distributions. Beta=1 means a doubling of flow, not a proposed patient challenge. Mathematical examples use r up to 3; this is not a physiological feasibility recommendation. Endpoint zero denotes an infinite-capacity limit. The midpoint may correspond to different structures at different flows; it is a pointwise predictor, not one identified mechanistic trajectory.

With second-residual uncertainty +/-0.005, the supplied noisy envelope is conservative, not sharp. Its midpoint guarantees its reported radius but is not asserted to be minimax for the smaller true moment class. Baseline and other inputs are exact in this demonstration. General measurement uncertainty needs propagation before using any guarantee.

Files: adequacy_grid.csv (1616 declared states); adequacy_examples.csv; tolerance_thresholds.csv; second_state_bounds.csv; direction_conditions.csv; empirical_model_choices.csv; output_scale_example.csv. The last file demonstrates why minimax in content is not generally minimax in oxygen pressure.
