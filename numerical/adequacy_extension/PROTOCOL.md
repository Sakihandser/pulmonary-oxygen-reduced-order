# Output-specific model adequacy

Date: 2026-11-14. Exploratory mathematical extension designed after reading the prior observations, counterexamples and review. This is not clinical preregistration.

Question: which reduced predictor is sufficient for a declared output error, and what additional measurement reduces structural uncertainty?

1. Reuse the parallel content-relaxation class. Derive attainable output intervals, the exact maximum error of any point prediction, and the absolute-error minimax prediction. The tolerance is user-specified, not a clinical cutoff.
2. Allow current end-capillary content and consumption to differ from baseline. They are supplied inputs, not inferred from PVR. Keep path weights fixed. Unknown common exponent is handled by a union over its declared interval, not by selecting its best value retrospectively.
3. Derive a sharp two-state special case: known exponent and a second state with exponent factor one half. With y=sqrt(z), the two measured residuals are E[y] and E[y^2]. Use extremal two-point distributions to bound all subsequent or intervening states. Endpoint zero represents an infinite-capacity limit. No arbitrary finite path distribution is claimed to have been identified.
4. For an uncertain second residual use a conservative Jensen envelope valid over the full stated interval; do not label the noisy envelope sharp. Do not confuse the residual's measurement interval with a confidence interval.
5. Derive the general oxygen-content response identity with varying demand and end-capillary content and a sufficient direction condition for heterogeneous path exponents. This is an application of established balances, not a new physiological law.
6. Numerical cases: baseline residuals .05, .10, .20, .40; common beta 0, .25, .5, 1; r from 1 to 2; supplied Q0=5, Ce=18.548851571 and M=250 unless stated otherwise. Illustrative tolerances .1 and .2 mL/dL. Two-state demonstrations use the two existing finite models without refitting. Inspect sensitivity at +/- .005 residual and beta interval [0,1].
7. Check analytic bounds against independent finite-grid linear programming, explicit attaining distributions, random mixtures, and independently solved Fick balance. These are implementation checks, not new observations. Retain all source observations and previous outputs unchanged.
