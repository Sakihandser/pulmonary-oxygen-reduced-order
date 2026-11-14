# Flow-responsive closure

This exploratory extension was specified after the matched-map comparison. It adds no observations, parameter fits or independent validation. The question is whether the fixed-admixture compensation thresholds remain valid when the capacity–flow rule used in the prediction comparison is also imposed on the closed balance.

Use the content map K(Q)=Kstar (Q/Q0)^(1-beta), x=K/Q, s=exp(-x), with beta in [0,1]. Hold end-capillary content and metabolic demand fixed when differentiating with respect to Q. Derive the sign of the arterial-content derivative and of the delivery derivative, subject to nonnegative venous content. Retain the original constant-s results as the beta=0 special case.

Reference: Q0=5 L/min, M=250 mL/min, Hb=14 g/dL, PA=100 mmHg, s0=0.15. A separate imposed exchange deterioration sets sstar=0.20 at Q0, so Kstar=-Q0 log(sstar). Subsequent changes of Q follow beta=0, 0.25, 0.5, 0.75 or 1. These are scenarios, not calibrated physiological responses. Locate delivery and arterial-content recovery roots relative to the original reference. Search Q/Q0>=1; report infinite arterial recovery for beta=1 by the analytic monotonicity result, and distinguish roots above the displayed Q/Q0<=2 range.

Numerical checks: independent closed-balance roots; finite-difference derivatives; randomized feasible states; response reversal across the analytic sign boundary; fixed-s recovery; zero-demand and invalid-domain handling; complete flow-range feasibility. Plot arterial content, delivery and the beta sign boundary. Check analytic capacity-scaling sensitivity of the conditional predictor against finite differences at the existing source states; report individual values and flow-method strata without confidence intervals.

Claims must distinguish a re-expression of classical shunt–Fick conservation from the new comparison and its consequences. In particular, the deficit-factor ratio applies to thresholds evaluated with the same prescribed changed-state s, not to two thresholds at which s(Q) takes different values.
