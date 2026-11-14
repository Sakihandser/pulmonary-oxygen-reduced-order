# Testing the model class against observed PH states

Run `python code/reproduce_applicability.py` from this directory, or the root `python reproduce.py`.

The nominal reconstruction places 10 of 11 nonhyperoxic paired changes outside the fixed-weight, common-endpoint, beta_i in [0,1] two-state envelope. This is a diagnostic at published aggregate values and two historical individuals, not a patient-level rejection test. A shared arterial-pressure halfwidth of 0.151 to 2.848 mmHg can restore compatibility in the ten cases when other inputs are exact. The alternative reported-saturation convention also fails nominal compatibility for all seven modern nonhyperoxic changes. Neither convention establishes the true patient residual.

Nine of the eleven nonhyperoxic pairs remain incompatible in all 108 shared-input scenarios, one becomes compatible in 3/108, and the historical M02 pair remains compatible in all 108. Scenarios are not probabilities. Changed-state curve shifts and flow-ratio uncertainty are reported separately. All seven CSVs are included, rather than only favorable sensitivity results.

For the four independently measured-flow nonhyperoxic changes, fixed versus proportional content-capacity predictions differ by 0, 1.0299, 0 and 0 mmHg. These conditions cannot establish a robust capacity-response ranking. The half-separation is an idealized discrimination budget, not demonstrated measurement precision.

Practical use: audit compatibility before applying an internal error certificate; avoid inferring pulmonary recruitment from a fitted effective coefficient; require independent repeated-flow gas measurements to test transportability. Failing the audit does not identify whether gas conversion, aggregation, redistribution, endpoint changes, or the capacity law is responsible. Passing it does not validate a physiological model. The underlying interval and minimax mathematics are elementary or established; the contribution is the reproducible application audit, with limited external validity.

The older second-state example remains illustrative. A special flow ratio of 2^(1/beta) and an assumed residual error do not constitute a feasible clinical design. The new compatibility diagnostic accepts any positive flow ratio but does not solve general multi-state optimal experimental design.
