# Ranking stability and observation design

Date: 2026-11-10. Retrospective diagnostic extension after model choice. No new participant data. No preregistration claim.

## Frozen empirical ranking

Keep the Carta air coefficients and all test predictions fixed. Resample entire paired model-error vectors within the seven-row pulmonary vascular sample, the four historical idiopathic-PH rows, and the three recurrent-emboli rows. Enumerate every count vector c summing to n and weight it by n!/(product c_i! n^n); this exactly represents all n^n ordered empirical bootstrap samples. Calculate each model's RMSE, the admixture-minus-comparator RMSE difference, and weighted quantiles at 0.025 and 0.975. Split best-rank weight equally on numerical ties (absolute tolerance 1e-10). Leave-one-row deletion also freezes all coefficients and differs from the existing refitted LOPO procedure.

These descriptive distributions condition on the observed, already-used sample, its disease mixture, calibration, blood relation and inputs. They do not correct selection bias, quantify uncertainty in the Carta aggregate fit, assert population confidence coverage, or represent probabilities of physiological truth. The tiny subgroups particularly limit inference. Aggregate condition means are not resampled as patients. Original and gas-equation results are not alternative independent validation cohorts.

## Outcome-free alveolar computation

Use FiO2(PB-PH2O)-PaCO2[FiO2+(1-FiO2)/R] for each Dantzker row. This computational path never reads test PaO2 or A-a gradient; PB and water pressure remain assumed source conventions, not independently measured alveolar oxygen. Compare against the original PA reconstructed from Pa+A-a without overwriting it. Cross P50=24,26.8,30 mmHg, Hill n=2.5,2.7,3, and R=.75,.8,.9. Refit only Carta air in each scenario, then freeze all model coefficients. Report group-specific RMSE and MAE. The calculation removes an outcome-bearing arithmetic path, but does not prove the original reconstruction was invalid: the original A-a gradient may itself have been calculated using another gas convention.

## Observation map

For fixed known Ce=C(PA), Ca=Ce-sM/[10Q(1-s)], Cv=Ce-M/[10Q(1-s)]. The unknowns are Q, M and s; use log Q, log M and logit s for the interior local Jacobian. Scale content rows by Ce and use log Q as the independently observed flow row. At positive M,Q and 0<s<1, arterial-only rank is one; arterial plus Q or arterial plus venous content has rank two; adding independent Q to both contents gives rank three. This is structural information, not estimation precision. The result requires fixed known blood/gas quantities; unknown Ce or blood-curve parameters introduce further ambiguity.

Exact inversion is s=(Ce-Ca)/(Ce-Cv), M=10Q(Ca-Cv). An independent flow observation is essential to determine the scale of M and delivery. Rescaling Q and M together leaves both contents unchanged. A derived DO2=10QCa is not an additional independent observation. Q derived from the same oxygen observations and an assumed M does not add independent rank.

At the published model reference, vary Ca,Cv,Ce independently by bounded errors of ±.05, ±.1, ±.2 mL/dL and Q by ±10%. Enumerate all 16 corners for each width. Under 0<=Cv<=Ca<=Ce these functions are coordinatewise monotone, so corner extrema bound each output over the specified box. Check all corners for feasibility rather than discard them. Bounds are illustrative metrology scenarios; no probability, observed device accuracy or population coverage is assigned. Unknown covariance would change attainable bounds.

## Future independent validation

Task A: predict Pa using contemporaneous observed Pv, Hb, gas inputs and independent flow; this tests only the conditional exchange map. Task B: predict both Pa and Pv using independent Q, independently measured oxygen consumption, Hb and gas inputs; this tests the closed oxygen balance. Do not compute consumption from the same held-out Pa/Pv and then label Task B a prediction. Test the cardiac Q(PVR) closure separately with held-out measured flow.

Freeze model coefficients, source conventions, endpoints, challenge strata and exclusion rules before an unused modern PAH cohort is accessed. Retain CTEPH as a distinct stratum. Keep repeated states paired by participant. Report participant-level RMSE/MAE/bias, paired model contrasts, calibration and feasible-state failures; bootstrap participants with all their states, not individual conditions. Report normoxia/hypoxia/exercise separately from hyperoxia. Do not tune R, blood curves or s on the held-out outcomes. Compare subject-specific baseline calibration only as a separately declared protocol, fitted before each participant's held-out states.

First establish repeatability and input-error covariance in a pilot, including pH, temperature, co-oximetry, measured gas conditions, independent flow and directly measured consumption. Determine sample size by participant-level simulation for a prospectively chosen meaningful error contrast and interval width. Neither a sample-size claim nor treatment-effect endpoint is justified by the seven historical cases. Ethics, consent/data access and clinical feasibility must be settled by investigators before collection; this package does not claim those steps occurred.

## Files and verification

Run code/review_extension.py and code/review_figures.py, or the root reproduce.py. Seven new CSVs record bootstrap distributions (summary), deletions, both gas-boundary grids, nominal gas inputs, equivalent-output witnesses, observation ranks and error bounds. review_summary.json records assumptions. Tests independently enumerate n=3 ordered samples, finite-difference the observation Jacobian, test exact inversion/scale symmetry, enforce input-domain errors, and verify that gas calculation ignores the arterial outcome.
