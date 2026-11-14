# Exchange maps, flow scaling and compensation

The matched comparison separates four exchange maps from two scaling rules. Every cell uses one baseline coefficient. Proportional content capacity is fixed effective admixture; proportional saturation capacity is a separate control. In the 11 nonhyperoxic changes their equal-family Pa RMSEs are 2.963 and 2.845 mmHg. The original admixture advantage does not uniquely identify the content map. The scaling difference is concentrated in historical high-flow changes; modern thermodilution states slightly favor fixed capacity.

At fixed end-capillary content Ce, the closed content balance gives QD/Q0 - 1 = epsilon (QA/Q0 - 1), where epsilon=(Ce-Ca0)/Ce. This holds for simultaneous changes in demand and admixture, subject to feasibility. The unchanged-Ce rule does not apply if Ce changes; use the full thresholds then.

Run `python review_extension/code/reproduce_review.py` from the archive root after installing `requirements-lock.txt`. This checks 10 saved result CSVs, nine tests and redraws final main Figures 1–2. The root `python reproduce.py` also executes all supporting analyses. The historical selection figures are retained as executable benchmark outputs; final submitted Figures 1–2 come from this extension, and Figures 3–5 from selection_extension.

`PROTOCOL.md` records the exploratory design; `FLOW_PROVENANCE.md` distinguishes thermodilution, direct Fick and unresolved historical methods. `results/matched_predictions.csv` contains all condition residuals; `matched_metrics.csv` includes modern, historical and measurement-method strata; `matched_sensitivity.csv` and `matched_winners.csv` retain all 108 common assumptions. Flow, curve and boundary stress outputs are deterministic scenarios, not confidence intervals or patient data. `dimensionless_compensation.csv` verifies the reduced compensation identity across 81 cases.

No independent cohort is added by these computations. The current flow cancels in all beta=0 controls; models consuming direct-Fick flow retain a dependency on arterial oxygen. Numerical tests establish implementation identities, not clinical validity.
