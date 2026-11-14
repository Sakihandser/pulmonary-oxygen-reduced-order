# Primary effective admixture model

The primary application is pulmonary hypertension mechanism analysis. No prospective patient-level or treatment-response accuracy is claimed. `equilibrium_admixture` is the stable code identifier. The four remaining laws are comparison/sensitivity models; the old controlled benchmark has a separate calibration and closure.

## State and units

Q is total blood flow in L/min; M is assigned systemic oxygen use in mL O2/min; Hb is g/dL; PA, Pa and Pv are mmHg. Oxygen contents Ce, Ca and Cv are mL O2/dL. Ce=C(PA), with the shared Hill curve. Effective admixture s_eff is a fraction in [0,1), not percent, not the regional capacity share and not a measured anatomical shunt. DO2=10 Q Ca is systemic oxygen delivery, distinct from M.

The two constraints are Ca=(1-s_eff)Ce+s_eff Cv and M=10 Q(Ca-Cv). They give:

    Cv = Ce - M/[10 Q(1-s_eff)]
    Ca = Ce - s_eff M/[10 Q(1-s_eff)]
    DO2 = 10 Q Ce - s_eff M/(1-s_eff)

Feasibility requires Q >= M/[10(1-s_eff)Ce]. Equality permits Cv=0 mathematically; it does not define a tolerable clinical state. At s_eff=0 the arterial endpoint is equilibrium regardless of Q, provided demand is feasible. At M=0, Ca=Cv=Ce. s_eff=1 is excluded because positive-demand steady balance is impossible and zero demand is underdetermined.

At fixed remaining inputs, dDO2/dQ=10Ce, dDO2/dM=-s_eff/(1-s_eff), dDO2/ds_eff=-M/(1-s_eff)^2, dDO2/dHb=10Q*1.34*S(PA), and dDO2/dPA=10Q*C'(PA). These are partial derivatives, not total treatment responses. There is no explicit PVR term. Pressure influences the plotted trajectory only through the assumed Q(PVR) relation. Common PA and global mixing also eliminate regional-flow sensitivity at fixed Q.

## Parameter provenance

| Input | Value or construction | Role and source |
| --- | --- | --- |
| Blood curve | P50 26.8 mmHg; n 2.7; binding 1.34 mL/g; dissolved 0.0031 mL/dL/mmHg | Shared approximation, not patient-specific estimation; configured centrally |
| Calibration | Carta 2022 air observation, DOI 10.3389/fmed.2022.791423 | Attributed source means; model coefficient fitted once from that baseline |
| s_eff | 0.14363059964245567 at nominal assumptions | Recomputed from calibration using (Ce-Ca)/(Ce-Cv); not derived from test hyperoxia |
| PA reference | 100.40705107350374 mmHg | Calculated from Carta air gases and PB/PH2O, R=0.8; not directly observed alveolar PO2 |
| Q and mPAP | Existing resistance network and cardiac closure; PAWP=8 mmHg | Hydraulic model outputs; closure not clinically validated |
| Hb and M reference | 14 g/dL and 250 mL/min | Declared scenario values; M is not inferred from an outcome |
| Perturbations | model_config.json | Chosen stress tests; not observed interventions or probability distributions |

## Output dictionary and reproducibility

`primary_flow_controls.csv`: PVR_WU, closure, Q_L_min and mPAP_mmHg identify each hydraulic counterfactual. `primary_counterfactuals.csv` records all five independent oxygen inputs for each named change. Shared output fields are feasible; Pa_mmHg; Sa_percent; Ca_ml_dl; Cv_ml_dl; Sv_percent; DO2_ml_min; mass_balance_error (mL/min); Q_min_L_min; equilibrium_content_ml_dl. Delta_Sa is in percentage points, while DO2_change_percent is relative percent. The local sensitivity CSV records units and perturbation size explicitly.

`primary_subgroup_metrics.csv` uses the same frozen predictions as predictions.csv, split by Table II of Dantzker and Bower 1979, DOI 10.1172/JCI109542. It contains no new fit or patient. RMSE, MAE and signed prediction-minus-observation bias are mmHg. Small subgroups are descriptive; the four idiopathic cases slightly favor saturation by RMSE.

The source data, source registry and model_config.json are protected by input_sha256.json. The root release manifest covers scripts, results and documentation. Work in a clean copy: reproduction compares saved CSVs, rewrites generated files and emits a new run summary. No web access is used after dependency installation.


## Observation limits

At known Ce, the arterial output identifies sM/[Q(1-s)]. A venous output additionally identifies M/Q but not the absolute flow–demand scale. The inverse with independent Q is s=(Ce-Ca)/(Ce-Cv), M=10Q(Ca-Cv). This requires a nondegenerate physical domain, a known blood relation and known Ce. Rank alone does not establish estimation precision; bounded-error examples and equal-output witnesses are in REVIEW_PROTOCOL.md and review_*.csv. Neither derived DO2 nor flow computed from the same oxygen observations adds an independent observation.
