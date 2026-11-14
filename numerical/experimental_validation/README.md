# Published observations and conditional oxygen predictions

This exploratory extension compares five one-parameter effective oxygenation models. It adds **22 quantitative records from five studies**, plus six contextual conditions from two other studies. A record is a published patient row or a condition mean, not an independent participant. No participant total is obtained by summing the repeated cohort sizes. These data cannot identify a universally accurate microscopic oxygen-exchange law.

Run from a clean copy of the numerical supplement with the root requirements installed:

```text
python experimental_validation/code/reproduce_empirical.py
```

The root `python reproduce.py` also executes this extension after the original benchmark. Neither reproduction entry point downloads data. `input_sha256.json` protects the transcriptions; `make_data.py` exposes the original units, conversions, assumed values and post-treatment arithmetic. `data/source_registry.csv` identifies the primary sources and tables. Full texts are not redistributed.

| Result | Output | Calculation |
| --- | --- | --- |
| Frozen baseline, excluded conditions, individual cross-validation and external transfer | results/predictions.csv, metrics.csv | code/validate.py |
| Blood curve, alveolar and missing-Hb sensitivities | results/assumption_sensitivity.csv, external_sensitivity.csv | code/validate.py |
| Baseline rounding and measured saturation discrepancy | results/baseline_rounding.csv, blood_curve_audit.csv | code/validate.py |
| Fixed-D recruitment and measured mixing controls | results/recruitment_sensitivity.csv, measured_mixing_controls.csv | code/couple_network.py |
| Closed oxygen balance in the existing vascular network | results/network_scenarios.csv, network_summary.json | code/couple_network.py |
| Mixing uncertainty | results/network_admixture_sensitivity.csv | code/couple_network.py |
| Figures 8–11 | figures/*.png, *.tif, *.svg | code/figures.py |

`PROTOCOL.md` specifies calibration/test separation. Published mixed-venous oxygen, flow and alveolar inputs are supplied at each tested condition: this is a conditional arterial endpoint prediction, not a prediction of the entire response from baseline alone. Candidate selection and performance assessment use the same small external cohort; an additional untouched cohort is still required to estimate the selected model's generalization without selection bias.

For the pulmonary vascular application, effective admixture is the adopted primary whole-lung closure. Its parameter also absorbs unresolved V/Q mismatch and blood-curve error; it is not an anatomical shunt measurement. The independent seven-person transfer includes four idiopathic-PH and three embolic cases. Modern idiopathic-PAH-specific individual data remain needed. COPD and fibrosis results are not PAH treatment evidence. The inherited cardiac pressure-flow closure, common alveolar pressure and fixed mixing are assumptions, not validated disease trajectories.

No posterior uncertainty, p-value, patient-level confidence interval or clinical decision rule is estimated. The numerical precision of CSVs is computational precision, not measurement precision. The original physiological model and all historical results remain available separately.


## Primary pulmonary hypertension model

`model_config.json` is the shared configuration for the empirical blood curve, calibration source and primary scenarios. `code/primary_model.py` supplies the exact oxygen-balance solution; `code/couple_network.py` uses that same implementation for the network. `code/analyze_primary.py` adds four result CSVs and an explanatory summary. The hydraulic solver and original benchmark remain unchanged. The empirical admixture parameter is fitted from Carta air on each run; it is never copied from the reported hyperoxic shunt.

| Result | Reproducible source | Manuscript location |
| --- | --- | --- |
| Frozen seven-person transfer split by published diagnosis | primary_subgroup_metrics.csv | Results 3.8, Table S30 |
| Inherited flow response versus fixed-flow control | primary_flow_controls.csv | Results 3.10, Figure 11A–B |
| Separate changes to flow, Hb, alveolar PO2, admixture and demand | primary_counterfactuals.csv | Results 3.10, Figure 11C–D, Table S31 |
| Exact local oxygen-delivery derivatives | primary_local_sensitivity.csv | Supplement K2 |
| Parameter roles and provenance | MODEL_SPECIFICATION.md, model_config.json | Supplement K3, Table S32 |

All four CSVs are in `results/`. Subgroup descriptions and one-factor scenarios were added after the initial comparison. They add no new patients and do not provide independent confirmation of the selected model. On the four historical idiopathic-PH rows, saturation and mixing RMSE are 13.79 and 13.96 mmHg; mixing is therefore **not** established as the most accurate pure-PAH model. It remains the adopted primary model for mechanistic analysis, based on the seven-person aggregate result, interpretability and exact coupling. Do not merge the two subgroup sample sizes with the seven-person total.

Run all 19 extension tests with `python -m unittest discover -s experimental_validation/tests -v`. The six added tests independently check the exact balance against a numerical root, the feasibility boundary, limiting cases, input domains, analytical derivatives and fixed-flow invariance.


## Ranking stability and measurement design

The seven-person nominal advantage is not stable to all changes in case composition: effective admixture ranks first in 76.1% of exact fixed-fit empirical resamples; deleting D06 changes the nominal winner. The four-case historical PH ranking also changes between source-reconstructed and gas-equation PA. These post-selection diagnostics are not additional validation cohorts.

`REVIEW_PROTOCOL.md` specifies exact count-vector enumeration, a gas-equation path that does not read held-out arterial oxygen, exact same-output witnesses, dimensionless observation ranks and bounded-error recovery. `code/review_extension.py` produces seven CSVs and `code/review_figures.py` produces Figure 12. Six further tests bring the empirical suite to 19 tests. `reproduce_empirical.py` compares 23 CSVs. Unknown Hb, PA or blood-curve parameters are outside the full-rank assertion and would require additional information.
