# Response model selection and compensation analysis

Run `python selection_extension/code/reproduce_selection.py` from the numerical archive root after installing the locked dependencies. The root `python reproduce.py` includes this module and the earlier supporting calculations. Reproduction is offline.

The primary task predicts a changed Pa after one permitted baseline calibration, using current Pv, flow and gas inputs. The full evidence inventory is 37 quantitative records from eight studies; 11 nonhyperoxic changes from six baselines and three cohort families enter the primary response score. Newly extracted observations are 15 records from three reports, including seven states in two historical people. They are not 15 new patients or three independent modern cohorts.

| Result | File |
| --- | --- |
| Inputs and baseline flags | results/analysis_inputs.csv |
| Frozen predictions and inclusive exponent predictions | results/response_predictions.csv |
| All source and scope scores | results/response_metrics.csv |
| 108 assumption scenarios | results/response_sensitivity.csv |
| Scenario winners | results/response_sensitivity_winners.csv |
| Complete source and baseline deletion | results/group_deletion.csv |
| Nonhyperoxic held-family exponent | results/nonhyperoxia_exponent_predictions.csv and mechanism_summary.json |
| Inclusive exponent and fit profiles | results/flow_exponent_folds.csv and flow_exponent_profile.csv |
| Alternative loss scales | results/endpoint_scale_scores.csv |
| Exact input attribution and residuals | results/response_decomposition.csv |
| Outcome-only parameter diagnostics | results/state_diagnostics.csv |
| Flow and pressure compensation boundaries | results/compensation_boundaries.csv and pressure_resistance_thresholds.csv |
| Two-output map and examples | results/compensation_map.csv and compensation_examples.csv |
| Main article figures | figures/Fig1–Fig5, PNG/TIFF/SVG |

Do not confuse `flow_adaptive` in response_metrics.csv (trained using the inclusive regime) with `flow_adaptive_nonhyperoxia` in the separate nonhyperoxic prediction CSV. The article's 3.08 mmHg result uses the latter. All-data fitted profiles are not validation scores. Assumption-scenario winner counts are not posterior model probabilities. The selected model is a conditional predictor, and the compensation maps are model scenarios rather than clinical treatment cutoffs.

The five main figures have distinct roles: evidence hierarchy, complete residual comparison, regime and capacity scaling, input attribution, and exact compensation maps. Plot generation checks that Figure 1 text bounding boxes remain inside their shapes.
