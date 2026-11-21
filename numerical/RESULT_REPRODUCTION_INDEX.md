# Arterial oxygenation in pulmonary hypertension

Run `python reproduce.py` from a fresh extraction with requirements-lock.txt installed. The runner checks 146 legacy result CSVs and nine oxygen-analysis CSVs, runs 137 legacy unit tests plus the oxygen-decision regression checks, and exports 26 submitted figures. Tables 1–12 are supplied in Figure_Sources.zip. Figure 12 and Supplement V use oxygen_decisions/verify.py and the nine CSVs documented in oxygen_decisions/README.md. Figures 10–11 use certainty_extension. The distributed manifest identifies the exact analysis.

## Archived module descriptions and historical numbering

# Current observed-state audit

Sections 2.8 and 3.8, Figure 9, Table 10 and Supplement S use applicability_extension/results: pair_compatibility.csv, prediction_separation.csv, reported_saturation_audit.csv, input_sensitivity.csv, state_curve_sensitivity.csv, sensitivity_summary.csv and measurement_separation.csv. Root reproduction checks 142 result CSVs and 131 tests. The audit does not add observations or independently validate the structural model.

The following retained sections describe individual earlier modules and their historical run counts.

# Output adequacy

Sections 2.6.3–2.6.4, 3.7, Figure 8, Tables 8–9 and Supplement R use adequacy_extension/code/analyze_adequacy.py and figures_adequacy.py. Seven CSVs and twelve tests extend the full reproduction to 126 CSVs and 131 tests. The older numbering below refers to the earlier modules.

# Current structural analysis

Main Section 2.6.2, Section 3.6, Figure 7, Table 7 and Supplement Q use `heterogeneity_extension/code/analyze_heterogeneity.py` and `figures_heterogeneity.py`. The five full-precision outputs are matched_models.csv, heterogeneous_response.csv, bounded_recovery.csv, sharp_bounds.csv and local_slope_bounds.csv. Ten tests independently check balance roots, derivatives, homogeneous limits, bound attainers, matching constraints, feasibility and recovery decisions. Run reproduce_heterogeneity.py for this extension or root reproduce.py for all 119 result CSVs and 96 tests.

The documentation below retains earlier benchmark numbering. Final submitted figure names are mapped in README_REPRODUCE.md and Figure_Sources.zip.

# Capacity scaling and oxygen compensation

The main article separates four exchange maps from two capacity–flow rules and derives the closed response for flow-dependent capacity. Final Figures 1–2 come from review_extension; Figures 3–5 come from selection_extension; Figure 6 comes from coupling_extension. Main Tables 1–6 are supplied as article_table_1.csv through article_table_6.csv in Figure_Sources.zip. Supporting benchmark figures use the B prefix in the submission; their stable executable filenames are retained below.

Install requirements-lock.txt, then run `python reproduce.py` in a disposable copy. Seven suites contain 86 tests. The complete offline run checks 114 result CSV files and 361,880 numeric cells, including four CSVs and eight tests in coupling_extension. Final submitted figures total 20: six main, twelve benchmark and two supplementary. Two additional selection-extension figure outputs are retained as comparison diagnostics, but are not the submitted Figures 1–2.

`review_extension/PROTOCOL.md` describes the exploratory matched comparison; `review_extension/FLOW_PROVENANCE.md` records flow methods. Carta 2022 and Lichtblau 2021 used thermodilution; Häfliger 2025 and the secondary Lichtblau 2025 study used direct Fick. Current Q cancels in every proportional-capacity control. No new participants or independent observations are created by this extension.

The principal results are the matched scaling contrast and the closed arterial-response boundary beta_crit=(1−s)/(−ln s). At fixed consumption and end-capillary content, increasing flow lowers arterial content above that boundary, while delivery rises in the feasible domain for 0≤beta≤1. The prescribed-fraction identity QD/Q0−1 = epsilon(QA/Q0−1) remains valid when the same fixed s is used at both thresholds. Coupled recovery roots generally use different s values and do not obey that identity. See coupling_extension/PROTOCOL.md and Supplement P for the derivation, counterexample and scope. These model consequences have not been independently validated as patient responses.

The release manifest covers the supplied files. Timing and environment in saved run summaries describe the recorded run; do not interpret their temporary paths as required locations. No network is used after dependency installation.

## Coupled response files

| Article result | Full-precision source | Calculation |
| --- | --- | --- |
| Equations (9)–(10), Figure 6 | coupling_extension/results/slope_boundary.csv; flow_response_curves.csv | coupling_extension/code/analyze_coupling.py |
| Table 6 and Supplement Table S43 | coupling_extension/results/recovery_thresholds.csv | Same script; independent recovery-root checks in tests/test_coupling.py |
| Supplement Equation (P4) | coupling_extension/results/conditional_scaling_sensitivity.csv | Same script, using the frozen conditional-response inputs |

## Supporting benchmark documentation

The following retained documentation uses the original benchmark table and figure numbers; these are not the current main-article numbers.

# Result reproduction index

All paths below are relative to S2. Install the exact packages in requirements-lock.txt, then run `python reproduce.py`. That command regenerates all result CSVs, verifies printed numerical tables and exports Figures 1–12, S1 and S2. S3 contains the same reporting inputs and exporters.

| Main result | Full-precision CSV | Calculation | Export |
| --- | --- | --- | --- |
| Table 4; Figure 2A–B, shared calibration | common_calibration/results/calibration.csv; transfer.csv | python common_calibration/code/run.py | python sensitivity/reporting/tools/plot_calibration_maps.py |
| Figure 2C–D, state maps and thresholds | common_calibration/results/signed_state_map.csv; discrimination_summary.csv; ranked_grid_states.csv | python common_calibration/code/mechanism.py (after run.py) | python sensitivity/reporting/tools/plot_calibration_maps.py |
| S1 G; Figure S2, dissolved-oxygen control | common_calibration/results/dissolved_grid.csv; dissolved_summary.csv; dissolved_calibration.csv | python common_calibration/code/mechanism.py | python sensitivity/reporting/tools/plot_calibration_maps.py |
| Figure 3, equivalent conductance | controlled_exchange/results/equivalent_conductance.csv | python controlled_exchange/code/reproduce_diagnostic.py | python sensitivity/reporting/tools/plot_main_figures.py |
| Table 5; Figure 4, convention-based comparisons | comparison_data/results/common_boundary_comparison.csv | python controlled_exchange/code/reproduce_diagnostic.py | python sensitivity/reporting/tools/plot_main_figures.py |
| Table 6; Figures 5–6, allocation and intermediate concentration | sensitivity/results/k_scan.csv; hydraulic_scan.csv; D_scan.csv; common_calibration/results/allocation_scan.csv | python sensitivity/code/run_experiments.py; python common_calibration/code/run.py | python sensitivity/reporting/tools/plot_calibration_maps.py; python sensitivity/reporting/tools/plot_sensitivity_figures.py |
| Table 7; Figure 7, ordered diagnostic path | controlled_exchange/results/sequential_decomposition.csv | python controlled_exchange/code/reproduce_diagnostic.py | python sensitivity/reporting/tools/plot_main_figures.py |
| Table 8, closure sensitivity | sensitivity/results/feedback_grid.csv; total_demand_closure.csv | python sensitivity/code/run_experiments.py | numerical table only |
| Section 3.6; Figure S1A, external numerical qualification | sensitivity/results/external_alvin_benchmark.csv; external_alvin_trajectories.csv | python sensitivity/code/run_experiments.py | python sensitivity/reporting/tools/plot_sensitivity_figures.py |

Figure state-map contours interpolate the sampled grid for display only. The CSV boolean columns define the exact qualifying grid states. The ranked states maximize a deterministic model difference on this grid; they do not optimize a physiological experiment or supply selection probabilities.


The bounded-calibration and grid-resolution analyses in S1 Section I are generated by `python robustness/code/run.py` in S2. The five CSVs in `robustness/results/` contain the seven target-offset fits, exact monotone marginal envelopes, their separation counts, all 882 baseline grid verification states, and 21/41/81 grid-resolution comparisons. The one-command S2 run executes these calculations and compares all five CSVs against the distributed values. The numerical allowance remains 1e-7 mL/dL per output; target halfwidth c and output bound e are separate assumed inputs. These are deterministic scenarios, not confidence intervals.

The public commit 7460853b82a8e72bf4ded4b1c94d14a3baaaa676 records the original computational core. The additional robustness analysis is distributed in the present S2 archive, whose release_sha256.json identifies the complete supplied files. This archive is the reproduction source for Section I; the earlier public commit alone does not contain that extension.


Sections 2.8–2.9, 3.7–3.9, 4.6 and S1 J, Tables S26–S29 and Figures 8–11 are mapped in experimental_validation/README.md. The 12 new numerical CSVs are regenerated and compared by experimental_validation/code/reproduce_empirical.py, which is also invoked by the root entry point.


Primary PH refinement: `experimental_validation/code/analyze_primary.py` generates `primary_flow_controls.csv`, `primary_counterfactuals.csv`, `primary_local_sensitivity.csv` and `primary_subgroup_metrics.csv` in its results directory. They support Figure 11, Results 3.10 and Supplement K (Tables S30–S32). `primary_model.py` is the single exact oxygen solver used by this analysis and the network. `model_config.json` and `MODEL_SPECIFICATION.md` supply configuration, units and parameter provenance. Full reproduction now compares 58 original plus 23 empirical CSVs and exports 12 main plus two supplementary figures.


## Ranking stability and observation requirements

Sections 2.10, 3.11, 4.7 and Supplement L (Tables S33–S36) use `experimental_validation/code/review_extension.py`. The seven `review_*.csv` files report exact empirical resampling, case deletion, two gas-input conventions, nominal PA values, equivalent-output witnesses, observation ranks and bounded measurement errors. `review_figures.py` exports Figure 12. `REVIEW_PROTOCOL.md` states the post-selection scope and future two-task validation design. Root reproduction compares 81 result CSVs and runs 55 tests across four suites. The numerical/data evidence has not gained additional participants.

The applicability_extension adds seven diagnostic CSVs, 10 tests and Figure 9. It is an exploratory audit of observed PH states, not new patients or independent clinical validation. See its PROTOCOL.md and README.md.

The certainty_extension adds nine CSVs, 13 tests and Figures 10–11. It evaluates common multi-state model constraints, finite countermodels and capacity exclusion. Older section counts refer to their individual modules. No new patients or independent clinical validation are claimed.

Measurement extension: run `python measurement_extension/code/reproduce_measurement.py`. Four additional CSV outputs and six tests are included in the full runner. The preceding 142-result-CSV and 131-test totals become 146 and 137. Measurement error boxes are stipulated and do not establish clinical coverage.
