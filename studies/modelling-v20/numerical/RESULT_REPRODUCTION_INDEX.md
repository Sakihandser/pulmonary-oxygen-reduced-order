# Result reproduction index

All paths below are relative to S2. Install the exact packages in requirements-lock.txt, then run `python reproduce.py`. That command regenerates all result CSVs, verifies printed numerical tables and exports Figures 1–7, S1 and S2. S3 contains the same reporting inputs and exporters.

| Main result | Full-precision CSV | Calculation | Export |
| --- | --- | --- | --- |
| Table 4; Figure 2A–B, shared calibration | common_calibration/results/calibration.csv; transfer.csv | python common_calibration/code/run.py | python sensitivity/reporting/tools/plot_v20.py |
| Figure 2C–D, state maps and thresholds | common_calibration/results/signed_state_map.csv; discrimination_summary.csv; ranked_grid_states.csv | python common_calibration/code/mechanism.py (after run.py) | python sensitivity/reporting/tools/plot_v20.py |
| S1 G; Figure S2, dissolved-oxygen control | common_calibration/results/dissolved_grid.csv; dissolved_summary.csv; dissolved_calibration.csv | python common_calibration/code/mechanism.py | python sensitivity/reporting/tools/plot_v20.py |
| Figure 3, equivalent conductance | controlled_exchange/results/equivalent_conductance.csv | python controlled_exchange/code/reproduce_diagnostic.py | python sensitivity/reporting/tools/plot_main_figures.py |
| Table 5; Figure 4, convention-based comparisons | comparison_data/results/common_boundary_comparison.csv | python controlled_exchange/code/reproduce_diagnostic.py | python sensitivity/reporting/tools/plot_main_figures.py |
| Table 6; Figures 5–6, allocation and intermediate concentration | sensitivity/results/k_scan.csv; hydraulic_scan.csv; D_scan.csv; common_calibration/results/allocation_scan.csv | python sensitivity/code/run_experiments.py; python common_calibration/code/run.py | python sensitivity/reporting/tools/plot_v20.py; python sensitivity/reporting/tools/plot_sensitivity_figures.py |
| Table 7; Figure 7, ordered diagnostic path | controlled_exchange/results/sequential_decomposition.csv | python controlled_exchange/code/reproduce_diagnostic.py | python sensitivity/reporting/tools/plot_main_figures.py |
| Table 8, closure sensitivity | sensitivity/results/feedback_grid.csv; total_demand_closure.csv | python sensitivity/code/run_experiments.py | numerical table only |
| Section 3.6; Figure S1A, external numerical qualification | sensitivity/results/external_alvin_benchmark.csv; external_alvin_trajectories.csv | python sensitivity/code/run_experiments.py | python sensitivity/reporting/tools/plot_sensitivity_figures.py |

Figure state-map contours interpolate the sampled grid for display only. The CSV boolean columns define the exact qualifying grid states. The ranked states maximize a deterministic model difference on this grid; they do not optimize a physiological experiment or supply selection probabilities.
