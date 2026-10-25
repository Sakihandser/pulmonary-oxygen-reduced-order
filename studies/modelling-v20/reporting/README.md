# Supplementary Material S3 Figure and table sources

Use Python 3.12 with the exact dependencies in S2/requirements-lock.txt. From a clean extraction run `python reproduce_figures.py`. This exports Figures 1–7, S1 and S2 as PNG, editable SVG and 600 dpi RGB TIFF without network access. Editable SVGs are included; submission TIFFs are also supplied separately.

main_table1.csv through main_table10.csv match the current article. table_source_map.json and RESULT_REPRODUCTION_INDEX.md identify the current numerical sources and exact calculation/export commands. Table 4 now uses the endpoint-matched rows of common_transfer.csv; cross_state_transfer.csv retains the older single-law fitting-state sensitivity reported in S1 Table S14.

Figure 1 is the computational framework. Figure 2 uses common_transfer.csv and signed_state_map.csv; Figure 3 uses equivalent_conductance.csv; Figure 4 uses common_boundary_comparison.csv; Figure 5 uses allocation_scan.csv; Figure 6 uses k_scan.csv, D_scan.csv, demand_path.csv and hydraulic_scan.csv; Figure 7 uses sequential_decomposition.csv and ablation_matrix.csv. Figure S1 uses external_alvin_benchmark.csv and ratio_scan.csv; Figure S2 uses dissolved_grid.csv. dissolved_summary.csv, dissolved_calibration.csv and ranked_grid_states.csv support the mechanism and state-selection text.

The four convention-based constructions use three transfer laws. Dmatch and conditional Dgas differ in calibration within one constant-D law. Ordered scenarios are path dependent. Numerical separation is conditional on the shared blood relation, frozen fit and declared error bound. Neither grid counts nor the external pressure submodule establish clinical discrimination or physiological accuracy. All text is UTF-8.
