# Supplementary Material S2 Numerical data and code

This archive supplies the complete numerical materials for the article. Use Python 3.12 and run from a clean extraction:

```text
python -m pip install -r requirements-lock.txt
python reproduce.py
```

After dependencies are installed, reproduction uses no network. The command recomputes aggregate observation analyses, the spatial identifiability calculation, all 128 controlled exchange conditions, independent endpoint checks and the sensitivity experiments. It compares generated numerical CSV values with the supplied full-precision outputs using an absolute tolerance of 2e-6 plus a relative tolerance of 1e-8. It verifies numerical main Tables 4–8 at printed precision and exports Figures 1–7, S1 and S2 into sensitivity/reporting/figures. Runtime depends on hardware. The generated reproduction_summary.json reports the actual checks and execution environment.

Directories and roles:

- hydraulic_network: steady resistance-network implementation.
- network_inputs: reference hydraulic states and flow vectors.
- aggregate_analysis: attributed published summary cells, flow-response fits and pressure-partition examples.
- observation_analysis: conditional gas inversion, coefficient alternatives and observation maps.
- identifiability: regional allocation witnesses at unchanged whole-lung gas observations.
- reference_matrix: full-precision reference for the 128-condition comparison.
- comparison_data: common-boundary comparisons, endpoint checks and output definitions.
- controlled_exchange: relaxation and content-state transfer implementations, the 128-condition calculation and independent quadrature audit.
- sensitivity: protocol.json, pressure-state solver, scans, restricted external pressure-kernel qualification and main-table sources.

The ten main-table files are sensitivity/reporting/main_table1.csv through main_table10.csv. table_source_map.json identifies numerical sources or the methodological basis of each table. Tables 1–3, 9 and 10 describe methods, evidence or future empirical design; they are not additional fitted observations. S3 supplies the same main-table CSVs and figure sources. CALIBRATION_CONVENTIONS.md explains the distinct reference conductance, scan setting and display rounding.

The optional native source-kernel audit on Windows is `powershell -File sensitivity/code/run_native_alvin.ps1`. The upstream C# source, pinned commit, GPL license, compilation harness and 250 native nodes are supplied. The portable Python command checks the supplied native node records; it does not compile C# on every platform. The source stores 50 nodes with 49 updates. The qualification covers this pressure submodule with a harmonized blood relation and does not evaluate the Unity application or establish physiological validity.

Published observations are attributed aggregate transcriptions. Publisher full text and participant-level records are unnecessary for these calculations and are not included. Source Hb normalization remains unresolved. Source provenance and third-party licenses are provided in LICENSE_DATA.md and sensitivity/external_alvin/PROVENANCE.md. release_sha256.json records the distributed file hashes; generated outputs change when reproduction is run.


## Common calibration and error accounting

The common_calibration directory supplies two shared targets for all three laws, six withheld boundary states, a 441-state inlet/flow grid per target, and a 1925-state allocation scan. Run `python common_calibration/code/run.py` to regenerate these calculations. `python reproduce.py` also runs this module and checks the updated Table 4, then exports seven main figures and supplementary Figures S1–S2. Scalar and semicolon-separated numeric cells use absolute tolerance 2e-6 plus relative tolerance 1e-8; nonnumeric text must match exactly. The observation-error scenarios are deterministic bounds, not measured error distributions.


See RESULT_REPRODUCTION_INDEX.md for the current numbered table/figure mapping. The ablation independently refits the same reference target for b=0 and b=0.0031; all other inputs remain fixed. Run `python -m unittest discover -s tests` for the supplied regression checks.
