# External pressure-kernel provenance

Schmid et al. 2022, DOI 10.3389/fbinf.2021.774300.
Source repository: https://github.com/scfischer/schmid-et-al-2022
Pinned commit: 4ef62ff6ecdc8b78c527bcb6fabe3d0defcdbf7a
Files: code/Assets/Simulation/Systems/GasExchange/PartialPressure.cs
Retrieved 16 October 2026; upstream GPL-3.0 LICENSE retained.

run_native_alvin.ps1 extracts the original CharacterizeSections and CalculatePressures methods unchanged, supplies their fields, and compiles them with Add-Type. It executes 50 stored nodes for five common-boundary states. The Python adaptation in run_experiments.py reproduces its pressure recurrence, using double precision and the common Hill relation afterwards. The section count, 49 updates to the last stored node, and original gain cap are explicit. beta is fitted once at reference, then frozen. No Unity application, native Dash saturation, native CO2 prediction or complete external physiological model validation is claimed. The harness and adapted recurrence are distributed with source under GPL-3.0; source measurements retain their original status.
