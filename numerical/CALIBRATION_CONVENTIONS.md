# Calibration constants and reporting conventions

The reference Dmatch is 5.47329138964624 mL/min/mmHg. The pressure-state implementation loads it from reference_matrix/results/ablation_matrix.csv. An independent pressure-state refit differs by approximately 3.2e-8; that numerical refit is reported separately and does not replace the fixed reference constant.

sensitivity/protocol.json specifies 5.473290260830361 as a deterministic D-scan setting. It differs from Dmatch by approximately 1.13e-6 and is not a second reference calibration. The 28.09197 scan setting is rounded, while the gas construction computes its coefficient mapping at full precision. The main tables round their source values; numerical result CSVs use full precision.

Dgas values 28.09197 at Hb=15 and 27.67245 at Hb=13.8 g/dL describe different blood-coefficient conditions. The thetaCO sensitivity slice holds Hb=13.8, gas signals and other coefficients fixed. It is not a confidence interval.
