# Measurement uncertainty and decision margins

Run `python measurement_extension/code/reproduce_measurement.py` from the numerical root.

This extension maps stipulated oxygen-content error boxes to exact residual intervals. It separates common affine calibration errors, which cancel in the residual, from differential errors, which do not. These are analytic measurement models, not empirical estimates of analyzer accuracy. They cannot be substituted for joint patient uncertainty or used to assign clinical probabilities.

For strictly ordered contents throughout the box, monotonicity gives sharp corner bounds on s=(Ce-Ca)/(Ce-Cv). Equal independent bounded content errors of halfwidth e produce residual halfwidth 2e/(Ce-Cv). Here independent means independently variable in a Cartesian set, not statistically independent. Common positive scale and common offset cancel exactly. Hb error is not exactly a common content scale when dissolved oxygen is retained.

The existing finite structures A and B generate three synthetic states at r=1,1.25,1.5. The fixed protocol evaluates e=0.001,0.003,0.005,0.01,0.02 mL/dL using 32768 exposure cells and the prior support [0.000001,32], query r=2, fixed Ce, M and target. No endpoints or target are retuned. All finite witnesses are retained. Query-input uncertainty is evaluated separately with an explicit sufficient margin frontier; it must not be added to observation intervals without considering dependencies.

The code uses the original floating-point dual guard, not directed-rounding certification. The analytic box bounds are exact real-arithmetic results. The numerical bounds and witnesses are checked within declared tolerances. No new participants, measurements or external validation are introduced.
