# Decision consistency under common multi state constraints

Frozen before extension calculations on 14 November 2026, after inspecting the earlier analyses. This is exploratory methodological research, not a clinical preregistration.

The decision concerns a specified oxygen-content threshold at a specified positive flow. All models in a nonempty class must meet the same observation intervals. Distinguish a conclusion proved by an outer bound, opposite conclusions exhibited by feasible finite models, incompatibility proved by a positive relaxation margin, and unresolved numerical cases. Failure of one candidate or one numerical search is not class-wide evidence.

For the computational example the fixed-capacity parallel class is s(r)=E[exp(-X/r)] with fixed weights, common path endpoints and nonnegative exposure X in [0.000001,32]. X=K_i/(w_i Q0). This support is a declared benchmark assumption, not measured PAH anatomy. All observation states use the same distribution. Optional E[X] information represents independently measured total capacity divided by reference flow; no gas-derived capacity is silently treated as this quantity.

Use cellwise monotonic bounds for a continuous-support outer linear program, and finite grid-node mixtures for inner witnesses. Refine 128, 512 and 2048 cells, with additional refinement only if a classification is unresolved. Reconstruct dual lower bounds from multipliers rather than accepting objective values as certificates. Floating calculations have a declared 1e-8 residual decision guard; these are numerically verified enclosures, not machine-checked interval arithmetic.

Synthetic examples reuse the earlier matched-capacity structures. Use arbitrary measurement ratios 1.25 and 1.5 and residual halfwidths 0.001, 0.005 and 0.01, without fitting them to patients. Evaluate r=2 and the prior illustrative arterial recovery target. Separately allow current Ce +/-0.1 mL/dL and M +/-10 percent to expose sensitivity to input assumptions. Retain all cases.

Apply simultaneous-state fitting to each existing primary nonhyperoxic source unit; report minimum common residual perturbation brackets, not patient rejection probabilities. Do not pool individuals, treat study means as individual trajectories, or claim independent validation.

A second experiment constructs finite fixed-capacity mixtures approximating the homogeneous beta=1/2 capacity-response curve over r in [0.5,2]. Use the explicit Levy density and a censored tail. Vary the maximum exposure 8, 16 and 32; report approximation bounds, maximum errors and total capacity. This tests whether an oxygen-response curve uniquely implies recruitment. The exact uncensored representation has infinite mean exposure and is not a finite-capacity lung. Independent capacity restrictions can rule out the constructed alternatives.

The subsequent minimum-capacity exclusion calculation includes an unbounded exposure tail cell [32,infinity), so its lower bound does not assume a physiological upper exposure of 32. A finite grid supplies an attainable upper bound. The two bounds are reported separately; their gap is not hidden. The fixed-support decision demonstrations remain explicitly conditional on the stated support.

Classical set membership, convex optimization, Laplace transforms, log convexity and minimax rules are attributed as existing methods. Novelty is not inferred from their application. No clinical observations are added.
