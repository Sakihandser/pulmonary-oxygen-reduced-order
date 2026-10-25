# Output and calibration dictionary

| Field or constant | Unit | Definition and source |
| --- | --- | --- |
| k | dimensionless | Distal arterial multiplier contrast; 1 and 25 are specified scenarios, not PAH severity grades |
| allocation fixed | dimensionless | CT air-volume weights; not measured capillary blood-volume shares |
| allocation flow | dimensionless | Capacity allocated by current flow shares; a prescribed benchmark, not a measured recruitment law |
| hb | g/dL | 15 or 13.8; the latter is a source-group mean, not a participant value paired with the network |
| cv_ml_dl | mL/dL | Inlet content from fixed Pv=40 through the Hill relation, or from systemic feedback |
| mixed_end_ml_dl | mL/dL | Cm; terminal exchange content mixed by regional flow |
| ca_ml_dl | mL/dL | Ca=(1-s)Cm+sCv, with s=0 or 0.03; prescribed post-exchange mixing |
| sa_fraction | fraction 0-1 | Inverted from Ca using the common Hill relation; not measured clinical co-oximetry |
| do2_ml_min | mL/min | 10QCa, systemic oxygen delivery; distinct from oxygen exchange conductance D |
| uptake_ml_min | mL/min | 10Q(Ca-Cv), actual uptake consistent with the specified mixing boundary |
| systemic_used_ml_min | mL/min | Actual systemic use under feedback, limited by extraction ceilings; need not reach 250 |
| redistribution_content_loss | mL/dL | Ca(k=1)-Ca(k=25) within the specified condition; Ca=Cm in the common no-mixing comparison |
| redistribution_sa_loss_pp | percentage points | 100[Sa(k=1)-Sa(k=25)]; not a relative percentage reduction |
| d_total | mL/min/mmHg | Total D for constant-conductance models; blank for relaxation models without constant D |
| Vc relaxation | mL | 86, prescribed relaxation capacity; distinct from gas-inverted Vc |
| tau | s | 0.30, prescribed relaxation time; not identified from independent oxygenation data across states |
| Dmatch | mL/min/mmHg | 5.47329; model-to-model endpoint match in one uniform bed at Q=3.94644, Pv=40, PA=100 and Hb=15 |
| Dgas Hb15 | mL/min/mmHg | 28.09197; conditional mapping at DLCO=16.4, DLNO=48, thetaCO=.45, thetaNO=4.5 and alpha=1.97 |
| Dgas Hb13.8 | mL/min/mmHg | 27.67245; exploratory scaling thetaO2=1.8*13.8/15; original DLCO not recorrected |
| thetaO2 | mL O2/mL blood/min/mmHg | Normoxic operational reference 1.8; 0.6 mmol/mL/min/kPa*22.414/7.50062 is approximately 1.793; not a universal hypoxic constant |
| thetaCO | mL CO/mL blood/min/mmHg | .45, a diagnostic assumption; matching to source Hb, alveolar oxygen and the correction chain is unresolved |
| lambda O2/CO | dimensionless | 1.23, membrane-conductance conversion; distinct from Hb coefficients, total oxygen conductance and blood-side reaction terms |

The previously recorded work checked inputs, conversions and implementation. Physiological applicability of reaction coefficients, original Hb normalization, source covariance and measured regional allocation remain unresolved. CSV decimal precision is computational precision, not experimental measurement precision.

In separate exploratory coupling CSV files, venous_fraction_fixture acts only on the LUL: Rv_LUL equals this fraction times reference-bed total resistance. Other lobes retain Rv=0.2 times reference-bed resistance. arterial_scale is a common upstream resistance ratio determined by an outer root solve to restore prescribed total PVR and Q; it is not a fitted PAH remodeling coefficient. oxygen_scale_fixture scales oxygen conductance only, without hydraulic-resistance scaling, and does not represent observed pathological microvascular-bed loss. These coupling files are not external validation for this diagnostic article.
