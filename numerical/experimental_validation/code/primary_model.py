"""Primary PH oxygen model: exact equilibrium/admixture and systemic balance.

Units: Q L/min, M mL O2/min, Hb g/dL, PA mmHg, C mL O2/dL.
This is an effective whole-lung closure, not a regional V/Q reconstruction.
"""
import math
from validate import CONFIG, content, invcontent, sat


def solve(q, hb, pa, demand, admixture, p50=None, hill=None):
    p50 = CONFIG['blood']['p50_mmHg'] if p50 is None else p50
    hill = CONFIG['blood']['hill_exponent'] if hill is None else hill
    values = (q, hb, pa, demand, admixture, p50, hill)
    if not all(math.isfinite(v) for v in values):
        raise ValueError('All model inputs must be finite')
    if q <= 0 or hb <= 0 or pa <= 0 or demand < 0 or not 0 <= admixture < 1 or p50 <= 0 or hill <= 0:
        raise ValueError('Require Q, Hb, PA, P50, n > 0; M >= 0; 0 <= s_eff < 1')
    ce = content(pa, hb, p50, hill)
    minimum_flow = demand / (10 * (1-admixture) * ce)
    cv = ce - demand / (10 * q * (1-admixture))
    empty = dict(feasible=False, Pa_mmHg='', Sa_percent='', Ca_ml_dl='', Cv_ml_dl='',
                 Sv_percent='', DO2_ml_min='', mass_balance_error='')
    # Admit the mathematical Cv=0 boundary; rounding tolerance is not an oxygen floor.
    if cv < -1e-12:
        return dict(empty, Q_min_L_min=minimum_flow, equilibrium_content_ml_dl=ce)
    cv = max(0., cv)
    ca = (1-admixture)*ce + admixture*cv
    p = invcontent(ca, hb, p50, hill)
    pv = invcontent(cv, hb, p50, hill)
    return dict(feasible=True, Pa_mmHg=p, Sa_percent=100*sat(p,p50,hill), Ca_ml_dl=ca,
                Cv_ml_dl=cv, Sv_percent=100*sat(pv,p50,hill), DO2_ml_min=10*q*ca,
                mass_balance_error=10*q*(ca-cv)-demand,
                Q_min_L_min=minimum_flow, equilibrium_content_ml_dl=ce)


def delivery_derivatives(q, hb, pa, demand, admixture, p50=None, hill=None):
    """Local derivatives at a feasible interior state; others held constant."""
    p50 = CONFIG['blood']['p50_mmHg'] if p50 is None else p50
    hill = CONFIG['blood']['hill_exponent'] if hill is None else hill
    state = solve(q,hb,pa,demand,admixture,p50,hill)
    if not state['feasible']:
        raise ValueError('Derivative requires a feasible state')
    s = sat(pa,p50,hill)
    a,b = CONFIG['blood']['binding_ml_g'],CONFIG['blood']['dissolved_ml_dl_mmHg']
    return dict(Q=10*state['equilibrium_content_ml_dl'],
                Hb=10*q*a*s,
                PA=10*q*(a*hb*hill*s*(1-s)/pa+b),
                M=-admixture/(1-admixture),
                s_eff=-demand/(1-admixture)**2)
