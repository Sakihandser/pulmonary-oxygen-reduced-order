"""Explicit gas observation and conditional oxygen-transfer maps.

This is an observation-model analysis, not a physiologically validated model.
Units: D in mL gas/min/mmHg, Vc in mL blood, theta in mL gas/mL blood/min/mmHg.
Alpha relates NO to CO membrane conductance; lambda relates O2 to CO membrane.
"""
import math
import numpy as np


def positive(*values):
    if not all(math.isfinite(x) and x > 0 for x in values):
        raise ValueError('Finite positive parameters required')


def gas_conductance(dm_co, volume, theta, membrane_ratio=1.):
    positive(dm_co,volume,membrane_ratio)
    if theta == math.inf:
        return membrane_ratio*dm_co
    positive(theta)
    return 1/(1/(membrane_ratio*dm_co)+1/(theta*volume))


def invert_two_gases(dlco,dlno,theta_co,theta_no=4.5,alpha=1.97):
    """Exact inversion conditional on known theta/alpha; no error model implied."""
    positive(dlco,dlno,theta_co,alpha)
    if theta_no != math.inf:
        positive(theta_no)
    matrix=np.array([[1.,1/theta_co],
                     [1/alpha,0. if theta_no==math.inf else 1/theta_no]])
    if abs(np.linalg.det(matrix)) < 1e-12:
        raise ValueError('Singular gas observation matrix')
    u,v=np.linalg.solve(matrix,[1/dlco,1/dlno])
    if u <= 0 or v <= 0:
        raise ValueError('Nonpositive inverse membrane or volume estimate')
    return dict(dm=1/u,vc=1/v)


def oxygen_estimate(dm_co,volume,theta_o2=1.8,lambda_o2_co=1.23):
    """RF conditional estimate. theta=1.8 is a normoxic operational reference,
    NOT a universal hypoxaemic or disease parameter. This does not predict SaO2.
    """
    positive(dm_co,volume,theta_o2,lambda_o2_co)
    membrane=lambda_o2_co*dm_co
    blood=theta_o2*volume
    total=membrane*blood/(membrane+blood)
    return dict(effective=total,membrane=membrane,blood=blood,
                membrane_resistance_fraction=blood/(membrane+blood),
                blood_resistance_fraction=membrane/(membrane+blood))


def log_gas_jacobian(dm_co,volume,theta_co,theta_no=4.5,alpha=1.97):
    """d(log DLCO,log DLNO)/d(log Dm,log Vc)."""
    co=gas_conductance(dm_co,volume,theta_co)
    no=gas_conductance(dm_co,volume,theta_no,alpha)
    co_m=co/dm_co; no_m=no/(alpha*dm_co)
    return np.array([[co_m,1-co_m],[no_m,1-no_m]])


def gain_threshold(membrane_ratio,blood_ratio):
    """Critical baseline z=B/M for gain when membrane grows and blood shrinks.
    No coefficients are calibrated to an oxygen outcome.
    """
    positive(membrane_ratio,blood_ratio)
    if not membrane_ratio > 1 or not blood_ratio < 1:
        raise ValueError('Threshold case requires membrane gain and blood loss')
    return membrane_ratio*(1-blood_ratio)/(blood_ratio*(membrane_ratio-1))


VOLUME_ROW=np.array([1.,2.,0.])
MEMBRANE_ROW=np.array([1.,1.,1.])
HYDRAULIC_ROW=np.array([1.,4.,0.])
AREA_ROW=np.array([1.,1.,0.])
GEOMETRY_NULL=np.array([-2.,1.,1.])


def geometry_observation_rows(membrane_fraction):
    if not 0 <= membrane_fraction <= 1:
        raise ValueError('Fraction outside [0,1]')
    base=np.vstack([VOLUME_ROW,MEMBRANE_ROW])
    oxygen_row=membrane_fraction*MEMBRANE_ROW+(1-membrane_fraction)*VOLUME_ROW
    return dict(base=base,oxygen=oxygen_row,
                plus_oxygen=np.vstack([base,oxygen_row]),
                plus_transit=np.vstack([base,VOLUME_ROW]),
                plus_hydraulics=np.vstack([base,HYDRAULIC_ROW]),
                plus_area=np.vstack([base,AREA_ROW]))


def multifio2_reparameterize(raw_intercept,raw_slope,a,b):
    """For fixed physiological state, 1/DLCO = intercept+slope*PAO2.

    theta_CO^-1 = a*PAO2+b. The observed line only identifies slope and
    intercept, not a,b,Dm,V independently.
    """
    positive(raw_intercept,raw_slope,a,b)
    vc=a/raw_slope
    inverse_dm=raw_intercept-b*raw_slope/a
    if inverse_dm<=0:
        raise ValueError('Nonphysical membrane estimate')
    return dict(dm=1/inverse_dm,vc=vc)


def pressure_driver_jacobian(q,ra,rc,rv):
    """Derivatives of PAP and midpoint Pc w.r.t. Pw,Q,Ra,Rc,Rv.

    Passive series-bed identity, evaluated at fixed other variables. It does
    not assume that a real intervention changes only one of these quantities.
    """
    positive(q,ra,rc,rv)
    return np.array([[1.,ra+rc+rv,q,q,q],
                     [1.,rv+rc/2,0.,q/2,q]])
