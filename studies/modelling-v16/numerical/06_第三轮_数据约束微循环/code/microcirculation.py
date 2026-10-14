"""Candidate mechanisms, not an experimentally validated physiological model.

Units: pressure mmHg, Q L/min, volume mL, conductance mL/min/mmHg,
oxygen content mL/dL, resistance Wood units. No systemic demand feedback.
"""
from dataclasses import dataclass
import math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


@dataclass(frozen=True)
class Coupling:
    reference_pressure: float
    reference_volume: float
    reference_dm: float
    reference_capillary_resistance: float
    recruitment_exponent: float
    radius_exponent: float
    contact_exponent: float
    resting_recruited_fraction: float

    def __post_init__(self):
        positives = (self.reference_pressure, self.reference_volume,
                     self.reference_dm, self.reference_capillary_resistance)
        if not all(math.isfinite(v) and v > 0 for v in positives):
            raise ValueError("Reference quantities must be finite and positive")
        if not 0 < self.resting_recruited_fraction <= 1:
            raise ValueError("Resting recruited fraction must be in (0,1]")
        if not all(math.isfinite(v) for v in (self.recruitment_exponent,
                                             self.radius_exponent, self.contact_exponent)):
            raise ValueError("Exponents must be finite")

    def equilibrium(self, pressure):
        """Local power law with finite recruitment reserve; radius/contact unbounded.

        Outside the fitted pressure interval this is only a hypothesis. Saturation
        fraction and all dynamic constants need additional experimental data.
        """
        if not math.isfinite(pressure) or pressure <= 0:
            raise ValueError("Positive transmural pressure required; collapse is unsupported")
        ratio = pressure / self.reference_pressure
        n = min(ratio ** self.recruitment_exponent,
                1 / self.resting_recruited_fraction)
        return np.array([n, ratio ** self.radius_exponent,
                         ratio ** self.contact_exponent])

    def observables(self, state):
        n, radius, contact = np.asarray(state, dtype=float)
        if not all(math.isfinite(v) and v > 0 for v in (n, radius, contact)):
            raise ValueError("All states must be finite and positive")
        return dict(volume=self.reference_volume * n * radius ** 2,
                    dm=self.reference_dm * n * radius * contact,
                    capillary_resistance=self.reference_capillary_resistance /
                    (n * radius ** 4))

    def step(self, state, pressure, duration, time_constants):
        """Exact relaxation for constant prescribed pressure, not fitted dynamics."""
        time_constants = np.asarray(time_constants, dtype=float)
        if duration < 0 or not math.isfinite(duration):
            raise ValueError("Nonnegative finite duration required")
        if time_constants.shape != (3,) or np.any(time_constants <= 0):
            raise ValueError("Three positive time constants required")
        if not np.isfinite(time_constants).all():
            raise ValueError("Finite time constants required")
        self.observables(state)
        target = self.equilibrium(pressure)
        return target + (np.asarray(state) - target) * np.exp(-duration / time_constants)


def solve_bed(coupling, arterial_pressure, venous_pressure,
              arterial_resistance, venous_resistance, external_pressure=0.0):
    """Self-consistent steady, single-bed Zone-3 candidate.

    The midpoint capillary pressure approximation presumes a uniform pressure
    drop within the bed. It is not an identification of actual capillary pressure.
    Pressure and structural parameters must be provided explicitly by the caller.
    """
    if not (arterial_pressure > venous_pressure > external_pressure):
        raise ValueError("Only positive-transmural-pressure Zone-3 beds supported")
    if arterial_resistance < 0 or venous_resistance < 0:
        raise ValueError("Negative segment resistance")

    def evaluate(pc):
        state = coupling.equilibrium(pc - external_pressure)
        out = coupling.observables(state)
        rc = out['capillary_resistance']
        q = (arterial_pressure - venous_pressure) / (arterial_resistance + rc + venous_resistance)
        return pc - venous_pressure - q * (venous_resistance + rc / 2), q, out, state

    pc = brentq(lambda v: evaluate(v)[0], venous_pressure, arterial_pressure,
                xtol=1e-12)
    residual, q, out, state = evaluate(pc)
    return dict(out, flow=q, capillary_pressure=pc, state=state.tolist(),
                pressure_residual=residual,
                total_resistance=(arterial_pressure - venous_pressure) / q)


def oxygen_content(po2, hb=15., p50=26.8, hill_exponent=2.7):
    if hb <= 0 or p50 <= 0 or hill_exponent <= 0 or po2 < 0:
        raise ValueError("Invalid oxygen dissociation parameters")
    fraction = po2 ** hill_exponent / (po2 ** hill_exponent + p50 ** hill_exponent)
    return 1.34 * hb * fraction + 0.0031 * po2


def pressure_from_content(content, hb=15., p50=26.8, hill_exponent=2.7):
    if content < 0:
        raise ValueError("Negative oxygen content")
    upper = max(200., content / 0.0031 + 1.)
    return brentq(lambda p: oxygen_content(p, hb, p50, hill_exponent) - content,
                  0., upper, xtol=1e-10)


def end_capillary_content(flow, oxygen_conductance, alveolar_po2, inlet_content,
                          hb=15., p50=26.8, hill_exponent=2.7, rtol=1e-9):
    """Plug flow, integrating content rather than an imposed transit-time tau.

    dC/dx = D/(10Q) * (PA-P(C)); x ranges 0..1. D must be independently
    supplied. DLCO or DmCO cannot be substituted without a transfer model.
    This implementation uses a Hill curve, not a measured dissociation curve.
    """
    if flow <= 0 or oxygen_conductance < 0 or alveolar_po2 < 0:
        raise ValueError("Invalid exchange boundary")
    ceq = oxygen_content(alveolar_po2, hb, p50, hill_exponent)
    if not 0 <= inlet_content <= ceq:
        raise ValueError("This oxygenation solver requires inlet content <= alveolar equilibrium")
    if oxygen_conductance == 0 or inlet_content == ceq:
        return float(inlet_content)

    def rhs(_x, y):
        content = min(max(float(y[0]), 0.), ceq)
        return [oxygen_conductance / (10 * flow) *
                (alveolar_po2 - pressure_from_content(content, hb, p50, hill_exponent))]

    result = solve_ivp(rhs, (0., 1.), [inlet_content], rtol=rtol, atol=rtol * 0.01)
    if not result.success:
        raise RuntimeError(result.message)
    return float(np.clip(result.y[0, -1], inlet_content, ceq))


def mix_beds(flows, contents, inlet_content, bypass_flow=0.):
    flows, contents = np.asarray(flows), np.asarray(contents)
    if flows.shape != contents.shape or np.any(flows <= 0) or bypass_flow < 0:
        raise ValueError("Invalid mixing flows")
    total = float(flows.sum() + bypass_flow)
    mixed = float((flows @ contents + bypass_flow * inlet_content) / total)
    uptake = float(10 * (flows @ (contents - inlet_content)))
    return dict(content=mixed, uptake=uptake, total_flow=total,
                balance_residual=10 * total * (mixed - inlet_content) - uptake)


OBSERVATION_MATRIX = np.array([[1., 2., 0.], [1., 1., 1.]])
HYDRAULIC_ROW = np.array([1., 4., 0.])
NULL_DIRECTION = np.array([-2., 1., 1.])
