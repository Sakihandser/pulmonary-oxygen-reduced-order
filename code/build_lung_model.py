#!/usr/bin/env python3
"""Steady reduced pulmonary resistance-network model.

This is a lumped, rigid-tube, steady Poiseuille network. It is not a
classical one-dimensional pulse-wave model because it has no time-dependent
mass/momentum PDEs or compliant-wall equations.

The primary PAH-like perturbation is distal arterial resistance. Venous
resistance is unchanged unless an explicit sensitivity analysis varies it.
The healthy PVR target is a calibration target, not independent validation.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
MU = 3.5e-3
WU = 133.322 * 60.0 / 1e-3
RB = 3.0

# Huang et al. 1996 Table 2; mean diameter and length in mm.
HUANG_ART = {
    15: (14.80, 25.30), 14: (7.34, 35.69), 13: (4.16, 25.97),
    12: (2.71, 18.07), 11: (1.75, 12.35), 10: (1.16, 6.58),
    9: (0.77, 3.73), 8: (0.51, 2.81), 7: (0.34, 1.92),
    6: (0.22, 1.08), 5: (0.15, 0.68), 4: (0.097, 0.45),
    3: (0.056, 0.36), 2: (0.036, 0.26), 1: (0.020, 0.22),
}
HUANG_VEIN = {
    15: (12.97, 35.68), 14: (8.65, 34.99), 13: (5.86, 19.49),
    12: (4.00, 26.49), 11: (2.88, 17.90), 10: (1.99, 14.78),
    9: (1.42, 11.24), 8: (0.90, 6.78), 7: (0.62, 4.79),
    6: (0.38, 2.92), 5: (0.23, 1.50), 4: (0.13, 1.06),
    3: (0.067, 0.38), 2: (0.031, 0.21), 1: (0.018, 0.13),
}

# Qureshi et al. 2014 Table 1; proximal diameter, distal diameter and length in cm.
QURESHI = {
    "MPA": (2.7, 2.6, 4.50), "RPA": (1.8, 1.2, 5.75),
    "LPA": (2.2, 2.2, 2.50), "RIA": (1.1, 1.1, 1.25),
    "RTA": (0.9, 0.9, 1.00), "LIA": (2.1, 1.8, 2.25),
    "LTA": (1.2, 1.2, 1.00),
}

# Male inspiratory lobar CT medians; used only as relative weights.
LOBE_VOL = {
    "RUL": 1038.02, "RML": 486.49, "RLL": 1204.04,
    "LUL": 1278.15, "LLL": 1053.67,
}
LOBE_W = {key: value / sum(LOBE_VOL.values()) for key, value in LOBE_VOL.items()}
LOBE_ORDER = list(LOBE_W)

# (cardiac-output share, baseline oxygen-demand share, maximum extraction)
DEFAULT_ORGANS = {
    "脑": (0.15, 0.20, 0.45),
    "右心室心肌": (0.05, 0.11, 0.70),
    "肾": (0.20, 0.07, 0.20),
    "肠": (0.18, 0.16, 0.35),
    "骨骼肌": (0.22, 0.20, 0.40),
    "其他": (0.20, 0.26, 0.45),
}


def poiseuille(length_m, radius_m, mu=MU):
    return 8.0 * mu * length_m / (math.pi * radius_m**4)


def qureshi_length_cm(radius_cm, kind="artery"):
    if kind == "artery":
        return 15.75 * radius_cm**1.10 if radius_cm >= 0.005 else 1.79 * radius_cm**0.47
    return 14.54 * radius_cm


def structured_resistance(
    kind, order_from, order_to, n_root, open_fraction=1.0,
    resistance_multiplier=1.0, rb=RB, mu=MU,
):
    """Equivalent resistance of a symmetric Strahler tree.

    open_fraction changes the number of perfused paths at orders 8 and below.
    resistance_multiplier changes the resistance of individual distal vessels.
    """
    if not 0.0 < open_fraction <= 1.0:
        raise ValueError("open_fraction must be in (0, 1]")
    if resistance_multiplier <= 0.0:
        raise ValueError("resistance_multiplier must be positive")
    table = HUANG_ART if kind == "artery" else HUANG_VEIN
    total = 0.0
    rows = []
    for order in range(order_from, order_to - 1, -1):
        diameter_mm, length_mm = table[order]
        count = n_root * rb ** (order_from - order)
        multiplier = 1.0
        if order <= 8:
            count *= open_fraction
            multiplier = resistance_multiplier
        single = poiseuille(length_mm / 1000.0, diameter_mm / 2000.0, mu) * multiplier
        parallel = single / count
        total += parallel
        rows.append({
            "kind": kind, "order": order, "D_mm": diameter_mm,
            "L_mm": length_mm, "N": count, "R_one_WU": single / WU,
            "R_parallel_WU": parallel / WU,
        })
    return total, rows


def explicit_large_resistance(mu=MU):
    vessels = {}
    for name, (proximal_diameter, distal_diameter, length_cm) in QURESHI.items():
        # Mean diameter divided by two gives the radius.
        radius_m = ((proximal_diameter + distal_diameter) / 4.0) / 100.0
        vessels[name] = poiseuille(length_cm / 100.0, radius_m, mu)
    for name, order in {"RUL": 14, "RML": 13, "RLL": 14, "LUL": 14, "LLL": 14}.items():
        diameter_mm, length_mm = HUANG_ART[order]
        vessels[name] = poiseuille(length_mm / 1000.0, diameter_mm / 2000.0, mu)
    series_guess = vessels["MPA"] + 0.5 * (vessels["RPA"] + vessels["LPA"])
    return {"vessels": vessels, "series_guess": series_guess}


def build_cases(rb=RB, mu=MU, n_root_override=None):
    large = explicit_large_resistance(mu)
    cap_all = 0.15 * WU * (mu / MU)
    if n_root_override is None:
        artery_unit = structured_resistance("artery", 12, 1, 1.0, rb=rb, mu=mu)[0]
        vein_unit = structured_resistance("vein", 12, 1, 1.0, rb=rb, mu=mu)[0]
        target_micro = 0.85 * WU * (mu / MU) - large["series_guess"]
        n_root = (artery_unit + vein_unit) / target_micro
    else:
        n_root = n_root_override
    artery_rows = structured_resistance("artery", 12, 1, n_root, rb=rb, mu=mu)[1]
    return {
        "n_root": n_root, "art_rows": artery_rows, "large": large,
        "cap_all": cap_all, "rb": rb, "mu": mu,
    }


def solve(
    case, arterial_multiplier=1.0, venous_multiplier=1.0,
    bed_open_fraction=1.0, focal=None, co_max=5.0, reserve=50.0,
    scale_cap=False, lobe_arterial_multiplier=None,
):
    """Solve the steady resistance network.

    arterial_multiplier is the primary PAH-like operator. Bed opening is
    reserved for a separate spatial loss-of-perfused-bed experiment.
    """
    large = case["large"]["vessels"]
    lobe_resistance, lobe_open = {}, {}
    lobe_arterial_multiplier = lobe_arterial_multiplier or {}
    for lobe, weight in LOBE_W.items():
        open_fraction = focal.get(lobe, bed_open_fraction) if focal else bed_open_fraction
        lobe_open[lobe] = open_fraction
        local_arterial = arterial_multiplier * lobe_arterial_multiplier.get(lobe, 1.0)
        artery = structured_resistance(
            "artery", 12, 1, weight * case["n_root"], open_fraction,
            local_arterial, case["rb"], case["mu"],
        )[0]
        vein = structured_resistance(
            "vein", 12, 1, weight * case["n_root"], open_fraction,
            venous_multiplier, case["rb"], case["mu"],
        )[0]
        capillary = case["cap_all"] / weight
        if scale_cap:
            capillary /= open_fraction
        lobe_resistance[lobe] = large[lobe] + artery + vein + capillary

    r_inter = 1.0 / (1.0 / lobe_resistance["RML"] + 1.0 / lobe_resistance["RLL"])
    r_ria = large["RIA"] + r_inter
    r_rul = large["RTA"] + lobe_resistance["RUL"]
    r_right = large["RPA"] + 1.0 / (1.0 / r_ria + 1.0 / r_rul)
    r_left = large["LPA"] + 1.0 / (
        1.0 / lobe_resistance["LUL"] +
        1.0 / (large["LTA"] + lobe_resistance["LLL"])
    )
    pvr = large["MPA"] + 1.0 / (1.0 / r_right + 1.0 / r_left)
    pvr_wu = pvr / WU

    mpap_ref, pawp = 14.0, 8.0
    co = co_max * (1.0 - (pawp - mpap_ref) / reserve)
    co /= 1.0 + co_max * pvr_wu / reserve
    co = max(1.2, co)
    mpap = pawp + co * pvr_wu

    pressure_drop_pa = (mpap - pawp) * 133.322
    q_total = pressure_drop_pa / pvr
    p_after_mpa = pressure_drop_pa - q_total * large["MPA"]
    q_right, q_left = p_after_mpa / r_right, p_after_mpa / r_left
    p_rpa = p_after_mpa - q_right * large["RPA"]
    p_lpa = p_after_mpa - q_left * large["LPA"]
    q_ria = p_rpa / r_ria
    p_ria = p_rpa - q_ria * large["RIA"]
    flows_m3_s = {
        "RUL": p_rpa / r_rul,
        "RML": p_ria / lobe_resistance["RML"],
        "RLL": p_ria / lobe_resistance["RLL"],
        "LUL": p_lpa / lobe_resistance["LUL"],
        "LLL": p_lpa / (large["LTA"] + lobe_resistance["LLL"]),
    }
    flows = {key: value * 60.0 * 1000.0 for key, value in flows_m3_s.items()}
    flow_sum = sum(flows.values())
    return {
        "CO_L_min": co, "mPAP_mmHg": mpap, "PAWP_mmHg": pawp,
        "PVR_WU": pvr_wu, "flows_L_min": flows,
        "flow_fraction": {key: value / flow_sum for key, value in flows.items()},
        "lobe_R_WU": {key: value / WU for key, value in lobe_resistance.items()},
        "open_fraction": lobe_open,
        "arterial_multiplier": arterial_multiplier,
        "venous_multiplier": venous_multiplier,
        "lobe_arterial_multiplier": {
            lobe: arterial_multiplier * lobe_arterial_multiplier.get(lobe, 1.0)
            for lobe in LOBE_ORDER
        },
        "scale_cap": scale_cap,
        "flow_conservation_error_L_min": flow_sum - co,
        "pvr_identity_error_WU": (mpap - pawp) / co - pvr_wu,
    }


def hill(po2, p50=26.8, exponent=2.7):
    po2 = np.asarray(po2, dtype=float)
    return po2**exponent / (po2**exponent + p50**exponent)


def oxygen_content(po2, hb=15.0):
    return 0.0031 * np.asarray(po2, dtype=float) + 1.34 * hb * hill(po2)


PO2_GRID = np.linspace(0.01, 650.0, 65000)
SATURATION_GRID = hill(PO2_GRID)


@lru_cache(maxsize=16)
def _content_grid(hb):
    return oxygen_content(PO2_GRID, hb)


def po2_from_content(value, hb=15.0):
    content_grid = _content_grid(float(hb))
    value = float(np.clip(value, content_grid[0], content_grid[-1]))
    return float(np.interp(value, content_grid, PO2_GRID))


def po2_from_saturation(value):
    value = float(np.clip(value, SATURATION_GRID[0], SATURATION_GRID[-1]))
    return float(np.interp(value, SATURATION_GRID, PO2_GRID))


def oxygen(
    result, vc_ml=86.0, tau_eq=0.30, vo2=250.0, shunt=0.03,
    hb=15.0, pao2=100.0, lobe_pao2=None, organs=None,
    rv_power_scale=False, ref_power=None, vc_mode="flow",
    tolerance=1e-9, max_iterations=200, exchange_mode="saturation",
    dl_o2=25.0, lobe_va_L_min=None, pio2=150.0,
):
    """Solve a closed oxygen-content balance for lung and organ compartments.

    exchange_mode is saturation (baseline), content, or diffusion.
    lobe_va_L_min turns on a lobe-wise alveolar mass balance.
    vc_mode fixed keeps anatomical volume shares and ignores bed opening.
    """
    if vc_mode not in {"flow", "anatomy", "fixed"}:
        raise ValueError("vc_mode must be flow, anatomy, or fixed")
    if exchange_mode not in {"saturation", "content", "diffusion"}:
        raise ValueError("exchange_mode must be saturation, content, or diffusion")
    if not 0.0 <= shunt < 1.0:
        raise ValueError("shunt must be in [0, 1)")
    organs = organs or DEFAULT_ORGANS
    if not math.isclose(sum(values[0] for values in organs.values()), 1.0, abs_tol=1e-9):
        raise ValueError("organ flow fractions must sum to 1")
    flows, co = result["flows_L_min"], result["CO_L_min"]

    demands = {}
    for name, (_, demand_share, _) in organs.items():
        demand = demand_share * vo2
        if name == "右心室心肌" and rv_power_scale:
            if not ref_power:
                raise ValueError("ref_power is required when rv_power_scale is true")
            demand *= max(result["mPAP_mmHg"] * co / ref_power, 1.0)
        demands[name] = demand

    def volume_share(lobe):
        if vc_mode == "flow":
            return result["flow_fraction"][lobe]
        if vc_mode == "fixed":
            return LOBE_W[lobe]
        return LOBE_W[lobe] * result["open_fraction"][lobe]

    def equilibrate(mixed_cv, pv, local_pao2, transit, local_vc):
        if exchange_mode == "saturation":
            s_alv = float(hill(local_pao2))
            s_in = float(hill(pv))
            s_end = s_alv - (s_alv - s_in) * math.exp(-transit / tau_eq)
            po2_end = po2_from_saturation(s_end)
            c_end = float(oxygen_content(po2_end, hb))
            return s_end, po2_end, c_end
        if exchange_mode == "content":
            c_alv = float(oxygen_content(local_pao2, hb))
            c_end = c_alv - (c_alv - mixed_cv) * math.exp(-transit / tau_eq)
            c_end = min(max(c_end, 0.0), c_alv)
            po2_end = po2_from_content(c_end, hb)
            return float(hill(po2_end)), po2_end, c_end
        c = mixed_cv
        pc = pv
        steps = 16
        dt = transit / steps
        c_alv = float(oxygen_content(local_pao2, hb))
        local_dl = dl_o2 * local_vc / vc_ml
        for _ in range(steps):
            dc = local_dl * (local_pao2 - pc) / max(local_vc, 1e-9) * (100.0 / 60.0) * dt
            c = min(max(c + dc, 0.0), c_alv)
            pc = po2_from_content(c, hb)
        return float(hill(pc)), pc, c

    def evaluate(ca):
        organ_output, mixed_cv, total_demand, total_used = {}, 0.0, 0.0, 0.0
        for name, (flow_share, _, emax) in organs.items():
            flow = flow_share * co
            supply = flow * ca * 10.0
            demand = demands[name]
            used = min(demand, emax * supply)
            cv_organ = max(ca - used / (flow * 10.0), 0.0)
            mixed_cv += flow_share * cv_organ
            total_demand += demand
            total_used += used
            organ_output[name] = {
                "Q_L_min": flow, "DO2_mL_min": supply,
                "VO2_demand_mL_min": demand, "VO2_mL_min": used,
                "deficit_mL_min": demand - used,
                "extraction": used / supply if supply else 1.0,
                "emax": emax, "CvO2_mL_dL": cv_organ,
            }

        pv = po2_from_content(mixed_cv, hb)
        local_pao2 = {
            lobe: (lobe_pao2.get(lobe, pao2) if lobe_pao2 else pao2)
            for lobe in flows
        }
        if lobe_va_L_min:
            for _ in range(18):
                trial = {}
                for lobe, flow in flows.items():
                    local_vc = vc_ml * volume_share(lobe)
                    transit = (local_vc / 1000.0) / (flow / 60.0)
                    _, _, c_end = equilibrate(
                        mixed_cv, pv, local_pao2[lobe], transit, local_vc,
                    )
                    uptake = flow * 10.0 * (c_end - mixed_cv)
                    updated = pio2 - 0.863 * uptake / lobe_va_L_min[lobe]
                    trial[lobe] = float(np.clip(updated, 15.0, pio2))
                if max(abs(trial[lobe] - local_pao2[lobe]) for lobe in trial) < 0.05:
                    local_pao2 = trial
                    break
                local_pao2 = {
                    lobe: 0.55 * local_pao2[lobe] + 0.45 * trial[lobe]
                    for lobe in trial
                }

        ends, endcap_content = {}, 0.0
        for lobe, flow in flows.items():
            local_vc = vc_ml * volume_share(lobe)
            transit = (local_vc / 1000.0) / (flow / 60.0)
            s_end, po2_end, c_end = equilibrate(
                mixed_cv, pv, local_pao2[lobe], transit, local_vc,
            )
            endcap_content += result["flow_fraction"][lobe] * c_end
            ends[lobe] = {
                "tau_s": transit, "tau_ratio": transit / tau_eq,
                "S_end": s_end, "PO2_end": po2_end,
                "C_end_mL_dL": c_end, "PAO2": local_pao2[lobe],
                "vc_share": volume_share(lobe), "vc_ml": local_vc,
                "VA_L_min": None if not lobe_va_L_min else lobe_va_L_min[lobe],
            }
        ca_new = (1.0 - shunt) * endcap_content + shunt * mixed_cv
        return ca_new, mixed_cv, organ_output, ends, total_demand, total_used, endcap_content

    ca = float(oxygen_content(min(pao2, 100.0), hb))
    converged = False
    iterations = 0
    for iterations in range(1, max_iterations + 1):
        ca_new = evaluate(ca)[0]
        if abs(ca_new - ca) < tolerance:
            ca, converged = ca_new, True
            break
        ca = 0.5 * (ca + ca_new)

    ca, cv, organ_output, ends, total_demand, total_used, endcap_content = evaluate(ca)
    pa, pv = po2_from_content(ca, hb), po2_from_content(cv, hb)
    s_values = [ends[lobe]["S_end"] for lobe in LOBE_ORDER]
    systemic_extraction = co * 10.0 * (ca - cv)
    return {
        "PaO2": pa, "SaO2": float(hill(pa)),
        "PvO2": pv, "SvO2": float(hill(pv)),
        "CaO2": ca, "CvO2": cv, "endcap_CaO2": endcap_content,
        "lobes": ends, "organs": organ_output,
        "S_end_var": float(np.var(s_values)),
        "S_end_range": float(max(s_values) - min(s_values)),
        "total_VO2_demand_mL_min": total_demand,
        "total_VO2_used_mL_min": total_used,
        "total_deficit_mL_min": total_demand - total_used,
        "oxygen_mass_balance_error_mL_min": systemic_extraction - total_used,
        "converged": converged, "iterations": iterations,
        "effective_vc_ml": sum(item["vc_ml"] for item in ends.values()),
        "mixing_basis": "oxygen_content",
    }


def write_geometry(case):
    rows = []
    for name, (proximal, distal, length_cm) in QURESHI.items():
        rows.append({
            "segment": name, "source": "Qureshi 2014 Table 1",
            "D_prox_mm": proximal * 10.0, "D_dist_mm": distal * 10.0,
            "L_mm": length_cm * 10.0,
            "R_WU": case["large"]["vessels"][name] / WU,
        })
    for name, order in {"RUL": 14, "RML": 13, "RLL": 14, "LUL": 14, "LLL": 14}.items():
        diameter, length = HUANG_ART[order]
        rows.append({
            "segment": name + "_artery_order" + str(order),
            "source": "Huang 1996 Table 2", "D_prox_mm": diameter,
            "D_dist_mm": diameter, "L_mm": length,
            "R_WU": case["large"]["vessels"][name] / WU,
        })
    pd.DataFrame(rows).to_csv(OUT / "explicit_vessels.csv", index=False)
    pd.DataFrame(case["art_rows"]).to_csv(OUT / "huang_arterial_orders.csv", index=False)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-plots", action="store_true", help="retained for command compatibility")
    parser.parse_args(argv)
    case = build_cases()
    scenarios = {
        "健康": solve(case),
        "远端动脉阻力x5": solve(case, arterial_multiplier=5.0),
        "全肺灌注床开放50%": solve(case, bed_open_fraction=0.5, scale_cap=True),
        "仅右下叶灌注床开放40%": solve(case, focal={"RLL": 0.4}, scale_cap=True),
    }
    oxygen_results = {name: oxygen(result) for name, result in scenarios.items()}
    write_geometry(case)
    rows = []
    for name, result in scenarios.items():
        ox = oxygen_results[name]
        rows.append({
            "scenario": name, "CO_L_min": result["CO_L_min"],
            "mPAP_mmHg": result["mPAP_mmHg"], "PVR_WU": result["PVR_WU"],
            "arterial_multiplier": result["arterial_multiplier"],
            "venous_multiplier": result["venous_multiplier"],
            "SaO2": ox["SaO2"], "SvO2": ox["SvO2"],
            "PaO2": ox["PaO2"], "PvO2": ox["PvO2"],
            "total_VO2_demand_mL_min": ox["total_VO2_demand_mL_min"],
            "total_VO2_used_mL_min": ox["total_VO2_used_mL_min"],
            "total_deficit_mL_min": ox["total_deficit_mL_min"],
            "oxygen_mass_balance_error_mL_min": ox["oxygen_mass_balance_error_mL_min"],
            **{"Q_" + key: value for key, value in result["flows_L_min"].items()},
            **{"S_" + key: ox["lobes"][key]["S_end"] for key in result["flows_L_min"]},
        })
    pd.DataFrame(rows).to_csv(OUT / "scenario_summary.csv", index=False)
    payload = {
        "model_type": "steady_reduced_resistance_network",
        "primary_disease_operator": "distal_arterial_resistance_multiplier",
        "scenarios": scenarios, "oxygen": oxygen_results,
        "lobe_weight": LOBE_W, "cap_all_WU": case["cap_all"] / WU,
        "n_root_order12": case["n_root"],
    }
    (OUT / "results.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    metadata = {
        "python": platform.python_version(), "numpy": np.__version__,
        "pandas": pd.__version__, "model_type": payload["model_type"],
        "mixing_basis": "oxygen_content",
    }
    (OUT / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(rows, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
