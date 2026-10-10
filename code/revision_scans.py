#!/usr/bin/env python3
"""Analyses added in response to the v04 manuscript review.

The script generates data tables for oxygen conservation, contrast definition,
regional-perfusion sensitivity, exchange-form calibration, numerical
convergence, sampling stability, parameter provenance, and simple independent
benchmarks.  It also creates the main regional-sensitivity figure without a
matplotlib dependency.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
DATA.mkdir(exist_ok=True)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


lung = load_module("lung_revision", HERE / "build_lung_model.py")
phase3 = load_module("phase3_revision", HERE / "phase3_scans.py")

CONTRASTS = (1.0, 2.0, 4.0, 10.0, 25.0, 60.0)
TAU_VALUES = (0.15, 0.30, 0.45, 0.60)
VOLUME_RULES = ("fixed", "hybrid", "flow")
TARGET_PVR = 4.19


def target_uniform(case):
    multiplier = phase3.arterial_multiplier_for(case, TARGET_PVR)
    return lung.solve(case, arterial_multiplier=multiplier)


def contrast_result(case, contrast):
    if contrast == 1.0:
        return target_uniform(case)
    return phase3.match_contrast(case, contrast, target_pvr=TARGET_PVR)


def oxygen_conservation(case):
    result = target_uniform(case)
    ventilation = {
        lobe: 5.0 * lung.LOBE_W[lobe] for lobe in lung.LOBE_ORDER
    }
    ventilation["RLL"] *= 0.30
    rows = []
    for boundary, va in (("fixed_alveolar_PO2", None), ("low_RLL_ventilation", ventilation)):
        for fraction in (0.0, 0.03, 0.10, 0.20):
            ox = lung.oxygen(
                result,
                venous_admixture=fraction,
                lobe_va_L_min=va,
                vc_mode="fixed",
            )
            rows.append({
                "boundary": boundary,
                "venous_admixture_fraction": fraction,
                "systemic_Fick_VO2_mL_min": result["CO_L_min"] * 10.0 * (ox["CaO2"] - ox["CvO2"]),
                "pulmonary_effective_uptake_mL_min": ox["pulmonary_oxygen_uptake_mL_min"],
                "organ_actual_VO2_mL_min": ox["total_VO2_used_mL_min"],
                "pulmonary_systemic_residual_mL_min": ox["pulmonary_systemic_balance_error_mL_min"],
                "organ_systemic_residual_mL_min": ox["oxygen_mass_balance_error_mL_min"],
                "SaO2": ox["SaO2"],
                "converged": ox["converged"],
                "outer_iterations": ox["iterations"],
                "ventilation_converged": ox["ventilation_converged"],
                "ventilation_iterations": ox["ventilation_iterations"],
                "ventilation_residual_mmHg": ox["ventilation_residual_mmHg"],
                "alveolar_bound_hits": ox["ventilation_bound_hits"],
            })
    frame = pd.DataFrame(rows)
    frame.to_csv(DATA / "revision_oxygen_conservation.csv", index=False)
    return frame


def contrast_definition(case):
    rows = []
    lul_reference = None
    for contrast in CONTRASTS:
        result = contrast_result(case, contrast)
        multipliers = result["lobe_arterial_multiplier"]
        if contrast == 1.0:
            scale = multipliers["LUL"]
        else:
            scale = multipliers["LUL"]
        row = {
            "distal_multiplier_contrast_k": contrast,
            "whole_lung_scale_s": scale,
            "PVR_WU": result["PVR_WU"],
            "CO_L_min": result["CO_L_min"],
        }
        for lobe in lung.LOBE_ORDER:
            row[f"multiplier_{lobe}"] = multipliers[lobe]
            row[f"lobe_R_WU_{lobe}"] = result["lobe_R_WU"][lobe]
            row[f"flow_fraction_{lobe}"] = result["flow_fraction"][lobe]
        if contrast == 25.0:
            lul_reference = result["lobe_R_WU"]["LUL"] * lung.LOBE_W["LUL"]
            for lobe in lung.LOBE_ORDER:
                row[f"weight_standardized_R_ratio_{lobe}_to_LUL"] = (
                    (result["lobe_R_WU"][lobe] * lung.LOBE_W[lobe]) / lul_reference
                )
        rows.append(row)
    frame = pd.DataFrame(rows)
    frame.to_csv(DATA / "revision_contrast_definition.csv", index=False)
    return frame


def core_sensitivity(case):
    results = {contrast: contrast_result(case, contrast) for contrast in CONTRASTS}
    rows = []
    for volume_rule in VOLUME_RULES:
        for tau in TAU_VALUES:
            baseline = lung.oxygen(
                results[1.0], vc_ml=86.0, tau_eq=tau,
                vc_mode=volume_rule, venous_admixture=0.03,
            )
            for contrast in CONTRASTS:
                result = results[contrast]
                ox = lung.oxygen(
                    result, vc_ml=86.0, tau_eq=tau,
                    vc_mode=volume_rule, venous_admixture=0.03,
                )
                taus = [ox["lobes"][lobe]["tau_s"] for lobe in lung.LOBE_ORDER]
                rows.append({
                    "volume_rule": volume_rule,
                    "tau_eq_s": tau,
                    "distal_multiplier_contrast_k": contrast,
                    "PVR_WU": result["PVR_WU"],
                    "LUL_flow_fraction": result["flow_fraction"]["LUL"],
                    "SaO2": ox["SaO2"],
                    "delta_SaO2_vs_uniform": ox["SaO2"] - baseline["SaO2"],
                    "tau_min_s": min(taus),
                    "tau_max_s": max(taus),
                    "chi_min": min(taus) / tau,
                    "endcap_saturation_range": ox["S_end_range"],
                })
    frame = pd.DataFrame(rows)
    frame.to_csv(DATA / "revision_core_sensitivity.csv", index=False)
    draw_core_figure(frame, DATA / "fig_revision_core_sensitivity.png")
    return frame


def saturation_target_content(transit, tau=0.30, pv=40.0, pa=100.0, hb=15.0):
    s_in = float(lung.hill(pv))
    s_alv = float(lung.hill(pa))
    s_end = s_alv - (s_alv - s_in) * math.exp(-transit / tau)
    po2 = lung.po2_from_saturation(s_end)
    return float(lung.oxygen_content(po2, hb))


def content_curve(transit, tau, pv=40.0, pa=100.0, hb=15.0):
    cv = float(lung.oxygen_content(pv, hb))
    ca = float(lung.oxygen_content(pa, hb))
    return ca - (ca - cv) * math.exp(-transit / tau)


def diffusion_curve(transit, dl_o2, pv=40.0, pa=100.0, hb=15.0, vc_ml=86.0, steps=2048):
    content = float(lung.oxygen_content(pv, hb))
    cap = float(lung.oxygen_content(pa, hb))
    dt = transit / steps
    pressure = pv
    for _ in range(steps):
        increment = dl_o2 * (pa - pressure) / vc_ml * (100.0 / 60.0) * dt
        content = min(max(content + increment, 0.0), cap)
        pressure = lung.po2_from_content(content, hb)
    return content


def exchange_calibration(case):
    reference_times = np.array([0.30, 0.60, 1.20])
    targets = np.array([saturation_target_content(value) for value in reference_times])
    tau_grid = np.linspace(0.05, 1.20, 2301)
    tau_errors = [
        np.mean([(content_curve(t, tau) - y) ** 2 for t, y in zip(reference_times, targets)])
        for tau in tau_grid
    ]
    tau_content = float(tau_grid[int(np.argmin(tau_errors))])
    dl_grid = np.linspace(0.5, 80.0, 796)
    dl_errors = [
        np.mean([(diffusion_curve(t, dl, steps=512) - y) ** 2 for t, y in zip(reference_times, targets)])
        for dl in dl_grid
    ]
    dl_fit = float(dl_grid[int(np.argmin(dl_errors))])

    rows = []
    for transit, target in zip(reference_times, targets):
        for mode, value, parameter, parameter_value in (
            ("saturation_exponential", target, "tau_eq_s", 0.30),
            ("content_exponential", content_curve(transit, tau_content), "tau_eq_s", tau_content),
            ("fixed_capacity", diffusion_curve(transit, dl_fit), "DL_mL_min_mmHg", dl_fit),
        ):
            rows.append({
                "transit_s": transit,
                "exchange_mode": mode,
                "calibrated_parameter": parameter,
                "parameter_value": parameter_value,
                "end_content_mL_dL": value,
                "target_end_content_mL_dL": target,
                "absolute_error_mL_dL": abs(value - target),
            })
    calibration = pd.DataFrame(rows)
    calibration.to_csv(DATA / "revision_exchange_calibration.csv", index=False)

    scenarios = []
    for contrast in (1.0, 4.0, 10.0, 25.0, 60.0):
        result = contrast_result(case, contrast)
        for mode, tau, dl in (
            ("saturation_exponential", 0.30, 25.0),
            ("content_exponential", tau_content, 25.0),
            ("fixed_capacity", 0.30, dl_fit),
        ):
            exchange_mode = {
                "saturation_exponential": "saturation",
                "content_exponential": "content",
                "fixed_capacity": "diffusion",
            }[mode]
            ox = lung.oxygen(
                result,
                vc_ml=86.0,
                vc_mode="fixed",
                tau_eq=tau,
                exchange_mode=exchange_mode,
                dl_o2=dl,
                diffusion_steps=1024,
            )
            scenarios.append({
                "distal_multiplier_contrast_k": contrast,
                "exchange_mode": mode,
                "tau_eq_s": tau,
                "DL_mL_min_mmHg": dl,
                "SaO2": ox["SaO2"],
                "endcap_saturation_range": ox["S_end_range"],
            })
    pd.DataFrame(scenarios).to_csv(DATA / "revision_exchange_model_sensitivity.csv", index=False)
    return tau_content, dl_fit


def numerical_convergence(case):
    result = contrast_result(case, 25.0)
    rows = []
    previous = None
    for steps in (16, 64, 256, 1024, 2048):
        ox = lung.oxygen(
            result,
            vc_ml=86.0,
            vc_mode="fixed",
            exchange_mode="diffusion",
            dl_o2=25.0,
            diffusion_steps=steps,
        )
        rows.append({
            "integration_steps": steps,
            "SaO2": ox["SaO2"],
            "endcap_saturation_range": ox["S_end_range"],
            "change_from_previous": np.nan if previous is None else ox["SaO2"] - previous,
            "pulmonary_systemic_residual_mL_min": ox["pulmonary_systemic_balance_error_mL_min"],
        })
        previous = ox["SaO2"]
    frame = pd.DataFrame(rows)
    frame.to_csv(DATA / "revision_numerical_convergence.csv", index=False)
    return frame


def latin_hypercube(count, dimensions, seed):
    rng = np.random.default_rng(seed)
    design = np.empty((count, dimensions))
    for column in range(dimensions):
        strata = (np.arange(count) + rng.random(count)) / count
        design[:, column] = strata[rng.permutation(count)]
    return design


def sampling_stability():
    lows = np.array([12.0, 43.0, 0.15, 80.0, 0.02, 40.0, 2.5])
    highs = np.array([18.0, 172.0, 0.60, 110.0, 0.08, 60.0, 4.5])
    rows = []
    for seed in (20261008, 20261009, 20261010):
        for count in (200, 500, 1000):
            samples = lows + (highs - lows) * latin_hypercube(count, 7, seed)
            aggregates = {target: [] for target in (3.0, 4.19, 6.0, 10.0)}
            for hb, vc, tau, pao2, admixture, reserve, viscosity in samples:
                local = lung.build_cases(mu=viscosity / 1000.0)
                healthy = lung.solve(local, reserve=reserve)
                healthy_ox = lung.oxygen(
                    healthy, hb=hb, vc_ml=vc, tau_eq=tau, pao2=pao2,
                    venous_admixture=admixture,
                )
                for target in aggregates:
                    multiplier = phase3.arterial_multiplier_for(local, target)
                    result = lung.solve(local, arterial_multiplier=multiplier, reserve=reserve)
                    ox = lung.oxygen(
                        result, hb=hb, vc_ml=vc, tau_eq=tau, pao2=pao2,
                        venous_admixture=admixture,
                    )
                    aggregates[target].append((
                        ox["SaO2"] - healthy_ox["SaO2"],
                        ox["SvO2"] - healthy_ox["SvO2"],
                    ))
            for target, values in aggregates.items():
                array = np.asarray(values)
                rows.append({
                    "seed": seed,
                    "sample_count": count,
                    "target_PVR_WU": target,
                    "SaO2_not_lower_fraction_tol_1e-6": np.mean(array[:, 0] >= -1e-6),
                    "SaO2_not_lower_fraction_tol_1e-5": np.mean(array[:, 0] >= -1e-5),
                    "SvO2_lower_fraction_tol_1e-4": np.mean(array[:, 1] < -1e-4),
                    "SvO2_lower_fraction_tol_1e-3": np.mean(array[:, 1] < -1e-3),
                    "median_delta_SaO2": np.median(array[:, 0]),
                    "median_delta_SvO2": np.median(array[:, 1]),
                })
    frame = pd.DataFrame(rows)
    frame.to_csv(DATA / "revision_sampling_stability.csv", index=False)
    return frame


def parameter_tables(case):
    parameters = [
        ("Hemoglobin", 15.0, 12.0, 18.0, "g/dL", "literature_direct", "Oxygen-content equation"),
        ("Total capillary blood volume", 86.0, 43.0, 172.0, "mL", "order_of_magnitude", "Scenario value based on resting estimates"),
        ("Saturation exchange time constant", 0.30, 0.15, 0.60, "s", "investigator_set", "Mechanism parameter, not measured transit time"),
        ("Fixed alveolar oxygen tension", 100.0, 80.0, 110.0, "mmHg", "investigator_set", "Prescribed boundary"),
        ("Effective venous admixture fraction", 0.03, 0.02, 0.08, "fraction", "investigator_set", "Post-exchange content mixing, not resolved anatomical shunt"),
        ("Pressure reserve width", 50.0, 40.0, 60.0, "mmHg", "investigator_set", "Prescribed pressure-flow relation"),
        ("Dynamic viscosity", 3.5, 2.5, 4.5, "mPa s", "literature_range", "Newtonian model"),
        ("Whole-body oxygen demand", 250.0, np.nan, np.nan, "mL/min", "literature_order", "Not varied in sampling"),
        ("Capacity-control parameter", 25.0, np.nan, np.nan, "mL/min/mmHg", "order_of_magnitude", "Numerical structure control"),
        ("Total alveolar ventilation", 5.0, np.nan, np.nan, "L/min BTPS", "investigator_set", "Ventilation scenarios"),
        ("Inspired oxygen tension", 150.0, np.nan, np.nan, "mmHg", "investigator_set", "Ventilation scenarios"),
        ("Reference cardiac-output scale", 5.0, np.nan, np.nan, "L/min", "investigator_set", "Scale parameter, not a hard maximum"),
        ("Minimum cardiac output", 1.2, np.nan, np.nan, "L/min", "investigator_set", "Lower clamp"),
        ("Branching ratio", 3.0, np.nan, np.nan, "dimensionless", "calibration", "Symmetric structured tree"),
    ]
    pd.DataFrame(parameters, columns=[
        "parameter", "baseline", "range_low", "range_high", "unit",
        "provenance_class", "interpretation",
    ]).to_csv(DATA / "revision_model_parameters.csv", index=False)

    organ_rows = []
    for organ, (flow, demand, extraction) in lung.DEFAULT_ORGANS.items():
        organ_rows.append({
            "organ": organ,
            "cardiac_output_fraction": flow,
            "oxygen_demand_fraction": demand,
            "maximum_extraction_fraction": extraction,
            "provenance_class": "investigator_set",
        })
    pd.DataFrame(organ_rows).to_csv(DATA / "revision_organ_parameters.csv", index=False)

    healthy = lung.solve(case)
    diagnostics = pd.DataFrame([{
        "approximate_large_plus_tree_target_WU": 0.85,
        "capillary_term_WU": case["cap_all"] / lung.WU,
        "complete_network_PVR_WU": healthy["PVR_WU"],
        "equivalent_order12_root_count": case["n_root"],
        "CO_reference_scale_L_min": healthy["CO_reference_L_min"],
        "solved_healthy_CO_L_min": healthy["CO_L_min"],
    }])
    diagnostics.to_csv(DATA / "revision_calibration_diagnostics.csv", index=False)


def validation_benchmarks(case):
    uniform = target_uniform(case)
    no_admixture = lung.oxygen(
        uniform, tau_eq=0.03, venous_admixture=0.0, vc_mode="fixed"
    )
    alveolar_saturation = float(lung.hill(100.0))
    flow_volume = lung.oxygen(
        uniform, tau_eq=0.30, venous_admixture=0.0, vc_mode="flow"
    )
    taus = [flow_volume["lobes"][lobe]["tau_s"] for lobe in lung.LOBE_ORDER]
    rows = [
        {
            "benchmark": "no_admixture_fast_exchange_upper_bound",
            "observed": no_admixture["SaO2"],
            "expected": alveolar_saturation,
            "absolute_error": abs(no_admixture["SaO2"] - alveolar_saturation),
        },
        {
            "benchmark": "flow_proportional_volume_equal_transit",
            "observed": max(taus) - min(taus),
            "expected": 0.0,
            "absolute_error": max(taus) - min(taus),
        },
        {
            "benchmark": "pulmonary_systemic_oxygen_conservation",
            "observed": abs(no_admixture["pulmonary_systemic_balance_error_mL_min"]),
            "expected": 0.0,
            "absolute_error": abs(no_admixture["pulmonary_systemic_balance_error_mL_min"]),
        },
    ]
    pd.DataFrame(rows).to_csv(DATA / "revision_validation_benchmarks.csv", index=False)


def draw_core_figure(frame: pd.DataFrame, path: Path):
    width, height = 2100, 700
    margin_left, margin_top, margin_bottom = 120, 95, 110
    panel_gap = 70
    panel_width = (width - margin_left - 80 - 2 * panel_gap) // 3
    plot_height = height - margin_top - margin_bottom
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    font_path = Path("C:/Windows/Fonts/arial.ttf")
    bold_path = Path("C:/Windows/Fonts/arialbd.ttf")
    font = ImageFont.truetype(str(font_path), 24)
    small = ImageFont.truetype(str(font_path), 20)
    bold = ImageFont.truetype(str(bold_path), 28)
    colors = {0.15: "#1f77b4", 0.30: "#d62728", 0.45: "#2ca02c", 0.60: "#9467bd"}
    y_min, y_max = -24.0, 2.0
    labels = {"fixed": "Fixed anatomical volume", "hybrid": "Hybrid allocation", "flow": "Flow-proportional volume"}

    for panel_index, rule in enumerate(VOLUME_RULES):
        x0 = margin_left + panel_index * (panel_width + panel_gap)
        y0 = margin_top
        x1 = x0 + panel_width
        y1 = y0 + plot_height
        draw.rectangle((x0, y0, x1, y1), outline="#444444", width=2)
        for y_tick in (-20, -15, -10, -5, 0):
            y = y1 - (y_tick - y_min) / (y_max - y_min) * plot_height
            draw.line((x0, y, x1, y), fill="#E5E7EB", width=1)
            if panel_index == 0:
                draw.text((x0 - 65, y - 12), str(y_tick), fill="#333333", font=small)
        x_positions = {}
        for index, contrast in enumerate(CONTRASTS):
            x = x0 + 35 + index * (panel_width - 70) / (len(CONTRASTS) - 1)
            x_positions[contrast] = x
            draw.line((x, y1, x, y1 + 8), fill="#333333", width=2)
            label = f"{contrast:g}"
            box = draw.textbbox((0, 0), label, font=small)
            draw.text((x - (box[2] - box[0]) / 2, y1 + 15), label, fill="#333333", font=small)
        draw.text((x0 + 10, 30), labels[rule], fill="#111111", font=bold)
        for tau in TAU_VALUES:
            subset = frame[(frame["volume_rule"] == rule) & (frame["tau_eq_s"] == tau)].sort_values("distal_multiplier_contrast_k")
            points = []
            for _, row in subset.iterrows():
                x = x_positions[row["distal_multiplier_contrast_k"]]
                value = 100.0 * row["delta_SaO2_vs_uniform"]
                y = y1 - (value - y_min) / (y_max - y_min) * plot_height
                points.append((x, y))
            draw.line(points, fill=colors[tau], width=4)
            for point in points:
                draw.ellipse((point[0] - 6, point[1] - 6, point[0] + 6, point[1] + 6), fill=colors[tau], outline="white", width=2)
    draw.text((width // 2 - 230, height - 48), "Distal arterial multiplier contrast k", fill="#111111", font=font)
    axis_label = Image.new("RGBA", (520, 45), (255, 255, 255, 0))
    axis_draw = ImageDraw.Draw(axis_label)
    axis_draw.text((0, 5), "Delta SaO2 vs uniform (percentage points)", fill="#111111", font=font)
    axis_label = axis_label.rotate(90, expand=True)
    image.paste(axis_label, (10, 110), axis_label)
    legend_x = width - 560
    for index, tau in enumerate(TAU_VALUES):
        y = 250 + index * 34
        draw.line((legend_x, y + 10, legend_x + 45, y + 10), fill=colors[tau], width=4)
        draw.text((legend_x + 55, y), f"tau_eq = {tau:.2f} s", fill="#111111", font=small)
    image.save(path, dpi=(300, 300))


def draw_vq_figure(frame: pd.DataFrame, path: Path):
    labels = list("ABCDEFGH")
    width, height = 1800, 720
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    font_path = Path("C:/Windows/Fonts/arial.ttf")
    bold_path = Path("C:/Windows/Fonts/arialbd.ttf")
    font = ImageFont.truetype(str(font_path), 24)
    small = ImageFont.truetype(str(font_path), 20)
    bold = ImageFont.truetype(str(bold_path), 28)
    panels = [
        ("SaO2", 0.90, 0.97, "Arterial oxygen saturation"),
        ("systemic_DO2_mL_min", 700.0, 1020.0, "Systemic oxygen delivery (mL/min)"),
    ]
    for panel_index, (column, y_min, y_max, title) in enumerate(panels):
        x0 = 105 + panel_index * 860
        y0, x1, y1 = 95, x0 + 730, 575
        draw.rectangle((x0, y0, x1, y1), outline="#444444", width=2)
        draw.text((x0 + 15, 32), title, fill="#111111", font=bold)
        for tick_index in range(6):
            value = y_min + tick_index * (y_max - y_min) / 5
            y = y1 - tick_index * (y1 - y0) / 5
            draw.line((x0, y, x1, y), fill="#E5E7EB", width=1)
            label = f"{value:.2f}" if column == "SaO2" else f"{value:.0f}"
            draw.text((x0 - 72, y - 12), label, fill="#333333", font=small)
        bar_space = (x1 - x0 - 60) / len(frame)
        for index, (_, row) in enumerate(frame.iterrows()):
            value = float(row[column])
            x_left = x0 + 35 + index * bar_space
            x_right = x_left + bar_space * 0.62
            y = y1 - (value - y_min) / (y_max - y_min) * (y1 - y0)
            color = "#4C78A8" if index < 4 else "#E45756"
            draw.rectangle((x_left, y, x_right, y1), fill=color)
            draw.text((x_left + 4, y1 + 15), labels[index], fill="#111111", font=font)
            value_label = f"{value:.3f}" if column == "SaO2" else f"{value:.0f}"
            draw.text((x_left - 2, y - 30), value_label, fill="#333333", font=small)
    draw.text((135, 640), "Blue: healthy-resistance scenarios    Red: PVR 4.19 WU background", fill="#333333", font=font)
    image.save(path, dpi=(300, 300))


def draw_sampling_figure(frame: pd.DataFrame, path: Path):
    focus = frame[frame["sample_count"] == 1000].copy()
    targets = sorted(focus["target_PVR_WU"].unique())
    width, height = 1800, 700
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    font_path = Path("C:/Windows/Fonts/arial.ttf")
    bold_path = Path("C:/Windows/Fonts/arialbd.ttf")
    font = ImageFont.truetype(str(font_path), 24)
    small = ImageFont.truetype(str(font_path), 20)
    bold = ImageFont.truetype(str(bold_path), 28)
    panels = [
        ("SaO2_not_lower_fraction_tol_1e-6", 0.50, 0.72, "SaO2 did not decrease"),
        ("SvO2_lower_fraction_tol_1e-4", 0.95, 1.005, "SvO2 decreased"),
    ]
    seed_colors = {20261008: "#1f77b4", 20261009: "#d62728", 20261010: "#2ca02c"}
    for panel_index, (column, y_min, y_max, title) in enumerate(panels):
        x0 = 105 + panel_index * 860
        y0, x1, y1 = 95, x0 + 730, 565
        draw.rectangle((x0, y0, x1, y1), outline="#444444", width=2)
        draw.text((x0 + 15, 32), title, fill="#111111", font=bold)
        for tick_index in range(6):
            value = y_min + tick_index * (y_max - y_min) / 5
            y = y1 - tick_index * (y1 - y0) / 5
            draw.line((x0, y, x1, y), fill="#E5E7EB", width=1)
            draw.text((x0 - 70, y - 12), f"{value:.2f}", fill="#333333", font=small)
        x_positions = {
            target: x0 + 55 + index * (x1 - x0 - 110) / (len(targets) - 1)
            for index, target in enumerate(targets)
        }
        for target, x in x_positions.items():
            draw.line((x, y1, x, y1 + 8), fill="#333333", width=2)
            label = f"{target:g}"
            draw.text((x - 14, y1 + 15), label, fill="#333333", font=small)
        for seed, color in seed_colors.items():
            subset = focus[focus["seed"] == seed].sort_values("target_PVR_WU")
            points = []
            for _, row in subset.iterrows():
                x = x_positions[row["target_PVR_WU"]]
                value = float(row[column])
                y = y1 - (value - y_min) / (y_max - y_min) * (y1 - y0)
                points.append((x, y))
            draw.line(points, fill=color, width=4)
            for point in points:
                draw.ellipse((point[0] - 6, point[1] - 6, point[0] + 6, point[1] + 6), fill=color, outline="white", width=2)
        draw.text((x0 + 265, height - 55), "Target PVR (WU)", fill="#111111", font=font)
    legend_x = 680
    for index, (seed, color) in enumerate(seed_colors.items()):
        x = legend_x + index * 235
        draw.line((x, 625, x + 45, 625), fill=color, width=4)
        draw.text((x + 55, 612), str(seed), fill="#111111", font=small)
    image.save(path, dpi=(300, 300))


def main():
    case = lung.build_cases()
    oxygen_conservation(case)
    contrast_definition(case)
    core_sensitivity(case)
    exchange_calibration(case)
    numerical_convergence(case)
    sampling = sampling_stability()
    parameter_tables(case)
    validation_benchmarks(case)
    draw_vq_figure(pd.read_csv(DATA / "goal10_vq.csv"), DATA / "fig_revision_vq.png")
    draw_sampling_figure(sampling, DATA / "fig_revision_sampling_stability.png")
    print("Revision analyses completed")


if __name__ == "__main__":
    main()

