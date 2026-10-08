#!/usr/bin/env python3
"""Phase-2 analyses for the steady reduced pulmonary resistance network.

Primary disease cases increase distal arterial resistance only. Loss of
perfused vascular bed is treated separately as a spatial limiting experiment.
All blood streams are mixed by oxygen content and the organ-to-lung oxygen
cycle is solved to mass balance.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

SRC = Path(__file__).resolve().parent / "build_lung_model.py"
SPEC = importlib.util.spec_from_file_location("lung", SRC)
lung = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lung)

OUT = Path(__file__).resolve().parent
LOBE_ORDER = ["RUL", "RML", "RLL", "LUL", "LLL"]
LOBE_CN = {
    "RUL": "右上叶", "RML": "右中叶", "RLL": "右下叶",
    "LUL": "左上叶", "LLL": "左下叶",
}
ORGANS = lung.DEFAULT_ORGANS


def network(case, **kwargs):
    return lung.solve(case, **kwargs)


def oxygen(result, **kwargs):
    return lung.oxygen(result, **kwargs)


def latin_hypercube(sample_count, dimension_count, seed):
    """Generate a reproducible Latin-hypercube design on the unit interval."""
    rng = np.random.default_rng(seed)
    design = np.empty((sample_count, dimension_count), dtype=float)
    for column in range(dimension_count):
        strata = (np.arange(sample_count) + rng.random(sample_count)) / sample_count
        design[:, column] = strata[rng.permutation(sample_count)]
    return design


def pack(name, result, ox):
    row = {
        "scenario": name,
        "CO_L_min": result["CO_L_min"],
        "mPAP_mmHg": result["mPAP_mmHg"],
        "PVR_WU": result["PVR_WU"],
        "arterial_multiplier": result["arterial_multiplier"],
        "venous_multiplier": result["venous_multiplier"],
        "SaO2": ox["SaO2"],
        "SvO2": ox["SvO2"],
        "PaO2": ox["PaO2"],
        "PvO2": ox["PvO2"],
        "CaO2_mL_dL": ox["CaO2"],
        "CvO2_mL_dL": ox["CvO2"],
        "S_end_var": ox["S_end_var"],
        "S_end_range": ox["S_end_range"],
        "effective_vc_ml": ox["effective_vc_ml"],
        "total_VO2_demand_mL_min": ox["total_VO2_demand_mL_min"],
        "total_VO2_used_mL_min": ox["total_VO2_used_mL_min"],
        "total_deficit_mL_min": ox["total_deficit_mL_min"],
        "oxygen_mass_balance_error_mL_min": ox["oxygen_mass_balance_error_mL_min"],
        "oxygen_converged": ox["converged"],
        "oxygen_iterations": ox["iterations"],
    }
    for lobe in LOBE_ORDER:
        row["f_" + lobe] = result["flow_fraction"][lobe]
        row["S_" + lobe] = ox["lobes"][lobe]["S_end"]
        row["tau_" + lobe] = ox["lobes"][lobe]["tau_s"]
        row["Vc_" + lobe] = ox["lobes"][lobe]["vc_ml"]
    for organ, values in ox["organs"].items():
        row["E_" + organ] = values["extraction"]
        row["def_" + organ] = values["deficit_mL_min"]
        row["DO2_" + organ] = values["DO2_mL_min"]
        row["dem_" + organ] = values["VO2_demand_mL_min"]
        row["used_" + organ] = values["VO2_mL_min"]
    return row


def arterial_multiplier_for(case, target_pvr):
    if target_pvr < network(case)["PVR_WU"]:
        raise ValueError("target PVR is below the calibrated healthy PVR")
    low, high = 1.0, 2.0
    while network(case, arterial_multiplier=high)["PVR_WU"] < target_pvr:
        high *= 2.0
        if high > 10000.0:
            raise RuntimeError("could not bracket target PVR")
    for _ in range(60):
        middle = 0.5 * (low + high)
        if network(case, arterial_multiplier=middle)["PVR_WU"] > target_pvr:
            high = middle
        else:
            low = middle
    return 0.5 * (low + high)


def match_open(case, lobes, target_pvr, scale_cap, low=0.001, high=1.0):
    def evaluate(open_fraction):
        focal = {lobe: open_fraction for lobe in lobes}
        return network(case, focal=focal, scale_cap=scale_cap)

    low_result, high_result = evaluate(low), evaluate(high)
    if low_result["PVR_WU"] < target_pvr:
        return {
            "reached": False, "open_fraction": low,
            "PVR_WU": low_result["PVR_WU"], "result": low_result,
        }
    if high_result["PVR_WU"] > target_pvr:
        return {
            "reached": False, "open_fraction": high,
            "PVR_WU": high_result["PVR_WU"], "result": high_result,
        }
    for _ in range(60):
        middle = 0.5 * (low + high)
        result = evaluate(middle)
        if result["PVR_WU"] > target_pvr:
            low = middle
        else:
            high = middle
    open_fraction = 0.5 * (low + high)
    result = evaluate(open_fraction)
    return {
        "reached": abs(result["PVR_WU"] - target_pvr) < 1e-8,
        "open_fraction": open_fraction,
        "PVR_WU": result["PVR_WU"],
        "result": result,
    }


def formal_focal_match(case, target_pvr):
    """Extreme four-lobe bed-loss pattern used only as a limiting experiment."""
    low, high = 0.05, 1.0

    def evaluate(spared_fraction):
        return network(
            case,
            focal={
                "RUL": 0.005, "RML": 0.005, "RLL": 0.005,
                "LLL": 0.005, "LUL": spared_fraction,
            },
            scale_cap=True,
        )

    if evaluate(low)["PVR_WU"] < target_pvr or evaluate(high)["PVR_WU"] > target_pvr:
        raise RuntimeError("formal focal pattern does not bracket the target PVR")
    for _ in range(60):
        middle = 0.5 * (low + high)
        if evaluate(middle)["PVR_WU"] > target_pvr:
            low = middle
        else:
            high = middle
    spared_fraction = 0.5 * (low + high)
    result = evaluate(spared_fraction)
    return {
        "reached": abs(result["PVR_WU"] - target_pvr) < 1e-8,
        "open_fraction": spared_fraction,
        "PVR_WU": result["PVR_WU"],
        "result": result,
    }


def verification_rows(case, scenarios):
    rows = []
    for name, result, ox in scenarios:
        checks = [
            (
                "flow_conservation",
                abs(result["flow_conservation_error_L_min"]),
                1e-9,
                "sum of five lobar flows minus cardiac output, L/min",
            ),
            (
                "pvr_identity",
                abs(result["pvr_identity_error_WU"]),
                1e-12,
                "(mPAP-PAWP)/CO minus network PVR, WU",
            ),
            (
                "oxygen_mass_balance",
                abs(ox["oxygen_mass_balance_error_mL_min"]),
                1e-6,
                "systemic extraction minus summed organ oxygen use, mL/min",
            ),
            (
                "oxygen_convergence",
                0.0 if ox["converged"] else 1.0,
                0.0,
                "closed oxygen-content iteration converged",
            ),
        ]
        for check, value, tolerance, definition in checks:
            rows.append({
                "scope": name, "check": check, "value": value,
                "tolerance": tolerance,
                "status": "PASS" if value <= tolerance else "FAIL",
                "definition": definition,
            })

    model_order1 = case["n_root"] * case["rb"] ** 11
    reference_order1 = 102.41e6
    rows.append({
        "scope": "形态学外部交叉核对",
        "check": "order1_artery_count",
        "value": model_order1,
        "tolerance": "",
        "status": "REFERENCE_ONLY",
        "definition": (
            "模型全肺order-1动脉数量；对照值为Huang 1996 Table 5"
            "左肺校正数量的两倍约102.41 million"
        ),
        "reference_value": reference_order1,
        "relative_difference": (model_order1 - reference_order1) / reference_order1,
    })
    return rows


def configure_matplotlib():
    try:
        import matplotlib.pyplot as plt
        from matplotlib import font_manager
    except ModuleNotFoundError:
        return None
    candidates = [
        "Microsoft YaHei", "SimHei", "Noto Sans CJK SC",
        "WenQuanYi Zen Hei", "Arial Unicode MS",
    ]
    installed = {font.name for font in font_manager.fontManager.ttflist}
    for candidate in candidates:
        if candidate in installed:
            plt.rcParams["font.sans-serif"] = [candidate]
            break
    plt.rcParams["axes.unicode_minus"] = False
    return plt


def make_figures(goal1, controlled, pao2_rows, shunt_rows, ladder, goal6):
    plt = configure_matplotlib()
    if plt is None:
        print("matplotlib is not installed; CSV/JSON outputs were written and figures were skipped.")
        return

    lookup = {row["scenario"]: row for row in goal1}
    keys = ["健康", "均匀动脉增阻_PVR4.19", "极限局灶床丢失_PVR4.19"]
    labels = ["健康", "均匀动脉增阻", "极限局灶床丢失"]
    x = np.arange(len(LOBE_ORDER))
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    width = 0.22
    for index, (key, label) in enumerate(zip(keys, labels)):
        ax.bar(
            x + (index - 1) * width,
            [lookup[key]["f_" + lobe] for lobe in LOBE_ORDER],
            width,
            label=label,
        )
    ax.set_xticks(x)
    ax.set_xticklabels([LOBE_CN[lobe] for lobe in LOBE_ORDER])
    ax.set_ylabel("流量份额")
    ax.set_title("相同PVR下的肺叶灌注分布")
    ax.set_ylim(0.0, 1.0)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_goal1_flow.png", dpi=180)
    plt.close(fig)

    primary_controlled = [row for row in controlled if row["analysis_role"] == "primary"]
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.0))
    labels_c = ["均匀", "极限局灶"]
    axes[0].bar(labels_c, [row["S_end_range"] for row in primary_controlled], color=["#4C78A8", "#F58518"])
    axes[0].set_ylabel("五叶末端饱和度极差")
    axes[0].set_title("等PVR且等总毛细血管容积")
    axes[1].bar(labels_c, [row["SaO2"] for row in primary_controlled], color=["#4C78A8", "#F58518"])
    axes[1].set_ylabel("SaO2")
    axes[1].set_ylim(
        min(row["SaO2"] for row in primary_controlled) - 0.01,
        max(row["SaO2"] for row in primary_controlled) + 0.01,
    )
    axes[1].set_title("全身动脉氧合")
    fig.tight_layout()
    fig.savefig(OUT / "fig_goalB_controlled.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.0))
    axes[0].plot(
        [row["RLL_PAO2"] for row in pao2_rows],
        [row["SaO2"] for row in pao2_rows], "o-",
    )
    axes[0].set_xlabel("右下叶肺泡氧分压 / mmHg")
    axes[0].set_ylabel("SaO2")
    axes[0].set_title("单叶低肺泡氧")
    axes[0].grid(True, alpha=0.3)
    axes[1].plot(
        [row["shunt"] for row in shunt_rows],
        [row["SaO2"] for row in shunt_rows], "s-",
    )
    axes[1].set_xlabel("假设性右向左分流比例")
    axes[1].set_ylabel("SaO2")
    axes[1].set_title("分流敏感性")
    axes[1].grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "fig_goal3_hypoxemia.png", dpi=180)
    plt.close(fig)

    ladder_sorted = sorted(ladder, key=lambda item: item["PVR_WU"])
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.7))
    x_pvr = [row["PVR_WU"] for row in ladder_sorted]
    series = [
        ("CO_L_min", "心输出量 / L min-1"),
        ("SvO2", "混合静脉氧饱和度"),
        ("SaO2", "动脉氧饱和度"),
    ]
    for axis, (key, label) in zip(axes, series):
        axis.plot(x_pvr, [row[key] for row in ladder_sorted], "o-")
        axis.set_xlabel("PVR / WU")
        axis.set_ylabel(label)
        axis.grid(True, alpha=0.3)
    fig.suptitle("远端动脉阻力阶梯")
    fig.tight_layout()
    fig.savefig(OUT / "fig_goal4_ladder.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    labels_6 = [row["scenario"].split("|")[0] for row in goal6]
    ax.plot(labels_6, [row["rv_demand"] for row in goal6], "o-", label="右室氧耗需求")
    ax.plot(
        labels_6,
        [row["used_右心室心肌"] for row in goal6],
        "s-",
        label="闭环实际氧耗",
    )
    ax.set_ylabel("mL/min")
    ax.set_title("预设器官分配下的右室供需情景")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "fig_goal6_rv.png", dpi=180)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args(argv)

    case = lung.build_cases()
    healthy = network(case)
    target_pvr = 4.19
    target_multiplier = arterial_multiplier_for(case, target_pvr)
    uniform_target = network(case, arterial_multiplier=target_multiplier)
    uniform_x5 = network(case, arterial_multiplier=5.0)
    focal = formal_focal_match(case, target_pvr)
    focal_result = focal["result"]
    reference_power = healthy["mPAP_mmHg"] * healthy["CO_L_min"]

    patterns = {
        "仅右下叶": ["RLL"],
        "双下叶": ["RLL", "LLL"],
        "右肺三叶": ["RUL", "RML", "RLL"],
        "除左上叶外四叶": ["RUL", "RML", "RLL", "LLL"],
    }
    ceilings = {
        name: match_open(case, lobes, target_pvr, False)
        for name, lobes in patterns.items()
    }
    matched = {
        name: match_open(case, lobes, target_pvr, True)
        for name, lobes in patterns.items()
    }

    named_goal1 = {
        "健康": healthy,
        "远端动脉阻力x5": uniform_x5,
        "均匀动脉增阻_PVR4.19": uniform_target,
        "仅右下叶_床不撤": ceilings["仅右下叶"]["result"],
        "四叶_床不撤": ceilings["除左上叶外四叶"]["result"],
        "极限局灶床丢失_PVR4.19": focal_result,
    }
    goal1 = [
        pack(name, result, oxygen(result, vc_mode="anatomy"))
        for name, result in named_goal1.items()
    ]
    pd.DataFrame(goal1).to_csv(OUT / "goal1_pvr_match.csv", index=False)

    tau_rows, vc_rows = [], []
    for tau_eq in (0.30, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0):
        ox = oxygen(healthy, tau_eq=tau_eq, vc_mode="anatomy")
        tau_rows.append({
            "tau_eq_s": tau_eq, "SaO2": ox["SaO2"], "SvO2": ox["SvO2"],
            "PaO2": ox["PaO2"], "S_end_range": ox["S_end_range"],
            "effective_vc_ml": ox["effective_vc_ml"],
        })
    for vc_ml in (86, 60, 40, 25, 15, 8, 5):
        ox = oxygen(healthy, vc_ml=vc_ml, vc_mode="anatomy")
        vc_rows.append({
            "vc_ml": vc_ml, "SaO2": ox["SaO2"], "SvO2": ox["SvO2"],
            "PaO2": ox["PaO2"], "S_end_range": ox["S_end_range"],
            "effective_vc_ml": ox["effective_vc_ml"],
        })
    pd.DataFrame(tau_rows).to_csv(OUT / "goal2_tau_eq.csv", index=False)
    pd.DataFrame(vc_rows).to_csv(OUT / "goal2_vc.csv", index=False)

    focal_effective_vc = 86.0 * sum(
        lung.LOBE_W[lobe] * focal_result["open_fraction"][lobe]
        for lobe in LOBE_ORDER
    )
    controlled_specs = [
        (
            "均匀动脉增阻_等PVR等总Vc_基线交换",
            uniform_target,
            0.30,
            "primary",
        ),
        (
            "极限局灶床丢失_等PVR等总Vc_基线交换",
            focal_result,
            0.30,
            "primary",
        ),
        (
            "均匀动脉增阻_等PVR等总Vc_慢交换压力测试",
            uniform_target,
            1.0,
            "stress_test",
        ),
        (
            "极限局灶床丢失_等PVR等总Vc_慢交换压力测试",
            focal_result,
            1.0,
            "stress_test",
        ),
    ]
    controlled = []
    for name, result, tau_eq, analysis_role in controlled_specs:
        vc_ml = focal_effective_vc if result is uniform_target else 86.0
        ox = oxygen(result, vc_ml=vc_ml, tau_eq=tau_eq, vc_mode="anatomy")
        row = pack(name, result, ox)
        row["target_PVR_WU"] = target_pvr
        row["controlled_total_vc_ml"] = focal_effective_vc
        row["tau_eq_s"] = tau_eq
        row["analysis_role"] = analysis_role
        row["comparison_purpose"] = "isolate spatial distribution at equal PVR and total capillary volume"
        controlled.append(row)
    pd.DataFrame(controlled).to_csv(OUT / "goalB_controlled.csv", index=False)

    b_specs = [
        ("均匀_PVR4.19_流量容积_0.3s", uniform_target, 86.0, 0.30, "flow"),
        ("局灶_PVR4.19_流量容积_0.3s", focal_result, 86.0, 0.30, "flow"),
        ("均匀_PVR4.19_解剖容积_0.3s", uniform_target, 86.0, 0.30, "anatomy"),
        ("局灶_PVR4.19_解剖容积_0.3s", focal_result, 86.0, 0.30, "anatomy"),
        ("均匀_PVR4.19_解剖容积_1s", uniform_target, 86.0, 1.0, "anatomy"),
        ("局灶_PVR4.19_解剖容积_1s", focal_result, 86.0, 1.0, "anatomy"),
    ]
    b_rows = []
    for name, result, vc_ml, tau_eq, mode in b_specs:
        b_rows.append(pack(name, result, oxygen(
            result, vc_ml=vc_ml, tau_eq=tau_eq, vc_mode=mode
        )))
    b_rows.extend(controlled)
    pd.DataFrame(b_rows).to_csv(OUT / "goalB_joint.csv", index=False)

    pao2_rows, shunt_rows = [], []
    for value in (100, 80, 60, 40):
        ox = oxygen(healthy, lobe_pao2={"RLL": value}, vc_mode="anatomy")
        pao2_rows.append({
            "RLL_PAO2": value, "SaO2": ox["SaO2"], "PaO2": ox["PaO2"],
            "SvO2": ox["SvO2"], "S_RLL": ox["lobes"]["RLL"]["S_end"],
            "S_LUL": ox["lobes"]["LUL"]["S_end"],
            "S_end_range": ox["S_end_range"],
        })
    for value in (0.03, 0.05, 0.10, 0.15, 0.20, 0.30):
        ox = oxygen(healthy, shunt=value, vc_mode="anatomy")
        shunt_rows.append({
            "shunt": value, "SaO2": ox["SaO2"], "PaO2": ox["PaO2"],
            "SvO2": ox["SvO2"], "S_end_range": ox["S_end_range"],
        })
    pd.DataFrame(pao2_rows).to_csv(OUT / "goal3_pao2.csv", index=False)
    pd.DataFrame(shunt_rows).to_csv(OUT / "goal3_shunt.csv", index=False)

    ladder = []
    for label, target in (
        ("健康", healthy["PVR_WU"]), ("PVR3", 3.0),
        ("PVR4.19", 4.19), ("PVR6", 6.0), ("PVR10", 10.0),
    ):
        if label == "健康":
            result = healthy
        else:
            result = network(case, arterial_multiplier=arterial_multiplier_for(case, target))
        ladder.append(pack(label, result, oxygen(result)))
    x5_row = pack("动脉阻力x5", uniform_x5, oxygen(uniform_x5))
    ladder.append(x5_row)
    pd.DataFrame(ladder).to_csv(OUT / "goal4_ladder.csv", index=False)

    maximum_multiplier = arterial_multiplier_for(case, 10.0)
    curve = []
    for multiplier in np.linspace(1.0, maximum_multiplier, 61):
        result = network(case, arterial_multiplier=float(multiplier))
        ox = oxygen(result)
        curve.append({
            "arterial_multiplier": float(multiplier),
            "PVR_WU": result["PVR_WU"], "mPAP_mmHg": result["mPAP_mmHg"],
            "CO_L_min": result["CO_L_min"], "SaO2": ox["SaO2"],
            "SvO2": ox["SvO2"],
            "E_RV": ox["organs"]["右心室心肌"]["extraction"],
            "def_RV": ox["organs"]["右心室心肌"]["deficit_mL_min"],
        })
    pd.DataFrame(curve).to_csv(OUT / "goal4_curve.csv", index=False)

    goal5 = []
    scenarios5 = {
        "健康": healthy,
        "均匀动脉增阻_PVR4.19": uniform_target,
        "全肺灌注床开放50%": network(case, bed_open_fraction=0.5, scale_cap=True),
        "仅右下叶灌注床开放40%": network(
            case, focal={"RLL": 0.4}, scale_cap=True
        ),
    }
    for name, result in scenarios5.items():
        goal5.append(pack(name, result, oxygen(result, organs=ORGANS)))
    pd.DataFrame(goal5).to_csv(OUT / "goal5_organs.csv", index=False)

    goal6 = []
    for label, target in (
        ("健康", healthy["PVR_WU"]), ("PVR3", 3.0),
        ("PVR4.19", 4.19), ("PVR6", 6.0), ("PVR10", 10.0),
    ):
        result = healthy if label == "健康" else network(
            case, arterial_multiplier=arterial_multiplier_for(case, target)
        )
        ox = oxygen(
            result, organs=ORGANS, rv_power_scale=True,
            ref_power=reference_power,
        )
        row = pack(label + "|右室功率情景", result, ox)
        row["rv_demand"] = ox["organs"]["右心室心肌"]["VO2_demand_mL_min"]
        row["power_scale"] = result["mPAP_mmHg"] * result["CO_L_min"] / reference_power
        row["interpretation_scope"] = "scenario-specific; depends on preset organ shares and extraction caps"
        goal6.append(row)
    pd.DataFrame(goal6).to_csv(OUT / "goal6_rv_power.csv", index=False)

    sensitivity = []

    def add_sensitivity(factor, level, result, ox):
        sensitivity.append({
            "factor": factor, "level": level,
            "PVR_WU": result["PVR_WU"], "CO_L_min": result["CO_L_min"],
            "SaO2": ox["SaO2"], "SvO2": ox["SvO2"],
            "E_RV": ox["organs"]["右心室心肌"]["extraction"],
            "def_RV": ox["organs"]["右心室心肌"]["deficit_mL_min"],
            "total_deficit_mL_min": ox["total_deficit_mL_min"],
        })

    add_sensitivity("基线", "健康", healthy, oxygen(healthy))
    for mu in (2.5e-3, 4.5e-3):
        case_mu = lung.build_cases(mu=mu, n_root_override=case["n_root"])
        result = network(case_mu, arterial_multiplier=5.0)
        add_sensitivity("黏度", str(mu * 1000.0) + " mPa s; 动脉x5", result, oxygen(result))
    for rb in (2.7, 3.36):
        case_rb = lung.build_cases(rb=rb, n_root_override=case["n_root"])
        result = network(case_rb, arterial_multiplier=5.0)
        add_sensitivity("分支比", "Rb=" + str(rb) + "; 动脉x5", result, oxygen(result))
    for hb in (12.0, 18.0):
        add_sensitivity(
            "血红蛋白", str(hb) + " g/dL; PVR4.19",
            uniform_target, oxygen(uniform_target, hb=hb),
        )
    for value in (80.0, 150.0):
        add_sensitivity(
            "肺泡氧分压", str(value) + " mmHg; PVR4.19",
            uniform_target, oxygen(uniform_target, pao2=value),
        )
    for value in (43.0, 172.0):
        add_sensitivity(
            "毛细血管容积", str(value) + " mL; PVR4.19",
            uniform_target, oxygen(uniform_target, vc_ml=value),
        )
    for value in (0.15, 0.30, 0.60, 1.0):
        add_sensitivity(
            "交换平衡时间", str(value) + " s; PVR4.19",
            uniform_target, oxygen(uniform_target, tau_eq=value),
        )
    for value in (40.0, 60.0):
        result = network(case, arterial_multiplier=target_multiplier, reserve=value)
        add_sensitivity("右心储备", str(value) + " mmHg", result, oxygen(result))
    venous_inclusive = network(
        case, arterial_multiplier=target_multiplier,
        venous_multiplier=target_multiplier,
    )
    add_sensitivity(
        "病变边界", "相同倍数同时作用于动脉和静脉",
        venous_inclusive, oxygen(venous_inclusive),
    )
    pd.DataFrame(sensitivity).to_csv(OUT / "goal7_sensitivity.csv", index=False)

    # Joint exploratory uncertainty analysis. These are deliberately broad
    # scenario ranges, not fitted probability distributions. PVR is held at
    # 4.19 WU in every draw so the analysis concerns oxygen-model uncertainty
    # at a common haemodynamic severity.
    uncertainty_seed = 20261008
    uncertainty_sample_count = 1000
    uncertainty_ranges = [
        ("viscosity_mPa_s", 2.5, 4.5, "mPa s"),
        ("branching_ratio", 2.70, 3.36, "dimensionless"),
        ("hemoglobin_g_dL", 12.0, 18.0, "g/dL"),
        ("alveolar_PO2_mmHg", 80.0, 150.0, "mmHg"),
        ("capillary_volume_mL", 43.0, 172.0, "mL"),
        ("exchange_time_s", 0.15, 0.60, "s"),
        ("right_to_left_shunt", 0.00, 0.10, "fraction"),
        ("rv_flow_fraction", 0.03, 0.07, "fraction of CO"),
        ("rv_max_extraction", 0.55, 0.85, "fraction"),
    ]
    pd.DataFrame(
        [
            {
                "parameter": name,
                "low": low,
                "high": high,
                "unit": unit,
                "distribution": "uniform exploratory range",
                "interpretation": "assumed scenario range; not an empirical prior",
            }
            for name, low, high, unit in uncertainty_ranges
        ]
    ).to_csv(OUT / "goal8_uncertainty_parameters.csv", index=False)

    unit_design = latin_hypercube(
        uncertainty_sample_count, len(uncertainty_ranges), uncertainty_seed
    )
    uncertainty_rows = []
    input_columns = [item[0] for item in uncertainty_ranges]
    for sample_id, unit_row in enumerate(unit_design, start=1):
        sampled = {
            name: low + float(unit_row[index]) * (high - low)
            for index, (name, low, high, _) in enumerate(uncertainty_ranges)
        }
        sample_case = lung.build_cases(
            mu=sampled["viscosity_mPa_s"] / 1000.0,
            rb=sampled["branching_ratio"],
            n_root_override=case["n_root"],
        )
        sample_healthy = network(sample_case)
        sample_multiplier = arterial_multiplier_for(sample_case, target_pvr)
        sample_result = network(
            sample_case, arterial_multiplier=sample_multiplier
        )
        sample_organs = dict(ORGANS)
        rv_flow = sampled["rv_flow_fraction"]
        other_flow = (
            ORGANS["右心室心肌"][0] + ORGANS["其他"][0] - rv_flow
        )
        sample_organs["右心室心肌"] = (
            rv_flow,
            ORGANS["右心室心肌"][1],
            sampled["rv_max_extraction"],
        )
        sample_organs["其他"] = (
            other_flow,
            ORGANS["其他"][1],
            ORGANS["其他"][2],
        )
        sample_ref_power = (
            sample_healthy["mPAP_mmHg"] * sample_healthy["CO_L_min"]
        )
        sample_ox = oxygen(
            sample_result,
            vc_ml=sampled["capillary_volume_mL"],
            tau_eq=sampled["exchange_time_s"],
            hb=sampled["hemoglobin_g_dL"],
            pao2=sampled["alveolar_PO2_mmHg"],
            shunt=sampled["right_to_left_shunt"],
            organs=sample_organs,
            rv_power_scale=True,
            ref_power=sample_ref_power,
        )
        uncertainty_rows.append({
            "sample_id": sample_id,
            "seed": uncertainty_seed,
            **sampled,
            "healthy_PVR_WU": sample_healthy["PVR_WU"],
            "arterial_multiplier_for_PVR4.19": sample_multiplier,
            "PVR_WU": sample_result["PVR_WU"],
            "CO_L_min": sample_result["CO_L_min"],
            "mPAP_mmHg": sample_result["mPAP_mmHg"],
            "SaO2": sample_ox["SaO2"],
            "SvO2": sample_ox["SvO2"],
            "PaO2_mmHg": sample_ox["PaO2"],
            "PvO2_mmHg": sample_ox["PvO2"],
            "rv_extraction": sample_ox["organs"]["右心室心肌"]["extraction"],
            "rv_deficit_mL_min": sample_ox["organs"]["右心室心肌"]["deficit_mL_min"],
            "total_deficit_mL_min": sample_ox["total_deficit_mL_min"],
            "oxygen_mass_balance_error_mL_min": sample_ox["oxygen_mass_balance_error_mL_min"],
            "oxygen_converged": sample_ox["converged"],
        })
    uncertainty_frame = pd.DataFrame(uncertainty_rows)
    uncertainty_frame.to_csv(OUT / "goal8_joint_uncertainty.csv", index=False)

    uncertainty_outputs = [
        "SaO2", "SvO2", "PaO2_mmHg", "PvO2_mmHg", "rv_extraction",
        "rv_deficit_mL_min", "total_deficit_mL_min",
    ]
    uncertainty_summary = []
    for output_name in uncertainty_outputs:
        values = uncertainty_frame[output_name]
        uncertainty_summary.append({
            "output": output_name,
            "q2.5": values.quantile(0.025),
            "median": values.quantile(0.5),
            "q97.5": values.quantile(0.975),
            "minimum": values.min(),
            "maximum": values.max(),
            "sample_count": uncertainty_sample_count,
            "interpretation": "exploratory range propagation; not a statistical confidence interval",
        })
    pd.DataFrame(uncertainty_summary).to_csv(
        OUT / "goal8_uncertainty_summary.csv", index=False
    )

    rank_matrix = uncertainty_frame[input_columns + uncertainty_outputs].corr(
        method="spearman"
    )
    uncertainty_correlations = []
    for output_name in uncertainty_outputs:
        for input_name in input_columns:
            coefficient = float(rank_matrix.loc[input_name, output_name])
            uncertainty_correlations.append({
                "output": output_name,
                "input": input_name,
                "spearman_rho": coefficient,
                "absolute_rho": abs(coefficient),
                "sample_count": uncertainty_sample_count,
            })
    pd.DataFrame(uncertainty_correlations).sort_values(
        ["output", "absolute_rho"], ascending=[True, False]
    ).to_csv(OUT / "goal8_rank_correlations.csv", index=False)

    verification = verification_rows(
        case,
        [
            ("健康", healthy, oxygen(healthy)),
            ("均匀动脉增阻_PVR4.19", uniform_target, oxygen(uniform_target)),
            (
                "极限局灶床丢失_PVR4.19",
                focal_result,
                oxygen(focal_result, vc_mode="anatomy"),
            ),
        ],
    )
    pd.DataFrame(verification).to_csv(OUT / "verification_validation.csv", index=False)

    parameters = [
        ("model_type", "steady reduced resistance network", "not applicable", "definition"),
        ("blood_viscosity", case["mu"] * 1000.0, "mPa s", "assumption"),
        ("branching_ratio", case["rb"], "dimensionless", "assumption"),
        ("healthy_capillary_sheet_resistance", case["cap_all"] / lung.WU, "WU", "calibrated residual"),
        ("order12_root_count", case["n_root"], "count", "calibrated"),
        ("cardiac_output_max", 5.0, "L/min", "assumption"),
        ("right_heart_reserve_width", 50.0, "mmHg", "assumption"),
        ("baseline_total_oxygen_demand", 250.0, "mL/min", "assumption"),
        ("hemoglobin", 15.0, "g/dL", "assumption"),
        ("capillary_blood_volume", 86.0, "mL", "assumption"),
        ("exchange_equilibration_time", 0.30, "s", "assumption"),
        ("right_to_left_shunt", 0.03, "fraction", "assumption"),
        ("target_comparison_PVR", target_pvr, "WU", "analysis design"),
        ("arterial_multiplier_for_target", target_multiplier, "dimensionless", "solved"),
        ("joint_uncertainty_sample_count", uncertainty_sample_count, "count", "analysis design"),
        ("joint_uncertainty_seed", uncertainty_seed, "integer", "analysis design"),
    ]
    pd.DataFrame(
        parameters, columns=["parameter", "value", "unit", "status"]
    ).to_csv(OUT / "model_parameters.csv", index=False)

    payload = {
        "model_type": "steady_reduced_resistance_network",
        "primary_disease_operator": "distal_arterial_resistance_only",
        "oxygen_mixing": "oxygen_content",
        "systemic_oxygen_cycle": "closed_mass_balance",
        "target_PVR_WU": target_pvr,
        "arterial_multiplier_for_target": target_multiplier,
        "focal_limiting_pattern": {
            "four_lobes_open_fraction": 0.005,
            "LUL_open_fraction": focal["open_fraction"],
            "effective_capillary_volume_mL": focal_effective_vc,
        },
        "controlled_comparison": controlled,
        "ladder": ladder,
        "organ_power_scenarios": goal6,
        "sensitivity": sensitivity,
        "joint_uncertainty": {
            "sample_count": uncertainty_sample_count,
            "seed": uncertainty_seed,
            "ranges_are_empirical_priors": False,
            "summary": uncertainty_summary,
            "rank_correlations": uncertainty_correlations,
        },
        "verification": verification,
    }
    (OUT / "phase2_results.json").write_text(
        json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    if not args.no_plots:
        make_figures(goal1, controlled, pao2_rows, shunt_rows, ladder, goal6)

    print(json.dumps({
        "target_PVR_WU": target_pvr,
        "arterial_multiplier_for_target": target_multiplier,
        "focal_LUL_open_fraction": focal["open_fraction"],
        "controlled_total_vc_ml": focal_effective_vc,
        "controlled": [
            {
                "scenario": row["scenario"],
                "SaO2": row["SaO2"],
                "S_end_range": row["S_end_range"],
            }
            for row in controlled
        ],
        "goal6": [
            {
                "scenario": row["scenario"],
                "rv_demand": row["rv_demand"],
                "rv_used": row["used_右心室心肌"],
                "rv_deficit": row["def_右心室心肌"],
                "total_deficit": row["total_deficit_mL_min"],
            }
            for row in goal6
        ],
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
