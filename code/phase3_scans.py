#!/usr/bin/env python3
"""Phase-3 checks that separate spatial position from transit time, ventilation,
coronary driving pressure, exchange form, and the PVR-ladder ordering.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

MODEL = Path(__file__).resolve().parent / "build_lung_model.py"
SPEC = importlib.util.spec_from_file_location("lung", MODEL)
lung = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(lung)

OUT = Path(__file__).resolve().parent
LOBE = lung.LOBE_ORDER
LOBE_CN = {
    "RUL": "右上叶", "RML": "右中叶", "RLL": "右下叶",
    "LUL": "左上叶", "LLL": "左下叶",
}
SPARED = "LUL"


def arterial_multiplier_for(case, target_pvr, **kwargs):
    low, high = 1.0, 2.0
    while lung.solve(case, arterial_multiplier=high, **kwargs)["PVR_WU"] < target_pvr:
        high *= 2.0
        if high > 20000.0:
            raise RuntimeError("could not bracket target PVR")
    for _ in range(50):
        middle = 0.5 * (low + high)
        if lung.solve(case, arterial_multiplier=middle, **kwargs)["PVR_WU"] > target_pvr:
            high = middle
        else:
            low = middle
    return 0.5 * (low + high)


def contrast_multiplier(scale, contrast):
    return {
        lobe: scale if lobe == SPARED else scale * contrast
        for lobe in LOBE
    }


def match_contrast(case, contrast, target_pvr=4.19):
    low, high = 0.05, 2.0
    while lung.solve(
        case, lobe_arterial_multiplier=contrast_multiplier(high, contrast)
    )["PVR_WU"] < target_pvr:
        high *= 2.0
        if high > 200.0:
            raise RuntimeError("contrast pattern does not reach target PVR")
    for _ in range(46):
        middle = 0.5 * (low + high)
        result = lung.solve(
            case, lobe_arterial_multiplier=contrast_multiplier(middle, contrast)
        )
        if result["PVR_WU"] > target_pvr:
            high = middle
        else:
            low = middle
    scale = 0.5 * (low + high)
    return lung.solve(
        case, lobe_arterial_multiplier=contrast_multiplier(scale, contrast)
    )


def pack_lobe(name, result, ox, **extra):
    taus = [ox["lobes"][lobe]["tau_s"] for lobe in LOBE]
    ratios = [ox["lobes"][lobe]["tau_ratio"] for lobe in LOBE]
    row = {
        "scenario": name,
        "PVR_WU": result["PVR_WU"],
        "mPAP_mmHg": result["mPAP_mmHg"],
        "CO_L_min": result["CO_L_min"],
        "SaO2": ox["SaO2"],
        "SvO2": ox["SvO2"],
        "PaO2": ox["PaO2"],
        "S_end_range": ox["S_end_range"],
        "tau_min_s": min(taus),
        "tau_max_s": max(taus),
        "tau_ratio_min": min(ratios),
        "tau_ratio_max": max(ratios),
        "effective_vc_ml": ox["effective_vc_ml"],
        "converged": ox["converged"],
        "mass_error": ox["oxygen_mass_balance_error_mL_min"],
    }
    row.update(extra)
    for lobe in LOBE:
        row["f_" + lobe] = result["flow_fraction"][lobe]
        row["S_" + lobe] = ox["lobes"][lobe]["S_end"]
        row["tau_" + lobe] = ox["lobes"][lobe]["tau_s"]
        row["PAO2_" + lobe] = ox["lobes"][lobe]["PAO2"]
        row["Vc_" + lobe] = ox["lobes"][lobe]["vc_ml"]
    return row


def transit_rows(case):
    rows = []
    uniform = lung.solve(
        case, arterial_multiplier=arterial_multiplier_for(case, 4.19)
    )
    for vc, mode, label in (
        (86.0, "fixed", "解剖容积86"),
        (21.31, "fixed", "固定总容积21.31"),
    ):
        ox = lung.oxygen(uniform, vc_ml=vc, vc_mode=mode)
        rows.append(pack_lobe(
            "均匀_" + label, uniform, ox, contrast=1.0, volume_rule=label,
        ))
    for contrast in (4.0, 10.0, 25.0, 60.0):
        result = match_contrast(case, contrast)
        for vc, label in ((86.0, "解剖容积86"), (21.31, "固定总容积21.31")):
            ox = lung.oxygen(result, vc_ml=vc, vc_mode="fixed")
            rows.append(pack_lobe(
                f"左上叶保留_对比{contrast:g}_{label}", result, ox,
                contrast=contrast, volume_rule=label,
            ))
    # Old proportional-volume operator, reported so the null is traceable.
    focal = {
        "RUL": 0.005, "RML": 0.005, "RLL": 0.005, "LLL": 0.005, "LUL": 0.9661,
    }
    old = lung.solve(case, focal=focal, scale_cap=True)
    old_ox = lung.oxygen(old, vc_ml=86.0, vc_mode="anatomy")
    rows.append(pack_lobe(
        "旧极限局灶_容积随开放比例", old, old_ox,
        contrast=np.nan, volume_rule="开放比例缩放",
    ))
    return rows


def ventilation_map(scale_rll=1.0, total=5.0, redistribute=False):
    weights = dict(lung.LOBE_W)
    raw = {lobe: total * weights[lobe] for lobe in LOBE}
    raw["RLL"] *= scale_rll
    if redistribute:
        lost = total - sum(raw.values())
        others = [lobe for lobe in LOBE if lobe != "RLL"]
        other_weight = sum(weights[lobe] for lobe in others)
        for lobe in others:
            raw[lobe] += lost * weights[lobe] / other_weight
    return raw


def solve_vq(case, va, hpv_gain, arterial_multiplier, rematch=False):
    lobe_mult = {lobe: 1.0 for lobe in LOBE}
    result, ox = None, None
    for _ in range(14):
        result = lung.solve(
            case, arterial_multiplier=arterial_multiplier,
            lobe_arterial_multiplier=lobe_mult,
        )
        ox = lung.oxygen(result, lobe_va_L_min=va, pio2=150.0, vc_mode="fixed")
        updated = {}
        for lobe in LOBE:
            stimulus = max(0.0, (100.0 - ox["lobes"][lobe]["PAO2"]) / 60.0)
            updated[lobe] = 1.0 + hpv_gain * stimulus
        if max(abs(updated[lobe] - lobe_mult[lobe]) for lobe in LOBE) < 0.03:
            lobe_mult = updated
            break
        lobe_mult = {
            lobe: 0.5 * lobe_mult[lobe] + 0.5 * updated[lobe] for lobe in LOBE
        }
    if rematch and hpv_gain > 0.0:
        scale = arterial_multiplier_for(
            case, result["PVR_WU"] if False else 4.19,
            lobe_arterial_multiplier=lobe_mult,
        ) if arterial_multiplier > 1.0 else None
        if scale is not None:
            result = lung.solve(
                case, arterial_multiplier=scale,
                lobe_arterial_multiplier=lobe_mult,
            )
            ox = lung.oxygen(result, lobe_va_L_min=va, pio2=150.0, vc_mode="fixed")
    return result, ox, lobe_mult


def vq_rows(case):
    rows = []
    healthy_m = 1.0
    disease_m = arterial_multiplier_for(case, 4.19)
    specs = [
        ("健康_通气按容积", 1.0, False, 0.0, healthy_m, False),
        ("健康_右下叶通气x0.3_不重分配", 0.3, False, 0.0, healthy_m, False),
        ("健康_右下叶通气x0.3_重分配", 0.3, True, 0.0, healthy_m, False),
        ("健康_右下叶通气x0.3_HPV", 0.3, False, 2.0, healthy_m, False),
        ("PVR4.19_通气按容积", 1.0, False, 0.0, disease_m, False),
        ("PVR4.19_右下叶通气x0.3", 0.3, False, 0.0, disease_m, False),
        ("PVR4.19_右下叶通气x0.3_HPV", 0.3, False, 2.0, disease_m, False),
        ("PVR4.19_右下叶通气x0.3_HPV后再配平", 0.3, False, 2.0, disease_m, True),
    ]
    for name, scale, redistribute, gain, multiplier, rematch in specs:
        va = ventilation_map(scale, redistribute=redistribute)
        result, ox, lobe_mult = solve_vq(
            case, va, gain, multiplier, rematch=rematch,
        )
        rows.append(pack_lobe(
            name, result, ox,
            hpv_gain=gain,
            rll_vent_scale=scale,
            ventilation_redistributed=redistribute,
            hpv_RLL=lobe_mult["RLL"],
            VA_RLL=va["RLL"],
        ))
    return rows


def coronary_organs(co, mpap, pao=90.0, mpap_ref=13.353, q_ref=0.05, co_ref=5.065):
    drive = max(pao - mpap, 5.0)
    drive_ref = pao - mpap_ref
    q_rv = q_ref * co_ref * drive / drive_ref
    q_rv = min(q_rv, 0.30 * co)
    rest = co - q_rv
    organs = {}
    for name, (flow_share, demand_share, emax) in lung.DEFAULT_ORGANS.items():
        if name == "右心室心肌":
            organs[name] = (q_rv / co, demand_share, emax)
        else:
            organs[name] = (flow_share / 0.95 * rest / co, demand_share, emax)
    return organs


def coronary_rows(case):
    rows = []
    targets = [None, 3.0, 4.19, 6.0, 10.0]
    healthy = lung.solve(case)
    ref_power = healthy["mPAP_mmHg"] * healthy["CO_L_min"]
    for target in targets:
        result = healthy if target is None else lung.solve(
            case, arterial_multiplier=arterial_multiplier_for(case, target)
        )
        variants = {
            "固定血流份额": lung.oxygen(result),
            "冠脉驱动压": lung.oxygen(
                result, organs=coronary_organs(result["CO_L_min"], result["mPAP_mmHg"])
            ),
            "冠脉驱动压加功率需求": lung.oxygen(
                result,
                organs=coronary_organs(result["CO_L_min"], result["mPAP_mmHg"]),
                rv_power_scale=True, ref_power=ref_power,
            ),
        }
        for name, ox in variants.items():
            rv = ox["organs"]["右心室心肌"]
            rows.append({
                "scenario": name,
                "PVR_WU": result["PVR_WU"],
                "mPAP_mmHg": result["mPAP_mmHg"],
                "CO_L_min": result["CO_L_min"],
                "SaO2": ox["SaO2"],
                "SvO2": ox["SvO2"],
                "Q_RV_L_min": rv["Q_L_min"],
                "E_RV": rv["extraction"],
                "def_RV": rv["deficit_mL_min"],
                "dem_RV": rv["VO2_demand_mL_min"],
                "used_RV": rv["VO2_mL_min"],
                "DO2_RV": rv["DO2_mL_min"],
                "total_deficit": ox["total_deficit_mL_min"],
                "drive_mmHg": max(90.0 - result["mPAP_mmHg"], 5.0),
            })
    return rows


def exchange_rows(case):
    rows = []
    uniform = lung.solve(
        case, arterial_multiplier=arterial_multiplier_for(case, 4.19)
    )
    focal = match_contrast(case, 25.0)
    for label, result in (("均匀", uniform), ("对比25", focal)):
        for vc in (86.0, 21.31):
            for mode in ("saturation", "content", "diffusion"):
                ox = lung.oxygen(
                    result, vc_ml=vc, vc_mode="fixed", exchange_mode=mode,
                )
                rows.append(pack_lobe(
                    f"{label}_{mode}_Vc{vc:g}", result, ox,
                    exchange_mode=mode, vc_ml=vc,
                ))
    return rows


def ordering_rows(case):
    rng = np.random.default_rng(20261008)
    sample_count = 200
    lows = np.array([12.0, 43.0, 0.15, 80.0, 0.02, 40.0, 2.5])
    highs = np.array([18.0, 172.0, 0.60, 110.0, 0.08, 60.0, 4.5])
    design = np.empty((sample_count, 7))
    for column in range(7):
        strata = (np.arange(sample_count) + rng.random(sample_count)) / sample_count
        design[:, column] = strata[rng.permutation(sample_count)]
    samples = lows + (highs - lows) * design
    targets = (3.0, 4.19, 6.0, 10.0)
    rows = []
    for index, (hb, vc, tau, pao2, shunt, reserve, viscosity) in enumerate(samples):
        local = lung.build_cases(mu=viscosity / 1000.0)
        healthy = lung.solve(local, reserve=reserve)
        healthy_ox = lung.oxygen(
            healthy, hb=hb, vc_ml=vc, tau_eq=tau, pao2=pao2, shunt=shunt,
        )
        for target in targets:
            multiplier = arterial_multiplier_for(local, target)
            result = lung.solve(local, arterial_multiplier=multiplier, reserve=reserve)
            ox = lung.oxygen(
                result, hb=hb, vc_ml=vc, tau_eq=tau, pao2=pao2, shunt=shunt,
            )
            rows.append({
                "sample": index,
                "target_PVR": target,
                "PVR_WU": result["PVR_WU"],
                "hb": hb,
                "vc_ml": vc,
                "tau_eq": tau,
                "pao2": pao2,
                "shunt": shunt,
                "reserve": reserve,
                "viscosity": viscosity,
                "SaO2": ox["SaO2"],
                "SvO2": ox["SvO2"],
                "SaO2_healthy": healthy_ox["SaO2"],
                "SvO2_healthy": healthy_ox["SvO2"],
                "SaO2_not_lower": ox["SaO2"] >= healthy_ox["SaO2"] - 1e-6,
                "SvO2_lower": ox["SvO2"] < healthy_ox["SvO2"] - 1e-4,
                "dSaO2": ox["SaO2"] - healthy_ox["SaO2"],
                "dSvO2": ox["SvO2"] - healthy_ox["SvO2"],
            })
    return rows


def configure_matplotlib():
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    installed = {font.name for font in font_manager.fontManager.ttflist}
    for candidate in ("WenQuanYi Zen Hei", "Noto Sans CJK SC", "SimHei"):
        if candidate in installed:
            plt.rcParams["font.sans-serif"] = [candidate]
            break
    plt.rcParams["axes.unicode_minus"] = False
    return plt


def make_figures(transit, vq, coronary, exchange, ordering):
    plt = configure_matplotlib()
    transit_df = pd.DataFrame(transit)
    focus = transit_df[transit_df["volume_rule"] == "解剖容积86"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    for ax, column, title in (
        (axes[0], "tau_", "通过时间"),
        (axes[1], "S_", "末端饱和度"),
    ):
        for contrast, label in ((1.0, "均匀"), (25.0, "左上叶保留，其余叶阻力×25")):
            row = focus[np.isclose(focus["contrast"], contrast)].iloc[0]
            ax.plot(
                [LOBE_CN[lobe] for lobe in LOBE],
                [row[column + lobe] for lobe in LOBE],
                marker="o", label=label,
            )
        ax.set_title(title)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("秒")
    fig.tight_layout()
    fig.savefig(OUT / "fig_goal9_transit.png", dpi=180)
    plt.close(fig)

    vq_df = pd.DataFrame(vq)
    fig, ax = plt.subplots(figsize=(8.4, 4.2))
    labels, sao2, rll = [], [], []
    for _, row in vq_df.iterrows():
        labels.append(row["scenario"].replace("_", "\n"))
        sao2.append(row["SaO2"])
        rll.append(row["PAO2_RLL"])
    x = np.arange(len(labels))
    ax.bar(x - 0.18, sao2, 0.36, label="SaO2")
    ax2 = ax.twinx()
    ax2.plot(x, rll, color="#b85c38", marker="o", label="右下叶肺泡氧")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("SaO2")
    ax2.set_ylabel("右下叶 PAO2, mmHg")
    ax.set_title("最小通气血流层：低通气与缺氧收缩")
    fig.tight_layout()
    fig.savefig(OUT / "fig_goal10_vq.png", dpi=180)
    plt.close(fig)

    cor = pd.DataFrame(coronary)
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0))
    for name in ("固定血流份额", "冠脉驱动压", "冠脉驱动压加功率需求"):
        part = cor[cor["scenario"] == name]
        axes[0].plot(part["PVR_WU"], part["Q_RV_L_min"], marker="o", label=name)
        axes[1].plot(part["PVR_WU"], part["def_RV"], marker="o", label=name)
    axes[0].set_title("右室血流")
    axes[1].set_title("右室氧耗缺口")
    axes[0].set_ylabel("L/min")
    axes[1].set_ylabel("mL/min")
    for ax in axes:
        ax.set_xlabel("PVR, Wood 单位")
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "fig_goal11_coronary.png", dpi=180)
    plt.close(fig)

    order = pd.DataFrame(ordering)
    summary = order.groupby("target_PVR").agg(
        sao2_not_lower=("SaO2_not_lower", "mean"),
        svo2_lower=("SvO2_lower", "mean"),
        dsao2=("dSaO2", "median"),
        dsvo2=("dSvO2", "median"),
    ).reset_index()
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(summary["target_PVR"], summary["sao2_not_lower"], marker="o", label="SaO2 不降的样本比例")
    ax.plot(summary["target_PVR"], summary["svo2_lower"], marker="o", label="SvO2 下降的样本比例")
    ax.set_ylim(0.0, 1.05)
    ax.set_xlabel("目标 PVR, Wood 单位")
    ax.set_ylabel("比例")
    ax.set_title("200 组探索性样本中的主排序")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "fig_goal12_ordering.png", dpi=180)
    plt.close(fig)
    return summary


def main():
    case = lung.build_cases()
    transit = transit_rows(case)
    vq = vq_rows(case)
    coronary = coronary_rows(case)
    exchange = exchange_rows(case)
    ordering = ordering_rows(case)
    pd.DataFrame(transit).to_csv(OUT / "goal9_transit.csv", index=False)
    pd.DataFrame(vq).to_csv(OUT / "goal10_vq.csv", index=False)
    pd.DataFrame(coronary).to_csv(OUT / "goal11_coronary.csv", index=False)
    pd.DataFrame(exchange).to_csv(OUT / "goal12_exchange.csv", index=False)
    pd.DataFrame(ordering).to_csv(OUT / "goal12_ordering.csv", index=False)
    summary = make_figures(transit, vq, coronary, exchange, ordering)
    summary.to_csv(OUT / "goal12_ordering_summary.csv", index=False)
    print(json.dumps({
        "transit": transit,
        "vq": vq,
        "coronary": coronary,
        "exchange": exchange,
        "ordering_summary": summary.to_dict(orient="records"),
    }, ensure_ascii=False, indent=2, default=float))


if __name__ == "__main__":
    main()
