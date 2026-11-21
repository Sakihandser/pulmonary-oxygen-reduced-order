"""Pairwise compatibility, query-input margins and finite counterexamples."""
from pathlib import Path
import csv
import json
import math
import sys

import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "certainty_extension/code"))
from analyze_certainty import structures
from consistent import content

def read(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))

def write(name, rows):
    with (HERE / name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)

def distance(s0, s1, ratio, epsilon, flow_error):
    # Extremes over the baseline residual and the flow ratio occur at corners.
    lo = max(0., s0 - epsilon)
    hi = min(1., s0 + epsilon)
    lower = min(lo, lo ** (1 / (ratio * (1 - flow_error))))
    upper = max(hi, hi ** (1 / (ratio * (1 + flow_error))))
    return max(0., s1 - epsilon - upper, lower - s1 - epsilon)

def observed():
    rows = []
    for filename in ("pair_compatibility.csv", "reported_saturation_audit.csv"):
        for row in read(ROOT / "applicability_extension/results" / filename):
            if row["primary"] != "1" or float(row["FiO2"]) > .30:
                continue
            for epsilon in (0., .005, .01, .02):
                for flow_error in (0., .05, .10):
                    gap = distance(float(row["s0"]), float(row["s1"]), float(row["r"]), epsilon, flow_error)
                    rows.append(dict(unit=row["unit"], condition=row["condition"], method=row["flow_method"],
                                     conversion=row["convention"], residual_halfwidth=epsilon,
                                     relative_flow_ratio_halfwidth=flow_error, residual_gap=gap))
    write("observed.csv", rows)

def synthetic():
    saved = json.loads((ROOT / "certainty_extension/results/summary.json").read_text())
    ce, target = saved["Ce"], saved["target_content"]
    models, _ = structures()
    x, w = models["A"]
    capacities = 5 * w * x
    def redistributed(delta):
        shares = w + np.array([-delta, delta])
        residual = shares @ np.exp(-capacities / (10 * shares))
        return content(residual, 10, ce, 250)
    crossing = brentq(lambda delta: redistributed(delta) - target, 0, .02)
    write("redistribution.csv", [dict(shift=d, Ca=redistributed(d), target=target)
                                 for d in (0, .01, crossing, .02)])
    obs = read(ROOT / "measurement_extension/results/measurement_precision.csv")
    rows = []
    original = read(ROOT / "measurement_extension/results/measurement_witnesses.csv")
    for direction in ("lower", "upper"):
        nodes = [row for row in original if row["model"] == "A" and
                 float(row["content_halfwidth"]) == .003 and row["direction"] == direction]
        xx = np.r_[[float(n["exposure"]) for n in nodes], x]
        ww = np.r_[.999 * np.array([float(n["weight"]) for n in nodes]), .001 * w]
        ww /= ww.sum()
        for row in (row for row in obs if row["model"] == "A"):
            s = ww @ np.exp(-xx / float(row["r"]))
            # Opposite common shifts of Ce/Cv versus Ca attain every residual
            # in the equal-width box with the denominator held fixed.
            shift = (s - float(row["s"])) * float(row["endpoint_gap"]) / 2
            assert abs(shift) < .003
            predicted = (float(row["Ce"]) + shift - float(row["Ca"]) + shift) / float(row["endpoint_gap"])
            assert abs(predicted - s) < 1e-12
        sq = ww @ np.exp(-xx / 2)
        for demand in (250., 252.5):
            rows.append(dict(direction=direction, demand=demand, query_residual=sq,
                             Ca=content(sq, 10, ce, demand), target=target))
    write("witnesses.csv", rows)
    margins = []
    for row in read(ROOT / "measurement_extension/results/propagated_decisions.csv"):
        for u in (0., .0025, .005, .01):
            for shift in (0., .001, .005):
                s = float(row["s_upper"])
                lower = content(s, 10, ce, 250 * (1 + u)) - target - shift
                margins.append(dict(model=row["model"], content_halfwidth=row["content_halfwidth"],
                                    demand_up_fraction=u, target_up=shift, lower_margin=lower))
    write("margin_grid.csv", margins)

if __name__ == "__main__":
    observed()
    synthetic()

