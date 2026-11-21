"""Retrospective external challenge of frozen oxygen laws on published summaries."""
import csv
from itertools import product
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "experimental_validation/code"))
from validate import content, invcontent


def read(name):
    with (HERE / name).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def save(name, rows):
    with (HERE / name).open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


def predict_pressure(rest_pa, rest_pv, rest_co2, ratio, query_pv, query_co2, p50, hill, rq, hb, shift):
    def endpoint(co2, affinity):
        alveolar_pressure = .21 * (760 - 47) - co2 * (.21 + .79 / rq)
        return content(alveolar_pressure, hb, affinity, hill)

    query_p50 = p50 + shift
    end0, end1 = endpoint(rest_co2, p50), endpoint(query_co2, query_p50)
    residual = (end0 - content(rest_pa, hb, p50, hill)) / (end0 - content(rest_pv, hb, p50, hill))
    gap = end1 - content(query_pv, hb, query_p50, hill)
    predictions = tuple(invcontent(end1 - gap * s, hb, query_p50, hill) for s in (residual, residual ** (1 / ratio)))
    return predictions, residual, end1, gap


def score(rest, observed, fixed, flow_scaled):
    low, high = sorted((fixed, flow_scaled))
    return dict(observed=observed, frozen_residual=fixed, flow_scaled=flow_scaled, lower=low, upper=high,
                outside_distance=max(low - observed, observed - high, 0.0),
                fixed_error=fixed-observed, unchanged_error=rest-observed,
                fixed_change=fixed-rest, observed_change=observed-rest)


def main():
    pressure = []
    for row in read("external_pressure.csv"):
        rest, pv0, co20, q0, q1, pv1, co21, observed = (
            float(row[k]) for k in ("baseline_Pa","baseline_Pv","baseline_PaCO2","baseline_CI",
                                    "query_CI","query_Pv","query_PaCO2","observed_Pa"))
        for p50, hill, rq, hb, shift in product((24.,26.8,30.), (2.5,2.7,3.), (.75,.8,.9,1.), (12.,15.,18.), (-2.,0.,2.)):
            (fixed, scaled), s0, end1, gap = predict_pressure(rest, pv0, co20, q1/q0, pv1, co21, p50, hill, rq, hb, shift)
            s1 = (end1 - content(observed, hb, p50+shift, hill)) / gap
            ratio = q1/q0
            s_low, s_high = max(0., s0-.01), min(1., s0+.01)
            lower = min(s_low, s_low**(1/(.9*ratio)))
            upper = max(s_high, s_high**(1/(1.1*ratio)))
            residual_gap = max(lower-s1-.01, s1-.01-upper, 0.)
            pressure.append(dict(condition=row["condition"], p50=p50, hill=hill, rq=rq, hb=hb, query_p50_shift=shift,
                                 joint_residual_gap=residual_gap, venous_delta_error=rest+pv1-pv0-observed,
                                 **score(rest, observed, fixed, scaled)))
    save("external_pressure_scores.csv", pressure)

    saturation = []
    for row in read("external_saturation.csv"):
        rest, ven0, ven1, observed, q0, q1 = (float(row[k]) for k in
            ("baseline_Sp","baseline_Sv","query_Sv","observed_Sp","baseline_flow","query_flow"))
        for endpoint in (98.,99.,100.):
            residual = (endpoint - rest) / (endpoint - ven0)
            fixed = endpoint - (endpoint - ven1)*residual
            scaled = endpoint - (endpoint - ven1)*residual**(q0/q1)
            saturation.append(dict(study=row["study"], endcap_percent=endpoint,
                                   **score(rest, observed, fixed, scaled)))
    save("external_saturation_scores.csv", saturation)

    partition = read("external_flow_partition.csv")
    base = partition[0]
    results = []
    for row in partition:
        q = float(row["CO_L_min"])
        predicted = q * float(base["shunt_L_min"]) / float(base["CO_L_min"])
        measured = float(row["shunt_L_min"])
        results.append(dict(workload_percent=row["workload_percent"], frozen_share_shunt_L_min=predicted,
                            reported_shunt_L_min=measured, difference_L_min=predicted-measured,
                            reported_share_change_pp=float(row["shunt_percent_reported"])-float(base["shunt_percent_reported"]),
                            saturation_change_pp=float(row["saturation_percent"])-float(base["saturation_percent"])))
    save("external_partition_scores.csv", results)


if __name__ == "__main__":
    main()

