"""Separate blood-curve, cross-population and within-person checks (Supplement V4)."""
from decimal import Decimal
from itertools import product

from scipy.optimize import brentq

from external_validation import content, predict_pressure, read, save
from validate import BINDING, DISSOLVED


def severinghaus_content(pressure, hb):
    cubic = pressure ** 3 + 150 * pressure
    return BINDING * hb * cubic / (cubic + 23400.) + DISSOLVED * pressure


def predict(rest_pa, rest_pv, rest_co2, pv, co2, curve0, curve1, rq, ratio=None):
    """Fit the resting residual and mix at measured follow-up Pv and PaCO2."""
    def alveolar(pco2):
        return .21 * (760. - 47.) - pco2 * (.21 + .79 / rq)

    end0, end1 = alveolar(rest_co2), alveolar(co2)
    if not rest_pv < rest_pa < end0 or pv >= end1:
        return None
    s0 = (curve0(end0) - curve0(rest_pa)) / (curve0(end0) - curve0(rest_pv))
    def pressure(residual):
        target = curve1(end1) - residual * (curve1(end1) - curve1(pv))
        return brentq(lambda p: curve1(p) - target, pv, end1)
    fixed = pressure(s0)
    other = pressure(s0 ** (1 / ratio)) if ratio is not None else None
    return fixed, other


def gas_inputs(row):
    return [float(row[key]) for key in ("baseline_Pa", "baseline_Pv", "baseline_PaCO2",
                                         "query_Pv", "query_PaCO2")]


def checks():
    scores = []
    for row in read("external_pressure.csv") + read("transfer_gases.csv"):
        rest, pv0, co20, pv1, co21 = gas_inputs(row)
        study = row.get("study", "AbdelKafi1998")
        ratio = float(row["query_CI"]) / float(row["baseline_CI"]) if "query_CI" in row else None
        observed = float(row["observed_Pa"])
        for name, curve in (("Hill", lambda p: content(p, 15., 26.8, 2.7)),
                            ("Severinghaus", lambda p: severinghaus_content(p, 15.))):
            fixed, other = predict(rest, pv0, co20, pv1, co21, curve, curve, .8, ratio)
            if ratio is not None and name == "Hill":
                reference = predict_pressure(rest, pv0, co20, ratio, pv1, co21, 26.8, 2.7, .8, 15., 0.)[0]
                assert max(abs(a-b) for a,b in zip((fixed, other), reference)) < 1e-8
            venous = rest + pv1 - pv0
            scores.append(dict(study=study, condition=row["condition"], curve=name, rest=rest,
                               observed=observed, fixed=fixed, fixed_capacity=other, venous=venous,
                               fixed_error=fixed-observed, venous_error=venous-observed,
                               unchanged_error=rest-observed))
    return scores


def sensitivity():
    results = []
    for row in read("transfer_gases.csv"):
        args = gas_inputs(row)
        observed = float(row["observed_Pa"])
        for p50, n, rq, hb, shift in product((24., 26.8, 30.), (2.5, 2.7, 3.), (.75, .8, .9, 1.),
                                              (12., 15., 18.), (-2., 0., 2.)):
            answer = predict(*args, lambda p: content(p, hb, p50, n),
                             lambda p: content(p, hb, p50+shift, n), rq)
            results.append(dict(study=row["study"], condition=row["condition"], p50=p50, hill=n,
                                rq=rq, hb=hb, shift=shift, feasible=int(answer is not None),
                                predicted=answer[0] if answer else None,
                                error=answer[0]-observed if answer else None))
    return results


def rounding_bounds(scores):
    """Exact separable range for fixed-minus-venous MAE at rounded outcome means."""
    results = []
    for curve in ("Hill", "Severinghaus"):
        rows = [r for r in scores if r["study"] == "AbdelKafi1998" and r["curve"] == curve]
        nominals, minima, maxima = [], [], []
        for row in rows:
            fixed, venous, target = (row[k] for k in ("fixed", "venous", "observed"))
            lo, hi = target - .5, target + .5
            knots = [lo, hi] + [p for p in (fixed, venous) if lo <= p <= hi]
            values = [abs(fixed-y)-abs(venous-y) for y in knots]
            nominals.append(abs(fixed-target)-abs(venous-target))
            minima.append(min(values))
            maxima.append(max(values))
        results.append(dict(curve=curve, settings=len(rows), nominal=sum(nominals)/len(rows),
                            lower=sum(minima)/len(rows), upper=sum(maxima)/len(rows)))
    return results


def saturation_changes():
    rows = read("exercise_saturation_pairs.csv")
    first = {}
    for row in rows:
        case = int(row["case"])
        if row["exercise_Sa_percent"]:
            key = (int(row["visit"]), int(row["work_kpm_min"]))
            first[case] = min(first.get(case, key), key)
    result = []
    for row in rows:
        case = int(row["case"])
        paired = bool(row["exercise_Sa_percent"])
        primary = paired and (int(row["visit"]), int(row["work_kpm_min"])) == first[case]
        delta = Decimal(row["exercise_Sa_percent"])-Decimal(row["rest_Sa_percent"]) if paired else None
        result.append(dict(**row, first_visit_lowest_workload=int(primary), delta_pp=delta))
    assert sum(r["first_visit_lowest_workload"] for r in result) == 8
    return result


if __name__ == "__main__":
    scores = checks()
    save("transfer_pressure_scores.csv", scores)
    save("transfer_rounding_bounds.csv", rounding_bounds(scores))
    save("transfer_hill_grid.csv", sensitivity())
    save("exercise_saturation_changes.csv", saturation_changes())
