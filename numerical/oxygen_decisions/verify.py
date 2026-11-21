"""Recompute reported products and exercise the allocation bound independently."""
import math
from pathlib import Path
import json
import random
import subprocess
import sys

from plot_decisions import plot
from oxygen_analysis import HERE, ROOT, observed, synthetic, read, write, distance, content

def robust_bound():
    saved = json.loads((ROOT / "certainty_extension/results/summary.json").read_text())
    ce = saved["Ce"]
    target = saved["target_content"]
    obs = next(row for row in read(ROOT / "measurement_extension/results/measurement_precision.csv")
               if row["model"] == "B" and float(row["r"]) == 1.5)
    residual_floor = float(obs["s"]) - .04 / float(obs["endpoint_gap"])
    rows = []
    for tv in (0., .01, .02, .03, .05):
        lower = max(0., residual_floor-tv)
        upper = content(lower, 10.5, ce+.02, 245.)
        rows.append(dict(allocation_TV=tv, residual_lower=lower, query_upper=upper,
                         target_lower=target-.01, exclusion_margin=target-.01-upper))
    write("allocation_bound.csv", rows)
    assert rows[2]["exclusion_margin"] > .01

def tests():
    rng = random.Random(20261121)
    for _ in range(2000):
        k = [rng.expovariate(.4) for _ in range(4)]
        def shares():
            v = [rng.random() for _ in k]
            return [x/sum(v) for x in v]
        def residual(w):
            return sum(a*math.exp(-b/a) for a,b in zip(w,k))
        before, after = shares(), shares()
        tv = sum(abs(a-b) for a,b in zip(before,after))/2
        assert abs(residual(before)-residual(after)) <= tv+1e-14
    assert abs(distance(.3,.2,1.5,.01,.1)-.08)<1e-12
    assert distance(.2,.2,1.,0.,0.)==0
    rows = read(HERE/'observed.csv')
    primary = [r for r in rows if r['conversion']=='hill' and float(r['residual_halfwidth'])==.01
               and float(r['relative_flow_ratio_halfwidth'])==.1]
    assert len(primary)==11 and sum(float(r['residual_gap'])>0 for r in primary)==6
    witnesses=read(HERE/'witnesses.csv')
    values=[float(r['Ca'])-float(r['target']) for r in witnesses if float(r['demand'])==252.5]
    assert min(values)<0<max(values)

if __name__=='__main__':
    # Compare output rows using the full archived schema and numeric tolerances.
    files=['observed.csv','witnesses.csv','redistribution.csv','margin_grid.csv','allocation_bound.csv','pah_exercise_budget.csv','external_pressure_scores.csv','external_saturation_scores.csv','external_partition_scores.csv']
    expected={f:read(HERE/f) for f in files if (HERE/f).exists()}
    observed()
    synthetic()
    robust_bound()
    exercise = read(HERE / "pah_exercise.csv")
    baseline = float(exercise[0]["Cav_ml_dL"])
    write("pah_exercise_budget.csv", [dict(stage=r["stage"], fixed_endpoint_odds_ratio_ceiling=baseline/float(r["Cav_ml_dL"])) for r in exercise])
    subprocess.run([sys.executable, str(HERE / "external_validation.py")], check=True)
    tests()
    for filename,old in expected.items():
        new=read(HERE/filename)
        assert len(new)==len(old)
        for a,b in zip(old,new):
            assert a.keys()==b.keys()
            for k in a:
                try:
                    ok=math.isclose(float(a[k]),float(b[k]),rel_tol=1e-9,abs_tol=1e-10)
                except ValueError:
                    ok=a[k]==b[k]
                assert ok,(filename,k)
    plot()
    print('Oxygen decision checks passed:', ', '.join(files))

