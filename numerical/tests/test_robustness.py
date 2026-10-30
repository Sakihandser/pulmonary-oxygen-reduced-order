import unittest, importlib.util, csv
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('robustness',ROOT/'robustness/code/run.py')
r=importlib.util.module_from_spec(spec);spec.loader.exec_module(r)
class RobustnessChecks(unittest.TestCase):
    def test_targets_outside_physical_range_rejected(self):
        y=r.m.end('saturation',r.m.Q0,40,17.2)
        for target in [r.m.C(40),r.m.C(100),y+.1]:
            with self.assertRaises(ValueError):r.fit(target)
    def test_independent_calibration_reaches_common_target(self):
        y=r.m.end('saturation',r.m.Q0,40,17.2)
        for offset in [-.05,0,.05]:
            for law,p in r.fit(y+offset).items():
                self.assertLess(abs(r.endpoint(law,r.m.Q0,40,p)-y-offset),1e-8)
    def test_envelopes_nest_and_robust_sets_shrink(self):
        with (ROOT/'robustness/results/calibration_envelopes.csv').open(newline='') as f:rows=list(csv.DictReader(f))
        groups={}
        for x in rows:groups.setdefault((x['Pv'],x['Q'],x['law1'],x['law2']),[]).append(x)
        for group in groups.values():
            group.sort(key=lambda x:float(x['calibration_halfwidth']))
            for a,b in zip(group,group[1:]):
                for j in [1,2]:
                    self.assertGreaterEqual(float(a['lower'+str(j)])+1e-9,float(b['lower'+str(j)]))
                    self.assertLessEqual(float(a['upper'+str(j)])-1e-9,float(b['upper'+str(j)]))
                self.assertGreaterEqual(int(a['robust_e_0.05']),int(b['robust_e_0.05']))
if __name__=='__main__':unittest.main()
