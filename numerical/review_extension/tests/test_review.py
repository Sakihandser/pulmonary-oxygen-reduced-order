from pathlib import Path
import sys,unittest,math
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import analyze_review as r

class StructuralControls(unittest.TestCase):
 def setUp(self):
  self.rows=r.a.load()
 def test_fixed_capacity_recovers_original_predictions(self):
  for _,rows in r.a.unit_groups(self.rows).items():
   b=next(x for x in rows if x['baseline'])
   for t in rows:
    for law in r.LAWS:
     self.assertAlmostEqual(r.matched(law,b,t,1),r.v.predict(law,r.v.fit_one(law,b),r.a.allowed_inputs(t)),places=8)
 def test_matched_scaling_is_explicit_capacity_rescaling(self):
  for _,rows in r.a.unit_groups(self.rows).items():
   b=next(x for x in rows if x['baseline'])
   for t in rows:
    for law in r.LAWS:
     expected=r.v.predict(law,r.v.fit_one(law,b)*t['flow']/b['flow'],r.a.allowed_inputs(t))
     self.assertAlmostEqual(r.matched(law,b,t,0),expected,places=8)
 def test_content_control_equals_mass_mixing(self):
  for _,rows in r.a.unit_groups(self.rows).items():
   b=next(x for x in rows if x['baseline']);s=r.v.fit_one('equilibrium_admixture',b)
   for t in rows:
    self.assertAlmostEqual(r.matched('content',b,t,0),r.v.predict('equilibrium_admixture',s,r.a.allowed_inputs(t)),places=8)
 def test_beta_zero_excludes_current_Q_and_outcome(self):
  for _,rows in r.a.unit_groups(self.rows).items():
   b=next(x for x in rows if x['baseline'])
   for t in rows:
    corrupt=dict(t,Pa_mmHg=-999,flow=math.nan)
    for law in r.LAWS:self.assertAlmostEqual(r.matched(law,b,t,0),r.matched(law,b,corrupt,0),places=9)
 def test_saturation_residual_direct_formula(self):
  for _,rows in r.a.unit_groups(self.rows).items():
   b=next(x for x in rows if x['baseline'])
   sat=lambda p:r.v.sat(p,b['p50'],b['hill'])
   frac=(sat(b['PA_mmHg'])-sat(b['Pa_mmHg']))/(sat(b['PA_mmHg'])-sat(b['Pv_mmHg']))
   for t in rows:
    y=r.v.invsat((1-frac)*sat(t['PA_mmHg'])+frac*sat(t['Pv_mmHg']),b['p50'],b['hill'])
    self.assertAlmostEqual(r.matched('saturation',b,t,0),y,places=8)
 def test_compensation_roots_random_feasible_references(self):
  rng=np.random.default_rng(261114)
  for _ in range(100):
   s0,s1=rng.uniform(.05,.35,2);ce=20.;q0=5.;M=250.
   qa,qd,qf=r.m.flow_thresholds(q0,M,s0,ce,M,s1,ce)
   baseline=r.m.closed(q0,M,s0,ce)
   actual=brentq(lambda q:r.m.closed(q,M,s1,ce)['Ca']-baseline['Ca'],.01,100.)
   delivery=brentq(lambda q:r.m.closed(q,M,s1,ce)['DO2']-baseline['DO2'],.01,100.)
   self.assertAlmostEqual(actual,qa,places=8);self.assertAlmostEqual(delivery,qd,places=8)
   self.assertAlmostEqual(r.m.closed(qf,M,s1,ce)['Cv'],0,places=8)
 def test_compensation_degenerate_targets(self):
  ce=20.;ca=r.m.closed(5,250,.15,ce)['Ca']
  self.assertTrue(math.isinf(r.m.flow_thresholds(5,250,.15,ce,250,.2,ca)[0]))
  self.assertEqual(r.m.flow_thresholds(5,250,.15,ce,0,.2,ca)[0],0)
  self.assertEqual(r.m.flow_thresholds(5,250,.15,ce,250,0,ce)[0],0)
 def test_dimensionless_rule_with_varying_demand(self):
  rng=np.random.default_rng(14011)
  for _ in range(100):
   s0,s1=rng.uniform(.05,.35,2);M0,M1=rng.uniform(100,400,2);ce=20.;q0=5.
   eps=r.m.odds(s0)*M0/(10*q0*ce);qa,qd,_=r.m.flow_thresholds(q0,M0,s0,ce,M1,s1,ce)
   self.assertAlmostEqual(qd/q0-1,eps*(qa/q0-1),places=12)
 def test_family_weighting_is_not_row_pooling(self):
  pp=[dict(model='m',family='a',error=1.,feasible=1),*[dict(model='m',family='b',error=3.,feasible=1) for _ in range(9)]]
  self.assertAlmostEqual(r.a.metric('x',pp,'m')['equal_family_RMSE'],math.sqrt(5))
if __name__=='__main__':unittest.main()
