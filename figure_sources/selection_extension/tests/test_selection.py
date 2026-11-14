from pathlib import Path
import unittest,sys,math
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import analyze as a
import mechanisms as m
v=a.v
class SelectionTests(unittest.TestCase):
 def test_source_structure(self):
  rows=a.read_csv(a.R/'data/additional_observations.csv')
  self.assertEqual(len(rows),15);self.assertEqual(sum(r['Pv_status']=='reported' for r in rows),11)
  self.assertEqual({r['cohort_group'] for r in rows if r['study']!='Melot1983'},{'Zurich2021_2023'})
  self.assertEqual(sum(r['baseline']=='1' for r in rows if r['study']=='Melot1983'),3)
 def test_venous_input_conversion(self):
  rr=a.load();r=next(r for r in rr if r['study']=='Lichtblau2025')
  self.assertAlmostEqual(v.sat(r['Pv_mmHg']),.69,12);self.assertFalse(r['is_primary_PH'])
 def test_test_outcome_is_not_predictor(self):
  rr=a.load();base=next(r for r in rr if r['study']=='Haefliger2025' and r['baseline']);target=next(r for r in rr if r['study']=='Haefliger2025' and not r['baseline']);changed=dict(target,Pa_mmHg=9999,Sa_percent=-99)
  for law in a.MODELS:
   par=a.base_fit(law,base);self.assertEqual(a.pred(law,par,target),a.pred(law,par,changed))
  self.assertEqual(a.beta_pred(.3,base,target),a.beta_pred(.3,base,changed))
 def test_gas_equation_no_outcome(self):
  for r in a.load():
   expected=r['FiO2']*(r['PB_mmHg']-r['PH2O_mmHg'])-r['PaCO2_mmHg']*(r['FiO2']+(1-r['FiO2'])/.8)
   self.assertEqual(r['PA_mmHg'],expected)
 def test_beta_endpoints_exact(self):
  for rr in a.unit_groups(a.load()).values():
   b=next(r for r in rr if r['baseline'])
   for r in rr:
    self.assertAlmostEqual(a.beta_pred(0,b,r),a.pred('equilibrium_admixture',v.fit_one('equilibrium_admixture',b),r),9)
    self.assertAlmostEqual(a.beta_pred(1,b,r),a.pred('content',v.fit_one('content',b),r),9)
 def test_scale_unit_invariance(self):
  rr=a.unit_groups(a.load())['M01_visit1'];b=rr[0];r=rr[1]
  for beta in [-.7,0,.3,1,2]:self.assertAlmostEqual(a.beta_pred(beta,b,r),a.beta_pred(beta,dict(b,flow=b['flow']*1.8),dict(r,flow=r['flow']*1.8)),12)
 def test_equal_family_weight(self):
  rows=[dict(model='x',family='A',error=2.,feasible=1)]*7+[dict(model='x',family='B',error=4.,feasible=1)]
  self.assertAlmostEqual(a.metric('test',rows,'x')['equal_family_RMSE'],math.sqrt(10),12)
 def test_flow_boundary_against_numerical_root(self):
  for s0,s1 in [(.1,.2),(.15,.25),(.2,.1)]:
   q0,m0,ce0,m1,ce1=5.,250.,19.,300.,18.5
   qa,qd,qf=m.flow_thresholds(q0,m0,s0,ce0,m1,s1,ce1);ref=m.closed(q0,m0,s0,ce0)
   self.assertAlmostEqual(qa,brentq(lambda q:m.closed(q,m1,s1,ce1)['Ca']-ref['Ca'],.01,1000),9)
   self.assertAlmostEqual(qd,brentq(lambda q:m.closed(q,m1,s1,ce1)['DO2']-ref['DO2'],.01,1000),9)
   self.assertAlmostEqual(m.closed(qf,m1,s1,ce1)['Cv'],0,10)
 def test_unreachable_arterial_target(self):
  self.assertTrue(math.isinf(m.flow_thresholds(5,250,.15,19,250,.2,10)[0]))
 def test_mass_balance(self):
  for q,s in [(3,.1),(5,.25),(8,.4)]:
   st=m.closed(q,250,s,19);self.assertAlmostEqual(10*q*(st['Ca']-st['Cv']),250,10);self.assertAlmostEqual(st['Ca'],(1-s)*19+s*st['Cv'],12)
 def test_compensation_nontrivial_discordance(self):
  ref=m.closed(5,250,.15,19);changed=m.closed(6,250,.25,19)
  self.assertLess(changed['Ca'],ref['Ca']);self.assertGreater(changed['DO2'],ref['DO2'])
 def test_shapley_adds_exactly(self):
  for rr in a.unit_groups(a.load()).values():
   b=next(r for r in rr if r['baseline']);s=v.fit_one('equilibrium_admixture',b)
   for r in rr:
    self.assertAlmostEqual(sum(m.shapley(b,r,s).values()),a.pred('equilibrium_admixture',s,r)-b['Pa_mmHg'],8)
 def test_positive_domain(self):
  with self.assertRaises(ValueError):m.closed(-1,250,.2,19)
  with self.assertRaises(ValueError):m.odds(1.)
 def test_group_train_test_disjoint(self):
  rr=[r for r in a.load() if r['is_primary_PH']]
  for f in a.FAMILIES:
   train={r['cohort_group'] for r in rr if r['cohort_group']!=f};test={r['cohort_group'] for r in rr if r['cohort_group']==f}
   self.assertFalse(train&test);self.assertEqual(len(train),2)
if __name__=='__main__':unittest.main()
