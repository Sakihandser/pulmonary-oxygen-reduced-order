import unittest,sys,math
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import analyze_applicability as x

class ApplicabilityTests(unittest.TestCase):
 def test_heterogeneous_exponents_and_both_flow_directions(self):
  rng=np.random.default_rng(606)
  for _ in range(600):
   z=rng.uniform(.0001,.9999,8);w=rng.dirichlet(np.ones(8));beta=rng.uniform(0,1,8);r=rng.uniform(.4,2.5)
   s0=w@z;s1=w@(z**(r**(-beta)));lo,hi=x.envelope(s0,r)
   self.assertLessEqual(lo-1e-12,s1);self.assertLessEqual(s1,hi+1e-12)
 def test_unknown_beta_interval_attainment(self):
  for s0 in [.01,.2,.9]:
   for r in [.6,1.2,2.]:
    lo,hi=x.envelope(s0,r)
    for s1 in np.linspace(lo,hi,20):
     beta=math.log(math.log(s1)/math.log(s0))/(-math.log(r))
     self.assertGreaterEqual(beta,-1e-12);self.assertLessEqual(beta,1+1e-12)
 def test_compatibility_at_unit_flow_and_radius(self):
  self.assertFalse(x.interval_compatible(.2,.2,.21,.21,1,1))
  self.assertAlmostEqual(x.residual_radius(.2,.21,1),.005,10)
  self.assertEqual(x.residual_radius(.2,.2,1),0)
 def test_interval_union_against_dense_direct_grid(self):
  svals=np.linspace(.13,.21,41);rvals=np.linspace(.85,1.4,51)
  betas=np.linspace(0,1,101)
  vals=[s**(r**(-b)) for s in svals for r in rvals for b in betas]
  self.assertTrue(x.interval_compatible(.13,.21,min(vals),min(vals),.85,1.4))
  self.assertFalse(x.interval_compatible(.13,.21,max(vals)+.001,max(vals)+.002,.85,1.4))
 def test_moment_lp_is_inside_union(self):
  z=np.linspace(0,1,2001);s0=.2
  for r in [.7,1.4,2.]:
   for beta in [.2,.6,1.]:
    f=z**(r**(-beta));bounds=x.envelope(s0,r)
    for sign in [-1,1]:
     lp=linprog(sign*f,A_eq=np.array([np.ones(len(z)),z]),b_eq=[1,s0],bounds=(0,None),method='highs')
     self.assertTrue(lp.success);val=lp.fun*sign
     self.assertLessEqual(bounds[0]-1e-9,val);self.assertLessEqual(val,bounds[1]+1e-9)
 def test_pressure_margin_crosses_boundary(self):
  for _,b,t in x.pairs(x.a.load()):
   r=t['flow']/b['flow'];e=x.pa_radius(b,t,r)
   self.assertTrue(x.interval_compatible(*x.pa_interval(b,e+1e-6),*x.pa_interval(t,e+1e-6),r,r))
   if e>1e-5:self.assertFalse(x.interval_compatible(*x.pa_interval(b,e-1e-5),*x.pa_interval(t,e-1e-5),r,r))
 def test_residual_margin_minimal(self):
  for s0,s1,r in [(.2,.3,1),(.2,.15,1.5),(.1,.5,.8),(.3,.02,.8)]:
   e=x.residual_radius(s0,s1,r)
   f=lambda d:x.interval_compatible(max(0,s0-d),min(1,s0+d),max(0,s1-d),min(1,s1+d),r,r)
   self.assertTrue(f(e+1e-6));self.assertFalse(f(e-1e-6))
 def test_common_content_scale_invariance(self):
  r=next(x.pairs(x.a.load()))[1];ce,ca,cv,s=x.contents(r)
  for scale in [.5,2,10]:self.assertAlmostEqual((scale*ce-scale*ca)/(scale*ce-scale*cv),s,12)
 def test_prediction_separation_matches_existing_predictions(self):
  import csv
  with (x.ROOT/'review_extension/results/matched_predictions.csv').open(encoding='utf-8') as f:old=list(csv.DictReader(f))
  for unit,b,t in x.pairs(x.a.load()):
   ce,_,cv,_=x.contents(t);s0=x.contents(b)[3];r=t['flow']/b['flow']
   for beta in [0,1]:
    pred=x.v.invcontent(ce-(ce-cv)*s0**(r**(-beta)),t['Hb_g_dl'],t['p50'],t['hill'])
    ref=next(a for a in old if a['unit']==unit and a['condition']==t['condition'] and a['model']==f'content_beta{beta}')
    self.assertAlmostEqual(pred,float(ref['predicted_Pa']),8)
 def test_invalid_inputs(self):
  for s,r in [(-.1,1),(.2,0),(1.1,2),(.2,float('inf'))]:
   with self.assertRaises(ValueError):x.envelope(s,r)
  with self.assertRaises(ValueError):x.interval_compatible(.2,.1,.2,.3,1,2)

if __name__=='__main__':unittest.main()
