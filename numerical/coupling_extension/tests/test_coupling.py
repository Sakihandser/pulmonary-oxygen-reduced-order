import unittest,sys,math,random
from pathlib import Path
from scipy.optimize import brentq
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import analyze_coupling as c

class CouplingChecks(unittest.TestCase):
 def test_independent_venous_balance(self):
  rng=random.Random(2711)
  for _ in range(100):
   q=rng.uniform(3,10);ce=rng.uniform(15,23);beta=rng.random();sstar=rng.uniform(.05,.4)
   preliminary=c.state(q,0,ce,5,sstar,beta);m=rng.uniform(.1,.8)*10*q*(1-preliminary['s'])*ce
   z=c.state(q,m,ce,5,sstar,beta);s=z['s']
   cv=brentq(lambda cv:10*q*((1-s)*ce+s*cv-cv)-m,0,ce)
   self.assertAlmostEqual(z['Cv'],cv,places=10);self.assertAlmostEqual(z['Ca'],(1-s)*ce+s*cv,places=10)
 def test_derivatives_and_feasible_delivery_sign(self):
  rng=random.Random(2712)
  for _ in range(150):
   q=rng.uniform(3,10);beta=rng.random();s=rng.uniform(.03,.5);ce=19.
   zz=c.state(q,0,ce,5,s,beta);m=.7*10*q*(1-zz['s'])*ce;z=c.state(q,m,ce,5,s,beta);h=q*1e-5
   lo=c.state(q-h,m,ce,5,s,beta);hi=c.state(q+h,m,ce,5,s,beta)
   self.assertAlmostEqual(z['dCa_dQ'],(hi['Ca']-lo['Ca'])/(2*h),delta=2e-7)
   self.assertAlmostEqual(z['dDO2_dQ'],(hi['DO2']-lo['DO2'])/(2*h),delta=2e-6)
   self.assertGreater(z['dDO2_dQ'],0)
 def test_critical_beta_reverses_content_slope(self):
  for s in [.03,.15,.2,.6]:
   crit=(1-s)/-math.log(s)
   self.assertGreater(c.state(5,100,19,5,s,crit-.01)['dCa_dQ'],0)
   self.assertAlmostEqual(c.state(5,100,19,5,s,crit)['dCa_dQ'],0,places=12)
   self.assertLess(c.state(5,100,19,5,s,crit+.01)['dCa_dQ'],0)
 def test_fixed_capacity_content_declines_delivery_increases(self):
  zz=[c.state(q,250,19,5,.2,1) for q in [5,6,8,10,20,100]]
  self.assertTrue(all(z['feasible'] for z in zz))
  for a,b in zip(zz,zz[1:]):self.assertLess(b['Ca'],a['Ca']);self.assertGreater(b['DO2'],a['DO2'])
  self.assertTrue(math.isinf(c.recovery(1,19)['Qa_ratio']))
 def test_recovery_roots_and_ratio_counterexample(self):
  ce=19.;ref=c.state(5,250,ce,5,.15,0)
  for beta in [0,.25,.5,.75,1]:
   r=c.recovery(beta,ce)
   self.assertAlmostEqual(c.state(5*r['Qd_ratio'],250,ce,5,.2,beta)['DO2'],ref['DO2'],places=8)
   if beta<1:self.assertAlmostEqual(c.state(5*r['Qa_ratio'],250,ce,5,.2,beta)['Ca'],ref['Ca'],places=9)
  self.assertAlmostEqual(c.recovery(0,ce)['Qa_ratio'],17/12,places=10)
  self.assertAlmostEqual(c.recovery(0,ce)['fixed_s_ratio_identity_residual'],0,places=10)
  self.assertGreater(abs(c.recovery(.25,ce)['fixed_s_ratio_identity_residual']),.01)
 def test_zero_demand_and_domain(self):
  z=c.state(5,0,19,5,.2,.7);self.assertEqual(z['Ca'],19);self.assertEqual(z['Cv'],19);self.assertEqual(z['dCa_dQ'],0);self.assertEqual(z['dDO2_dQ'],190)
  for args in [(0,250,19,5,.2,0),(5,-1,19,5,.2,0),(5,250,0,5,.2,0),(5,250,19,0,.2,0),(5,250,19,5,0,0),(5,250,19,5,1,0),(5,250,19,5,.2,-.1),(5,250,19,5,.2,1.1)]:
   with self.assertRaises(ValueError):c.state(*args)
 def test_conditional_beta_sensitivity(self):
  for unit,rr in c.a.unit_groups(c.a.load()).items():
   b=next(r for r in rr if r['baseline'])
   for r in rr:
    if r['baseline'] or not r['is_primary_PH']:continue
    p=c.a.beta_pred(0,b,r);s=c.v.fit_one('equilibrium_admixture',b);sat=c.v.sat(p,r['p50'],r['hill']);cp=1.34*r['Hb_g_dl']*r['hill']*sat*(1-sat)/p+.0031
    g=-(c.v.content(r['PA_mmHg'],r['Hb_g_dl'],r['p50'],r['hill'])-c.v.content(r['Pv_mmHg'],r['Hb_g_dl'],r['p50'],r['hill']))*s*(-math.log(s))*math.log(r['flow']/b['flow'])/cp
    h=1e-4;fd=(c.a.beta_pred(h,b,r)-c.a.beta_pred(-h,b,r))/(2*h)
    self.assertAlmostEqual(g,fd,delta=2e-6)
 def test_equal_flow_has_no_conditional_beta_information(self):
  for unit,rr in c.a.unit_groups(c.a.load()).items():
   b=next(r for r in rr if r['baseline'])
   for r in rr:
    t=dict(r,flow=b['flow']);ys=[c.a.beta_pred(beta,b,t) for beta in [0,.25,.5,.75,1]]
    self.assertLess(max(ys)-min(ys),1e-10)
if __name__=='__main__':unittest.main()
