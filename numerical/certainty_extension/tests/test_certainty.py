import unittest,sys,math
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from consistent import bounds,margin,content,classify,minimum_mean,matrices,solve,dual_bound
from analyze_certainty import structures,levy_discrete,levy_mean,levy_bound

class CertaintyTests(unittest.TestCase):
 def test_prior_structures_and_matched_capacity(self):
  models,mu=structures()
  for x,w in models.values():
   self.assertAlmostEqual(float(np.exp(-x)@w),.2,places=11)
   self.assertAlmostEqual(float(x@w),mu,places=11)
  self.assertAlmostEqual(mu,26.214385069/5,places=8)
 def test_levy_transform_independent_quadrature(self):
  a=-math.log(.2)
  # Substitute u=y^2 to regularize the infinite-mean density.
  for t in [.25,.5,1,2,4]:
   fn=lambda y:a/math.sqrt(math.pi)/y**2*math.exp(-a*a/(4*y*y)-t*y*y) if y else 0
   v,e=quad(fn,0,np.inf,epsabs=1e-11)
   self.assertAlmostEqual(v,math.exp(-a*math.sqrt(t)),places=9)
 def test_censored_mean_quadrature(self):
  from scipy.special import erfc
  a=-math.log(.2)
  for U in [8,16,32]:
   numerical=quad(lambda u:1-erfc(a/(2*math.sqrt(u))) if u else 1,0,U,epsabs=1e-10)[0]
   self.assertAlmostEqual(numerical,levy_mean(a,U),places=8)
 def test_alias_uniform_bound_and_finite_capacity(self):
  rng=np.random.default_rng(7131);a=-math.log(.2)
  for U,n in [(8,64),(16,128),(32,256)]:
   x,w,_=levy_discrete(a,U,n)
   self.assertTrue(np.all(x>0) and np.all(x<=U) and np.all(w>0))
   self.assertAlmostEqual(w.sum(),1,places=13)
   for r in np.r_[.5,2,rng.uniform(.5,2,300)]:
    self.assertLessEqual(abs(np.exp(-x/r)@w-math.exp(-a/math.sqrt(r))),levy_bound(a,U,n)+1e-12)
 def test_random_true_models_enclosed(self):
  rng=np.random.default_rng(9102)
  for k in range(24):
   x=rng.uniform(.01,20,5);w=rng.dirichlet(np.ones(5));r=np.array([.8,1,1.23,1.7])
   s=np.exp(-x[None,:]/r[:,None])@w;q=float(rng.uniform(.7,2.2))
   o=bounds(r,np.maximum(0,s-.005),np.minimum(1,s+.005),q,n=256)
   self.assertIn(o['status'],['bounded','unresolved'])
   val=float(np.exp(-x/q)@w)
   self.assertLessEqual(o['outer_lower'],val+1e-10);self.assertGreaterEqual(o['outer_upper'],val-1e-10)
 def test_dual_bound_against_independent_simplex_candidates(self):
  x,A,I,b=matrices([1],[.19],[.21],n=32);c=np.exp(-x[1:]/1.5)
  res=solve(c,A,b);lower,pad=dual_bound(c,A,b,res)
  self.assertLessEqual(lower,res.fun)
  rng=np.random.default_rng(6)
  for _ in range(200):
   w=rng.dirichlet(np.ones(len(c)))
   if np.max(A@w-b)<0:self.assertLessEqual(lower,float(c@w))
 def test_refinement_nests_outer_bounds(self):
  x,w=structures()[0]['A'];r=np.array([1,1.25,1.5]);s=np.exp(-x[None,:]/r[:,None])@w
  a=bounds(r,s-.001,s+.001,2,n=128);b=bounds(r,s-.001,s+.001,2,n=512)
  self.assertLessEqual(a['outer_lower'],b['outer_lower']+1e-6)
  self.assertGreaterEqual(a['outer_upper'],b['outer_upper']-1e-6)
 def test_opposite_witnesses_meet_same_data(self):
  o=bounds([1],[.199],[.201],2,n=512)
  vals=[]
  for w in o['witnesses']:
   x=np.array(w['x']);p=np.array(w['weights']);s=float(np.exp(-x)@p)
   self.assertTrue(.199-1e-8<=s<=.201+1e-8)
   vals.append(content(float(np.exp(-x/2)@p),10,18.548852,250))
  self.assertGreater(max(vals),17.6665);self.assertLess(min(vals),17.6665)
 def test_logconvex_incompatibility_and_noise_margin(self):
  e=(.27**2-.2*.3)/1.04
  self.assertAlmostEqual((.27-e)**2,(.2+e)*(.3+e),places=14)
  m=margin([1,1.25,1/.6],[.2,.27,.3],n=2048)
  self.assertGreater(m['lower'],0);self.assertGreaterEqual(m['upper'],e)
  self.assertLess(m['upper']-e,1e-6)
 def test_empty_class_not_vacuously_decided(self):
  o=bounds([1,1],[.2,.3],[.21,.31],2,n=512)
  self.assertEqual(o['status'],'incompatible')
  self.assertEqual(classify(o,10,18,250,17),'incompatible')
 def test_fick_closure_independent_root(self):
  for s in [.1,.2,.4]:
   C=content(s,8,19,250)
   cv=brentq(lambda v:10*8*((1-s)*19+s*v-v)-250,0,19)
   self.assertAlmostEqual(C,(1-s)*19+s*cv,places=11)
 def test_capacity_lower_bound_independent_inequality(self):
  # Jensen at r=1 yields E[X]>=-log(s_upper), independently of LP.
  m=minimum_mean([1],[.199],[.201],n=512)
  self.assertLessEqual(m['mean_lower'],m['mean_upper'])
  self.assertGreaterEqual(m['mean_upper'],-math.log(.201)-1e-8)
  self.assertLess(abs(m['mean_lower']+math.log(.201)),.1)
 def test_invalid_inputs(self):
  for r,lo,hi in [([0],[.2],[.3]),([1],[.3],[.2]),([1],[-.1],[.2])]:
   with self.assertRaises(ValueError):bounds(r,lo,hi,2)
  with self.assertRaises(ValueError):content(.9,1,18,250)

if __name__=='__main__':unittest.main()
