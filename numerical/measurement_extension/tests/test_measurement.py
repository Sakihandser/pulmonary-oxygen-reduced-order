import unittest,sys,itertools
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from analyze_measurement import residual,box_interval
class MeasurementTests(unittest.TestCase):
 def test_corners(self):
  ce,ca,cv=19.,18.,14.;err=(.03,.02,.01)
  vals=[residual(ce+a*err[0],ca+b*err[1],cv+c*err[2]) for a,b,c in itertools.product([-1,1],repeat=3)]
  np.testing.assert_allclose(box_interval(ce,ca,cv,*err),[min(vals),max(vals)],atol=1e-14)
 def test_equal_box_identity(self):
  for e in [.001,.01,.1]:np.testing.assert_allclose(box_interval(19,18,14,e,e,e),[.2-2*e/5,.2+2*e/5])
 def test_common_affine(self):
  for a,b in itertools.product([.7,1.3,2.],[-2,0,3]):self.assertAlmostEqual(residual(a*19+b,a*18+b,a*14+b),.2)
 def test_domain(self):
  for e in [-.1,.6]:
   with self.assertRaises(ValueError):box_interval(19,18,14,e,e,e)
 def test_random_interior(self):
  rng=np.random.default_rng(20261118);lo,hi=box_interval(19,18,14,.03,.02,.01)
  for row in rng.uniform(-1,1,(2000,3)):
   s=residual(*(np.array([19,18,14])+row*np.array([.03,.02,.01])))
   self.assertTrue(lo<=s<=hi)
 def test_frontier(self):
  import json
  h=Path(__file__).resolve().parents[1];v=json.loads((h/'results/summary.json').read_text())
  margin=v['content_margin'];u=v['relative_consumption_margin'];self.assertGreater(margin,0);self.assertLess(u,.02)
if __name__=='__main__':unittest.main()
