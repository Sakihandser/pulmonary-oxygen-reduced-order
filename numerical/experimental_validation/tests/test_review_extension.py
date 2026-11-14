import sys, unittest, itertools
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
from review_extension import empirical_weights, weighted_quantile, algebraic_jacobian, recover, gas_boundary
from primary_model import solve

class ReviewTests(unittest.TestCase):
    def test_exact_weights_against_ordered_resampling(self):
        counts,w=empirical_weights(3)
        brute={}
        for seq in itertools.product(range(3),repeat=3):
            c=tuple(seq.count(i) for i in range(3));brute[c]=brute.get(c,0)+1
        for c,p in zip(counts,w):self.assertAlmostEqual(p,brute[tuple(c)]/27)
    def test_quantile_definition(self):
        self.assertEqual(weighted_quantile(np.array([3.,1.,2.]),np.array([.1,.7,.2]),.5),1.)
        self.assertEqual(weighted_quantile(np.array([3.,1.,2.]),np.array([.1,.7,.2]),.95),3.)
    def test_jacobian_against_finite_difference(self):
        q,m,s=4.,250.,.15;ce=solve(q,14,100,m,s)['equilibrium_content_ml_dl']
        theta=np.array([np.log(q),np.log(m),np.log(s/(1-s))])
        def f(t):
            out=solve(np.exp(t[0]),14,100,np.exp(t[1]),1/(1+np.exp(-t[2])))
            return np.array([out['Ca_ml_dl']/ce,out['Cv_ml_dl']/ce,t[0]])
        numeric=np.column_stack([(f(theta+np.eye(3)[i]*1e-5)-f(theta-np.eye(3)[i]*1e-5))/2e-5 for i in range(3)])
        np.testing.assert_allclose(numeric,algebraic_jacobian(q,m,s,ce),atol=2e-10,rtol=1e-7)
    def test_recovery_and_scale_invariance(self):
        a=solve(4,14,100,250,.15);b=solve(6,14,100,375,.15)
        self.assertAlmostEqual(a['Pa_mmHg'],b['Pa_mmHg'])
        self.assertAlmostEqual(a['Cv_ml_dl'],b['Cv_ml_dl'])
        s,m=recover(a['Ca_ml_dl'],a['Cv_ml_dl'],a['equilibrium_content_ml_dl'],4)
        self.assertAlmostEqual(s,.15);self.assertAlmostEqual(m,250)
        self.assertAlmostEqual(b['DO2_ml_min']/a['DO2_ml_min'],1.5)
    def test_gas_boundary_does_not_read_outcome(self):
        r=dict(FiO2=.21,PB_mmHg=760.,PH2O_mmHg=47.1,PaCO2_mmHg=30.)
        self.assertAlmostEqual(gas_boundary(r,.8),113.784)
        self.assertEqual(gas_boundary(r,.8),gas_boundary(dict(r,Pa_mmHg=-999,AA_gradient=9999),.8))
    def test_recovery_rejects_nonphysical(self):
        for values in [(20,10,19,4),(15,16,19,4),(19,19,19,4),(15,10,19,-1)]:
            with self.assertRaises(ValueError):recover(*values)

if __name__=='__main__':unittest.main()
