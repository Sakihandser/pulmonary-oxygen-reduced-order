import sys, unittest, math
from pathlib import Path
import numpy as np
from scipy.optimize import linprog,brentq
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import analyze_adequacy as a

class AdequacyTests(unittest.TestCase):
    def test_one_state_random_and_beta_union(self):
        rng=np.random.default_rng(1429)
        for _ in range(400):
            w=rng.dirichlet(np.ones(7));z=rng.uniform(.001,.999,7);s=float(w@z);r=rng.uniform(1,3);b=rng.uniform(0,1)
            lo,hi=a.residual_bounds(s,r,0,1)
            truth=float(w@z**r**(-b))
            self.assertLessEqual(lo,truth+1e-12);self.assertLessEqual(truth,hi+1e-12)

    def test_minimax_is_sharp(self):
        for low,high in [(1,5),(17.2,18.1),(3,3)]:
            center=(low+high)/2
            risk,ok=a.certificate(low,high,center,(high-low)/2)
            self.assertTrue(ok)
            for v in np.linspace(low-2,high+2,251):self.assertGreaterEqual(a.certificate(low,high,v,0)[0]+1e-12,risk)

    def test_closed_balance_independent_root(self):
        for ce,m,q,s in [(18.5,250,5,.2),(17,375,10,.4),(20,0,8,.1)]:
            cv=brentq(lambda cv:10*q*((1-s)*ce+s*cv-cv)-m,-100,ce)
            ca=(1-s)*ce+s*cv
            b=a.content_interval(s,s,q,ce,m)
            self.assertAlmostEqual(ca,b['midpoint'],places=11)

    def test_two_state_attainers(self):
        for s0 in [.05,.2,.5]:
            for s1 in np.linspace(s0,math.sqrt(s0),9):
                ds=a.two_state_distributions(s0,float(s1))
                for w,y in ds:
                    self.assertAlmostEqual(float(w@y),s1,places=12);self.assertAlmostEqual(float(w@(y*y)),s0,places=12)
                for r in [1,1.2,1.5,2,2.5,3]:
                    lo,hi=a.two_state_bounds(s0,float(s1),r)
                    vals=[float(w@y**(2/r)) for w,y in ds]
                    self.assertAlmostEqual(lo,min(vals),places=12);self.assertAlmostEqual(hi,max(vals),places=12)

    def test_two_state_random_mixtures(self):
        rng=np.random.default_rng(20261114)
        for _ in range(600):
            y=rng.uniform(0,1,8);w=rng.dirichlet(np.ones(8));s0=float(w@y**2);s1=float(w@y);r=rng.uniform(1,4)
            lo,hi=a.two_state_bounds(s0,s1,r);truth=float(w@y**(2/r))
            self.assertLessEqual(lo,truth+1e-12);self.assertLessEqual(truth,hi+1e-12)

    def test_independent_linear_program(self):
        # Dense grid is an independent numerical approximation, not a proof.
        y=np.linspace(0,1,10001);s0=.2;s1=.3
        A=np.array([np.ones_like(y),y,y*y]);b=[1,s1,s0]
        for r in [1.25,1.5,2.5,3.]:
            f=y**(2/r);lo,hi=a.two_state_bounds(s0,s1,r)
            lp=linprog(f,A_eq=A,b_eq=b,bounds=(0,None),method='highs')
            hp=linprog(-f,A_eq=A,b_eq=b,bounds=(0,None),method='highs')
            self.assertTrue(lp.success and hp.success)
            self.assertAlmostEqual(lo,lp.fun,delta=2e-6);self.assertAlmostEqual(hi,-hp.fun,delta=2e-6)

    def test_noisy_envelope(self):
        rng=np.random.default_rng(41)
        for _ in range(300):
            y=rng.uniform(0,1,5);w=rng.dirichlet(np.ones(5));s0=float(w@y**2);s1=float(w@y)
            for r in [1,1.5,2,3]:
                lo,hi=a.noisy_two_state_bounds(s0,s1-.005,s1+.005,r)
                exactlo,exacthi=a.two_state_bounds(s0,s1,r)
                self.assertLessEqual(lo,exactlo+1e-12);self.assertLessEqual(exacthi,hi+1e-12)

    def test_general_derivative_with_changing_inputs(self):
        q=5.;ce=18.5;m=250;s=.2;eta=.6;ge=-.3;k=.04
        def ca(t):return ce+ge*t-m*math.exp(eta*t)/(10*q*math.exp(t))*(s+k*t)/(1-s-k*t)
        fd=(ca(1e-5)-ca(-1e-5))/2e-5
        self.assertAlmostEqual(fd,a.general_slope(q,ce,m,s,eta,ge,k),places=8)

    def test_heterogeneous_exponent_direction_bound(self):
        rng=np.random.default_rng(23)
        for _ in range(300):
            w=rng.dirichlet(np.ones(5));z=rng.uniform(.01,.9,5);betas=rng.uniform(0,.7,5);s=float(w@z)
            k=float(w@(betas*(-z*np.log(z))))
            actual=a.general_slope(5,18,250,s,.4,.2,k)
            lower=a.guaranteed_slope(5,18,250,s,.4,.2,.7)
            self.assertGreaterEqual(actual+1e-12,lower)

    def test_feasibility_and_zero_demand(self):
        self.assertEqual(a.content_interval(.6,.8,1,10,100)['status'],'infeasible')
        b=a.content_interval(.1,.8,5,18,0);self.assertEqual(b['width'],0)
        b=a.content_interval(.1,.95,5,18,250);self.assertAlmostEqual(b['low'],5)

    def test_invalid_moments(self):
        for s1 in [.1,.5]:
            with self.assertRaises(ValueError):a.two_state_bounds(.2,s1,2)
        with self.assertRaises(ValueError):a.noisy_two_state_bounds(.2,.6,.7,2)

    def test_exact_and_zero_exponent_collapse(self):
        self.assertEqual(a.residual_bounds(.2,2,0,0),(.2,.2))
        for s0 in [.1,.2,.4]:
            for s1 in [s0,math.sqrt(s0)]:
                lo,hi=a.two_state_bounds(s0,s1,1.5);self.assertAlmostEqual(lo,hi,places=12)

if __name__=='__main__':unittest.main()
